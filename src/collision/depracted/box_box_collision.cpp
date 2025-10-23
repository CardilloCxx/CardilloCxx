#include "box_box_collision.hpp"
#include <algorithm>
#include <array>
#include <cmath>

namespace cardillo::collision::boxbox {

static inline void computeAxes(const ObbCollider& B, std::array<Vector3r,3>& axes) {
    axes[0] = B.R.col(0);
    axes[1] = B.R.col(1);
    axes[2] = B.R.col(2);
}

static inline void projectBox(const ObbCollider& B, const Vector3r& axis, real_t& minp, real_t& maxp) {
    // Center projection
    real_t c = axis.dot(B.center);
    // Radius along axis is sum_i he[i] * |axis dot e_i|
    real_t r = 0;
    for (int i=0;i<3;++i) r += B.halfExtents[i] * std::abs(axis.dot(B.R.col(i)));
    minp = c - r; maxp = c + r;
}

static inline bool satOverlap(const ObbCollider& A, const ObbCollider& B, Vector3r& outNormal, real_t& outDepth) {
    // Test 15 axes: A's 3, B's 3, and 9 cross products
    std::array<Vector3r,3> aAxes, bAxes;
    computeAxes(A, aAxes); computeAxes(B, bAxes);
    real_t minDepth = std::numeric_limits<real_t>::max();
    Vector3r bestAxis = Vector3r::UnitX();
    auto testAxis = [&](const Vector3r& ax){
        Vector3r n = ax.normalized();
        if (!std::isfinite(n.squaredNorm())) return true;
        real_t a0,a1,b0,b1; projectBox(A, n, a0,a1); projectBox(B, n, b0,b1);
        real_t d0 = a1 - b0; real_t d1 = b1 - a0; real_t pen = std::min(d0,d1);
        if (pen < (real_t)0) return false; // separating axis
        if (pen < minDepth) { minDepth = pen; bestAxis = (d0 < d1) ? n : -n; }
        return true;
    };
    // Face axes
    for (int i=0;i<3;++i) if (!testAxis(aAxes[i])) return false;
    for (int i=0;i<3;++i) if (!testAxis(bAxes[i])) return false;
    // Edge-edge axes (cross products)
    for (int i=0;i<3;++i) for (int j=0;j<3;++j) {
        Vector3r ax = aAxes[i].cross(bAxes[j]);
        if (ax.squaredNorm() < (real_t)1e-14) continue;
        if (!testAxis(ax)) return false;
    }
    outNormal = bestAxis; outDepth = minDepth; return true;
}

static inline std::vector<Vector3r> getBoxVertices(const ObbCollider& B) {
    std::array<Vector3r,3> ax; computeAxes(B, ax);
    const Vector3r& e = B.halfExtents;
    std::vector<Vector3r> v; v.reserve(8);
    for (int sx=-1;sx<=1;sx+=2)
    for (int sy=-1;sy<=1;sy+=2)
    for (int sz=-1;sz<=1;sz+=2) {
        Vector3r p = B.center + (real_t)sx * ax[0] * e.x() + (real_t)sy * ax[1] * e.y() + (real_t)sz * ax[2] * e.z();
        v.push_back(p);
    }
    return v;
}

// Support point of an oriented box in world space along direction d
static inline Vector3r supportPointOBB(const ObbCollider& B, const Vector3r& d) {
    Vector3r p = B.center;
    for (int i = 0; i < 3; ++i) {
        const Vector3r ai = B.R.col(i);
        real_t s = (d.dot(ai) >= (real_t)0) ? (real_t)1 : (real_t)-1;
        p += s * B.halfExtents[i] * ai;
    }
    return p;
}

static inline void clipPolyAgainstPlane(std::vector<Vector3r>& poly, const Vector3r& n, real_t o) {
    if (poly.empty()) return;
    std::vector<Vector3r> out; out.reserve(poly.size()+4);
    auto side = [&](const Vector3r& p){ return n.dot(p) - o; };
    for (size_t i=0;i<poly.size();++i) {
        Vector3r a = poly[i]; Vector3r b = poly[(i+1)%poly.size()];
        real_t da = side(a), db = side(b);
        bool ina = da <= 0, inb = db <= 0;
        if (ina && inb) { out.push_back(b); }
        else if (ina && !inb) { real_t t = da/(da-db); out.push_back(a + t*(b-a)); }
        else if (!ina && inb) { real_t t = da/(da-db); out.push_back(a + t*(b-a)); out.push_back(b); }
    }
    poly.swap(out);
}

static inline void buildFacePolygon(const ObbCollider& B, int faceAxis, real_t sgn, std::vector<Vector3r>& poly) {
    // Build the quad face of B with normal sgn*axis
    Matrix33r R = B.R; Vector3r n = R.col(faceAxis) * sgn;
    Vector3r u = R.col((faceAxis+1)%3);
    Vector3r v = R.col((faceAxis+2)%3);
    Vector3r hx = u * B.halfExtents[(faceAxis+1)%3];
    Vector3r hy = v * B.halfExtents[(faceAxis+2)%3];
    Vector3r c = B.center + n * B.halfExtents[faceAxis];
    poly = { c - hx - hy, c + hx - hy, c + hx + hy, c - hx + hy };
}

static inline void emitContact(const ObbCollider& A, const ObbCollider& B, const Vector3r& p, const Vector3r& n, real_t penetration, std::vector<Contact>& out) {
    // Provide translational part here; assembler augments torque terms based on body type.
    MatrixXXr wA(1,3); wA << -n[0], -n[1], -n[2];
    MatrixXXr wB(1,3); wB <<  n[0],  n[1],  n[2];
    Contact c{A.e, B.e, p, n,
              /*pointA_body*/ A.R.transpose() * (p - A.center),
              /*pointB_body*/ B.R.transpose() * (p - B.center),
              /*normalA_body*/ A.R.transpose() * n,
              /*normalB_body*/ B.R.transpose() * n,
              wA, wB, penetration};
    out.push_back(c);
}

void collideOBB(const ObbCollider& A, const ObbCollider& B, std::vector<Contact>& out) {
    Vector3r n; real_t depth;
    if (!satOverlap(A, B, n, depth)) return;

    // Orient SAT normal to point from A to B so that downstream logic can form Ref->Inc consistently
    if ((B.center - A.center).dot(n) < (real_t)0) n = -n;

    // Determine reference and incident boxes based on face-aligned normal
    std::array<Vector3r,3> aAxes, bAxes; computeAxes(A, aAxes); computeAxes(B, bAxes);
    real_t aMax = std::max({std::abs(n.dot(aAxes[0])), std::abs(n.dot(aAxes[1])), std::abs(n.dot(aAxes[2]))});
    real_t bMax = std::max({std::abs(n.dot(bAxes[0])), std::abs(n.dot(bAxes[1])), std::abs(n.dot(bAxes[2]))});
    bool refA = (aMax >= bMax);
    const ObbCollider& Ref = refA ? A : B;
    const ObbCollider& Inc = refA ? B : A;
    Vector3r N = refA ? n : -n; // normal from Ref to Inc

    // Reference face selection (axis most aligned with N)
    int refAxis = 0;
    {
        real_t d0 = std::abs(N.dot(Ref.R.col(0)));
        real_t d1 = std::abs(N.dot(Ref.R.col(1)));
        real_t d2 = std::abs(N.dot(Ref.R.col(2)));
        refAxis = (d0 > d1 && d0 > d2) ? 0 : (d1 > d2 ? 1 : 2);
    }
    real_t refSgn = (N.dot(Ref.R.col(refAxis)) > 0) ? (real_t)1 : (real_t)-1;

    // Incident face selection: choose axis with largest |N·Inc.axis|, then orient opposite to N
    int incAxis = 0; real_t maxAbs = (real_t)-1;
    for (int i=0;i<3;++i) {
        real_t d = N.dot(Inc.R.col(i));
        real_t ad = std::abs(d);
        if (ad > maxAbs) { maxAbs = ad; incAxis = i; }
    }
    real_t incSgn = (N.dot(Inc.R.col(incAxis)) > 0) ? (real_t)-1 : (real_t)1;

    // Build the incident face quad (facing towards Ref)
    std::vector<Vector3r> polyInc;
    buildFacePolygon(Inc, incAxis, incSgn, polyInc);

    // Clip the incident face polygon against the reference face side planes and the reference face plane
    Vector3r u = Ref.R.col((refAxis+1)%3);
    Vector3r v = Ref.R.col((refAxis+2)%3);
    Vector3r fc = Ref.center + Ref.R.col(refAxis) * refSgn * Ref.halfExtents[refAxis];
    real_t hu = Ref.halfExtents[(refAxis+1)%3];
    real_t hv = Ref.halfExtents[(refAxis+2)%3];

    auto clip = [&](std::vector<Vector3r>& poly){
        clipPolyAgainstPlane(poly,  u,  u.dot(fc) + hu);
        clipPolyAgainstPlane(poly, -u, -u.dot(fc) + hu);
        clipPolyAgainstPlane(poly,  v,  v.dot(fc) + hv);
        clipPolyAgainstPlane(poly, -v, -v.dot(fc) + hv);
        clipPolyAgainstPlane(poly,  N,  N.dot(fc)); // keep points behind/ref-side of the face
    };
    
    // Full-face containment check: are all incident face corners inside Ref rectangle and behind its plane?
    bool fullContainment = true;
    {
        std::array<Vector3r,4> incCorners;
        buildFacePolygon(Inc, incAxis, incSgn, polyInc); // polyInc currently holds these corners
        for (const auto& pi : polyInc) {
            Vector3r d = pi - fc;
            real_t du = d.dot(u);
            real_t dv = d.dot(v);
            real_t planeSide = N.dot(pi) - N.dot(fc);
            const real_t eps = (real_t)1e-7;
            if (std::abs(du) > hu + eps || std::abs(dv) > hv + eps || planeSide > eps) {
                fullContainment = false; break;
            }
        }
    }

    if (fullContainment) {
        // Emit the four projections of the incident face corners onto the reference plane.
        // This yields a stable 4-point manifold matching the intruding face footprint.
        for (const auto& pi : polyInc) {
            real_t distToPlane = N.dot(pi) - N.dot(fc);
            Vector3r cp = pi - distToPlane * N;           // project onto Ref face plane
            real_t pen = std::max<real_t>(0, -distToPlane); // penetration along N
            emitContact(refA?A:B, refA?B:A, cp, N, pen, out);
        }
        return;
    }

    clip(polyInc);

    if (polyInc.empty()) {
        // Edge/vertex fallback: project incident support point onto reference face rectangle.
        // Get incident support point most toward Ref along -N
        Vector3r pI = supportPointOBB(Inc, -N);
        // Project onto Ref face plane and clamp to rectangle bounds
        Vector3r d = pI - fc;
        real_t au = d.dot(u);
        real_t av = d.dot(v);
        au = std::clamp(au, -hu, hu);
        av = std::clamp(av, -hv, hv);
        Vector3r pR = fc + au * u + av * v; // closest point on Ref face rectangle
        // Penetration along N from Ref to Inc
        real_t pen = std::max((real_t)0, N.dot(pI - pR));
        // Use the reference-side point for contact position (lies on surface/edge)
        // Emit with a=Ref, b=Inc and normal N (from Ref to Inc)
        emitContact(refA?A:B, refA?B:A, pR, N, pen, out);
        return;
    }

    // Project resulting polygon vertices onto the reference face plane and emit contacts (up to 4)
    const int maxContacts = 4;
    for (size_t i = 0; i < polyInc.size() && (int)i < maxContacts; ++i) {
        Vector3r p = polyInc[i];
        // Project to plane: remove distance along N to the face plane at fc
        real_t distToPlane = N.dot(p) - N.dot(fc);
        Vector3r cp = p - distToPlane * N;
        real_t pen_i = std::max((real_t)0, -distToPlane);
        // Emit with a=Ref, b=Inc and normal N (from Ref to Inc)
        emitContact(refA?A:B, refA?B:A, cp, N, pen_i, out);
    }
}

}
