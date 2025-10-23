#include "CPG.h"
#include <iostream>

namespace cardillo::solver {

cardillo::VectorXr CPG::solve(const cardillo::VectorXr& v_pre, real_t tol)
{
    const auto& W = m_dyn.W();
    const auto& MinvDiag = m_dyn.MinvDiag();
    const int C = (int)W.rows();
    if (C == 0) return v_pre;

    // Scalars as vectors
    cardillo::VectorXr p = cardillo::VectorXr::Zero(C), p_old = p, u = cardillo::VectorXr::Zero(C),
                       q = cardillo::VectorXr::Zero(C), w = cardillo::VectorXr::Zero(C),
                       z = cardillo::VectorXr::Zero(C), b = cardillo::VectorXr::Zero(C),
                       Wq = cardillo::VectorXr::Zero(C), Wp = cardillo::VectorXr::Zero(C);

    // b = W * v_pre
    W_times_v(v_pre, b);
    // Initial u = -(b + G p) with p=0
    u = -b;
    // Initial w and direction q
    for (int cid = 0; cid < C; ++cid) w[cid] = std::max<real_t>((real_t)0, u[cid]);
    q = w;

    size_t k = 0;
    real_t err = std::numeric_limits<real_t>::max();
    while (err > tol && k < 10000) {
        ++k;
        // Wq = G * q = W * (Minv * (W^T * q))
        G_times_p(q, Wq);
        real_t qWq = x_dot_y(q, Wq);
        if (qWq <= (real_t)1e-20) break;
        real_t alpha = x_dot_y(u, q) / qWq;
        p_old = p;
        for (int cid = 0; cid < C; ++cid) p[cid] = std::max<real_t>((real_t)0, p[cid] + alpha * q[cid]);
        // u_{k+1} = -(b + G * p_{k+1})
        G_times_p(p, Wp);
        u = -(b + Wp);
        // w_{k+1}
        for (int cid = 0; cid < C; ++cid) w[cid] = (p[cid] > (real_t)1e-20) ? u[cid] : std::max<real_t>(0, u[cid]);
        // z_{k+1}
        z = q;
        // beta
        real_t beta = x_dot_y(w, Wq) / qWq;
        // q_{k+1}
        q = w + beta * z;
        // error
        err = (p - p_old).norm();
        if (k % 100 == 0) std::cout << "[CPG] iter " << k << ", err = " << err << std::endl;
    }

    // v = v_pre + Minv * W^T * p
    cardillo::VectorXr out = v_pre;
    cardillo::VectorXr tmp = W.transpose() * p;
    out.noalias() += MinvDiag.cwiseProduct(tmp);
    return out;
}

void CPG::G_times_p(const cardillo::VectorXr& p, cardillo::VectorXr& result)
{
    const auto& W = m_dyn.W();
    const auto& MinvDiag = m_dyn.MinvDiag();
    cardillo::VectorXr tmp = W.transpose() * p;
    tmp = MinvDiag.cwiseProduct(tmp);
    result = W * tmp;
}

void CPG::W_times_v(const cardillo::VectorXr& v, cardillo::VectorXr& result)
{
    const auto& W = m_dyn.W();
    result = W * v;
}

real_t CPG::x_dot_y(const cardillo::VectorXr& x, const cardillo::VectorXr& y)
{
    return x.dot(y);
}

CPG::~CPG() = default;

} // namespace cardillo::solver