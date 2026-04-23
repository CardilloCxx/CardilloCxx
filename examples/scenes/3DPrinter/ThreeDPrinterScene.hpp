#pragma once

#include "../SceneBase.hpp"

#include <Eigen/Geometry>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include "misc/spline.hpp"

using namespace cardillo;

class ThreeDPrinterScene : public SceneBase {
   public:
    const char* sceneName() const override { return "3DPrinter"; }

    void populate(cardillo::physics::PhysicsEngine& engine) override {
        using namespace cardillo;

        struct ObjectSpec {
            const char* name;
            real_t x;
            real_t y;
            real_t z;
            real_t rxDeg;
            real_t ryDeg;
            real_t rzDeg;
            real_t sx;
            real_t sy;
            real_t sz;
            bool isPin;
            bool isStatic;
            real_t density;
        };

        const std::vector<ObjectSpec> specs = {
            {"BeamLeft", (real_t)0.000000, (real_t)0.385000, (real_t)0.000000, (real_t)0.000000, (real_t)-0.000000, (real_t)0.000000, (real_t)0.400000, (real_t)0.015000, (real_t)0.015000, false,
             false, (real_t)780.0},
            {"BeamRight", (real_t)0.000000, (real_t)-0.385000, (real_t)0.000000, (real_t)0.000000, (real_t)-0.000000, (real_t)0.000000, (real_t)0.400000, (real_t)0.015000, (real_t)0.015000, false,
             false, (real_t)780.0},
            {"MovingBeam", (real_t)0.000000, (real_t)0.000291, (real_t)0.000000, (real_t)0.000000, (real_t)-0.000000, (real_t)0.000000, (real_t)0.015, (real_t)0.36700, (real_t)0.015000, false, false,
             (real_t)780.0},
            {"Gantry", (real_t)0.000000, (real_t)0.000000, (real_t)0.000000, (real_t)0.000000, (real_t)-0.000000, (real_t)0.000000, (real_t)0.050000, (real_t)0.050000, (real_t)0.050000, false, false,
             (real_t)780.0},
            {"MotorR", (real_t)-0.425106, (real_t)-0.384989, (real_t)-0.009735, (real_t)0.000000, (real_t)-0.000000, (real_t)0.000000, (real_t)0.025000, (real_t)0.025000, (real_t)0.025000, false,
             false, (real_t)2780.0},
            {"MotorPinR", (real_t)-0.426000, (real_t)-0.385000, (real_t)0.027423, (real_t)0.000000, (real_t)-0.000000, (real_t)0.000000, (real_t)0.010000, (real_t)0.010000, (real_t)0.016000, true,
             false, (real_t)780.0},
            {"BackPinR", (real_t)0.374000, (real_t)-0.385000, (real_t)0.027423, (real_t)0.000000, (real_t)-0.000000, (real_t)0.000000, (real_t)0.010000, (real_t)0.010000, (real_t)0.016000, true,
             false, (real_t)780.0},
            {"BackPinL", (real_t)0.374000, (real_t)0.385000, (real_t)0.027423, (real_t)0.000000, (real_t)-0.000000, (real_t)0.000000, (real_t)0.010000, (real_t)0.010000, (real_t)0.016000, true, false,
             (real_t)780.0},
            {"MiddlePinL", (real_t)-0.001000, (real_t)0.385000, (real_t)0.027423, (real_t)0.000000, (real_t)-0.000000, (real_t)0.000000, (real_t)0.010000, (real_t)0.010000, (real_t)0.016000, true,
             false, (real_t)780.0},
            {"MiddlePinR", (real_t)-0.001000, (real_t)-0.364, (real_t)0.027423, (real_t)0.000000, (real_t)-0.000000, (real_t)0.000000, (real_t)0.010000, (real_t)0.010000, (real_t)0.016000, true,
             false, (real_t)780.0},
            {"MotorPinR.001", (real_t)-0.426000, (real_t)-0.385000, (real_t)0.017423, (real_t)0.000000, (real_t)-0.000000, (real_t)0.000000, (real_t)0.015000, (real_t)0.015000, (real_t)0.004800, true,
             false, (real_t)780.0},
            {"MiddlePinR.001", (real_t)-0.001000, (real_t)-0.364, (real_t)0.017423, (real_t)0.000000, (real_t)-0.000000, (real_t)0.000000, (real_t)0.015000, (real_t)0.015000, (real_t)0.004800, true,
             false, (real_t)780.0},
            {"MiddlePinL.001", (real_t)-0.001000, (real_t)0.385000, (real_t)0.017423, (real_t)0.000000, (real_t)-0.000000, (real_t)0.000000, (real_t)0.015000, (real_t)0.015000, (real_t)0.004800, true,
             false, (real_t)780.0},
            {"BackPinR.001", (real_t)0.374000, (real_t)-0.385000, (real_t)0.017423, (real_t)0.000000, (real_t)-0.000000, (real_t)0.000000, (real_t)0.015000, (real_t)0.015000, (real_t)0.004800, true,
             false, (real_t)780.0},
            {"BackPinL.001", (real_t)0.374000, (real_t)0.385000, (real_t)0.017423, (real_t)0.000000, (real_t)-0.000000, (real_t)0.000000, (real_t)0.015000, (real_t)0.015000, (real_t)0.004800, true,
             false, (real_t)780.0},
            {"MotorPinR.002", (real_t)-0.426000, (real_t)-0.385000, (real_t)0.042423, (real_t)0.000000, (real_t)-0.000000, (real_t)0.000000, (real_t)0.015000, (real_t)0.015000, (real_t)0.004800, true,
             false, (real_t)780.0},
            {"MiddlePinR.002", (real_t)-0.001000, (real_t)-0.364, (real_t)0.042423, (real_t)0.000000, (real_t)-0.000000, (real_t)0.000000, (real_t)0.015000, (real_t)0.015000, (real_t)0.004800, true,
             false, (real_t)780.0},
            {"MiddlePinL.002", (real_t)-0.001000, (real_t)0.385000, (real_t)0.042423, (real_t)0.000000, (real_t)-0.000000, (real_t)0.000000, (real_t)0.015000, (real_t)0.015000, (real_t)0.004800, true,
             false, (real_t)780.0},
            {"BackPinR.002", (real_t)0.374000, (real_t)-0.385000, (real_t)0.042423, (real_t)0.000000, (real_t)-0.000000, (real_t)0.000000, (real_t)0.015000, (real_t)0.015000, (real_t)0.004800, true,
             false, (real_t)780.0},
            {"BackPinL.002", (real_t)0.374000, (real_t)0.385000, (real_t)0.042423, (real_t)0.000000, (real_t)-0.000000, (real_t)0.000000, (real_t)0.015000, (real_t)0.015000, (real_t)0.004800, true,
             false, (real_t)780.0},
        };

        auto degToRad = [](real_t deg) { return deg * (real_t)M_PI / (real_t)180.0; };
        std::unordered_map<std::string, entt::entity> eByName;
        std::unordered_map<std::string, Vector3r> pByName;
        std::unordered_map<std::string, real_t> rByName;

        for (const auto& o : specs) {
            physics::RigidState rs;
            rs.position = Vector3r(o.x, o.y, o.z);

            const Quaternion4r qx(Eigen::AngleAxis<real_t>(degToRad(o.rxDeg), Vector3r::UnitX()));
            const Quaternion4r qy(Eigen::AngleAxis<real_t>(degToRad(o.ryDeg), Vector3r::UnitY()));
            const Quaternion4r qz(Eigen::AngleAxis<real_t>(degToRad(o.rzDeg), Vector3r::UnitZ()));
            rs.orientation = qz * qy * qx;

            if (o.isPin) {
                const real_t radius = ((real_t)0.5) * (o.sx + o.sy);
                rByName[o.name] = radius;
                const physics::CylinderShape shape(radius, o.sz);
                if (o.isStatic) {
                    eByName[o.name] = engine.addStaticBody(shape, rs);
                } else {
                    physics::RigidProps props = physics::RigidProps::withDensity(o.density);
                    props.friction = (real_t)0.7;
                    eByName[o.name] = engine.addRigidBody(shape, rs, props);
                }
            } else {
                const physics::CubeShape shape(Vector3r(o.sx, o.sy, o.sz));
                if (o.isStatic) {
                    eByName[o.name] = engine.addStaticBody(shape, rs);
                } else {
                    physics::RigidProps props = physics::RigidProps::withDensity(o.density);
                    props.friction = (real_t)0.7;
                    eByName[o.name] = engine.addRigidBody(shape, rs, props);
                }
            }

            pByName[o.name] = rs.position;
            engine.track(eByName[o.name], o.name);
        }

        auto getE = [&](const char* name) { return eByName.at(std::string(name)); };

        // Rigid connectors between each base pin and its .001/.002 counterparts.
        const std::array<std::array<const char*, 2>, 10> rigidPairs = {{
            {{"MotorPinR.001", "MotorPinR"}},
            {{"MotorPinR.002", "MotorPinR"}},
            {{"MiddlePinR.001", "MiddlePinR"}},
            {{"MiddlePinR.002", "MiddlePinR"}},
            {{"MiddlePinL.001", "MiddlePinL"}},
            {{"MiddlePinL.002", "MiddlePinL"}},
            {{"BackPinR.001", "BackPinR"}},
            {{"BackPinR.002", "BackPinR"}},
            {{"BackPinL.001", "BackPinL"}},
            {{"BackPinL.002", "BackPinL"}},
        }};
        for (const auto& p : rigidPairs) engine.addRigidConstraint(getE(p[0]), getE(p[1]));

        // Requested hinge constraints with axis Z.
        m_motor_constraint = engine.addRigidConstraint(getE("MotorR"), getE("MotorPinR"));

        engine.addHingeConstraint(getE("BackPinR"), getE("BeamRight"), physics::JointFrame::fromAxis(Vector3r::Zero(), Vector3r::UnitZ(), getE("BackPinR")));
        engine.addHingeConstraint(getE("BackPinL"), getE("BeamRight"), physics::JointFrame::fromAxis(Vector3r::Zero(), Vector3r::UnitZ(), getE("BackPinL")));
        engine.addHingeConstraint(getE("MiddlePinL"), getE("MovingBeam"), physics::JointFrame::fromAxis(Vector3r::Zero(), Vector3r::UnitZ(), getE("MiddlePinL")));
        engine.addHingeConstraint(getE("MiddlePinR"), getE("MovingBeam"), physics::JointFrame::fromAxis(Vector3r::Zero(), Vector3r::UnitZ(), getE("MiddlePinR")));

        // Rigid attachments between bars/motor.
        engine.addRigidConstraint(getE("MotorR"), getE("BeamRight"));
        engine.addRigidConstraint(getE("BeamLeft"), getE("BeamRight"));
        engine.addRigidConstraint(getE("BeamLeft"));

        // Slider constraints using translation-rotation constraints:
        // - MovingBeam relative to both side beams: free only in X
        // - Gantry relative to MovingBeam: free only in Y
        const real_t inf = std::numeric_limits<real_t>::infinity();
        const Vector3r lockRot = Vector3r::Constant(inf);
        const Vector3r freeX((real_t)0.0, inf, inf);
        const Vector3r freeY(inf, (real_t)0.0, inf);

        engine.addTranslationRotationConstraint(getE("MovingBeam"), getE("BeamLeft"), physics::JointFrame(getE("MovingBeam")), freeX, Vector3r::Zero(), lockRot, Vector3r::Zero());
        engine.addTranslationRotationConstraint(getE("Gantry"), getE("MovingBeam"), physics::JointFrame(getE("Gantry")), freeY, Vector3r::Zero(), lockRot, Vector3r::Zero());

        // Belt path (piecewise linear + circular) in absolute world coordinates.
        const real_t beltZ = ((real_t)0.5) * (pByName.at("MotorPinR.001").z() + pByName.at("MotorPinR.002").z());
        const Vector3r cMotor = pByName.at("MotorPinR");
        const Vector3r cBackR = pByName.at("BackPinR");
        const Vector3r cBackL = pByName.at("BackPinL");
        const Vector3r cMidR = pByName.at("MiddlePinR");
        const Vector3r cMidL = pByName.at("MiddlePinL");
        const real_t pinR = rByName.at("MotorPinR");
        const real_t beltBufferR = (real_t)0.0005;
        const real_t beltR = pinR + beltBufferR;

        const Vector3r cMotorB(cMotor.x(), cMotor.y(), beltZ);
        const Vector3r cBackRB(cBackR.x(), cBackR.y(), beltZ);
        const Vector3r cBackLB(cBackL.x(), cBackL.y(), beltZ);
        const Vector3r cMidRB(cMidR.x(), cMidR.y(), beltZ);
        const Vector3r cMidLB(cMidL.x(), cMidL.y(), beltZ);

        const Vector3r gantryC = pByName.at("Gantry");
        const real_t gantryHalfY = (real_t)0.050000;

        // User convention: left/right along Y, top/bottom along X.
        auto edgeTop = [&](const Vector3r& c) { return Vector3r(c.x() + beltR, c.y(), c.z()); };
        auto edgeBottom = [&](const Vector3r& c) { return Vector3r(c.x() - beltR, c.y(), c.z()); };
        auto edgeLeft = [&](const Vector3r& c) { return Vector3r(c.x(), c.y() + beltR, c.z()); };
        auto edgeRight = [&](const Vector3r& c) { return Vector3r(c.x(), c.y() - beltR, c.z()); };

        const Vector3r mrTop = edgeTop(cMidRB);
        const Vector3r mrRight = edgeRight(cMidRB);

        const Vector3r motorLeft = edgeLeft(cMotorB);
        const Vector3r motorRight = edgeRight(cMotorB);

        const Vector3r backRRight = edgeRight(cBackRB);
        const Vector3r backRTop = edgeTop(cBackRB);

        const Vector3r backLTop = edgeTop(cBackLB);
        const Vector3r backLLeft = edgeLeft(cBackLB);

        const Vector3r midLLeft = edgeLeft(cMidLB);
        const Vector3r midLBottom = edgeBottom(cMidLB);

        // Start at gantry right edge with requested X anchor at moving-right top edge.
        const Vector3r gantryRightAnchor(mrTop.x(), gantryC.y() - gantryHalfY, beltZ);
        // End at gantry left edge; X differs after final winding around middle-left bottom.
        const Vector3r gantryLeftAnchor(midLBottom.x(), gantryC.y() + gantryHalfY, beltZ);

        std::vector<std::unique_ptr<misc::SplinePattern>> beltOwned;
        std::vector<const misc::SplinePattern*> beltSegs;
        auto addLine = [&](const Vector3r& a, const Vector3r& b) {
            if ((b - a).norm() > (real_t)1e-10) {
                beltOwned.push_back(std::make_unique<misc::LinearSpline>(a, b));
                beltSegs.push_back(beltOwned.back().get());
            }
        };
        auto addArc = [&](const Vector3r& c, real_t thetaStart, real_t thetaSpan, bool swapInsideOut) {
            Vector3r n = Vector3r::UnitZ();
            if (swapInsideOut) n = -n;
            // Negating angles while flipping the circle normal preserves geometric path
            // endpoints but inverts local frame orientation (inside-out).
            const real_t t0 = swapInsideOut ? -thetaStart : thetaStart;
            const real_t ts = swapInsideOut ? -thetaSpan : thetaSpan;
            beltOwned.push_back(std::make_unique<misc::CircleSpline>(c, beltR, n, Vector3r::UnitX(), t0, ts));
            beltSegs.push_back(beltOwned.back().get());
        };
        auto addAxisAligned = [&](const Vector3r& a, const Vector3r& b, bool xFirst) {
            if ((b - a).norm() <= (real_t)1e-10) return;
            if (xFirst) {
                const Vector3r mid(b.x(), a.y(), a.z());
                addLine(a, mid);
                addLine(mid, b);
            } else {
                const Vector3r mid(a.x(), b.y(), a.z());
                addLine(a, mid);
                addLine(mid, b);
            }
        };

        // MiddleRight: top -> right
        // MotorR: left -> right around bottom
        // BackR: right -> top
        // BackL: top -> left
        // MiddleL: left -> bottom
        addAxisAligned(gantryRightAnchor, mrTop, false);
        addArc(cMidRB, (real_t)0.0, (real_t)-M_PI_2, false);

        addAxisAligned(mrRight, motorLeft, true);
        addArc(cMotorB, (real_t)M_PI_2, (real_t)M_PI, true);

        addAxisAligned(motorRight, backRRight, false);
        addArc(cBackRB, (real_t)-M_PI_2, (real_t)M_PI_2, true);

        addAxisAligned(backRTop, backLTop, true);
        addArc(cBackLB, (real_t)0.0, (real_t)M_PI_2, true);

        addAxisAligned(backLLeft, midLLeft, false);
        addArc(cMidLB, (real_t)M_PI_2, (real_t)M_PI_2, true);

        addAxisAligned(midLBottom, gantryLeftAnchor, false);

        const physics::BeamCrossSection beltSection((real_t)0.0005, (real_t)0.01, physics::BeamBodyType::Cube);
        physics::BeamSpringParams beltSprings = physics::BeamSpringParams::fromMaterial((real_t)1.0e9, (real_t)0.30, (real_t)400, (real_t)40, (real_t)0.4, (real_t)0.4, (real_t)5.0);
        beltSprings.gamma0 = Vector3r::Zero();
        beltSprings.gamma0->x() = (real_t)-1.0e-4;
        beltSprings.kappa0 = Vector3r::Zero();
        const physics::RigidProps beltProps = physics::RigidProps::withDensity((real_t)1100.0);
        const physics::RigidState beltStateDefaults(Vector3r::Zero(), Vector3r::Zero(), Quaternion4r::Identity(), Vector3r::Zero());

        auto beltEnds = engine.createBeams(beltSegs, beltSection, beltSprings, beltStateDefaults, beltProps, (size_t)1040);
        if (beltEnds.first != entt::null) {
            engine.track(beltEnds.first, "BeltStart");
            engine.addRigidConstraint(beltEnds.first, getE("Gantry"));
            engine.disableCollisionBetween(beltEnds.first, getE("Gantry"));
        }
        if (beltEnds.second != entt::null) {
            engine.track(beltEnds.second, "BeltEnd");
            engine.addRigidConstraint(beltEnds.second, getE("Gantry"));
            engine.disableCollisionBetween(beltEnds.second, getE("Gantry"));
        }

        // Collision disable between all parts
        for (size_t i = 0; i < specs.size(); ++i) {
            for (size_t j = i + 1; j < specs.size(); ++j) {
                engine.disableCollisionBetween(getE(specs[i].name), getE(specs[j].name));
            }
        }
    }

    void updateScene(cardillo::physics::PhysicsEngine& engine, real_t t, real_t /*dt*/) override {
        const real_t motorSpeed = M_2_PI * 100.0 * std::sin(M_2_PI * (real_t)10.0 * t);
        engine.setConstraintAngularVelocity(m_motor_constraint, Vector3r(0, 0, motorSpeed));
    }

   private:
    index_t m_motor_constraint{0};
};
