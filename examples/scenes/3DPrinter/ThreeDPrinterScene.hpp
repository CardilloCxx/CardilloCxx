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

    void populate(physics::PhysicsEngine& engine) override {
        using namespace cardillo;

        m_motor_constraints.clear();

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
        auto rotX180Pos = [](const Vector3r& p) { return Vector3r(p.x(), -p.y(), -p.z()); };

        std::unordered_map<std::string, entt::entity> eByName;
        std::unordered_map<std::string, Vector3r> pByName;
        std::unordered_map<std::string, real_t> rByName;
        std::unordered_map<std::string, const ObjectSpec*> specByName;
        for (const auto& o : specs) specByName[o.name] = &o;

        auto addObject = [&](const ObjectSpec& o, const std::string& finalName, bool applyX180) {
            physics::RigidState rs;
            rs.position = Vector3r(o.x, o.y, o.z);

            const Quaternion4r qx(Eigen::AngleAxis<real_t>(degToRad(o.rxDeg), Vector3r::UnitX()));
            const Quaternion4r qy(Eigen::AngleAxis<real_t>(degToRad(o.ryDeg), Vector3r::UnitY()));
            const Quaternion4r qz(Eigen::AngleAxis<real_t>(degToRad(o.rzDeg), Vector3r::UnitZ()));
            rs.orientation = qz * qy * qx;

            if (applyX180) {
                const Quaternion4r qRotX180(Eigen::AngleAxis<real_t>((real_t)M_PI, Vector3r::UnitX()));
                rs.position = rotX180Pos(rs.position);
                rs.orientation = qRotX180 * rs.orientation;
            }

            if (o.isPin) {
                const real_t radius = ((real_t)0.5) * (o.sx + o.sy);
                rByName[finalName] = radius;
                const physics::CylinderShape shape(radius, o.sz);
                if (o.isStatic) {
                    eByName[finalName] = engine.addStaticBody(shape, rs);
                } else {
                    physics::RigidProps props = physics::RigidProps::withDensity(o.density);
                    props.friction = (real_t)0.7;
                    eByName[finalName] = engine.addRigidBody(shape, rs, props);
                }
            } else {
                const physics::CubeShape shape(Vector3r(o.sx, o.sy, o.sz));
                if (o.isStatic) {
                    eByName[finalName] = engine.addStaticBody(shape, rs);
                } else {
                    physics::RigidProps props = physics::RigidProps::withDensity(o.density);
                    props.friction = (real_t)0.7;
                    eByName[finalName] = engine.addRigidBody(shape, rs, props);
                }
            }

            pByName[finalName] = rs.position;
        };

        // Build base scene objects from specs.
        for (const auto& o : specs) addObject(o, o.name, false);

        // Add mirrored CoreXY side: motor + pins only, rotated around world X by 180deg.
        const std::vector<std::string> mirroredCoreXY = {"MotorR",         "MotorPinR",    "BackPinR",     "BackPinL",      "MiddlePinL",     "MiddlePinR",     "MotorPinR.001", "MiddlePinR.001",
                                                         "MiddlePinL.001", "BackPinR.001", "BackPinL.001", "MotorPinR.002", "MiddlePinR.002", "MiddlePinL.002", "BackPinR.002",  "BackPinL.002"};
        for (const auto& baseName : mirroredCoreXY) {
            addObject(*specByName.at(baseName), baseName + "_X180", true);
        }

        auto getE = [&](const std::string& name) { return eByName.at(name); };

        struct AssemblyNames {
            std::string motor;
            std::string motorPin;
            std::string backPinR;
            std::string backPinL;
            std::string middlePinR;
            std::string middlePinL;
            std::string motorPin001;
            std::string motorPin002;
            std::string middlePinR001;
            std::string middlePinR002;
            std::string middlePinL001;
            std::string middlePinL002;
            std::string backPinR001;
            std::string backPinR002;
            std::string backPinL001;
            std::string backPinL002;
            std::string beltStartTrack;
            std::string beltEndTrack;
        };

        auto makeAssembly = [](const std::string& suffix) {
            return AssemblyNames{"MotorR" + suffix,        "MotorPinR" + suffix,     "BackPinR" + suffix,       "BackPinL" + suffix,       "MiddlePinR" + suffix,     "MiddlePinL" + suffix,
                                 "MotorPinR.001" + suffix, "MotorPinR.002" + suffix, "MiddlePinR.001" + suffix, "MiddlePinR.002" + suffix, "MiddlePinL.001" + suffix, "MiddlePinL.002" + suffix,
                                 "BackPinR.001" + suffix,  "BackPinR.002" + suffix,  "BackPinL.001" + suffix,   "BackPinL.002" + suffix,   "BeltStart" + suffix,      "BeltEnd" + suffix};
        };

        const real_t inf = std::numeric_limits<real_t>::infinity();
        const Vector3r lockRot = Vector3r::Constant(inf);
        const Vector3r freeX((real_t)0.0, inf, inf);
        const Vector3r freeY(inf, (real_t)0.0, inf);

        engine.addTranslationRotationConstraint(getE("MovingBeam"), getE("BeamLeft"), physics::JointFrame(getE("MovingBeam")), freeX, Vector3r::Zero(), lockRot, Vector3r::Zero());
        engine.addTranslationRotationConstraint(getE("Gantry"), getE("MovingBeam"), physics::JointFrame(getE("Gantry")), freeY, Vector3r::Zero(), lockRot, Vector3r::Zero());
        engine.addRigidConstraint(getE("BeamLeft"), getE("BeamRight"));

        auto attachAssembly = [&](const AssemblyNames& a, bool flippedX180) {
            const std::array<std::array<std::string, 2>, 10> rigidPairs = {{
                {{a.motorPin001, a.motorPin}},
                {{a.motorPin002, a.motorPin}},
                {{a.middlePinR001, a.middlePinR}},
                {{a.middlePinR002, a.middlePinR}},
                {{a.middlePinL001, a.middlePinL}},
                {{a.middlePinL002, a.middlePinL}},
                {{a.backPinR001, a.backPinR}},
                {{a.backPinR002, a.backPinR}},
                {{a.backPinL001, a.backPinL}},
                {{a.backPinL002, a.backPinL}},
            }};
            for (const auto& p : rigidPairs) engine.addRigidConstraint(getE(p[0]), getE(p[1]));

            m_motor_constraints.push_back(engine.addRigidConstraint(getE(a.motorPin), getE(a.motor)));

            auto beamByY = [&](const std::string& pinName) {
                const bool worldLeft = pByName.at(pinName).y() >= (real_t)0.0;
                if (!flippedX180) return worldLeft ? std::string("BeamLeft") : std::string("BeamRight");
                // Under 180deg rotation around X, local left/right is mirrored in world Y.
                return worldLeft ? std::string("BeamRight") : std::string("BeamLeft");
            };

            engine.addHingeConstraint(getE(a.backPinR), getE(beamByY(a.backPinR)), physics::JointFrame::fromAxis(Vector3r::Zero(), Vector3r::UnitZ(), getE(a.backPinR)));
            engine.addHingeConstraint(getE(a.backPinL), getE(beamByY(a.backPinL)), physics::JointFrame::fromAxis(Vector3r::Zero(), Vector3r::UnitZ(), getE(a.backPinL)));
            engine.addHingeConstraint(getE(a.middlePinL), getE("MovingBeam"), physics::JointFrame::fromAxis(Vector3r::Zero(), Vector3r::UnitZ(), getE(a.middlePinL)));
            engine.addHingeConstraint(getE(a.middlePinR), getE("MovingBeam"), physics::JointFrame::fromAxis(Vector3r::Zero(), Vector3r::UnitZ(), getE(a.middlePinR)));

            engine.addRigidConstraint(getE(a.motor), getE(beamByY(a.motor)));

            const real_t beltZ = ((real_t)0.5) * (pByName.at(a.motorPin001).z() + pByName.at(a.motorPin002).z());
            const Vector3r cMotor = pByName.at(a.motorPin);
            const Vector3r cBackR = pByName.at(a.backPinR);
            const Vector3r cBackL = pByName.at(a.backPinL);
            const Vector3r cMidR = pByName.at(a.middlePinR);
            const Vector3r cMidL = pByName.at(a.middlePinL);
            const real_t pinR = rByName.at(a.motorPin);
            const real_t beltBufferR = 0.0005;  //(real_t)0.0005;
            const real_t beltR = pinR + beltBufferR;

            const Vector3r cMotorB(cMotor.x(), cMotor.y(), beltZ);
            const Vector3r cBackRB(cBackR.x(), cBackR.y(), beltZ);
            const Vector3r cBackLB(cBackL.x(), cBackL.y(), beltZ);
            const Vector3r cMidRB(cMidR.x(), cMidR.y(), beltZ);
            const Vector3r cMidLB(cMidL.x(), cMidL.y(), beltZ);

            const Vector3r gantryC = pByName.at("Gantry");
            const real_t gantryHalfY = (real_t)0.050000;

            auto edgeTop = [&](const Vector3r& c) { return Vector3r(c.x() + beltR, c.y(), c.z()); };
            auto edgeBottom = [&](const Vector3r& c) { return Vector3r(c.x() - beltR, c.y(), c.z()); };
            auto edgeLeft = [&](const Vector3r& c) { return !flippedX180 ? Vector3r(c.x(), c.y() + beltR, c.z()) : Vector3r(c.x(), c.y() - beltR, c.z()); };
            auto edgeRight = [&](const Vector3r& c) { return !flippedX180 ? Vector3r(c.x(), c.y() - beltR, c.z()) : Vector3r(c.x(), c.y() + beltR, c.z()); };

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

            const Vector3r gantryRightAnchor(mrTop.x(), !flippedX180 ? (gantryC.y() - gantryHalfY) : (gantryC.y() + gantryHalfY), beltZ);
            const Vector3r gantryLeftAnchor(midLBottom.x(), !flippedX180 ? (gantryC.y() + gantryHalfY) : (gantryC.y() - gantryHalfY), beltZ);

            std::vector<std::unique_ptr<misc::SplinePattern>> beltOwned;
            std::vector<const misc::SplinePattern*> beltSegs;
            auto addLine = [&](const Vector3r& a0, const Vector3r& b0) {
                if ((b0 - a0).norm() > (real_t)1e-10) {
                    beltOwned.push_back(std::make_unique<misc::LinearSpline>(a0, b0));
                    beltSegs.push_back(beltOwned.back().get());
                }
            };
            auto addArc = [&](const Vector3r& c, real_t thetaStart, real_t thetaSpan, bool swapInsideOut) {
                // Keep the rotated side consistent by flipping the base belt normal.
                Vector3r n = Vector3r::UnitZ();
                if (flippedX180) n = -n;
                if (swapInsideOut) n = -n;
                const real_t t0 = swapInsideOut ? -thetaStart : thetaStart;
                const real_t ts = swapInsideOut ? -thetaSpan : thetaSpan;
                beltOwned.push_back(std::make_unique<misc::CircleSpline>(c, beltR, n, Vector3r::UnitX(), t0, ts));
                beltSegs.push_back(beltOwned.back().get());
            };
            auto addAxisAligned = [&](const Vector3r& a0, const Vector3r& b0, bool xFirst) {
                if ((b0 - a0).norm() <= (real_t)1e-10) return;
                if (xFirst) {
                    const Vector3r mid(b0.x(), a0.y(), a0.z());
                    addLine(a0, mid);
                    addLine(mid, b0);
                } else {
                    const Vector3r mid(a0.x(), b0.y(), a0.z());
                    addLine(a0, mid);
                    addLine(mid, b0);
                }
            };

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
            physics::BeamSpringParams beltSprings = physics::BeamSpringParams::fromMaterial((real_t)1.0e9, (real_t)0.30, 1e6, (real_t)4000, (real_t)0.4, (real_t)0.4, (real_t)0.5);
            beltSprings.gamma0 = Vector3r::Zero();
            beltSprings.gamma0->x() = (real_t)-1.0e-4;
            beltSprings.kappa0 = Vector3r::Zero();
            const physics::RigidProps beltProps = physics::RigidProps::withDensity((real_t)1100.0);
            const physics::RigidState beltStateDefaults(Vector3r::Zero(), Vector3r::Zero(), Quaternion4r::Identity(), Vector3r::Zero());

            auto beltEnds = engine.createBeams(beltSegs, beltSection, beltSprings, beltStateDefaults, beltProps, (size_t)1040);
            if (beltEnds.first != entt::null) {
                engine.addRigidConstraint(beltEnds.first, getE("Gantry"));
                engine.disableCollisionBetween(beltEnds.first, getE("Gantry"));
            }
            if (beltEnds.second != entt::null) {
                engine.addRigidConstraint(beltEnds.second, getE("Gantry"));
                engine.disableCollisionBetween(beltEnds.second, getE("Gantry"));
            }
        };

        attachAssembly(makeAssembly(""), false);
        attachAssembly(makeAssembly("_X180"), true);

        // Attach Printer frame to the world so it doesnt fall.
        engine.addRigidConstraint(getE("BeamLeft"));

        // Collision disable between all created parts.
        std::vector<std::string> allNames;
        allNames.reserve(eByName.size());
        for (const auto& kv : eByName) allNames.push_back(kv.first);
        for (size_t i = 0; i < allNames.size(); ++i) {
            for (size_t j = i + 1; j < allNames.size(); ++j) {
                engine.disableCollisionBetween(getE(allNames[i]), getE(allNames[j]));
            }
        }

        engine.track(getE("Gantry"), "Gantry");
    }

    void updateScene(physics::PhysicsEngine& engine, real_t t, real_t /*dt*/) override {
        auto right_motor = m_motor_constraints[0];
        auto left_motor = m_motor_constraints[1];

        real_t start_ramp = std::min(1.0, t / (real_t)0.1);
        real_t speed = 10.0 * (real_t)M_PI * start_ramp;
        real_t vx = speed * std::sin(M_PI * 2.0 * t);
        real_t vy = speed * std::cos(M_PI * 2.0 * t);
        real_t targetR = vx + vy;
        real_t targetL = -vy + vx;

        engine.setConstraintAngularVelocity(right_motor, Vector3r(0, 0, targetR));
        engine.setConstraintAngularVelocity(left_motor, Vector3r(0, 0, targetL));
    }

   private:
    std::vector<index_t> m_motor_constraints;
};
