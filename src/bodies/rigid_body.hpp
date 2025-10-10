#pragma once

#include "cardillo.hpp"

#include <Eigen/Dense>

namespace cardillo {

template<typename Scalar_>
class RigidBody {
    private:
        Scalar_ m_;
        Eigen::Matrix<Scalar_, 3, 3> K_Theta_S_;
        Eigen::Matrix<Scalar_, 6, 6> M_;
        Eigen::Matrix<Scalar_, 6, 6> M_inv_;

        Eigen::Vector<Scalar_, 3, 1> position_;
        Eigen::Quaternion<Scalar_> orientation_;

        Eigen::Vector<Scalar_, 3, 1> linear_velocity_;
        Eigen::Vector<Scalar_, 3, 1> angular_velocity_;
};

}