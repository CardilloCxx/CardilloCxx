#ifndef CPG_H
#define CPG_H

#pragma once

#include "../misc/types.hpp"
#include "../physics/dynamics_assembler.hpp"

namespace cardillo::solver {

class CPG
{
public:
    explicit CPG(cardillo::physics::DynamicsAssembler& dyn)
        : m_dyn(dyn) {}

    std::vector<cardillo::VectorXr> solve(const std::vector<cardillo::VectorXr>& v_pre_blocks, real_t tol = (real_t)1e-5);
    void G_times_p(const std::vector<cardillo::VectorXr>& p, std::vector<cardillo::VectorXr>& temp, std::vector<cardillo::VectorXr>& result);
    void W_times_v(const std::vector<cardillo::VectorXr>& v, std::vector<cardillo::VectorXr>& result);
    real_t x_dot_y(const std::vector<cardillo::VectorXr>& x, const std::vector<cardillo::VectorXr>& y);
    ~CPG();

private:
    cardillo::physics::DynamicsAssembler& m_dyn;
};

} // namespace cardillo::solver

#endif