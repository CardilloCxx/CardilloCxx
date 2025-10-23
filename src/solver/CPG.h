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

    cardillo::VectorXr solve(const cardillo::VectorXr& v_pre, real_t tol = (real_t)1e-5);
    void G_times_p(const cardillo::VectorXr& p, cardillo::VectorXr& result);
    void W_times_v(const cardillo::VectorXr& v, cardillo::VectorXr& result);
    real_t x_dot_y(const cardillo::VectorXr& x, const cardillo::VectorXr& y);
    ~CPG();

private:
    cardillo::physics::DynamicsAssembler& m_dyn;
};

} // namespace cardillo::solver

#endif