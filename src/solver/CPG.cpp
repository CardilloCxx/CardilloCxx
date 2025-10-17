#include "CPG.h"
#include <iostream>

namespace cardillo::solver {

std::vector<cardillo::VectorXr> CPG::solve(const std::vector<cardillo::VectorXr>& v_pre_blocks, real_t tol)
{
    const auto& Wblocks = m_dyn.WBlocks();
    const auto& MinvBlocks = m_dyn.MinvBlocks();
    const auto& bodyToBlocks = m_dyn.WBlocksFromBodyAll();
    const auto& blockToBody = m_dyn.WBlockToBodyAll();
    const auto& blockToContact = m_dyn.WBlockToContactAll();
    const auto& contactToBlocks = m_dyn.WBlocksFromContactAll();
    const int C = (int)contactToBlocks.size();
    const int Nb = (int)bodyToBlocks.size();

    if (C == 0) return v_pre_blocks;

    // Original CPG variables (1D normal impulses)
    std::vector<cardillo::VectorXr> p(C), p_old(C), u(C), q(C), w(C), z(C), b(C), Wq(C), Wp(C);
    for (int cid = 0; cid < C; ++cid){
        p[(size_t)cid] = cardillo::VectorXr::Zero(1);
        p_old[(size_t)cid] = cardillo::VectorXr::Zero(1);
        u[(size_t)cid] = cardillo::VectorXr::Zero(1);
        q[(size_t)cid] = cardillo::VectorXr::Zero(1);
        b[(size_t)cid] = cardillo::VectorXr::Zero(1);
        w[(size_t)cid] = cardillo::VectorXr::Zero(1);
        z[(size_t)cid] = cardillo::VectorXr::Zero(1);
        Wq[(size_t)cid] = cardillo::VectorXr::Zero(1);
        Wp[(size_t)cid] = cardillo::VectorXr::Zero(1);
    }

    // Per-body temporary buffers
    std::vector<cardillo::VectorXr> temp((size_t)Nb);
    for (int bidx = 0; bidx < Nb; ++bidx) temp[(size_t)bidx] = cardillo::VectorXr::Zero(v_pre_blocks[(size_t)bidx].size());

    // b = W * v_pre
    W_times_v(v_pre_blocks, b);
    // Initial u = -(b + G p) with p=0 to measure violation via positive part
    for (int cid = 0; cid < C; ++cid) u[(size_t)cid][0] = -b[(size_t)cid][0];
    // Initial w and direction q
    for (int cid = 0; cid < C; ++cid) w[(size_t)cid][0] = std::max<real_t>((real_t)0, u[(size_t)cid][0]);
    q = w;

    size_t k = 0;
    real_t err = std::numeric_limits<real_t>::max();
    while (err > tol && k < 10000) {
        ++k;
        // Wq = G * q
        G_times_p(q, temp, Wq);
        real_t qWq = x_dot_y(q, Wq);
        if (qWq <= (real_t)1e-20) break; // safeguard
        real_t alpha = x_dot_y(u, q) / qWq;
        // p_{k+0.5} = p_k + alpha * q_k; project
        p_old = p;
        for (int cid = 0; cid < C; ++cid) p[(size_t)cid][0] = std::max<real_t>((real_t)0, p[(size_t)cid][0] + alpha * q[(size_t)cid][0]);

        // u_{k+1} = -(b + G * p_{k+1}) to keep u = -y
        G_times_p(p, temp, Wp);
        for (int cid = 0; cid < C; ++cid) u[(size_t)cid][0] = -(b[(size_t)cid][0] + Wp[(size_t)cid][0]);

        // w_{k+1} = Prox_T(u_{k+1}) with active-set rule based on p
        for (int cid = 0; cid < C; ++cid) w[(size_t)cid][0] = (p[(size_t)cid][0] > (real_t)1e-20) ? u[cid][0] : std::max<real_t>(0, u[cid][0]);
        //(p[(size_t)cid][0] == 0) ? std::max<real_t>(0, u[cid][0]) : u[cid][0]; // w_{k+1} = Prox_Transpose(u_{k+1})

        // z_{k+1} = Proj_{A(w_{k+1})}(q_k)
        for (int cid = 0; cid < C; ++cid) z[(size_t)cid] = q[(size_t) cid];
        //(w[(size_t)cid][0] > (real_t)1e-20) ? q[(size_t)cid] : cardillo::VectorXr::Zero(1);

        // beta_{k+1}
        real_t beta = x_dot_y(w, Wq) / qWq;
        // q_{k+1} = w_{k+1} + beta z_{k+1}
        for (int cid = 0; cid < C; ++cid) q[(size_t)cid][0] = w[(size_t)cid][0] + beta * z[(size_t)cid][0];

        // error = ||p - p_old||
        err = 0;
        for (int cid = 0; cid < C; ++cid) err += (p[(size_t) cid] - p_old[(size_t) cid]).squaredNorm();
        err = std::sqrt((double)err);
        if (k % 100 == 0) std::cout << "[CPG] iter " << k << ", err = " << err << std::endl;
    }

    // Assemble velocities: v = v_pre + Minv * W^T * p
    std::vector<cardillo::VectorXr> sumW((size_t)Nb);
    for (int bidx = 0; bidx < Nb; ++bidx) sumW[(size_t)bidx] = cardillo::VectorXr::Zero(v_pre_blocks[(size_t)bidx].size());
    for (int bidx = 0; bidx < Nb; ++bidx) {
        for (int idx : bodyToBlocks[(size_t)bidx]) sumW[(size_t)bidx].noalias() += Wblocks[(size_t)idx].transpose() * p[(size_t)blockToContact[(size_t)idx]];
    }
    std::vector<cardillo::VectorXr> out((size_t)Nb);
    for (int bidx = 0; bidx < Nb; ++bidx) out[(size_t)bidx] = v_pre_blocks[(size_t)bidx] + MinvBlocks[(size_t)bidx] * sumW[(size_t)bidx];
    return out;
}

void CPG::G_times_p(const std::vector<cardillo::VectorXr>& p, std::vector<cardillo::VectorXr>& temp, std::vector<cardillo::VectorXr>& result)
{
    const auto& Wblocks = m_dyn.WBlocks();
    const auto& MinvBlocks = m_dyn.MinvBlocks();
    const auto& bodyToBlocks = m_dyn.WBlocksFromBodyAll();
    const auto& blockToBody = m_dyn.WBlockToBodyAll();
    const auto& blockToContact = m_dyn.WBlockToContactAll();
    const auto& contactToBlocks = m_dyn.WBlocksFromContactAll();
    const int C = (int)contactToBlocks.size();
    const int Nb = (int)bodyToBlocks.size();

    // temp = W^T * p per body
    for (int b = 0; b < Nb; ++b) {
        temp[(size_t)b].setZero();
        for (int idx : bodyToBlocks[(size_t)b]) temp[(size_t)b].noalias() += Wblocks[(size_t)idx].transpose() * p[(size_t)blockToContact[(size_t)idx]];
    }
    // temp = Minv * temp
    for (int b = 0; b < Nb; ++b) temp[(size_t)b] = MinvBlocks[(size_t)b] * temp[(size_t)b];
    // result = W * temp per contact
    if ((int)result.size() != C) result.resize((size_t)C);
    for (int cid = 0; cid < C; ++cid) {
        auto [idxA, idxB] = contactToBlocks[(size_t)cid];
    real_t y = 0;
        if (idxA >= 0) y += Wblocks[(size_t)idxA].row(0) * temp[(size_t)blockToBody[(size_t)idxA]];
        if (idxB >= 0) y += Wblocks[(size_t)idxB].row(0) * temp[(size_t)blockToBody[(size_t)idxB]];
        if ((int)result[(size_t)cid].size() != 1) result[(size_t)cid] = cardillo::VectorXr::Zero(1);
        result[(size_t)cid][0] = y;
    }
}

void CPG::W_times_v(const std::vector<cardillo::VectorXr>& v, std::vector<cardillo::VectorXr>& result)
{
    const auto& Wblocks = m_dyn.WBlocks();
    const auto& blockToBody = m_dyn.WBlockToBodyAll();
    const auto& contactToBlocks = m_dyn.WBlocksFromContactAll();
    const int C = (int)contactToBlocks.size();
    if ((int)result.size() != C) result.resize((size_t)C);
    for (int cid = 0; cid < C; ++cid) {
        auto [idxA, idxB] = contactToBlocks[(size_t)cid];
    real_t y = 0;
        if (idxA >= 0) y += Wblocks[(size_t)idxA].row(0) * v[(size_t)blockToBody[(size_t)idxA]];
        if (idxB >= 0) y += Wblocks[(size_t)idxB].row(0) * v[(size_t)blockToBody[(size_t)idxB]];
        if ((int)result[(size_t)cid].size() != 1) result[(size_t)cid] = cardillo::VectorXr::Zero(1);
        result[(size_t)cid][0] = y;
    }
}

real_t CPG::x_dot_y(const std::vector<cardillo::VectorXr>& x, const std::vector<cardillo::VectorXr>& y)
{
    real_t result = 0;
    const int C = (int)x.size();
    for (int cid = 0; cid < C; ++cid) result += x[(size_t)cid].dot(y[(size_t)cid]);
    return result;
}

CPG::~CPG() = default;

} // namespace cardillo::solver