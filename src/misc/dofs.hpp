// Prevent multiple inclusion
#pragma once

#include "types.hpp"


template <index_t size_>
struct Coordinates
{
    static constexpr index_t size = size_;
    Eigen::Matrix<real_t, size, 1> value;
};

template <index_t size_>
struct Index
{
    index_t start;
    static constexpr index_t size = size_;
};