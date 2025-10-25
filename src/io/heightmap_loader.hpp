#pragma once

#include <string>
#include <vector>
#include <utility>
#include <Eigen/Core>
#include <coal/data_types.h>

namespace cardillo::io {

// Load a single-channel (or RGB) EXR heightmap into a row-major matrix of CoalScalar.
// Returns {rows, cols} in dims and fills heights (rows x cols) where (0,0) is top-left of the image.
// Returns true on success; false on failure (error message printed).
bool load_exr_heightmap(const std::string& path,
                        Eigen::Matrix<coal::CoalScalar, Eigen::Dynamic, Eigen::Dynamic>& heights_out,
                        int& rows,
                        int& cols);

}
