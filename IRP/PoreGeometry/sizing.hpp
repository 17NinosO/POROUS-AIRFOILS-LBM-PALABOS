// sizing.hpp
#pragma once
#include <vector>
#include "config.hpp"

namespace pore_geometry {

std::vector<double> assignPoreRadii(const GeometryConfig& config, std::size_t n_pores);

} // namespace pore_geometry