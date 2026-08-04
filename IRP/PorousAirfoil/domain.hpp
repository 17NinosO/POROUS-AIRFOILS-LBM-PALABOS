// domain.hpp
#pragma once
#include <cstdint>
#include <vector>
#include "config.hpp"

namespace pore_geometry {

class DomainMask {
public:
    DomainMask(int nx, int ny) : nx_(nx), ny_(ny), data_(std::size_t(nx) * ny, 0) {}

    bool at(int row, int col) const { return data_[std::size_t(row) * nx_ + col] != 0; }
    void set(int row, int col, bool value) { data_[std::size_t(row) * nx_ + col] = value ? 1 : 0; }

    int nx() const { return nx_; }
    int ny() const { return ny_; }

    double porosity() const {
        std::size_t count = 0;
        for (auto v : data_) count += (v != 0);
        return double(count) / double(data_.size());
    }

private:
    int nx_, ny_;
    std::vector<std::uint8_t> data_;
};

DomainMask buildDomainMask(const GeometryConfig& config);

}  // namespace pore_geometry