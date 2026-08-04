// pgm_writer.hpp
#pragma once
#include <fstream>
#include <stdexcept>
#include <string>
#include "domain.hpp"

namespace pore_geometry {

inline void writePGM(const DomainMask& mask, const std::string& path) {
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("Could not open " + path + " for writing");

    out << "P5\n" << mask.nx() << " " << mask.ny() << "\n255\n";

    // flip vertically: mask row 0 is y=0 (bottom), PGM row 0 is the top
    for (int row = mask.ny() - 1; row >= 0; --row) {
        for (int col = 0; col < mask.nx(); ++col) {
            unsigned char value = mask.at(row, col) ? 255 : 0;  // fluid=white, solid=black
            out.put(static_cast<char>(value));
        }
    }
}

}  // namespace pore_geometry