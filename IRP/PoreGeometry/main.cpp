// main.cpp
#include <iostream>
#include <random>
#include "config.hpp"
#include "pore_shape.hpp"

int main() {
    std::mt19937_64 rng(42);

    auto circle = pore_geometry::makeOrganicPore({0, 0}, 10.0, /*n_harmonics=*/0, /*roughness=*/3.0, rng);
    auto organic = pore_geometry::makeOrganicPore({0, 0}, 10.0, /*n_harmonics=*/4, /*roughness=*/3.0, rng);

    std::cout << "n_harmonics=0 (should be exactly 10 everywhere):\n";
    for (double theta = 0; theta < 6.29; theta += 1.57)
        std::cout << "  theta=" << theta << " r=" << circle.radiusAt(theta) << "\n";

    std::cout << "n_harmonics=4 (should vary around 10):\n";
    for (double theta = 0; theta < 6.29; theta += 1.57)
        std::cout << "  theta=" << theta << " r=" << organic.radiusAt(theta) << "\n";

    return 0;
}