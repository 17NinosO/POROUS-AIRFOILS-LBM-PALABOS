#pragma once
#include <cmath>

inline double naca0012Thickness(double x) {
    if (x < 0.0) x = 0.0;
    if (x > 1.0) x = 1.0;
    return 5.0 * 0.12 * (
        0.2969 * std::sqrt(x)
        - 0.1260 * x
        - 0.3516 * x * x
        + 0.2843 * x * x * x
        - 0.1015 * x * x * x * x
    );
}

inline bool isInsideAirfoil(double x_norm, double y_norm) {
    if (x_norm < 0.0 || x_norm > 1.0) return false;
    double yt = naca0012Thickness(x_norm);
    return std::fabs(y_norm) < yt;
}