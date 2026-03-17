/********************************************************************************
 * File: numericutils.hpp
 *
 * Utility functions for safe and consistent floating-point operations.
 *
 * Provides helper functions for:
 *   - Floating-point comparison with combined relative and absolute tolerance.
 *   - Portable numeric constants (epsilon and tolerance).
 *   - Special mathematical functions (e.g. digamma ψ).
 *
 * Intended for scientific computation where double precision stability matters.
 *
 * Author: Wenbo Lyu (Github: @SpatLyu)
 * License: GPL-3
 ********************************************************************************/

#ifndef NUMERICUTILS_HPP
#define NUMERICUTILS_HPP

#include <cmath>
#include <algorithm>
#include <limits>
#include <initializer_list>

namespace NumericUtils
{

    // ==============================
    // Common numeric constants
    // ==============================

    constexpr double DOUBLE_EPS     = std::numeric_limits<double>::epsilon(); // ≈ 2.22e-16
    constexpr double DOUBLE_TOL_ABS = 1.5e-16;  // Absolute tolerance
    constexpr double DOUBLE_TOL_REL = 1.5e-8;   // Relative tolerance

    // ==============================
    // Floating-point comparison
    // ==============================

    /**
     * @brief Compare two double values with combined relative and absolute tolerance.
     *
     * Implements a numerically stable test for near equality:
     *
     * |x - y| <= max(rel_tol * max(|x|, |y|, 1.0), abs_tol)
     *
     * @param x First value
     * @param y Second value
     * @param rel_tol Relative tolerance (default DOUBLE_TOL_REL)
     * @param abs_tol Absolute tolerance (default DOUBLE_TOL_ABS)
     * @return true if x and y are considered equal within tolerance
     */
    inline bool doubleNearlyEqual(double x,
                                  double y,
                                  double rel_tol = DOUBLE_TOL_REL,
                                  double abs_tol = DOUBLE_TOL_ABS) noexcept
    {
        double diff  = std::fabs(x - y);
        double scale = std::max({1.0, std::fabs(x), std::fabs(y)});
        return diff <= std::max(rel_tol * scale, abs_tol);
    }

    // ==============================
    // Special functions
    // ==============================

    /**
     * @brief Computes an approximation of the digamma function ψ(x).
     *
     * The digamma function is defined as the derivative of the logarithm
     * of the gamma function:
     *
     *      ψ(x) = d/dx log Γ(x)
     *
     * Implementation strategy:
     *
     * 1. Recurrence relation
     *
     *      ψ(x) = ψ(x + 1) - 1/x
     *
     *    Used to shift small x values upward to improve numerical stability.
     *
     * 2. Asymptotic expansion (for large x)
     *
     *      ψ(x) ≈ log(x) - 1/(2x)
     *              - 1/(12x²)
     *              + 1/(120x⁴)
     *              - 1/(252x⁶)
     *              + ...
     *
     * This combination provides a fast and sufficiently accurate
     * approximation for most statistical and scientific applications.
     *
     * @param x Input value
     * @return Approximation of ψ(x)
     */
    inline double Digamma(double x) noexcept
    {
        double result = 0.0;

        // Shift small x upward using recurrence
        while (x <= 5.0)
        {
            result -= 1.0 / x;
            x += 1.0;
        }

        double inv_x  = 1.0 / x;
        double inv_x2 = inv_x * inv_x;

        // Asymptotic expansion
        double series =
            inv_x2 * (-1.0 / 12.0 +
            inv_x2 * ( 1.0 / 120.0 +
            inv_x2 * (-1.0 / 252.0 +
            inv_x2 * ( 1.0 / 240.0 +
            inv_x2 * (-1.0 / 132.0 +
            inv_x2 * ( 691.0 / 32760.0 +
            inv_x2 * (-1.0 / 12.0 +
                      3617.0 / 8160.0)))))));

        return result + std::log(x) - 0.5 * inv_x + series;
    }

} // namespace NumericUtils

#endif // NUMERICUTILS_HPP
