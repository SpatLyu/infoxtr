/********************************************************************************************************
 * File: symdync.hpp
 *
 * High-performance symbolic dynamics utilities for pattern-based causal analysis.
 *
 * Provides lightweight functions for:
 *   - Transforming continuous state space data into discrete symbolic patterns.
 *   - Computing signature spaces via absolute or relative successive differences.
 *   - Encoding symbolic patterns with compact uint8 representation.
 *   - Quantifying sign agreement proportions between pattern sequences.
 *   - Constructing causal heatmaps and classifying interactions as positive, negative, dark, or null.
 *
 * Designed for large-scale, high-frequency time series and spatial cross-sectional data
 * where memory efficiency, deterministic indexing, and cache-friendly access patterns are critical.
 *
 * Author: Wenbo Lyu (Github: @SpatLyu)
 * License: GPL-3
 *******************************************************************************************************/

#ifndef INFOXTR_SYMDYNC_HPP
#define INFOXTR_SYMDYNC_HPP

#include <vector>
#include <cmath>
#include <limits>
#include <numeric>
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <cstdint>
#include <iterator>
#include <functional> 
#include <unordered_map>
#include <unordered_set>
#include "infoxtr/numericutils.hpp"

namespace infoxtr
{

namespace symdync
{
    /**
     * Generates a symbolic pattern representation from a continuous state space matrix.
     *
     * This function combines signature space computation and pattern encoding:
     *   1. Compute differences between successive columns for each row:
     *      - relative = true  : (x[i+1] - x[i]) / x[i]
     *      - relative = false : x[i+1] - x[i]
     *   2. Map the resulting values to discrete symbols (uint8_t):
     *        0 -> NA / undefined (NaN)
     *        1 -> negative change (value < 0)
     *        2 -> no change       (value == 0)
     *        3 -> positive change (value > 0)
     *   3. Optional na_rm behavior:
     *        - na_rm = true  : rows with any NaN are replaced by {0}
     *        - na_rm = false : NaN values encoded as 0, row otherwise kept
     *
     * @param mat       Input state space matrix [n_rows x n_cols].
     * @param relative  If true, compute relative changes; else absolute changes.
     * @param na_rm     Whether to remove rows containing NaN (default: true).
     * @return          Symbolic pattern matrix [n_rows x (n_cols - 1)], each row is a uint8_t vector.
     * @throws std::invalid_argument if input is empty or has fewer than 2 columns.
     */
    inline std::vector<std::vector<uint8_t>> sympat(
        const std::vector<std::vector<double>>& mat,
        bool relative = false,
        bool na_rm = true
    ) {
        if (mat.empty()) {
            throw std::invalid_argument("Input matrix must not be empty.");
        }

        const size_t n_rows = mat.size();
        const size_t n_cols = mat[0].size();

        if (n_cols < 2) {
            throw std::invalid_argument("State space matrix must have at least 2 columns.");
        }

        const size_t out_cols = n_cols - 1;
        std::vector<std::vector<uint8_t>> patterns;
        patterns.reserve(n_rows);

        for (size_t i = 0; i < n_rows; ++i) {
            const auto& row = mat[i];
            std::vector<uint8_t> pat;
            pat.reserve(out_cols);

            bool has_nan = false;

            for (size_t j = 0; j < out_cols; ++j) {
                double diff = row[j + 1] - row[j];

                // Compute relative change if requested
                if (relative && !std::isnan(diff) && !infoxtr::numericutils::doubleNearlyEqual(row[j], 0.0)) {
                    diff /= row[j];
                }

                if (std::isnan(diff)) {
                    pat.push_back(static_cast<uint8_t>(0));
                    has_nan = true;
                } 
                else if (infoxtr::numericutils::doubleNearlyEqual(diff, 0.0)) {
                    pat.push_back(static_cast<uint8_t>(2));
                } 
                else if (diff > 0.0) {
                    pat.push_back(static_cast<uint8_t>(3));
                } 
                else { // diff < 0
                    pat.push_back(static_cast<uint8_t>(1));
                }
            }

            // Handle row-level NA removal
            if (na_rm && has_nan) {
                patterns.emplace_back(std::vector<uint8_t>{0});
            } else {
                patterns.emplace_back(std::move(pat));
            }
        }

        return patterns;
    }

    #include <algorithm>

    /**
     * Encodes embeddings into compact identifiers using 
     * lexicographic ordering  and deduplication.
     *
     * Encoding scheme:
     *   - {0}                  -> 0   (reserved for NA / undefined patterns)
     *   - unique pattern #k    -> k   (k ∈ [1, n_unique], lexicographic order)
     */
    inline std::vector<uint64_t> symbolize(
        const std::vector<std::vector<double>>& mat,
        bool relative = false,
        bool na_rm = true
    ) {
        auto patterns = sympat(mat, relative, na_rm);

        std::vector<uint64_t> labels(patterns.size(), 0);

        // Step 1: collect non-{0} patterns
        std::vector<std::vector<uint8_t>> unique_patterns;
        unique_patterns.reserve(patterns.size());

        for (const auto& pat : patterns) {
            if (!(pat.size() == 1 && pat[0] == 0)) {
                unique_patterns.push_back(pat);
            }
        }

        // Step 2: sort
        std::sort(unique_patterns.begin(), unique_patterns.end());

        // Step 3: unique
        unique_patterns.erase(
            std::unique(unique_patterns.begin(), unique_patterns.end()),
            unique_patterns.end()
        );

        // Step 4: assign IDs via binary search
        for (size_t i = 0; i < patterns.size(); ++i) {
            const auto& pat = patterns[i];

            // special case
            if (pat.size() == 1 && pat[0] == 0) {
                labels[i] = 0;
                continue;
            }

            auto it = std::lower_bound(
                unique_patterns.begin(),
                unique_patterns.end(),
                pat
            );

            labels[i] = static_cast<uint64_t>(
                std::distance(unique_patterns.begin(), it) + 1
            );
        }

        return labels;
    }

    /**
     * Compute sign agreement proportions between two pattern spaces.
     *
     * This function compares two symbolic pattern spaces and evaluates the 
     * proportion of positive and negative sign agreements.
     *
     * Comparison logic:
     *
     *   Only positions where both symbols are non-zero are counted.
     *
     *   Let pat1[i][j] and pat2[i][j] be compared.
     *
     *   Valid symbols:
     *       1  -> negative
     *       2  -> stable
     *       3  -> positive
     *
     *   Positive agreement:
     *       (1,1), (2,2), (3,3)
     *
     *   Negative agreement:
     *       (1,3), (3,1)
     *
     *   Other combinations are ignored.
     *
     * Output:
     *
     *   Returns vector<double> of size 2:
     *       result[0] = positive_ratio
     *       result[1] = negative_ratio
     *
     *   Ratios are computed as:
     *
     *       pos_tot / pat_tot
     *       neg_tot / pat_tot
     *
     *   where pat_tot is the total number of valid comparisons.
     *
     *   If pat_tot == 0:
     *       both outputs are NaN.
     *
     * @param pat1 First pattern space.
     * @param pat2 Second pattern space.
     * @return     Vector containing positive and negative proportions.
     */
    inline std::vector<double> pairprop(
        const std::vector<std::vector<uint8_t>>& pat1,
        const std::vector<std::vector<uint8_t>>& pat2
    )
    {
        if (pat1.size() != pat2.size()) {
            throw std::invalid_argument("Pattern spaces must have same number of rows.");
        }

        size_t pos_tot = 0;
        size_t neg_tot = 0;
        size_t pat_tot = 0;

        const size_t n_rows = pat1.size();

        for (size_t i = 0; i < n_rows; ++i) {

            const auto& row1 = pat1[i];
            const auto& row2 = pat2[i];

            // Skip rows that represent invalid pattern {0}
            if ((row1.size() == 1 && row1[0] == 0) ||
                (row2.size() == 1 && row2[0] == 0)) {
                continue;
            }

            if (row1.size() != row2.size()) {
                throw std::invalid_argument("Pattern rows must have same length.");
            }

            const size_t n_cols = row1.size();

            for (size_t j = 0; j < n_cols; ++j) {

                uint8_t s1 = row1[j];
                uint8_t s2 = row2[j];

                // Only count non-zero pairs
                if (s1 != 0 && s2 != 0) {

                    ++pat_tot;

                    // Positive agreement
                    if (s1 == s2) {
                        ++pos_tot;
                    }
                    // Negative agreement
                    else if ((s1 == 1 && s2 == 3) ||
                            (s1 == 3 && s2 == 1)) {
                        ++neg_tot;
                    }
                }
            }
        }

        std::vector<double> result(2);

        if (pat_tot == 0) {
            double nan = std::numeric_limits<double>::quiet_NaN();
            result[0] = nan;
            result[1] = nan;
        } else {
            result[0] = static_cast<double>(pos_tot) / static_cast<double>(pat_tot);
            result[1] = static_cast<double>(neg_tot) / static_cast<double>(pat_tot);
        }

        return result;
    }
    
} // namespace symdync

}

#endif // INFOXTR_SYMDYNC_HPP
