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
    inline std::vector<std::vector<uint8_t>> gensympat(
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
    inline std::vector<double> countSignProp(
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

    /***************************************************************
     *  Pattern Causality Analysis Result
     ***************************************************************/
    struct PatternCausalityRes
    {
        // Per-observation causality strength
        std::vector<double> NoCausality;
        std::vector<double> PositiveCausality;
        std::vector<double> NegativeCausality;
        std::vector<double> DarkCausality;

        // Pattern classification per valid observation
        // 0: No causality
        // 1: Positive
        // 2: Negative
        // 3: Dark
        std::vector<size_t> PatternTypes;

        // Aggregated statistics (mean over valid entries)
        double TotalPos  = std::numeric_limits<double>::quiet_NaN();
        double TotalNeg  = std::numeric_limits<double>::quiet_NaN();
        double TotalDark = std::numeric_limits<double>::quiet_NaN();
    };

    /**
     * Compute pattern causality from Y to X.
     *
     *
     * Pipeline:
     *
     *  1. Collect unique patterns from X, Y_real and Y_pred.
     *  2. Remove patterns containing symbol 0.
     *  3. Add symmetric opposite patterns (1 <-> 3).
     *  4. Sort lexicographically to obtain deterministic indexing.
     *  5. Build K x K causal heatmap.
     *  6. Classify each observation:
     *
     *       0  No causality
     *       1  Positive     (main diagonal)
     *       2  Negative     (anti diagonal)
     *       3  Dark         (other off diagonal)
     *
     *  7. Optional strength weighting:
     *
     *       erf( ||pred_Y|| / (||Y|| + 1e-6) )
     */
    inline PatternCausalityRes computePatternCausality(
        const std::vector<std::vector<double>>& SMx,
        const std::vector<std::vector<double>>& SMy,
        const std::vector<std::vector<double>>& pred_SMy,
        bool weighted = true,
        bool save_detail = true
    )
    {   
        PatternCausalityRes res;

        const size_t n = SMx.size();
        if (n == 0) return res;
        
        /* ------------------------------------------------------------
        *  1. Generate symbolic pattern
        * ------------------------------------------------------------ */
        std::vector<std::vector<uint8_t>> PX = genPatternSpace(SMx, true);
        std::vector<std::vector<uint8_t>> PY_real = genPatternSpace(SMy, true);
        std::vector<std::vector<uint8_t>> PY_pred = genPatternSpace(pred_SMy, true);

        /* ------------------------------------------------------------
        *  2. Collect and filter pattern space
        * ------------------------------------------------------------ */
        std::vector<std::vector<uint8_t>> all_patterns;
        all_patterns.reserve(n * 3);

        auto contains_zero = [](const std::vector<uint8_t>& p)
        {
            for (uint8_t v : p) if (v == 0) return true;
            return false;
        };

        for (size_t i = 0; i < n; ++i)
        {
            if (!contains_zero(PX[i]))       all_patterns.push_back(PX[i]);
            if (!contains_zero(PY_real[i]))  all_patterns.push_back(PY_real[i]);
            if (!contains_zero(PY_pred[i]))  all_patterns.push_back(PY_pred[i]);
        }

        if (all_patterns.empty()) return res;

        std::sort(all_patterns.begin(), all_patterns.end());
        all_patterns.erase(
            std::unique(all_patterns.begin(), all_patterns.end()),
            all_patterns.end()
        );

        /* ------------------------------------------------------------
        *  3. Symmetric closure
        * ------------------------------------------------------------ */
        size_t original_size = all_patterns.size();

        for (size_t i = 0; i < original_size; ++i)
        {
            std::vector<uint8_t> opp = all_patterns[i];
            for (auto& v : opp)
            {
                if (v == 1) v = 3;
                else if (v == 3) v = 1;
            }
            all_patterns.push_back(std::move(opp));
        }

        std::sort(all_patterns.begin(), all_patterns.end());
        all_patterns.erase(
            std::unique(all_patterns.begin(), all_patterns.end()),
            all_patterns.end()
        );

        const size_t K = all_patterns.size();
        if (K == 0) return res;

        /* ------------------------------------------------------------
         * 4. Precompute opposite index
         * ------------------------------------------------------------ */
        std::vector<size_t> opposite_id(K);

        for (size_t i = 0; i < K; ++i)
        {
            auto opp = all_patterns[i];
            for (auto& v : opp)
            {
                if (v == 1) v = 3;
                else if (v == 3) v = 1;
            }

            auto it = std::lower_bound(all_patterns.begin(), all_patterns.end(), opp);
            opposite_id[i] = std::distance(all_patterns.begin(), it);
        }

        /* ------------------------------------------------------------
         *  5. Init result
         * ------------------------------------------------------------ */
        if (save_detail)
        {
            res.NoCausality.assign(n, 0.0);
            res.PositiveCausality.assign(n, 0.0);
            res.NegativeCausality.assign(n, 0.0);
            res.DarkCausality.assign(n, 0.0);
            res.PatternTypes.reserve(n);
        }

        std::vector<std::vector<double>> heatmap(
            K, std::vector<double>(K, std::numeric_limits<double>::quiet_NaN())
        );

        std::vector<std::vector<size_t>> counts(
            K, std::vector<size_t>(K, 0)
        );

        // Norm utility
        auto norm_ignore_nan = [](const std::vector<double>& v)
        {
            double sum = 0.0;
            for (double x : v)
                if (!std::isnan(x)) sum += x * x;
            return std::sqrt(sum);
        };

        /* ------------------------------------------------------------
         *  6. Main loop
         * ------------------------------------------------------------ */
        for (size_t t = 0; t < n; ++t)
        {
            if (contains_zero(PX[t]) ||
                contains_zero(PY_real[t]) ||
                contains_zero(PY_pred[t]))
                continue;

            /* --- causality existence --- */
            if (PY_pred[t] != PY_real[t])
            {   
                if (save_detail)
                {
                    res.NoCausality[t] = 1.0;
                    res.PatternTypes.push_back(0);
                }
                continue;
            }

            /* --- strength --- */
            double strength = weighted
                ? std::erf(
                    norm_ignore_nan(pred_SMy[t]) /
                    (norm_ignore_nan(SMy[t]) + 1e-6))
                : 1.0;

            /* --- index lookup --- */
            auto it_i = std::lower_bound(all_patterns.begin(), all_patterns.end(), PX[t]);
            auto it_j = std::lower_bound(all_patterns.begin(), all_patterns.end(), PY_pred[t]);

            if (it_i == all_patterns.end() || it_j == all_patterns.end())
                continue;

            size_t i = std::distance(all_patterns.begin(), it_i);
            size_t j = std::distance(all_patterns.begin(), it_j);

            /* --- classification --- */
            if (save_detail)
            {
                if (i == j)
                {
                    res.PositiveCausality[t] = strength;
                    res.PatternTypes.push_back(1);
                }
                else if (j == opposite_id[i])
                {
                    res.NegativeCausality[t] = strength;
                    res.PatternTypes.push_back(2);
                }
                else
                {
                    res.DarkCausality[t] = strength;
                    res.PatternTypes.push_back(3);
                }
            }

            /* --- heatmap --- */
            if (std::isnan(heatmap[i][j]))
            {
                heatmap[i][j] = strength;
                counts[i][j] = 1;
            }
            else
            {
                heatmap[i][j] += strength;
                counts[i][j] += 1;
            }
        }

        /* ------------------------------------------------------------
         * 7. Normalize + aggregate
         * ------------------------------------------------------------ */
        double diag_sum = 0.0, anti_sum = 0.0, other_sum = 0.0;
        size_t diag_cnt = 0, anti_cnt = 0, other_cnt = 0;

        for (size_t i = 0; i < K; ++i)
        {
            for (size_t j = 0; j < K; ++j)
            {
                if (counts[i][j] == 0) continue;

                double val = heatmap[i][j] / counts[i][j];

                if (i == j)
                {
                    diag_sum += val;
                    diag_cnt++;
                }
                else if (j == opposite_id[i])
                {
                    anti_sum += val;
                    anti_cnt++;
                }
                else
                {
                    other_sum += val;
                    other_cnt++;
                }
            }
        }

        res.TotalPos  = (diag_cnt  > 0) ? diag_sum  / diag_cnt  : 0.0;
        res.TotalNeg  = (anti_cnt  > 0) ? anti_sum  / anti_cnt  : 0.0;
        res.TotalDark = (other_cnt > 0) ? other_sum / other_cnt : 0.0;

        return res;
    }
    
} // namespace symdync

}

#endif // INFOXTR_SYMDYNC_HPP
