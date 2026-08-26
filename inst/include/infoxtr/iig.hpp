/******************************************************************************
 * File: iig.hpp
 *
 * Information Imbalance Gain (IIG)
 * -----------------------------------
 *
 * This module implements the information imbalance framework for quantifying
 * the degree to which neighbourhood structure in one representation is
 * preserved in another representation.
 *
 * Given two representations X and Y, the method identifies nearest neighbours
 * in a combined space
 *
 *      A = (alpha * X, Y)
 *
 * and evaluates their corresponding ranks in the Y space. The resulting
 * information imbalance measures how strongly the neighbourhood structure
 * induced by X and Y is associated.
 *
 * ---------------------------------------------------------------------------
 * Algorithm overview
 * ---------------------------------------------------------------------------
 *
 * Let X and Y denote two multivariate representations with the same number
 * of observations. For a prediction set P and a library set L, the algorithm
 * proceeds as follows:
 *
 * 1. Y-space rank construction
 *
 *      For each prediction point p in P, compute distances from Y_p to all
 *      library points Y_q, q in L.
 *
 *      The distances are ranked increasingly using average ranks for tied
 *      distances, following the convention of R's rank() with
 *      ties.method = "average".
 *
 *      Missing or non-finite distances are placed after valid distances.
 *      If p is also contained in the library, its self-distance is treated
 *      as missing and therefore cannot affect the nearest-neighbour search.
 *
 * 2. Combined-space neighbourhood construction
 *
 *      For each scaling parameter alpha, construct the combined representation
 *
 *          A = (alpha * X, Y).
 *
 *      For every prediction point p, identify its k nearest neighbours among
 *      the library points in A.
 *
 *      The prediction point itself is excluded when the prediction and library
 *      sets overlap.
 *
 *      Neighbours are selected using partial sorting. Distance is the primary
 *      criterion, while the original library index provides a deterministic
 *      tie-breaking rule.
 *
 * 3. Conditional Y-space ranks
 *
 *      For each prediction point, retrieve the Y-space ranks corresponding to
 *      its k nearest neighbours in A and compute their mean:
 *
 *          r_p = (1 / k) sum_j rank_Y(p, NN_j^A)
 *
 *      This quantity measures how close the neighbours selected in A are to
 *      the prediction point in the Y space.
 *
 * 4. Information imbalance
 *
 *      The information imbalance is computed as
 *
 *          II = 2 / N * mean_p(r_p)
 *
 *      where N is the number of prediction points.
 *
 *      Smaller values indicate stronger preservation of neighbourhood
 *      structure between the two representations, whereas larger values
 *      indicate greater neighbourhood mismatch.
 *
 * 5. Scaling analysis
 *
 *      When multiple values of alpha are supplied, the information imbalance
 *      is evaluated independently for each alpha. This allows the relative
 *      contribution of the X and Y components to the neighbourhood structure
 *      to be examined.
 *
 * ---------------------------------------------------------------------------
 * Missing values and ties
 * ---------------------------------------------------------------------------
 *
 * Distances that are NaN or otherwise non-finite are placed after valid
 * distances during ranking and neighbour selection, following the intended
 * na.last = TRUE behaviour of the reference R implementation.
 *
 * Tied Y-space distances receive average ranks. When selecting nearest
 * neighbours in the combined space, equal distances are ordered according
 * to the original library index to ensure deterministic results.
 *
 * ---------------------------------------------------------------------------
 * Parallel computation
 * ---------------------------------------------------------------------------
 *
 * Prediction points are processed independently using RcppThread. Both the
 * Y-space rank construction and the nearest-neighbour calculations are
 * parallelized across prediction points.
 *
 * Per-prediction intermediate results are stored independently and aggregated
 * after parallel execution to avoid data races.
 *
 * ---------------------------------------------------------------------------
 * Input data format
 * ---------------------------------------------------------------------------
 *
 *   Mx       : Matrix-like container of X representations.
 *   My       : Matrix-like container of Y representations.
 *   alpha    : Scaling parameters applied to X.
 *   lib      : Indices of observations eligible as library neighbours.
 *   pred     : Indices of prediction observations.
 *   k        : Number of nearest neighbours.
 *   threads  : Number of threads used for parallel computation.
 *   method   : Distance metric used for neighbourhood construction.
 *
 * Output:
 *
 *   A vector containing the information imbalance corresponding to each
 *   value of alpha.
 *
 * Author: Wenbo Lyu (Github: @SpatLyu)
 * License: GPL-3
 ******************************************************************************/

#ifndef INFOXTR_IIG_HPP
#define INFOXTR_IIG_HPP

#include <vector>
#include <cmath>
#include <limits>
#include <string>
#include <thread>
#include <cstdint>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include "infoxtr/distance.hpp"
#include "infoxtr/numericutils.hpp"
#include <RcppThread.h>

namespace infoxtr 
{

namespace iig
{
    
    using ContMat = std::vector<std::vector<double>>;

    inline std::vector<double> infoImbalance(
        const ContMat& Mx,
        const ContMat& My,
        const std::vector<double>& alpha,
        const std::vector<size_t>& lib,
        const std::vector<size_t>& pred,
        size_t k = 1,
        size_t threads = 1,
        const std::string& method = "euclidean")
    {
        const size_t Npred = pred.size();
        const size_t Nlib  = lib.size();

        if (Npred == 0 || Nlib == 0 || k == 0) {
            return std::vector<double>(alpha.size(),
                                       std::numeric_limits<double>::quiet_NaN());
        }

        if (Mx.empty() || My.empty()) {
            throw std::invalid_argument("Mx and My must not be empty.");
        }

        if (Mx.size() != My.size()) {
            throw std::invalid_argument(
                "Mx and My must contain the same number of observations.");
        }

        if (k > Nlib) {
            k = Nlib;
        }

        if (threads == 0) threads = 1;
        size_t hw = std::thread::hardware_concurrency();
        if (hw > 0) threads = std::min(threads, hw);

        const size_t dimX = Mx[0].size();
        const size_t dimY = My[0].size();

        // ------------------------------------------------------------------
        // Compute the ranks in the Y space.
        //
        // For each prediction point, distances to all library points are
        // ranked increasingly. Ties receive their average rank, following
        // R's rank(..., ties.method = "average").
        //
        // If the prediction point is also in the library, its own rank is
        // explicitly set to the largest rank, corresponding to the Inf
        // diagonal used in the R implementation.
        // ------------------------------------------------------------------
        std::vector<std::vector<double>> y_rank(
            Npred, std::vector<double>(Nlib));

        RcppThread::parallelFor(
            size_t(0), Npred, [&](size_t ip) {

            const size_t p = pred[ip];

            std::vector<double> distances(Nlib);

            for (size_t il = 0; il < Nlib; ++il) {

                const size_t q = lib[il];

                if (p == q) {
                    distances[il] = std::numeric_limits<double>::quiet_NaN();
                    continue;
                }

                distances[il] =
                    infoxtr::distance::distance(
                        My[p], My[q], method, true);
            }

            // Sort indices according to Y-space distance.
            // Finite distances come first; NA/NaN/Inf distances are moved
            // to the end, corresponding to R's na.last = TRUE.
            // For finite ties, library indices provide a deterministic ordering.
            std::vector<size_t> order(Nlib);
            std::iota(order.begin(), order.end(), size_t(0));

            std::sort(
                order.begin(),
                order.end(),
                [&](size_t a, size_t b) {

                    const bool a_finite = std::isfinite(distances[a]);
                    const bool b_finite = std::isfinite(distances[b]);

                    // Missing/non-finite distances are placed last.
                    if (a_finite != b_finite) {
                        return a_finite;
                    }

                    // Both are non-finite: keep a deterministic order.
                    if (!a_finite && !b_finite) {
                        return lib[a] < lib[b];
                    }

                    // Both are finite.
                    if (distances[a] < distances[b]) {
                        return true;
                    }

                    if (distances[a] > distances[b]) {
                        return false;
                    }

                    return lib[a] < lib[b];
                });

            // Assign average ranks for tied distances.
            // All non-finite distances form the final tie group, corresponding
            // to na.last = TRUE in R.
            size_t pos = 0;

            while (pos < Nlib) {

                size_t end = pos + 1;

                const bool current_finite =
                    std::isfinite(distances[order[pos]]);

                while (end < Nlib) {

                    const bool next_finite =
                        std::isfinite(distances[order[end]]);

                    // Non-finite distances are all treated as one final tie group.
                    if (!current_finite && !next_finite) {
                        ++end;
                        continue;
                    }

                    // Finite distances are tied only when their values are equal.
                    if (current_finite &&
                        next_finite &&
                        infoxtr::numericutils::doubleNearlyEqual(distances[order[end]], distances[order[pos]])) {
                        ++end;
                        continue;
                    }

                    break;
                }

                // R rank() uses 1-based ranks.
                const double rank =
                    0.5 * (
                        static_cast<double>(pos + 1) +
                        static_cast<double>(end)
                    );

                for (size_t j = pos; j < end; ++j) {
                    y_rank[ip][order[j]] = rank;
                }

                pos = end;
            }
            
            // // If the prediction point is also in the library, its own distance
            // // is treated as NA and therefore placed last. The self observation
            // // is excluded from nearest-neighbour selection in the A space.
            // for (size_t il = 0; il < Nlib; ++il) {
            //     if (lib[il] == p) {
            //         y_rank[ip][il] = static_cast<double>(Nlib);
            //     }
            // }

        }, threads);

        // ------------------------------------------------------------------
        // Compute information imbalance for each alpha.
        // ------------------------------------------------------------------
        std::vector<double> out(alpha.size());

        for (size_t ia = 0; ia < alpha.size(); ++ia) {

            const double a = alpha[ia];

            std::vector<double> rank_sums(Npred, 0.0);

            // --------------------------------------------------------------
            // For every prediction point, find its k nearest neighbours
            // in A = (alpha * X, Y).
            // --------------------------------------------------------------
            RcppThread::parallelFor(
                size_t(0), Npred, [&](size_t ip) {

                const size_t p = pred[ip];

                struct Candidate {
                    size_t lib_pos;
                    double distance;
                };

                std::vector<Candidate> candidates;
                candidates.reserve(Nlib);

                std::vector<double> vec_p(dimX + dimY);
                std::vector<double> vec_q(dimX + dimY);

                for (size_t d = 0; d < dimX; ++d) vec_p[d] = a * Mx[p][d];
                for (size_t d = 0; d < dimY; ++d) vec_p[dimX + d] = My[p][d];

                for (size_t il = 0; il < Nlib; ++il) {
                    const size_t q = lib[il];

                    // Exclude the prediction point itself when pred and lib overlap,
                    // consistent with the NA diagonal used in the R implementation.
                    if (q == p) continue;

                    for (size_t d = 0; d < dimX; ++d) vec_q[d] = a * Mx[q][d];
                    for (size_t d = 0; d < dimY; ++d) vec_q[dimX + d] = My[q][d];

                    const double d = infoxtr::distance::distance(vec_p, vec_q, method, true);
                    candidates.push_back({il, d});
                }

                // ----------------------------------------------------------
                // Select the k nearest neighbours.
                //
                // Distance is the primary ordering criterion. If several
                // library points have exactly the same distance, their
                // original library indices determine the ordering.
                // This reproduces a deterministic tie-breaking rule for
                // the partial selection.
                // ----------------------------------------------------------
                const size_t nk = std::min(k, candidates.size());

                if (nk == 0) {
                    continue;
                }

                std::partial_sort(
                    candidates.begin(),
                    candidates.begin() + nk,
                    candidates.end(),
                    [&](const Candidate& lhs,
                        const Candidate& rhs) {

                        const bool lhs_finite =
                            std::isfinite(lhs.distance);

                        const bool rhs_finite =
                            std::isfinite(rhs.distance);

                        // Non-finite distances are always placed after finite
                        // distances, corresponding to na.last = TRUE.
                        if (lhs_finite != rhs_finite) {
                            return lhs_finite;
                        }

                        // Both distances are non-finite.
                        // Use the library index as the deterministic tie breaker.
                        if (!lhs_finite && !rhs_finite) {
                            return lib[lhs.lib_pos] < lib[rhs.lib_pos];
                        }

                        // Both distances are finite.
                        if (lhs.distance < rhs.distance) {
                            return true;
                        }

                        if (lhs.distance > rhs.distance) {
                            return false;
                        }

                        // Equal distances are ordered by library index.
                        return lib[lhs.lib_pos] < lib[rhs.lib_pos];
                    });

                // ----------------------------------------------------------
                // Average the corresponding ranks in Y.
                // ----------------------------------------------------------
                double conditional_rank = 0.0;

                for (size_t j = 0; j < nk; ++j) {
                    conditional_rank +=
                        y_rank[ip][candidates[j].lib_pos];
                }

                conditional_rank /= static_cast<double>(nk);

                rank_sums[ip] = conditional_rank;
            }, threads);

            const double rank_sum = 
                std::accumulate(rank_sums.begin(), rank_sums.end(), 0.0);

            // --------------------------------------------------------------
            // Information imbalance:
            //
            //     II = 2 / N * mean(conditional ranks)
            //
            // where N corresponds to the number of prediction points.
            //
            // This is the direct generalization of the R implementation,
            // with ranks computed over the supplied library.
            // --------------------------------------------------------------
            const double mean_rank =
                rank_sum / static_cast<double>(Npred);

            out[ia] =
                2.0 / static_cast<double>(Npred) * mean_rank;
        }

        return out;
    }

    inline double infoImbalanceGain(
        const ContMat& Mx,
        const ContMat& My,
        const std::vector<double>& alpha,
        const std::vector<size_t>& lib,
        const std::vector<size_t>& pred,
        size_t k = 1,
        size_t threads = 1,
        const std::string& method = "euclidean")
    {
        std::vector<double> alpha_ext = alpha;
        alpha_ext.push_back(0.0);
        std::sort(alpha_ext.begin(), alpha_ext.end());

        std::vector<double> unique_alpha;
        if (!alpha_ext.empty()) {
            unique_alpha.push_back(alpha_ext[0]);
            for (size_t i = 1; i < alpha_ext.size(); ++i) {
                if (!infoxtr::numericutils::doubleNearlyEqual(alpha_ext[i], unique_alpha.back())) {
                    unique_alpha.push_back(alpha_ext[i]);
                }
            }
        }

        if (unique_alpha.empty()) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        std::vector<double> ii_vals = infoImbalance(
            Mx, My, unique_alpha, lib, pred, k, threads, method
        );

        if (ii_vals.empty()) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        double min_val = std::numeric_limits<double>::infinity();
        bool has_valid = false;
        
        for (double v : ii_vals) {
            if (std::isfinite(v)) {
                if (!has_valid || v < min_val) {
                    min_val = v;
                    has_valid = true;
                }
            }
        }

        if (!has_valid) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        double first_val = ii_vals[0];
        
        if (!std::isfinite(first_val) || infoxtr::numericutils::doubleNearlyEqual(first_val, 0.0)) {
            return std::numeric_limits<double>::quiet_NaN(); 
        }

        return (first_val - min_val) / first_val;
    }
    
} // namespace iig

}

#endif // INFOXTR_IIG_HPP
