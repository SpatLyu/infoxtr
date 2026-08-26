/******************************************************************************
 * File: iig.hpp
 *
 * Information Imbalance Gain (IIG)
 * -------------------------------------------------
 *
 * This module implements the SURD framework for decomposing the mutual
 * information between a target variable and multiple source variables into
 * interpretable components:
 *
 *   Redundant information
 *       Information simultaneously provided by multiple sources.
 *
 *   Unique information
 *       Information provided exclusively by a single source.
 *
 *   Synergistic information
 *       Information that only emerges when variables are considered jointly.
 *
 *   Information leak
 *       Remaining uncertainty in the target after conditioning on all sources.
 *
 * ---------------------------------------------------------------------------
 * Algorithm overview
 * ---------------------------------------------------------------------------
 *
 * Let Y denote the target variable and X = {X1, X2, ..., Xn} the source
 * variables. The algorithm proceeds in the following stages:
 *
 * 1. Subset enumeration
 *
 *      Generate all subsets of source variables up to order max_order.
 *      Each subset represents a candidate information channel.
 *
 * 2. Joint distribution construction
 *
 *      A joint state table is constructed for the full variable set using
 *      the joint entropy utilities in infotheo.hpp.
 *
 * 3. Conditional pointwise mutual information
 *
 *      For each target state s, compute pointwise mutual information
 *
 *          I_s(X_set ; Y)
 *
 *      for every subset X_set using grouped projections of the joint
 *      state table.
 *
 *      Mutual information is accumulated as:
 *
 *          I(X_set ; Y) = sum_s p(s) * I_s(X_set ; Y)
 *
 * 4. Monotonic SURD filtering
 *
 *      Subsets are sorted by their pointwise mutual information values.
 *      A monotonic constraint is enforced so that higher-order subsets
 *      cannot contain less information than the maximum of lower-order
 *      subsets. Violations are clipped to zero.
 *
 * 5. Information layer decomposition
 *
 *      A ladder-style decomposition is applied to the sorted values.
 *      Incremental differences between successive layers determine
 *      the information contributions.
 *
 *          delta_i = max(I_i - I_{i-1}, 0)
 *
 *      Contributions are classified as:
 *
 *          |subset| = 1   → redundant / unique layer
 *          |subset| > 1   → synergistic layer
 *
 * 6. Aggregation across target states
 *
 *      Contributions are weighted by p(s) and accumulated across
 *      target states.
 *
 * 7. Information leak
 *
 *      Remaining uncertainty in the target is measured as
 *
 *          H(Y | X_all) / H(Y)
 *
 * Input data format:
 *
 *   Row 0  : target variable
 *   Row 1+ : source variables
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

    inline std::vector<double> infoimbalance(
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

                    // The R implementation sets the diagonal distance to
                    // NA before ranking and subsequently replaces its rank
                    // with Inf. Thus, when pred and lib overlap, the point
                    // itself must not be selected as a neighbour.
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

    
    
} // namespace iig

}

#endif // INFOXTR_IIG_HPP
