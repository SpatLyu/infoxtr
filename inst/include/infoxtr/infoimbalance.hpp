/******************************************************************************
 * File: infoimbalance.hpp
 *
 * Information Imbalance and Imbalance Gain
 * ----------------------------------------
 *
 * This module implements the Information Imbalance framework for comparing
 * the information content of distance spaces through distance ranks.
 *
 * Three related functions are provided:
 *
 *   infoImbalance()
 *       Computes the Information Imbalance between two representations:
 *
 *           Delta(X -> Y)
 *
 *       For each prediction point, the k nearest neighbours are identified
 *       in X space, and their corresponding ranks in Y space are averaged.
 *
 *   imbalanceGain()
 *       Computes the Information Imbalance for a series of scaling
 *       parameters alpha using the combined space:
 *
 *           A(alpha) = (alpha * X, Y)
 *
 *       and returns Delta(alpha) for each supplied alpha.
 *
 *   imbalanceGainCausality()
 *       Converts the series of Delta(alpha) values into a normalized
 *       Information Imbalance Gain:
 *
 *           IIG(X -> Y)
 *             = [Delta(0) - min_alpha Delta(alpha)] / Delta(0).
 *
 *       The result measures the relative reduction in Information Imbalance
 *       obtained by incorporating X into the representation used to predict
 *       Y.
 *
 * ---------------------------------------------------------------------------
 * Information Imbalance
 * ---------------------------------------------------------------------------
 *
 * Let X and Y denote two representations of the same set of observations.
 * Let P denote the prediction set and L the library set from which neighbours
 * are selected.
 *
 * For each prediction point p in P, distances to all library points q in L
 * are computed in the Y space and ranked increasingly:
 *
 *      r^Y_{pq}.
 *
 * The k nearest neighbours of p are then identified in the X space. The
 * corresponding ranks in the Y space are averaged:
 *
 *      r^Y_{p|X}
 *          = (1 / k) sum_j r^Y_{p,NN_j^X}.
 *
 * The Information Imbalance from X to Y is
 *
 *      Delta(X -> Y)
 *          = 2 / N sum_p r^Y_{p|X},
 *
 * where N is the number of library points.
 *
 * A small Information Imbalance indicates that points that are close in
 * X space also tend to be close in Y space. A larger value indicates weaker
 * correspondence between the two neighbourhood structures.
 *
 * The parameter k specifies the number of neighbours used in the comparison.
 * k = 1 corresponds to the nearest-neighbour formulation.
 *
 * ---------------------------------------------------------------------------
 * Information Imbalance Gain
 * ---------------------------------------------------------------------------
 *
 * For directional information analysis, the current states of X and Y are
 * combined using a scaling parameter alpha:
 *
 *      A(alpha) = (alpha * X, Y).
 *
 * The Information Imbalance is then evaluated between this combined current
 * space and the future state of Y:
 *
 *      Delta(alpha)
 *          = Delta(A(alpha) -> Y_future).
 *
 * When alpha = 0, the X contribution vanishes and the baseline becomes
 *
 *      Delta(0) = Delta(Y_current -> Y_future).
 *
 * If incorporating X provides additional information about the future state
 * of Y, an appropriate value of alpha should reduce the Information
 * Imbalance.
 *
 * ---------------------------------------------------------------------------
 * Information Imbalance Gain and Causality
 * ---------------------------------------------------------------------------
 *
 * Given a series of Information Imbalance values Delta(alpha), the normalized
 * Information Imbalance Gain is defined as
 *
 *      IIG(X -> Y)
 *          = [Delta(0) - min_alpha Delta(alpha)]
 *            / Delta(0).
 *
 * Thus:
 *
 *   Delta(0)
 *       Baseline Information Imbalance obtained without X.
 *
 *   min_alpha Delta(alpha)
 *       Best Information Imbalance obtained after incorporating X.
 *
 *   Delta(0) - min_alpha Delta(alpha)
 *       Absolute reduction in Information Imbalance.
 *
 *   IIG(X -> Y)
 *       Relative reduction normalized by the baseline.
 *
 * A value of zero indicates that adding X does not reduce the Information
 * Imbalance. A positive value indicates that incorporating X improves the
 * correspondence between the current representation and the future state
 * of Y.
 *
 * When applied to time series with a positive prediction horizon h, this
 * provides a model free measure of directional predictive information from
 * X to Y.
 *
 * ---------------------------------------------------------------------------
 * Distance ranks
 * ---------------------------------------------------------------------------
 *
 * The framework operates on distance ranks rather than raw distance values.
 * For a prediction point p, the rank of library point q is obtained by
 * ordering distances from smallest to largest, with rank 1 assigned to the
 * closest point.
 *
 * Tied distances receive their average rank, following the convention of
 *
 *      R::rank(..., ties.method = "average").
 *
 * Non-finite distances are placed after finite distances, corresponding to
 * the intended na.last = TRUE behaviour of the reference R implementation.
 *
 * ---------------------------------------------------------------------------
 * Library, prediction and common parameters
 * ---------------------------------------------------------------------------
 *
 * The following parameters are shared by all three functions:
 *
 *   Mx
 *       Matrix-like container containing the X representations.
 *
 *   My
 *       Matrix-like container containing the Y representations.
 *
 *   lib
 *       Indices of observations eligible to serve as library neighbours.
 *
 *   pred
 *       Indices of observations for which the Information Imbalance is
 *       evaluated.
 *
 *   k
 *       Number of nearest neighbours used in the rank comparison.
 *
 *   threads
 *       Number of threads used for parallel computation.
 *
 *   method
 *       Distance metric used to construct the distance spaces.
 *
 * When prediction and library sets overlap, the prediction point itself is
 * excluded from the nearest-neighbour search.
 *
 * ---------------------------------------------------------------------------
 * Function-specific parameters
 * ---------------------------------------------------------------------------
 *
 * infoImbalance()
 *
 *   No additional parameters.
 *
 *   Returns a single Information Imbalance value:
 *
 *       Delta(X -> Y).
 *
 *
 * imbalanceGain()
 *
 *   alpha
 *       Scaling parameters applied to X.
 *
 *   h
 *       Prediction horizon. The current space uses t, while the target
 *       Y space uses t + h.
 *
 *   Returns one Information Imbalance value Delta(alpha) for each supplied
 *   alpha.
 *
 *
 * imbalanceGainCausality()
 *
 *   alpha
 *       Scaling parameters used to evaluate Delta(alpha).
 *
 *   h
 *       Prediction horizon.
 *
 *   Returns a single normalized Information Imbalance Gain:
 *
 *       [Delta(0) - min_alpha Delta(alpha)] / Delta(0).
 *
 * ---------------------------------------------------------------------------
 * Prediction horizon
 * ---------------------------------------------------------------------------
 *
 * infoImbalance() directly compares the representations supplied through
 * Mx and My and therefore has no explicit prediction-horizon parameter.
 *
 * For imbalanceGain() and imbalanceGainCausality(), h defines the prediction
 * horizon:
 *
 *      A_t(alpha) = (alpha * X_t, Y_t)
 *
 *      B_{t+h}    = Y_{t+h}.
 *
 * For h = 0, the comparison is synchronous. For h > 0, the framework
 * evaluates whether the current state contains information about the future
 * state of Y.
 *
 * ---------------------------------------------------------------------------
 * Parallel computation
 * ---------------------------------------------------------------------------
 *
 * Calculations for different prediction points are independent and are
 * parallelized using RcppThread.
 *
 * In infoImbalance(), the construction of Y-space distance ranks and the
 * nearest-neighbour search in X space are parallelized across prediction
 * points.
 *
 * In imbalanceGain(), the same procedure is repeated independently for each
 * alpha value.
 *
 * Per-prediction results are stored independently and aggregated after
 * parallel execution to avoid data races.
 *
 * ---------------------------------------------------------------------------
 * Author: Wenbo Lyu (Github: @SpatLyu)
 * License: GPL-3
 ******************************************************************************/

#ifndef INFOXTR_INFOIMBALANCE_HPP
#define INFOXTR_INFOIMBALANCE_HPP

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

namespace infoimbalance
{
    
    using ContMat = std::vector<std::vector<double>>;

    inline double infoImbalance(
        const ContMat& Mx,
        const ContMat& My,
        const std::vector<size_t>& lib,
        const std::vector<size_t>& pred,
        size_t k = 1,
        size_t threads = 1,
        const std::string& method = "euclidean")
    {
        const size_t Npred = pred.size();
        const size_t Nlib  = lib.size();

        if (Npred == 0 || Nlib == 0 || k == 0) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        if (Mx.empty() || My.empty()) {
            throw std::invalid_argument(
                "Mx and My must not be empty.");
        }

        if (Mx.size() != My.size()) {
            throw std::invalid_argument(
                "Mx and My must contain the same number of observations.");
        }

        if (k > Nlib) {
            k = Nlib;
        }

        if (threads == 0) {
            threads = 1;
        }

        const size_t hw = std::thread::hardware_concurrency();

        if (hw > 0) {
            threads = std::min(threads, hw);
        }

        // ------------------------------------------------------------------
        // Compute the ranks in the Y space.
        //
        // For each prediction point, distances to all library points are
        // ranked increasingly. Ties receive their average rank, following
        // R's rank(..., ties.method = "average").
        //
        // The prediction point itself is treated as missing when it is also
        // present in the library and is therefore placed at the end.
        // ------------------------------------------------------------------
        std::vector<std::vector<double>> y_rank(
            Npred,
            std::vector<double>(Nlib));

        RcppThread::parallelFor(
            size_t(0), Npred, [&](size_t ip) {
                const size_t p = pred[ip];

                std::vector<double> distances(Nlib);

                for (size_t il = 0; il < Nlib; ++il) {
                    const size_t q = lib[il];

                    if (p == q) {
                        distances[il] =
                            std::numeric_limits<double>::quiet_NaN();
                        continue;
                    }

                    distances[il] =
                        infoxtr::distance::distance(My[p], My[q], method, true);
                }

                // Sort library indices according to Y-space distance.
                // Finite distances come first and non-finite distances
                // are placed at the end.
                std::vector<size_t> order(Nlib);
                std::iota(order.begin(), order.end(), size_t(0));
                std::sort(
                    order.begin(),
                    order.end(),
                    [&](size_t a, size_t b) {
                        const bool a_finite =
                            std::isfinite(distances[a]);

                        const bool b_finite =
                            std::isfinite(distances[b]);

                        if (a_finite != b_finite) {
                            return a_finite;
                        }

                        if (!a_finite && !b_finite) {
                            return lib[a] < lib[b];
                        }

                        if (distances[a] < distances[b]) {
                            return true;
                        }

                        if (distances[a] > distances[b]) {
                            return false;
                        }

                        return lib[a] < lib[b];
                    });

                // Assign average ranks to tied distances.
                size_t pos = 0;

                while (pos < Nlib) {
                    size_t end = pos + 1;

                    const bool current_finite =
                        std::isfinite(distances[order[pos]]);

                    while (end < Nlib) {
                        const bool next_finite =
                            std::isfinite(distances[order[end]]);

                        // All non-finite distances form one final group.
                        if (!current_finite && !next_finite) {
                            ++end;
                            continue;
                        }

                        // Finite distances are tied when their values are
                        // numerically equal within the specified tolerance.
                        if (current_finite &&
                            next_finite &&
                            infoxtr::numericutils::doubleNearlyEqual(
                                distances[order[end]],
                                distances[order[pos]])) {
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

            }, threads);

        // ------------------------------------------------------------------
        // For each prediction point, find the k nearest neighbours in the
        // X space and calculate the mean of their corresponding Y ranks.
        // ------------------------------------------------------------------
        std::vector<double> rank_sums(Npred, 0.0);

        RcppThread::parallelFor(
            size_t(0), Npred, [&](size_t ip) {
                const size_t p = pred[ip];

                struct Candidate {
                    size_t lib_pos;
                    double distance;
                };

                std::vector<Candidate> candidates;
                candidates.reserve(Nlib);

                // Compute distances to all library points in X space.
                for (size_t il = 0; il < Nlib; ++il) {

                    const size_t q = lib[il];

                    // Exclude the prediction point itself.
                    if (q == p) {
                        continue;
                    }

                    const double d =
                        infoxtr::distance::distance(Mx[p], Mx[q], method, true);

                    candidates.push_back({il, d});
                }

                const size_t nk =
                    std::min(k, candidates.size());

                if (nk == 0) {
                    return;
                }

                // Select the k nearest neighbours in X space.
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

                        if (lhs_finite != rhs_finite) {
                            return lhs_finite;
                        }

                        if (!lhs_finite && !rhs_finite) {
                            return lib[lhs.lib_pos] <
                                lib[rhs.lib_pos];
                        }

                        if (lhs.distance < rhs.distance) {
                            return true;
                        }

                        if (lhs.distance > rhs.distance) {
                            return false;
                        }

                        return lib[lhs.lib_pos] <
                            lib[rhs.lib_pos];
                    });

                // Average the corresponding Y-space ranks.
                double conditional_rank = 0.0;

                for (size_t j = 0; j < nk; ++j) {
                    conditional_rank +=
                        y_rank[ip][candidates[j].lib_pos];
                }

                rank_sums[ip] =
                    conditional_rank / static_cast<double>(nk);

            }, threads);

        // ------------------------------------------------------------------
        // Information Imbalance:
        //
        //     II(X -> Y)
        //       = 2 / Npred
        //         * mean_i[
        //             mean_{j in NN_k^X(i)} r^Y_ij
        //           ]
        // ------------------------------------------------------------------
        const double rank_sum = 
            std::accumulate(rank_sums.begin(), rank_sums.end(), 0.0);

        const double mean_rank =
            rank_sum / static_cast<double>(Npred);

        return 2.0 / static_cast<double>(Nlib) * mean_rank;
    }
    
    inline std::vector<double> imbalanceGain(
        const ContMat& Mx,
        const ContMat& My,
        const std::vector<double>& alpha,
        const std::vector<size_t>& lib,
        const std::vector<size_t>& pred,
        size_t h = 1,
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
                        My[p + h], My[q + h], method, true);
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
                    return;
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
            // where N corresponds to the number of library points.
            //
            // This is the direct generalization of the R implementation,
            // with ranks computed over the supplied library.
            // --------------------------------------------------------------
            const double mean_rank =
                rank_sum / static_cast<double>(Npred);

            out[ia] =
                2.0 / static_cast<double>(Nlib) * mean_rank;
        }

        return out;
    }

    inline double imbalanceGainCausality(
        const ContMat& Mx,
        const ContMat& My,
        const std::vector<double>& alpha,
        const std::vector<size_t>& lib,
        const std::vector<size_t>& pred,
        size_t h = 1,
        size_t k = 1,
        size_t threads = 1,
        const std::string& method = "euclidean")
    {   
        std::vector<double> alpha_ext;
        alpha_ext.reserve(alpha.size() + 1);

        for (double a : alpha) {
            if (std::isfinite(a) && a >= 0.0) {
                alpha_ext.push_back(a);
            }
        }

        // Always include alpha = 0 as the baseline.
        alpha_ext.push_back(0.0);

        std::sort(alpha_ext.begin(), alpha_ext.end());

        std::vector<double> unique_alpha;
        if (!alpha_ext.empty()) {
            unique_alpha.push_back(alpha_ext[0]);

            for (size_t i = 1; i < alpha_ext.size(); ++i) {
                if (!infoxtr::numericutils::doubleNearlyEqual(
                        alpha_ext[i], unique_alpha.back())) {
                    unique_alpha.push_back(alpha_ext[i]);
                }
            }
        }

        if (unique_alpha.empty()) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        std::vector<double> ii_vals = imbalanceGain(
            Mx, My, unique_alpha, lib, pred, h, k, threads, method
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
    
} // namespace infoimbalance

}

#endif // INFOXTR_INFOIMBALANCE_HPP
