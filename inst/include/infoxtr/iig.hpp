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
        const std::vector<size_t>& lib,
        const std::vector<size_t>& pred,
        double alpha = 1.0,  
        size_t k = 1,  
        size_t h = 0,
        size_t threads = 1,
        const std::string& method = "euclidean")
    {
        
    }

    
    
} // namespace iig

}

#endif // INFOXTR_IIG_HPP
