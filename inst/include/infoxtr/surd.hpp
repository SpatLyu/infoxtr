

#ifndef INFOXTR_SURD_HPP
#define INFOXTR_SURD_HPP

#include <vector>
#include <cmath>
#include <limits>
#include <thread>
#include <cstdint>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include "infoxtr/combn.hpp"
#include "infoxtr/numericutils.hpp"
#include <RcppThread.h>

namespace infoxtr 
{

namespace surd 
{

    using DiscMat = std::vector<std::vector<uint64_t>>;
    using ContMat = std::vector<std::vector<double>>;

    /***********************************************************
     * Utilities
     ***********************************************************/

    

    /***********************************************************
     * Result structure
     ***********************************************************/
    struct SURDRes
    {
        std::vector<double> values;
        std::vector<uint8_t> types;
        std::vector<std::vector<size_t>> var_indices;

        size_t size() const noexcept { return values.size(); }
    };

    /***************************************************************
     * Synergistic-Unique-Redundant Decomposition for Discrete Data
     ***************************************************************/
    inline SURDRes surd(
        const DiscMat& mat,
        size_t max_order = std::numeric_limits<size_t>::max(),
        size_t threads = 1,
        double base = 2.0,
        bool normalize = false)
    {
        if (threads == 0) threads = 1;
        size_t hw = std::thread::hardware_concurrency();
        if (hw > 0) threads = std::min(threads, hw);

        SURDRes result;

        if (mat.size() < 2) return result;

        const size_t n_vars = mat.size();
        const size_t n_sources = mat.size() - 1;
        max_order = std::min(max_order, n_sources);

        // Construct variable combination vector
        std::vector<size_t> ag_idx(n_sources.size());
        std::iota(ag_idx.begin(), ag_idx.end(), 1);
        const std::vector<std::vector<size_t>> combs =
            infoxtr::combn::genSubsets(ag_idx, max_order);

        // Compute joint entropies
        std::vector<double> H_sources(combs.size(), 
                                      std::numeric_limits<double>::quiet_NaN());
        std::vector<double> H_joint(combs.size(), 
                                    std::numeric_limits<double>::quiet_NaN());
        double H_target = infoxtr::infotheo::je(mat, {0}, base, true);

        for (size_t i = 0; i < combs.size(); ++i)
        {   
            std::vector<size_t> joint_idx = {0};
            joint_idx.insert(joint_idx.end(), combs[i].begin(), combs[i].end());
            H_sources[i] = Infoxtr::infotheo::je(mat, combs[i], base, true);
            H_joint[i] = Infoxtr::infotheo::je(mat, joint_idx, base, true);
        }
            
    }

    // /*****************************************************************
    //  * Synergistic-Unique-Redundant Decomposition for Continuous Data
    //  *****************************************************************/
    // inline SURDRes surd(
    //     const ContMat& mat,
    //     size_t max_order = std::numeric_limits<size_t>::max(),
    //     size_t k = 3,
    //     size_t alg = 0,
    //     size_t threads = 1,
    //     double base = 2.0,
    //     bool normalize = false)
    
} // namespace surd

}

#endif // INFOXTR_SURD_HPP
