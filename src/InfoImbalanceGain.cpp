#include <vector>
#include <cmath>
#include <limits>
#include <string>
#include <numeric>
#include <cstdint>
#include <algorithm>
#include "infoxtr.h"

// Wrapper function to calculate information imbalance with all coordinates data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppInfoImbalance(
    const Rcpp::NumericMatrix& Mx,
    const Rcpp::NumericMatrix& My,
    const Rcpp::NumericVector& alpha,
    const Rcpp::IntegerVector& lib,
    const Rcpp::IntegerVector& pred,
    int k = 3,
    int threads = 1,
    const std::string& method = "euclidean")
{
    std::vector<std::vector<double>> mx = infoxtr::convert::mat_r2std(Mx, true);
    std::vector<std::vector<double>> my = infoxtr::convert::mat_r2std(My, true);

    std::vector<double> alpha_std = Rcpp::as<std::vector<double>>(alpha);

    const int n_obs = Mx.nrow(); 

    // Convert and check that lib and pred indices are within bounds & convert R based 1 index to C++ based 0 index
    std::vector<size_t> lib_std;
    lib_std.reserve(lib.size());
    for (int i = 0; i < lib.size(); ++i) {
        if (lib[i] < 1 || lib[i] > n_obs) {
            Rcpp::stop("lib contains out-of-bounds index at position %d (value: %d)", i + 1, lib[i]);
        }
        lib_std.push_back(static_cast<size_t>(lib[i] - 1));
    }

    std::vector<size_t> pred_std;
    pred_std.reserve(pred.size());
    for (int i = 0; i < pred.size(); ++i) {
        if (pred[i] < 1 || pred[i] > n_obs) {
            Rcpp::stop("pred contains out-of-bounds index at position %d (value: %d)", i + 1, pred[i]);
        }
        pred_std.push_back(static_cast<size_t>(pred[i] - 1));
    }
    
    return infoxtr::infoimbalance::infoImbalance(
                mx, my, alpha_std, lib_std, pred_std,
                static_cast<size_t>(std::abs(k)), 
                static_cast<size_t>(std::abs(threads)), method);
}

// Wrapper function to calculate information imbalance gain
// [[Rcpp::export(rng = false)]]
double RcppInfoImbalanceGain(
    const Rcpp::NumericMatrix& mat,
    const Rcpp::IntegerVector& target,
    const Rcpp::IntegerVector& agent,
    const Rcpp::NumericVector& alpha,
    const Rcpp::IntegerVector& lib,
    const Rcpp::IntegerVector& pred,
    const Rcpp::IntegerVector& E = Rcpp::IntegerVector::create(3),
    const Rcpp::IntegerVector& tau = Rcpp::IntegerVector::create(1),
    int style = 1,
    int h = 0,
    int k = 3,
    int threads = 1,
    const std::string& method = "euclidean",
    Rcpp::Nullable<Rcpp::List> nb = R_NilValue,
    Rcpp::Nullable<int> nrows = R_NilValue)
{
    std::vector<std::vector<double>> m = infoxtr::convert::mat_r2std(mat, false);
    const size_t n_cols = m.size();

    // ============================================================================
    // Convert target / agent and construct the corresponding E / tau.
    //
    // IMPORTANT:
    // E and tau are associated with the ORIGINAL target / agent positions.
    // Therefore E / tau are expanded BEFORE sorting and removing duplicates.
    //
    // Final order:
    //     target: sorted unique variable indices
    //     agent : sorted unique variable indices
    //
    // The corresponding E / tau entries are carried together with each variable.
    // ============================================================================

    // ----------------------------------------------------------------------------
    // 1. Convert E and tau
    // ----------------------------------------------------------------------------
    std::vector<int> E_std = Rcpp::as<std::vector<int>>(E);
    std::vector<int> tau_std = Rcpp::as<std::vector<int>>(tau);

    if (E_std.empty()) {
        Rcpp::stop(
            "E must contain at least one value."
        );
    }

    if (tau_std.empty()) {
        Rcpp::stop(
            "tau must contain at least one value."
        );
    }

    if (E_std.size() != tau_std.size()) {
        Rcpp::stop(
            "E and tau must have the same length."
        );
    }

    // E must be >= 2
    for (size_t i = 0; i < E_std.size(); ++i) {
        if (E_std[i] < 2) {
            Rcpp::stop(
                "All E values must be >= 2; "
                "E[%d] = %d.",
                static_cast<int>(i + 1),
                E_std[i]
            );
        }
    }

    // ----------------------------------------------------------------------------
    // 2. Raw target / agent
    // ----------------------------------------------------------------------------
    std::vector<size_t> target_raw = Rcpp::as<std::vector<size_t>>(target);
    std::vector<size_t> agent_raw = Rcpp::as<std::vector<size_t>>(agent);

    if (target_raw.empty()) {
        Rcpp::stop(
            "target must contain at least one variable."
        );
    }

    if (agent_raw.empty()) {
        Rcpp::stop(
            "agent must contain at least one variable."
        );
    }

    // ----------------------------------------------------------------------------
    // 3. Check indices and convert R 1-based -> C++ 0-based
    // ----------------------------------------------------------------------------
    for (size_t i = 0; i < target_raw.size(); ++i) {
        if (target_raw[i] < 1 || target_raw[i] > n_cols){
            Rcpp::stop(
                "Target index %d out of bounds [1, %d].",
                static_cast<int>(i + 1),
                static_cast<int>(n_cols)
            );
        }
        target_raw[i] -= 1;
    }

    for (size_t i = 0; i < agent_raw.size(); ++i) {
        if (agent_raw[i] < 1 || agent_raw[i] > n_cols)
        {
            Rcpp::stop(
                "Agent index %d out of bounds [1, %d].",
                static_cast<int>(i + 1),
                static_cast<int>(n_cols)
            );
        }
        agent_raw[i] -= 1;
    }

    // ----------------------------------------------------------------------------
    // 4. Expand E / tau according to the ORIGINAL target / agent lengths
    //
    // Cases:
    //
    //   length == 1:
    //       one value recycled to all target and agent variables.
    //
    //   length == 2:
    //       E[0] / tau[0] -> target
    //       E[1] / tau[1] -> agent
    //
    //   length >= 3:
    //       last value is reserved for agent;
    //       preceding values form the target sequence;
    //       target uses R-style recycling if necessary;
    //       unused preceding values are passed to agent;
    //       the reserved last value is appended to agent;
    //       agent uses R-style recycling if necessary.
    // ----------------------------------------------------------------------------
    const size_t nt_raw = target_raw.size();
    const size_t na_raw = agent_raw.size();
    const size_t nE = E_std.size();

    std::vector<int> E_target_raw(nt_raw);
    std::vector<int> tau_target_raw(nt_raw);
    std::vector<int> E_agent_raw(na_raw);
    std::vector<int> tau_agent_raw(na_raw);

    if (nE == 1) {
        // ----------------------------------------------------------------------------
        // Case 1: one E / tau value
        // ----------------------------------------------------------------------------
        for (size_t i = 0; i < nt_raw; ++i) {
            E_target_raw[i] = E_std[0];
            tau_target_raw[i] = tau_std[0];
        }

        for (size_t i = 0; i < na_raw; ++i) {
            E_agent_raw[i] = E_std[0];
            tau_agent_raw[i] = tau_std[0];
        }
    } else if (nE == 2) {
        // ----------------------------------------------------------------------------
        // Case 2: two E / tau values
        //
        //     first  -> target
        //     second -> agent
        // ----------------------------------------------------------------------------
        for (size_t i = 0; i < nt_raw; ++i) {
            E_target_raw[i] = E_std[0];
            tau_target_raw[i] = tau_std[0];
        }

        for (size_t i = 0; i < na_raw; ++i) {
            E_agent_raw[i] = E_std[1];
            tau_agent_raw[i] = tau_std[1];
        }
    } else {
        // ----------------------------------------------------------------------------
        // Case 3: three or more E / tau values
        //
        // Example:
        //
        //     E = c(2, 3, 4, 5)
        //
        //     last value:
        //         5 -> reserved for agent
        //
        //     front:
        //         2, 3, 4 -> target first
        //
        // If target has 2 variables:
        //
        //     target = 2, 3
        //     agent pool = 4, 5
        //
        // If target has 5 variables:
        //
        //     target = 2, 3, 4, 2, 3
        //     agent pool = 5
        //
        // If target consumes fewer front values, all remaining front values are
        // transferred to the agent pool, followed by the final reserved value.
        // ----------------------------------------------------------------------------
        const size_t n_front = nE - 1;

        // ------------------------------------------------------------------------
        // Target
        // ------------------------------------------------------------------------
        for (size_t i = 0; i < nt_raw; ++i) {
            const size_t j = i % n_front;
            E_target_raw[i] = E_std[j];
            tau_target_raw[i] = tau_std[j];
        }

        // ------------------------------------------------------------------------
        // Agent pool
        // ------------------------------------------------------------------------
        std::vector<int> E_agent_pool;
        std::vector<int> tau_agent_pool;

        // Remaining front values after target has consumed nt_raw values.
        if (n_front > nt_raw) {
            for (size_t i = nt_raw; i < n_front; ++i) {
                E_agent_pool.push_back(E_std[i]);
                tau_agent_pool.push_back(tau_std[i]);
            }
        }

        // Last value is always reserved for agent.
        E_agent_pool.push_back(E_std[nE - 1]);
        tau_agent_pool.push_back(tau_std[nE - 1]);

        // ------------------------------------------------------------------------
        // Agent with R-style recycling
        // ------------------------------------------------------------------------
        const size_t n_agent_pool = E_agent_pool.size();

        for (size_t i = 0; i < na_raw; ++i) {
            const size_t j = i % n_agent_pool;
            E_agent_raw[i] = E_agent_pool[j];
            tau_agent_raw[i] = tau_agent_pool[j];
        }
    }

    // ============================================================================
    // 5. Bind variable index with its E / tau BEFORE sorting and deduplication
    // ============================================================================
    struct VarEmbedParam
    {
        size_t index;
        int E;
        int tau;
        size_t original_position;
    };

    std::vector<VarEmbedParam> target_params;
    target_params.reserve(nt_raw);

    for (size_t i = 0; i < nt_raw; ++i) {
        target_params.push_back(
            {
                target_raw[i],
                E_target_raw[i],
                tau_target_raw[i],
                i
            }
        );
    }

    std::vector<VarEmbedParam> agent_params;
    agent_params.reserve(na_raw);

    for (size_t i = 0; i < na_raw;  ++i) {
        agent_params.push_back(
            {
                agent_raw[i],
                E_agent_raw[i],
                tau_agent_raw[i],
                i
            }
        );
    }

    // ============================================================================
    // 6. Sort target / agent by variable index
    //
    // original_position is used as a secondary key so that duplicated variables
    // retain their original input order.
    // ============================================================================
    std::sort(
        target_params.begin(),
        target_params.end(),
        [](const VarEmbedParam& a,
        const VarEmbedParam& b)
        {
            if (a.index != b.index) {
                return a.index < b.index;
            }

            return a.original_position <
                b.original_position;
        }
    );

    std::sort(
        agent_params.begin(),
        agent_params.end(),
        [](const VarEmbedParam& a,
        const VarEmbedParam& b)
        {
            if (a.index != b.index) {
                return a.index < b.index;
            }

            return a.original_position <
                b.original_position;
        }
    );

    // ============================================================================
    // 7. Remove duplicates while retaining the corresponding E / tau
    //
    // If the same variable occurs multiple times, the first occurrence in the
    // original input order is retained together with its E / tau.
    // ============================================================================
    std::vector<size_t> tg;
    std::vector<int> E_target;
    std::vector<int> tau_target;

    tg.reserve(target_params.size());
    E_target.reserve(target_params.size());
    tau_target.reserve(target_params.size());

    for (const auto& x : target_params) {
        if (!tg.empty() && tg.back() == x.index) continue;
        tg.push_back(x.index);
        E_target.push_back(x.E);
        tau_target.push_back(x.tau);
    }

    std::vector<size_t> ag;
    std::vector<int> E_agent;
    std::vector<int> tau_agent;

    ag.reserve(agent_params.size());
    E_agent.reserve( agent_params.size());
    tau_agent.reserve(agent_params.size());

    for (const auto& x : agent_params){
        if (!ag.empty() && ag.back() == x.index) continue;
        ag.push_back(x.index);
        E_agent.push_back(x.E);
        tau_agent.push_back(x.tau);
    }

    std::vector<size_t> tg = Rcpp::as<std::vector<size_t>>(target);
    std::vector<size_t> ag = Rcpp::as<std::vector<size_t>>(agent);

    for (auto& idx : tg) {
        if (idx < 1 || idx > n_cols) {
            Rcpp::stop("Target index %d out of bounds [1, %d]", 
                       static_cast<int>(idx), 
                       static_cast<int>(n_cols));
        }
        idx -= 1;  // to 0-based
    }
    std::sort(tg.begin(), tg.end());
    tg.erase(
        std::unique(tg.begin(), tg.end()),
        tg.end()
    );

    for (auto& idx : ag) {
        if (idx < 1 || idx > n_cols) {
            Rcpp::stop("Interact index %d out of bounds [1, %d]", 
                       static_cast<int>(idx), 
                       static_cast<int>(n_cols));
        }
        idx -= 1;  // to 0-based
    }
    std::sort(ag.begin(), ag.end());
    ag.erase(
        std::unique(ag.begin(), ag.end()),
        ag.end()
    );

    std::vector<double> alpha_std = Rcpp::as<std::vector<double>>(alpha);

    // Generate shadow manifolds for target/agent variables
    std::vector<std::vector<double>> mx;
    std::vector<std::vector<double>> my;

    if (nb.isNotNull()) 
    {
        // Convert Rcpp::List to std::vector<std::vector<size_t>>
        std::vector<std::vector<size_t>> nb_std = infoxtr::convert::nb2std(nb.get());
        lagged_values = infoxtr::lagg::lagg(
            cppMat, nb_std, static_cast<size_t>(std::abs(lag)), false);
    } 
    else if (nrows.isNotNull())
    {
        lagged_values = infoxtr::lagg::lagg(
            cppMat, 
            static_cast<size_t>(std::abs(Rcpp::as<int>(nrows))), 
            static_cast<size_t>(std::abs(lag)), false);
    }
    else  
    {
        lagged_values = infoxtr::lagg::lagg(
            cppMat, static_cast<size_t>(std::abs(lag)), false);
    }

    const int n_obs = mx.size(); 

    // Convert and check that lib and pred indices are within bounds & convert R based 1 index to C++ based 0 index
    std::vector<size_t> lib_std;
    lib_std.reserve(lib.size());
    for (int i = 0; i < lib.size(); ++i) {
        if (lib[i] < 1 || lib[i] > n_obs) {
            Rcpp::stop("lib contains out-of-bounds index at position %d (value: %d)", i + 1, lib[i]);
        }
        lib_std.push_back(static_cast<size_t>(lib[i] - 1));
    }

    std::vector<size_t> pred_std;
    pred_std.reserve(pred.size());
    for (int i = 0; i < pred.size(); ++i) {
        if (pred[i] < 1 || pred[i] > n_obs) {
            Rcpp::stop("pred contains out-of-bounds index at position %d (value: %d)", i + 1, pred[i]);
        }
        pred_std.push_back(static_cast<size_t>(pred[i] - 1));
    }
    
    return infoxtr::infoimbalance::infoImbalanceGain(
                mx, my, alpha_std, lib_std, pred_std,
                static_cast<size_t>(std::abs(k)), 
                static_cast<size_t>(std::abs(threads)), method);
}
