#include <vector>
#include <cmath>
#include <limits>
#include <string>
#include <numeric>
#include <cstdint>
#include <algorithm>
#include "infoxtr.h"

// Wrapper function to quantify interventional causality by knockoff operation
// [[Rcpp::export(rng = false)]]
Rcpp::List RcppKOCMI(const Rcpp::NumericMatrix& mat,
                     const Rcpp::IntegerVector& target,
                     const Rcpp::IntegerVector& agent,
                     const Rcpp::IntegerVector& conds,
                     const Rcpp::NumericMatrix& knockoff,
                     Rcpp::Nullable<Rcpp::NumericMatrix> null_knockoff = R_NilValue,
                     const std::string& type = "cont",
                     int nboots = 10000,
                     int k = 5,
                     int alg = 0,
                     int threads = 1,
                     int seed = 123456789,
                     double base = 2.0,
                     const std::string& method = "equal",
                     bool contain_null = true)
{   
    if (contain_null || !null_knockoff.isNotNull())
    {
        Rcpp::stop("[KOCMI] When `contain_null` is true, the `null_knockoff` matrix for the source variable must be supplied.")
    }

    const size_t n_cols = static_cast<size_t>(mat.ncol());
    const size_t n_obs = static_cast<size_t>(mat.nrow());
    
    size_t tg_idx = target[0];
    if (tg_idx < 1 || tg_idx > n_cols) {
        Rcpp::stop("Target index %d out of bounds [1, %d]", 
                   static_cast<int>(tg_idx), 
                   static_cast<int>(n_cols));
    }
    tg_idx -= 1; // to 0-based

    size_t ag_idx = agent[0];
    if (ag_idx < 1 || ag_idx > n_cols) {
        Rcpp::stop("Agent index %d out of bounds [1, %d]", 
                   static_cast<int>(ag_idx), 
                   static_cast<int>(n_cols));
    }
    ag_idx -= 1; // to 0-based
    
    std::vector<size_t> cg = Rcpp::as<std::vector<size_t>>(conds);
    for (auto& idx : cg) {
        if (idx < 1 || idx > n_cols) {
            Rcpp::stop("Conditioning index %d out of bounds [1, %d]", 
                       static_cast<int>(idx), 
                       static_cast<int>(n_cols));
        }
        idx -= 1;  // to 0-based
    }
    std::sort(cg.begin(), cg.end());
    cg.erase(
        std::unique(cg.begin(), cg.end()),
        cg.end()
    );
    if (cg.empty())
        Rcpp::stop("Conditioning indices should not be empty");

    // Convert original values from R to C++
    std::vector<double> tg_std(n_obs);
    std::vector<double> ag_std(n_obs);
    for (size_t r = 0; r < n_obs; ++r)
    {
        tg_std[r] = mat(r, tg_idx);
        ag_std[r] = mat(r, ag_idx);
    }
    std::vector<std::vector<double>> cg_mat(
        n_obs, std::vector<double>(cg.size())
    );
    for (size_t j = 0; j < cg.size(); ++j)
    {   
        for (size_t r = 0; r < n_obs; ++r)
        {
            cg_mat[j][r] = mat(r, cg[j]);
        }
    }
    
    std::vector<std::vector<double>> nkm;
    if (null_knockoff.isNotNull())
    {
        nkm = infoxtr::convert::mat_r2std(null_knockoff, false);
    }
    std::vector<std::vector<double>> km = infoxtr::convert::mat_r2std(knockoff, false);
    
    // Initialize result container
    infoxtr::kocmi::KOCMIRes res;
   
    if (type == "cont")
    {
        res = infoxtr::kocmi::kocmi(
            tg_std, ag_std, cg_mat, km, nkm,
            static_cast<size_t>(std::abs(nboots)),
            static_cast<size_t>(std::abs(k)),
            static_cast<size_t>(std::abs(alg)),
            static_cast<size_t>(std::abs(threads)),
            static_cast<uint64_t>(std::abs(seed)),
            contain_null
        );
    }
    else  
    {  
        std::vector<uint64_t> tg_vec = infoxtr::discretize::discretize(
            tg_std, method, static_cast<size_t>(std::abs(k))
        );

        std::vector<uint64_t> ag_vec = infoxtr::discretize::discretize(
            ag_std, method, static_cast<size_t>(std::abs(k))
        );

        std::vector<std::vector<uint64_t>> cg_discm(
            cg.size(), std::vector<uint64_t>(n_obs)
        );
        for (size_t j = 0; j < cg.size(); ++j)
        {   
            cg_discm[j] = infoxtr::discretize::discretize(
                cg_mat[j], method, static_cast<size_t>(std::abs(k))
            );    
        }

        std::vector<std::vector<uint64_t>> kdiscm(
            km.size(), std::vector<uint64_t>(n_obs)
        );
        for (size_t j = 0; j < km.size(); ++j)
        {   
            kdiscm[j] = infoxtr::discretize::discretize(
                km[j], method, static_cast<size_t>(std::abs(k))
            );    
        }

        std::vector<std::vector<uint64_t>> nkdiscm(
            nkm.size(), std::vector<uint64_t>(n_obs)
        );
        for (size_t j = 0; j < nkm.size(); ++j)
        {   
            nkdiscm[j] = infoxtr::discretize::discretize(
                nkm[j], method, static_cast<size_t>(std::abs(k))
            );    
        }

        res = infoxtr::kocmi::kocmi(
            tg_vec, ag_cev, cg_discm, kdiscm, nkdiscm,
            static_cast<size_t>(std::abs(nboots)),
            static_cast<size_t>(std::abs(threads)),
            static_cast<uint64_t>(std::abs(seed)),
            base, contain_null
        );
    }
    

    // Construct discrete data matrix
    std::vector<std::vector<uint64_t>> pm(
        ag.size() + 1, std::vector<uint64_t>(n_obs, 0)
    );
    
    // Preserve original values in target variable and discretize it
    std::vector<double> vec(n_obs);
    for (size_t r = 0; r < n_obs; ++r)
    {
        vec[r] = mat(r, tg_idx);
    }
    pm[0] = infoxtr::discretize::discretize(
        vec, method, static_cast<size_t>(std::abs(n))
    );


    // Discrete lagged values for agent variables
    for (size_t j = 0; j < lagged_values.size(); ++j)
    {   
        
    }

    infoxtr::surd::SURDRes res = infoxtr::surd::surd(
        pm, static_cast<size_t>(std::abs(max_order)),
        static_cast<size_t>(std::abs(threads)), base, normalize);

    std::vector<std::string> names;
    std::vector<std::string> types;
    std::vector<double> values;

    auto make_name = [](const std::vector<size_t>& vars)
    {
        std::string nm;

        for (size_t j = 0; j < vars.size(); ++j)
        {
            if (j > 0)
                nm += "_";

            nm += "V";
            nm += std::to_string(vars[j]);
        }

        return nm;
    };

    /**************************************************
     * Unique
     **************************************************/

    for (size_t i = 0; i < res.unique_vals.size(); ++i)
    {
        names.push_back(make_name(res.unique_vars[i]));
        types.push_back("U");
        values.push_back(res.unique_vals[i]);
    }

    /**************************************************
     * Redundant
     **************************************************/

    for (size_t i = 0; i < res.redundant_vals.size(); ++i)
    {
        names.push_back(make_name(res.redundant_vars[i]));
        types.push_back("R");
        values.push_back(res.redundant_vals[i]);
    }

    /**************************************************
     * Synergy
     **************************************************/

    for (size_t i = 0; i < res.synergy_vals.size(); ++i)
    {
        names.push_back(make_name(res.synergy_vars[i]));
        types.push_back("S");
        values.push_back(res.synergy_vals[i]);
    }

    /**************************************************
     * InfoLeak
     **************************************************/

    names.push_back("InfoLeak");
    types.push_back("InfoLeak");
    values.push_back(res.info_leak);

    Rcpp::CharacterVector names_r(names.begin(), names.end());
    Rcpp::CharacterVector types_r(types.begin(), types.end());
    Rcpp::NumericVector values_r(values.begin(), values.end());

    return Rcpp::List::create(
        Rcpp::Named("vars")  = names_r,
        Rcpp::Named("types")  = types_r,
        Rcpp::Named("values") = values_r
    );
}
