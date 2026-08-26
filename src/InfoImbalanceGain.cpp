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
    int k,
    int threads = 1,
    const std::string& method = "euclidean")
{
    std::vector<std::vector<double>> m = infoxtr::convert::mat_r2std(mat, false);

    std::vector<size_t> tg = Rcpp::as<std::vector<size_t>>(target);
    std::vector<size_t> ag = Rcpp::as<std::vector<size_t>>(agent);

    const size_t n_cols = m.size();
    for (auto& idx : tg) {
        if (idx < 1 || idx > n_cols) {
            Rcpp::stop("Target index %d out of bounds [1, %d]", 
                       static_cast<int>(idx), 
                       static_cast<int>(n_cols));
        }
        idx -= 1;  // to 0-based
    }
    for (auto& idx : ag) {
        if (idx < 1 || idx > n_cols) {
            Rcpp::stop("Interact index %d out of bounds [1, %d]", 
                       static_cast<int>(idx), 
                       static_cast<int>(n_cols));
        }
        idx -= 1;  // to 0-based
    }
    
    return infoxtr::transferentropy::transferEntropy(
                m, tg, ag, 
                static_cast<size_t>(std::abs(lag_p)), 
                static_cast<size_t>(std::abs(lag_q)), 
                static_cast<size_t>(std::abs(k)), 
                static_cast<size_t>(std::abs(alg)), 
                std::abs(base), normalize, lag_single);
}