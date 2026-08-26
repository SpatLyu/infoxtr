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
    
    return infoxtr::infoimbalance::Infoimbalance(
                mx, my, alpha_std, lib_std, pred_std,
                static_cast<size_t>(std::abs(k)), 
                static_cast<size_t>(std::abs(threads)), method);
}

// Wrapper function to calculate information imbalance gain
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppInfoImbalanceGain(
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
    
    return infoxtr::infoimbalance::Infoimbalance(
                mx, my, alpha_std, lib_std, pred_std,
                static_cast<size_t>(std::abs(k)), 
                static_cast<size_t>(std::abs(threads)), method);
}
