#include <vector>
#include <cmath>
#include <limits>
#include <numeric>
#include <algorithm>
#include "lag.hpp"
#include "DataTrans.h"

// Wrapper function to calculate spatial lag value for vector spatial data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppGenLatticeLag(const Rcpp::NumericMatrix& mat,
                                      const Rcpp::List& nb, 
                                      int lag = 1) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  int numRows = mat.nrow();
  int numCols = mat.ncol();
  std::vector<std::vector<double>> cppMat(numRows, std::vector<double>(numCols));
  for (int r = 0; r < numRows; ++r) {
    for (int c = 0; c < numCols; ++c) {
      cppMat[r][c] = mat(r, c);
    }
  }

  // Convert Rcpp::List to std::vector<std::vector<size_t>>
  std::vector<std::vector<size_t>> nb_std = nb2std(nb);

  // Calculate lagged values
  std::vector<std::vector<double>> lagged_values =
    Lag::GenLatticeLag(cppMat, nb_std, static_cast<size_t>(std::abs(lag)));

  return lagged_values;
}
