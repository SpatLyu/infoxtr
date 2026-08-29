info_imbalance = \(mx, my, lib = NULL, pred = NULL,
                   k = 3, threads = 1, method = "euclidean") {
  if (is.null(lib)) lib = seq_len(nrow(mx))
  if (is.null(pred)) pred = lib
  return(RcppInfoImbalance(as.matrix(mx), as.matrix(my), 
                           abs(lib), abs(pred),
                           abs(k), abs(threads), method))
}

#' Information Imbalance Gain
#'
#' @param mx Numeric matrix of hypothesised driving variable measurements.
#' @param my Numeric matrix of hypothesised response variable measurements.
#' @param alpha (optional) Scaling parameter weighting the putative driver measurements.
#' @param lib (optional) Library indices.
#' @param pred (optional) Prediction indices.
#' @param h (optional) Prediction horizon.
#' @param k (optional) Number of nearest neighbors when estimating ranks.
#' @param threads (optional) Number of parallel threads.
#' @param method (optional) Distance measure to be used: `"euclidean"`, `"manhattan"`, or `maximum"`.
#'
#' @returns A numeric vector.
#' @export
#' @references
#' Del Tatto, V., Fortunato, G., Bueti, D., Laio, A., 2024. Robust inference of causality in high-dimensional dynamical processes from the Information Imbalance of distance ranks. Proceedings of the National Academy of Sciences 121.
#'
#' @examples
#' set.seed(42)
#' mx = embed(rnorm(100), 3)
#' my = embed(rnorm(100), 3)
#' infoxtr::imbalance_gain(mx, my)
#' 
imbalance_gain = \(mx, my, alpha = seq(0,1,0.1), lib = NULL, pred = NULL,
                   h = 1, k = 3, threads = 1, method = "euclidean") {
  if (is.null(lib)) lib = seq_len(nrow(mx))
  if (is.null(pred)) pred = lib
  return(RcppImbalanceGain(as.matrix(mx), as.matrix(my), 
                           abs(alpha), abs(lib), abs(pred),
                           abs(h), abs(k), abs(threads), method))
}
