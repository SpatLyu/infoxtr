#' Mutual Information
#'
#' Estimate the mutual information between target and interacting variables.
#'
#' @inheritParams mi
#' @param agent 
#' @param lag_p (optional) 
#' @param lag_q (optional) 
#'
#' @returns A numerical value.
#' @export
#'
#' @examples
#' set.seed(42)
#' infocaus::te(matrix(stats::rnorm(100,1,10),ncol=2),1,2)
#'
te = \(data, target, agent, lag_p = 3, lag_q = 3, base = exp(1), 
       type = c("cont", "disc"), k = 3, normalize = FALSE) {
  type = match.arg(type)
  mat = as.matrix(data)
  if (type == "disc") {
    return(RcppDiscTE(mat, target, agent, lag_p, lag_q, base, TRUE, normalize))
  } else {
    return(RcppContTE(mat, target, agent, lag_p, lag_q, k, 0, base, normalize))
  }
}
