infoimbalance = \(mx, my, alpha = seq(0,1,0.1), lib = NULL, pred = NULL, 
                  h = 1, k = 3, threads = 1, method = "euclidean") {
  if (is.null(lib)) lib = seq_len(nrow(mx))
  if (is.null(pred)) pred = lib
  return(RcppInfoImbalance(as.matrix(mx), as.matrix(my), abs(alpha), abs(lib), abs(pred), 
                           abs(h), abs(k), abs(threads), method))
}
