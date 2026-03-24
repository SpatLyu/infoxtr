je = \(data, indices, base = exp(1), type = c("cont", "disc"), k = 3) {
  type = match.arg(type)
  mat = as.matrix(data)
  if (type == "disc") {
    return(RcppDiscJE(mat, indices, base, TRUE))
  } else {
    return(RcppContJE(mat, indices, k, 0, base))
  }
}

ce = \(data, target, conds, base = exp(1), type = c("cont", "disc"), k = 3) {
  type = match.arg(type)
  mat = as.matrix(data)
  if (type == "disc") {
    return(RcppDiscCE(mat, target, conds, base, TRUE))
  } else {
    return(RcppContCE(mat, target, conds, k, 0, base))
  }
}

mi = \(data, target, interact, base = exp(1), type = c("cont", "disc"), k = 3, normalize = FALSE) {
  type = match.arg(type)
  mat = as.matrix(data)
  if (type == "disc") {
    return(RcppDiscMI(mat, target, interact, base, TRUE, normalize))
  } else {
    return(RcppContMI(mat, target, interact, k, 0, base, normalize))
  }
}
