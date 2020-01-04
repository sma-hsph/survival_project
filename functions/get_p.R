get_p <- function(coef, stderr) {
  1 - pnorm(abs(coef / stderr))
}
