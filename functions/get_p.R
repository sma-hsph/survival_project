get_p <- function(coef, stderr) {
  pnorm(-abs(coef / stderr)) * 2
}
