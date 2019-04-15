evaluate <- function(truth, # a length-p vector of true coefficients
                     l_results # a list of fitted objects under the same simulation setup; each
                     # element should be a list of a coef and Sigma (can be NULL) component 
                     ) {
  p <- length(truth)
  (1:p) %>% 
    purrr::map_dfr(function(i) {
      i_truth <- truth[i]
      i_coef <- l_results %>% 
        purrr::map_dbl(~.x$coef[i])
      if(!is.null(l_results[[1]]$Sigma)) {
        i_Sigma <- l_results %>% 
          purrr::map_dbl(~.x$Sigma[i, i]) %>% 
          sqrt
      } else {
        i_Sigma <- NA
      }
      data.frame(variable = i,
                 mean = mean(i_coef),
                 sd = sd(i_coef),
                 bias = mean(i_coef) - i_truth,
                 MSE = sqrt(mean((i_coef - i_truth)^2)),
                 coverage = mean(i_coef + 1.96*i_Sigma > i_truth & 
                                   i_coef - 1.96*i_Sigma < i_truth)
                 )
    })
}
