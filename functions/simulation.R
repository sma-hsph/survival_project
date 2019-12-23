simulation <- function(b, a, # regression parameters.
                       study, # length-n vector specifying studies
                       missingness, # length-k vector indicating missingness
                       dist_eT, sd_eT, # length-k vectors for distribution of error terms
                       percCens, # length-k vectors for distribution of censoring variable
                       seed
) {
  set.seed(seed)
  # establishing dimensions of data
  b <- matrix(b, ncol=1)
  if(is.vector(a)) a <- matrix(a, ncol=1)
  pall <- nrow(b)
  pZ <- ncol(a)
  pX <- pall - pZ
  if(pX != nrow(a)) stop('a should be pX x pZ matrix!')
  k <- length(unique(study))
  if(k != length(missingness) | 
     k != length(dist_eT) | 
     k != length(sd_eT) |
     k != length(percCens))
    stop("Number of studies indicated by study, dist_eT, and sd_eT should agree!")
  
  dfData <- lapply((1:k), function(i) {
    iStudy <- unique(study)[i]
    n_tmp <- sum(study %in% iStudy)
    matX <- MASS::mvrnorm(n_tmp, mu=rep(0, pX), Sigma=diag(rep(1, pX))) # could induce correlation structure
    matZ <- matX %*% a + MASS::mvrnorm(n_tmp, mu=rep(0, pZ), Sigma=diag(rep(1, pZ))) # could induce correlation structure
    eT <- switch(dist_eT[i],
                 normal = rnorm(n_tmp),
                 evd = evd::rgev(n_tmp))
    eT <- (eT - mean(eT)) / sd(eT) * sd_eT[i] + mean(eT) # scale to specified standard deviation
    logT <- as.vector((cbind(matX, matZ) %*% b + eT))
    
    # generate censoring
    muCens <- 50
    sdCens <- 10
    obsCens <- 0
    while(obsCens < percCens[i]) {
      # Keep changing censoring dist'n until the censoring percentage gets close to the disired percentage.
      muCens <- muCens - 0.25
      logC <- rnorm(n_tmp, muCens, sdCens )
      delta <- 1 * (logT <= logC)
      obsCens <- 1 - sum(delta) / sum(n_tmp)
    }
    logY <- pmin(logT, logC)
    
    data.frame(matX, matZ,
               logY = logY, delta = delta, logT = logT,
               logC = logC, study = iStudy, missing = missingness[i])
  })
  dfData <-  Reduce('rbind', dfData)
  
  colnames(dfData)[1:pX] <- paste0("X", 1:pX)
  colnames(dfData)[(pX + 1):(pX + pZ)] <- paste0("Z", 1:pZ)
  return(dfData)
}

simulation_sanitycheck <- function(tb_sim) {
  tb_check <- tb_sim %>% 
    dplyr::mutate(check_k = purrr::map_dbl(missingness, length) == nStudies &
                    purrr::map_dbl(dist_eT, length) == nStudies &
                    purrr::map_dbl(sd_eT, length) == nStudies &
                    purrr::map_dbl(percCens, length) == nStudies)
  if(any(!tb_check$check_k))
    stop("Number of studies doesn't agree between nStudies, missingness, dist_eT, sd_eT, and percCens!")
  tb_check <- tb_sim %>% 
    dplyr::mutate(check_p = purrr::map_chr(a, class) == "matrix" &
                    purrr::map_dbl(a, ~ nrow(.x) + ncol(.x))  == purrr::map_dbl(b, length))
  if(any(!tb_check$check_p))
    stop("Covariate dimension doesn't agree between b and a!")
  print("Sanity check passed.")
}
