gehan.obj <- function(beta, y, delta, matX, wt, study) {
  if(is.null(wt)) wt <- rep(1, length(y))
  residual <- y - as.vector(matX %*% beta)
  sapply(unique(study), function(iStudy) {
    index <- study %in% iStudy
    wt_tmp <- wt[index]
    delta_tmp <- delta[index]
    residual_tmp <- residual[index]

    index_reorder <- order(residual_tmp)
    residual.order <- residual_tmp[index_reorder]
    delta.order <- delta_tmp[index_reorder]
    wt.order <- wt_tmp[index_reorder]

    s1 <- cumsum(delta.order*wt.order)
    s2 <- delta.order*rev(cumsum(rev(wt.order)))
    return(sum((s1-s2)*residual.order*wt.order))
  }) %>% sum %>% return
}

gehan.fit <- function(y, delta, matX, matZ,
                      study, missing, wt = rep(1, length(y)), beta.ini = NULL) {
  if(any(missing) & ncol(matZ) != 0) {
    # estimating alphahat and imputing
    mat_alphahat <- sapply(1:ncol(matZ), function(iCol) {
      lm( matZ[!missing, iCol, drop = F] ~ matX[!missing, , drop = F], weights=wt[!missing] )$coef[-1]
    })
    matZ[missing, ] <- matX[missing, , drop = F] %*% mat_alphahat
  }

  fit <- optim(beta.ini, gehan.obj,
               y = y,
               delta = delta,
               matX = cbind(matX, matZ),
               wt = wt,
               study = study
  )
  fit$par
}

perturbfn <- function(y, delta, matX, matZ, study, missing,
                      matWt, B = 500, beta.ini=NULL)
{
  p <- ncol(matX) + ncol(matZ)

  #perturb and store new beta B times
  matBetaPt <- matrix(NA, B, p) #each column is a beta sample
  for(b in 1:B) matBetaPt[b ,] <- gehan.fit(
    y = y,
    delta = delta,
    matX = matX,
    matZ = matZ,
    study = study,
    missing = missing,
    wt = matWt[,b],
    beta.ini = beta.ini
  )

  return(matBetaPt)
}
