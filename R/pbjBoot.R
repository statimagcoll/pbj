#' Performs (semi)Parametric Bootstrap Joint ((s)PBJ) Inference
#'
#' @param sqrtSigma list from statmap object
#' @param rboot a function that draws a bootstrapped sample. Should return an n vector. Defaults to Rademacher random variable.
#' @param null Is this a simulation under the null hypothesis?
#' @param method character, method to use for resampling procedure. Wild bootstrap, permutation, or nonparametric
#'
#' @return Returns vector of test statistics computed from the bootstrapped sample.
#' @export
#
pbjBoot = function(sqrtSigma, rboot=function(n){ (2*stats::rbinom(n, size=1, prob=0.5)-1)}, null=TRUE, method=c('wild', 'permutation', 'nonparametric')){
  method = tolower(method[1])
  eps = 0.001
  V = ncol(sqrtSigma$res)
  id = sqrtSigma$id
  n = sqrtSigma$n
  df = sqrtSigma$df
  rdf = sqrtSigma$rdf
  HC3 = sqrtSigma$HC3
  robust = sqrtSigma$robust
  transform = sqrtSigma$transform

  # Leverage or residual-forming matrix construction
  if (HC3 && !is.null(id)) {
    XtX_inv = solve(t(sqrtSigma$XW) %*% sqrtSigma$XW)
    Hmat = sqrtSigma$XW %*% XtX_inv %*% t(sqrtSigma$XW)
    Rmat = diag(n) - Hmat
    sqrtSigma$Rmat = Rmat
    h = rep(NA, n)  # not used in longitudinal case
  } else if (HC3) {
    h = rowSums(qr.Q(sqrtSigma$QR)^2)
    h = ifelse(h >= 1, 1 - eps, h)
  } else {
    h = rep(0, n)
  }


  if(!is.null(id)){ # For longitudinal data
    grouped_id = split(1:n, sort(id))

    if(method=='wild'){
      if (HC3) {
        # Adjust residuals using subject-specific R blocks
        for (i in seq_along(grouped_id)) {
          idx = grouped_id[[i]]
          Ri = sqrtSigma$Rmat[idx, idx, drop = FALSE]
          ei = sqrtSigma$res[idx, , drop = FALSE]
          xi = rboot(length(idx))
          sqrtSigma$res[idx, ] = solve(Ri + diag(eps, length(idx))) %*% (xi * ei)
        }
      } else {
        # Standard wild bootstrap
        sqrtSigma$res = sweep(sqrtSigma$res, 1, rboot(n), '*')
      }

    }else if(method=='permutation'){
      # Permute within each cluster (HC3 not used here)
      sampleIndex = unlist(sapply(grouped_id, function(x) {
        if(length(x) == 1) x else sample(x, size = length(x), replace = TRUE)
      }))
      sqrtSigma$res = sqrtSigma$res[sampleIndex, ]


    } else if (method=='nonparametric'){ # sample within and between clusters
      # Resample clusters and apply R-based correction if HC3
      sampleID = sample(unique(id), size = length(unique(id)), replace=TRUE)
      samp = unlist(lapply(sampleID, function(x) grouped_id[[as.character(x)]])) # collect observations by sampled ID

      sqrtSigma$res = sqrtSigma$res[samp, , drop = FALSE]
      sqrtSigma$X1W = sqrtSigma$X1W[samp, , drop = FALSE]
      sqrtSigma$XredW = sqrtSigma$XredW[samp, , drop = FALSE]
      sqrtSigma$X1res = qr.resid(qr(sqrtSigma$XredW), sqrtSigma$X1W)
      sqrtSigma$XW = sqrtSigma$XW[samp, , drop = FALSE]
      sqrtSigma$QR = qr(sqrtSigma$XW)

      if (HC3) {
        # Recompute R matrix for bootstrap sample
        n_boot = nrow(X)
        XtX_inv = solve(t(sqrtSigma$XW) %*% sqrtSigma$XW)
        Hmat = sqrtSigma$XW %*% XtX_inv %*% t(sqrtSigma$XW)
        Rmat = diag(n) - Hmat
        sqrtSigma$Rmat = Rmat

        grouped_id_boot = split(1:n_boot, rep(1:length(unique(id)), each = length(samp)/length(unique(id))))
        for (i in seq_along(grouped_id_boot)) {
          idx = grouped_id_boot[[i]]
          Ri = Rmat[idx, idx, drop = FALSE]
          ei = sqrtSigma$res[idx, , drop = FALSE]
          sqrtSigma$res[idx, ] = solve(Ri + diag(eps, length(idx))) %*% ei
        }
      }
    }

  }else{
    # Cross-sectional case
    if(method == 'wild'){
      sqrtSigma$res = sweep(sqrtSigma$res, 1, rboot(n)/(1-h), '*')
    } else if(method == 'permutation'){
      sqrtSigma$res = sqrtSigma$res[sample(n), ]
    } else if (method == 'nonparametric'){
      samp = sample(n, replace=TRUE)
      sqrtSigma$res = sweep(sqrtSigma$res[samp,,drop=FALSE], 1, 1-h[samp], '/')
      sqrtSigma$X1W = sqrtSigma$X1W[samp,,drop=FALSE]
      sqrtSigma$XredW = sqrtSigma$XredW[samp,,drop=FALSE]
      sqrtSigma$X1res = qr.resid(qr(sqrtSigma$XredW), sqrtSigma$X1W)
      sqrtSigma$XW = sqrtSigma$XW[samp,,drop=FALSE]
      sqrtSigma$QR = qr(sqrtSigma$XW)
    }
  }


  # for bootstrapping under the alternative
  if(!null) sqrtSigma$res = sqrtSigma$XW %*% sqrtSigma$coef + sqrtSigma$res

  # robust estimator or not
  if(robust){
    if (method == 'nonparametric') {
      h = h[samp]
      # Generate new id vector matching resampled structure
      id = rep(1:length(unique(id)), unlist(lapply(sampleID, function(x) table(id)[[as.character(x)]])))
    }
    statimg = .Call("pbj_pbjBootRobustX", sqrtSigma$QR, sqrtSigma$res, sqrtSigma$X1res, id, h, df)
  } else {
    sigmas = sqrt(colSums(qr.resid(sqrtSigma$QR, sqrtSigma$res)^2)/(rdf))
    sqrtSigma$res = sweep(sqrtSigma$res, 2, sigmas, FUN = '/')
    # this could be performed outside of the bootstrap function
    AsqrtInv = backsolve(r=qr.R(qr(sqrtSigma$X1res)), x=diag(df) )
    statimg = crossprod(AsqrtInv, matrix(sqrtSigma$X1res, nrow=df, ncol=nrow(sqrtSigma$X1res), byrow=TRUE))
    # used to compute chi-squared statistic
    statimg = statimg %*% sqrtSigma$res
  }


  statimg = switch(tolower(transform[1]),
                   none=statimg,
                   f=statimg,
                   t=qnorm(pt(statimg, df=rdf, log.p = TRUE), log.p=TRUE),
                   edgeworth={message('Computing edgeworth transform.')
                     matrix(qnorm(vpapx_edgeworth(stat=statimg, mu3=colSums(sqrtSigma$res^3, dims=1), mu4=colSums(sqrtSigma$res^4, dims=1) ) ), nrow=df)
                   })

  statimg = colSums(statimg^2)

  if(tolower(transform)=='f'){
    statimg = qchisq(pf(statimg/df, df1=df, df2=rdf, log.p = TRUE ), df=df, log.p=TRUE )
  }
  return(statimg)
}
