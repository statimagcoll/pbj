#' A vectorized version of papx_edgeworth
#'
#'
#'
#' @useDynLib pbj, .registration=TRUE
#' @param stat Vector of test statistics
#' @param mu3 The third moment of the test statistics
#' @param mu4 The fourth moment of the test statistics
#' @importFrom PDQutils papx_edgeworth
#' @importFrom PDQutils moment2cumulant
vpapx_edgeworth = Vectorize(function (stat, mu3, mu4) PDQutils::papx_edgeworth(stat, raw.cumulants=PDQutils::moment2cumulant(c(0,1, mu3, mu4) ) ))


#' Computes contiguous clusters from a statistical image given a threshold
#'
#' @param stat A statistical Nifti image as an RNifti image object.
#' @param mask A statistical Nifti image mask used in the analysis or a character path to one.
#' @param cft A vector of cluster forming thresholds (on the scale of the test statistic
#' image, which is usually chi-squared for pbj) for the test statistic image. Will compute cluster sizes or masses for each threshold.
#' @param method character string 'extent' or 'mass' indicating whether the cluster extent or cluster mass statistic should be used.
#' @param kernel The kernel type to compute connected components.
#' @param rois If TRUE, return image with voxel values having the indices of the clusters returned if rois=FALSE.
#' @return Returns list of tables of sizes of the connected components above cft.
#' @export
#' @importFrom mmand shapeKernel
#'
cluster = function(stat, mask, cft, method=c('extent', 'mass'), kernel='box', rois=FALSE){
  method = tolower(method[1])
  if(is.character(mask)) mask = readNifti(mask)
  ndims = length(dim(mask))
  tmp = mask
  k = mmand::shapeKernel(3, ndims, type=kernel)
  tmp = lapply(cft, function(th){ tmp[ mask!=0] = (stat[mask!=0]>th); tmp})
  if(rois){
    ccomps = lapply(tmp, function(tm){cc = mmand::components(tm, k); mask[mask!=0]=0; mask[!is.na(cc)] = cc[!is.na(cc)]; mask} )
  } else {
    ccomps = switch(method,
                    'extent'={
                      lapply(tmp, function(tm) c(0, table(c(mmand::components(tm, k))) ) )
                    },
                    'mass'={
                      lapply(tmp, function(tm) c(0, by(c(stat), c(mmand::components(tm, k)), sum) ))
                    })
    # modifies ccomps attribute in cluster function
    ccomps = lapply(1:length(cft), function(ind){ attributes(ccomps[[ind]]) <- list('cft'=cft[ind]); ccomps[[ind]]})

  }
  if(method=='extent'){
    names(ccomps) = paste('CEI', 1:length(ccomps))
  }else{
    names(ccomps) = paste('CMI', 1:length(ccomps))
  }
  return(ccomps)
}

#' Computes local maxima from an nifti image
#'
#' @param stat A statistical Nifti image as an RNifti image object.
#' @param kernel Type of kernel to use for max/dilation filter
#' @param width Width of kernel (assumes isotropic). If zero, then returns global maximum.
#' @param rois If TRUE, return image with voxel values having the indices of the local maxima.
#' @param roisWidth Width of kernel for computing local maxima used for table output later. Only used if width is zero
#' @return Default returns vector of local maxima in the image.
#' @export
#' @importFrom mmand shapeKernel
#'
maxima = function(stat, kernel='box', width=0, rois=FALSE, roisWidth=15){
  if(width==0 & !rois){
    list(maxima=max(stat))
  } else {
    if(width==0){width = roisWidth}
    ndims = length(dim(stat))
    dil = dilate(stat, shapeKernel(width, ndims, type=kernel))
    stat[which(stat<dil)] = 0
    imginds = which(stat!=0)
    if(rois){
      stat[imginds] = 1:length(imginds)
      list(maxima=stat)
    } else {
      list(maxima=stat[ imginds])
    }
  }
}

#' Computes empirical weighted cdf. Modified from ecdf
#'
#' @param x vector of values
#' @param w vector with length(w)=length(x) of weights
#' @return Returns list of tables of sizes of the connected components above thr.
#' @export
#' @importFrom stats approxfun
#'
wecdf = function (x, w=rep(1, length(x)))
{
  o = order(x)
  x <- x[o]
  w <- w[o]
  n <- length(x)
  sw <- sum(w)
  if (n < 1)
    stop("'x' must have 1 or more non-missing values")
  vals <- unique(x)
  # by command sorts x values (again)
  rval <- approxfun(vals, cumsum(c(by(w, x, sum) ))/sw,
                    method = "constant", yleft = 0, yright = 1-1/n, f = 0, ties = "ordered")
  class(rval) <- c("ecdf", "stepfun", class(rval))
  assign("nobs", n, envir = environment(rval))
  attr(rval, "call") <- sys.call()
  rval
}

#' Gets numeric index of observed stat corresponding to the method and CFT
#'
#' @param obsStat Field in statMap$pbj$obsStat
#' @param method What statistic to provide summary for? must have run that
#' analysis using the pbjInference and mmeStat functions.
#' @param cft cluster forming threshold to display. If NULL, just display the first.
#' @return Returns an index that corresponds to the method and CFT.
#' @seealso [mmeStat], [cluster], [maxima], [pbjInference]
#' @export
inferenceIndex = function(obsStat, method, cft=NULL){
  ind = grep(method, tolower(names(obsStat)))
  if(length(ind)==0){
    stop('Method ', method, ' was not run on this statMap.')
  }
  # get indices corresponding to this method
  cfts = sapply(obsStat[ind], attr, which='cft')
  # maxima don't have cft attribute
  if(!is.null(cft) & method!='maxima'){
    ind = ind[which(cfts==cft)]
    if(length(ind)==0){
      stop('Specified cft is ', cft, '. Existing cfts are ', paste(cfts, collapse=', '))
    }
  }
  ind[1]
}

#' Compute maxima, CMI, or CEI inference statistics
#'
#' @param stat statistical image
#' @param rois passed to maxima and cluster functions. Returns image with ROI indices.
#' @param mask Mask image.
#' @param maxima Compute local maxima?
#' @param CMI Compute cluster masses?
#' @param CEI Compute cluster extents?
#' @param cft A single threshold on the scale of the statistical image (chi-squared) for CEI or CMI.
#' @return Returns a list with the maxima and CEI for the given image.
#' This function is used  as the default `statistic` argument in [pbjInference()].
#' @export
#'
mmeStat = function(stat, rois=FALSE, mask, cft, maxima=FALSE, CMI=FALSE,
                   CEI=TRUE){
  res = c()
  if(maxima){
    res = maxima(stat, rois=rois)
  }
  if(CMI) {
    res = c(res, cluster(stat, mask=mask, cft=cft, rois=rois, method='mass'))
  }
  if(CEI){
    res = c(res, cluster(stat, mask=mask, cft=cft, rois=rois,
                         method='extent'))
  }
  res
}


#' Performs one bootstrap for RESI confidence set construction
#'
#' @param stat bootstrap value of the test statistic image. It is a nifti image array.
#' @param res The RESI estimate computed in the observed data.
#' @param n The sample size/ number of independent sampling units.
#' @param df Numerator degrees of freedom
#' @param rdf Denominator degrees of freedom
#' @param normMethod The value used to normalize the test statistic in each
#'  bootstrap. "param" means it's estimated assuming the test
#'  statistic if F-distributed. "1" means no scaling. "1" is the default.
#' @param mask Mask image identifying where measurements were taken in the image.
#' @export
#'
resics <- function(stat, res, n, df, rdf, normMethod=c( '1', 'param'), mask){
  # This is S hat
  #stat <- RESI::chisq2S(stat[mask!=0], df, N)
  # This is the test statistic in the mask
  lambda_b = stat[mask!=0]
  # This is S tilde squared
  stat = lambda_b/n
  w_max <- c()
  w_min <- c()

  if(normMethod == 1){
    SD <- 1
    boots <- stat - res
    w_max <- max(boots)
    w_min <- min(boots)
  }

  if(normMethod == "param"){
    # The estimated variance of the test statistic
    # SNV: Looks right to me
    varStildeSq <- (df/n)^2 *(2*((df + lambda_b)^2 + (df + 2*lambda_b)*(rdf - 2))*((rdf/df)^2) / ((rdf - 2)^2 * (rdf - 4)))
    # this is the standard error of the test statistic
    SD_b <- sqrt(varStildeSq)
    boots <- stat - res
    boots <- boots / SD_b
    w_max <- max(boots)
    w_min <- min(boots)
  }
  # returns max and min across the image
  return(c(w_min,w_max))
}


#' Uses bootstrap results to construct simultaneous confidence intervals.
#'
#' @param statmap The statMap object returned by lmPBJ using the `resics` function.
#' @param alpha The function returns 1-alpha simultaneous confidence intervals.
#'  @export
#'
sci <- function(statmap, alpha){
  normMethod = statmap$pbj$statArgs$normMethod
  a <- do.call(rbind, statmap$pbj$boots)
  w_min_vec <- a[,1]
  w_max_vec <- a[,2]
  n <- statmap$sqrtSigma$n
  # the observed test statistic image
  lambda <- statmap$stat
  res = lambda/n
  df <- statmap$sqrtSigma$df
  rdf <- statmap$sqrtSigma$rdf
  cu <- quantile(w_max_vec, 1 - alpha / 2)
  cl <- quantile(w_min_vec, alpha / 2)

  if(normMethod == 1){
    SD <- 1
  } else  if(normMethod == "param"){
    # The estimated standard error of the test statistic
    SD <- 2 * sqrt((df/n)^2 *((df + lambda)^2 + (df + 2*lambda)*(rdf - 2))*(rdf/df)^2 / (rdf - 2)^2 / (rdf - 4)  )
  }
  statmap$SCI <- data.frame('SLCI'=res - cu*SD, 'SUCI' = res - cl*SD)
  return(statmap)
}
