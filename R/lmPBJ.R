#' Computes Statistical Map for Neuroimaging Data
#'
#' This function computes a statistical map and residuals which are the
#'  objects necessary to perform the parametric bootstrap joint (PBJ) inference
#'  procedure.
#' @param images Character vector of subject images to be modeled as an outcome
#'  variable OR 4d array of imaging data OR 4d nifti object.
#' @param form formula or character that can be coerced to a formula or a design
#' matrix for the full model.
#' @param formred formula or character that can be coerced to a formula or a
#' design matrix for reduced model. If robust=TRUE then this must have one less
#' column than X.
#' @param mask File name for a binary mask file or niftiImage object.
#' @param id Vector to identify measurements from the same observation.
#' @param data R data frame containing variables in form. If form and formred
#' are matrices then this can be NULL.
#' @param W Numeric vector of weights for regression model. Can be used to deweight noisy
#'  observations. Same as what should be passed to lm.
#' @param Winv Inverse weights for regression model. Inverse of W.
#' @param W_structure The working correlation structure if id is provided. "independent" and "exchangeable" are accepted.
#' @param template Template image used for visualization.
#' @param formImages n X p matrix of images where n is the number of subjects and
#'  each column corresponds to an imaging covariate. Currently, not supported.
#' @param robust Logical, compute robust standard error estimates?
#' @param HC3 Logical, Uses HC3 SE estimates from Long and Ervin 2000? Defaults to TRUE.
#' @param transform character indicating type of transformation to use. "none", "t", "f", or "edgeworth" are currently accepted. Edgeworth is slow.
#' @param outdir If specified, output is saved as NIfTIs and statMap object is
#' saved as strings. This approach conserves memory, but has longer IO time.
#' Currently, not supported.
#' @param mc.cores Argument passed to mclapply for parallel things.
#' @param zeros Exclude voxels that have zeros? Zeros may exist due to differences in masking and
#' coverage or they may represent locations where the data took the value zero.
#' @keywords parametric bootstrap, statistical parametric map, semiparametric bootstrap
#' @return Returns a list with the following values:
#' \describe{
#'   \item{stat}{The statistical values where mask!=0. It is a chi^2-statistic map.}
#'   \item{coef}{A 4d niftiImage giving the parameter estimates for covariates only in the full model.}
#'   \item{sqrtSigma}{The covariance object used to sample from the joint distribution of the statistical image.}
#'   \item{mask}{The input mask.}
#'   \item{template}{The background template used for visualization.}
#'   \item{formulas}{A list containing the full and reduced models.}
#' }
#' @importFrom stats coefficients lm pf pt qnorm qchisq residuals
#' @importFrom RNifti writeNifti readNifti
#' @importFrom parallel mclapply
#' @importFrom PDQutils papx_edgeworth
#' @importFrom Matrix bdiag
#' @export
#' @example inst/examples/lmPBJ.R
lmPBJ = function(images, form, formred=~1, mask, id=NULL, data=NULL, W=NULL, W_structure="independent",
                 Winv=NULL, template=NULL, formImages=NULL, robust=TRUE, transform=c('none', 't', 'f', 'edgeworth'),
                 outdir=NULL, zeros=FALSE, HC3=TRUE, mc.cores = getOption("mc.cores", 2L)){

  # hard coded epsilon for rounding errors in computing hat values
  eps=0.001
  transform = tolower(transform[1])
  # N = length(unique(id)) # number of total subjects

  X = getDesign(form, formred, data=data)
  Xred = X[['Xred']]
  df = X[['df']]
  images = images[X[['nas']] ]
  X = X[['X']]
  if(nrow(X) < nrow(data)){
    message(nrow(data)-nrow(X), 'observations deleted due to missingness.')
  }

  n = nrow(X)
  if(class(images[[1]])[1] != 'niftiImage'){
    if(n!=length(images))
      stop('length(images) and nrow(X) must be the same.')
    images = as.character(images)
    # removes white space after images if there is any
    images = gsub(" +$", "", images)
    Y = simplify2array(RNifti::readNifti(images))
  } else {
    Y = simplify2array(images)
    images=NULL
  }
  dims = dim(Y)

  if(is.character(W)){
    stop('Image valued weights are not supported.')
  }

  # check if inverse weights are given
  if(is.null(Winv)){
    Winv = FALSE
    if(is.null(W)) W = rep(1,n)
  } else {
    W = Winv
    Winv = TRUE
  }

  # load mask
  if(class(mask)[1] !='niftiImage'){
    maskimg=as.character(mask)
    mask = RNifti::readNifti(maskimg)
  }

  # check that first input image and mask dimensions are the same
  ndims = length(dim(mask))
  if(any(dims[1:ndims] != dim(mask))  ){
    stop('images and mask dimensions do not match.')
  }

  # check that template and mask dimensions are the same
  if( !is.null(template)){
    if(class(template)[1]!='niftiImage'){
      temp = readNifti(as.character(template))
    } else {
      temp = template
    }
    dims = dim(temp)
    rm(temp)
    if(any(dims[1:ndims] != dim(mask))  ){
      stop('template image and mask dimensions (or pixel dimensions) do not match.')
    }
  }

  # load images
  if(zeros){
    # removes locations where there are any zeros
    mask = mask * c(apply(Y!=0, 1:ndims, all))
  }
  Y = t(apply(Y, (ndims+1), function(x) x[mask!=0]))
  V = ncol(Y)

  # assumes column names in X which aren't in Xred are of interest.
  peind = which(!colnames(X) %in% colnames(Xred))
  # rdf = N - ncol(X) # this is true unless X is rank deficient
  rdf = n - ncol(X)

  # inverse weights were passed
  if(Winv) W[W!=0] = 1/W[W!=0]

  # reorder everything by id
  if(!is.null(id)){
    X = X[order(id),, drop=FALSE]
    Xred = Xred[order(id),, drop=FALSE]
    Y = Y[order(id),, drop=FALSE]
    W = W[order(id)]
  }
  id = sort(id)

  # w is returned in output
  w = W
  X1 = X[,peind]
  W = sqrt(W)
  # this is a point wise matrix multiplication if W was passed as images
  # Y = Y * W

  # fit model to all image data
  QR = qr(X * W)
  # coef = qr.coef(QR, Y * W)[peind,,drop = FALSE]    # coef of X1
  coef = qr.coef(QR, Y * W)  # coef of X = [X0, X1]
  res = qr.resid(QR, Y * W) # residual
  X1res = qr.resid(qr(Xred * W), X1 * W)

  # model fitting for different W structure
  if(W_structure == "independent"){
    # coef0 = NULL # record the coefficients before estimating rho
    summary_rho = 0
    rho = NULL
    Y = Y * W
    W = diag(W)
  }else if(W_structure == "exchangeable"){
    # coef0 = coef

    # split observations by id
    grouped_id = split(seq(nrow(Y)), sort(id))
    # estimate variance of residual
    sigmas = sqrt(colSums(res^2)/rdf)
    # standardized residual
    stand_res = sweep(res, 2, sigmas, FUN ='/')
    # estimate rho
    ni_list = sapply(grouped_id, length) # number of observations for each id
    denominator = sum(ni_list * (ni_list - 1)/2) - ncol(X)  # denominator of the formula

    numerator_per_id = sapply(grouped_id, function(y) { # compute the difference in numerator for each sub-residual matrix
      x = stand_res[y, ]
      if (is.matrix(x) || is.data.frame(x)) {
        return(((colSums(x))^2 - colSums(x^2))/2)
      } else {
        return(0) # when there's only one observation for some id
      }
    }, simplify = "matrix")
    rho_numerator = Reduce('+', numerator_per_id) # summation across ids
    rho = rho_numerator / denominator
    summary_rho = mean(rho) # mean value as the summary

    # close form of square root of the inversed correlation matrix
    Cinv_sqrt = bdiag(lapply(ni_list, function(x) {
      if(x==1){
        inv_sqrt_matrix = 1
      }else{
        a = 1/(1-summary_rho)
        b = a * (-summary_rho/(1+(x-1)*summary_rho))
        Cinv = a *  diag(rep(1,x), ncol=x) +
          b * c(rep(1, x)) %*% t(c(rep(1, x)))
        D = diag(c(rep(a, x-1), a+x*b))
        Dsqrt = diag(sqrt(c(rep(a, x-1), a+x*b)))
        U = matrix(0, nrow = x, ncol = x)
        U[, x] = rep(-1/sqrt(x), x)

        for (i in 1:(x-1)) {
          U[1:i, i] = 1
          U[i+1, i] = -i
          U[, i] = U[, i] / sqrt(sum(U[, i]^2))
        }
        inv_sqrt_matrix = U %*% Dsqrt %*% t(U)
      }
      return(inv_sqrt_matrix)
    }))

    W = as.matrix(Cinv_sqrt %*% diag(as.numeric(W))) # new weight in the second iteration

    # refit model to all image data
    QR = qr(W %*% X)
    Y = as.matrix(W %*% Y)
    # coef = qr.coef(QR, Y)[peind,,drop = FALSE]
    coef = qr.coef(QR, Y)
    res = qr.resid(QR, Y)
    X1res = qr.resid(qr(W %*% Xred), W %*% X1)
  }else{
    stop('Error: W_structure must be "independent" or "exchangeable".')
  }

  # we use non-robust estimator under exchangeable structure
  if (W_structure == "exchangeable" && robust == FALSE) {
    stop('When W_structure is "exchangeable", robust has to be TRUE.')
  }

  if(HC3 & robust){
    h=rowSums(qr.Q(QR)^2); h = ifelse(h>=1, 1-eps, h)
    X1resQ = sweep(simplify2array(rep(list(res/(1-h)), df)),  c(1,3), X1res, '*')
  } else {
    HC3 = FALSE
    # returns nXVXm_1 array
    h=rep(0, n)
    X1resQ = sweep(simplify2array(rep(list(res), df)),  c(1,3), X1res, '*')
  }

  if(!is.null(id)){
    id = factor(id)
    IDmat = model.matrix(~-1+id)
    id = as.integer(id)
    X1resQ = array(apply(X1resQ, 3, function(mat) crossprod(IDmat, mat)), dim=c(ncol(IDmat), V, df))
  }

  # Things needed to resample the robust statistics
  sqrtSigma = list(res=Y, X1res=as.matrix(X1res), X1W=W %*% X1, XredW=W %*% Xred, QR=QR, XW=W %*% X, W=w, coef=coef, rho_avg = summary_rho, rho = rho,
                   n=n, df=df, rdf=rdf, robust=robust, HC3=HC3, transform=transform, id=id, nsub=length(unique(id)))

    # Computes test statistic
  if(robust){
    normedCoef = .Call("pbj_pbjBootRobustX", sqrtSigma$QR, sqrtSigma$res, sqrtSigma$X1res, id, h, df)
  } else {
    sigmas = sqrt(colSums(qr.resid(sqrtSigma$QR, sqrtSigma$res)^2)/(rdf))
    sqrtSigma$res = sweep(sqrtSigma$res, 2, sigmas, FUN = '/')
    # this could be performed outside of the bootstrap function
    AsqrtInv = backsolve(r=qr.R(qr(sqrtSigma$X1res)), x=diag(df) )
    normedCoef = crossprod(AsqrtInv, matrix(sqrtSigma$X1res, nrow=df, ncol=n, byrow=TRUE))
    # used to compute chi-squared statistic
    normedCoef = normedCoef %*% sqrtSigma$res
  }
  sqrtSigma$res = res
  # warnings due to BsqrtInv only existing under some settings
  # suppressWarnings(rm(BsqrtInv,X1resQ, AsqrtInv, Y, res, sigmas, X1res))

  # use transform to compute chi-squared statistic
  normedCoef = switch(transform,
                      none=normedCoef,
                      f=normedCoef,
                      t={ qnorm(pt(normedCoef, df=rdf, log.p = TRUE ), log.p=TRUE )},
                      edgeworth={message('Computing edgeworth transform.')
                        matrix(qnorm(vpapx_edgeworth(stat=normedCoef, mu3=colSums(sqrtSigma^3, dims=1), mu4=colSums(sqrtSigma^4, dims=1) ) ), nrow=df)
                      })
  stat = colSums(normedCoef^2)
  if(transform=='f'){
    stat = qchisq(pf(stat/df, df1=df, df2=rdf, log.p = TRUE), df=df, log.p=TRUE )
  }

  # used later to indicated t-statistic
  out = list(stat=stat, normedCoef=normedCoef, sqrtSigma=sqrtSigma, mask=mask, template=template,
             images=images, formulas=list(full=form, reduced=formred), data = get_all_vars(form, data = data))
  class(out) = c('statMap', 'list')

  # if outdir is specified the stat and sqrtSigma images are saved in outdir
  # and mask tries to get saved as a character.
  if(!is.null(outdir)){
    files = write.statMap(out, outdir)
    out$stat = files$stat
    out$sqrtSigma = files$sqrtSigma
    # if mask was a character then pass that forward instead if the niftiImage
    if(exists('maskimg'))
      out$mask = maskimg
  }
  return(out)
}
