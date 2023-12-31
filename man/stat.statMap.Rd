% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statmap.R
\name{stat.statMap}
\alias{stat.statMap}
\title{Returns a statistical niftiImage object from a statMap object}
\usage{
stat.statMap(x, method = c("p", "S", "chisq"))
}
\arguments{
\item{x}{the statMap object to extract a coefficient niftiImage from}

\item{method}{the scale of visualization, chi-squared statistic, effect size (S=RESI), p-value. If df is equal to 1, the maps are scaled by the sign of the coefficient.}
}
\value{
a niftiImage object of the chi-square statistical image.
}
\description{
Returns a statistical niftiImage object from a statMap object
}
\examples{
# loading example data
library(pain21)
pain = pain21()
pdata = pain$data

# fitting regression of images onto study sample size, weights proportional to study sample size
pbjModel2 = lmPBJ(images=pdata$images, form=~n, formred=~1, W = pdata$n, mask=pain$mask, data=pdata)
stat.statMap(pbjModel2)
stat.statMap(pbjModel2, method='chisq')
stat.statMap(pbjModel2, method='S')

}
