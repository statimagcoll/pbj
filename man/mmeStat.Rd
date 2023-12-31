% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbj.R
\name{mmeStat}
\alias{mmeStat}
\title{Compute maxima, CMI, or CEI inference statistics}
\usage{
mmeStat(stat, rois = FALSE, mask, cft, maxima = FALSE, CMI = FALSE, CEI = TRUE)
}
\arguments{
\item{stat}{statistical image}

\item{rois}{passed to maxima and cluster functions. Returns image with ROI indices.}

\item{mask}{Mask image.}

\item{cft}{A single threshold on the scale of the statistical image (chi-squared) for CEI or CMI.}

\item{maxima}{Compute local maxima?}

\item{CMI}{Compute cluster masses?}

\item{CEI}{Compute cluster extents?}
}
\value{
Returns a list with the maxima and CEI for the given image.
This function is used  as the default `statistic` argument in [pbjInference()].
}
\description{
Compute maxima, CMI, or CEI inference statistics
}
