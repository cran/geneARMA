\name{geneARMA-package}
\Rdversion{1.0}
\alias{geneARMA-package}
\alias{geneARMA}
\docType{package}
\title{
Simulate, model, and display data from a time-course microarray experiment with periodic gene expression
}
\description{
Fit models for periodic gene expression to time-course microarray data in a normal mixture model framework with mean approximated by a truncated Fourier series and covariance structure modeled by an ARMA(p,q) process.  Estimation is performed with the EM algorithm. 
}
\details{
\tabular{ll}{
Package: \tab geneARMA\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2009-10-07\cr
License: \tab GPL-3\cr
LazyLoad: \tab yes\cr
}

}
\author{
Timothy McMurry and Arthur Berg

Maintainer: Timothy McMurry <tmcmurry@depaul.edu>
}
\references{
	Ning Li, et al.  Functional clustering of periodic transcriptional profiles through ARMA(p,q)
}
\examples{
set.seed(100)
Data <- geneARMAsim(400, ars=c(.5, .1))
f1 <- geneARMAfit(Data$Y, Data$tm, 2, 2, 2, 0, eps.conv = .001, max.iter = 15, tau.init=c(.25, .45))
plot(f1, y=NULL, "all.means")
plot(f1, y=NULL, "single.cluster", j=2)
}
