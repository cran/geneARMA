\name{geneARMAsim}
\Rdversion{1.1}
\alias{geneARMAsim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate a periodic model of time-course gene expression data.}
\description{
	A function which generates data to test the geneARMAfit results.
}
\usage{
geneARMAsim(n = 100, J = 2, K = 2, sigsq = 0.01, ars = numeric(0), mas = numeric(0), length.out = 21, fourier.coefs, tau)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
	simulated sample size.
}
  \item{J}{
	number of clusters.
}
  \item{K}{
	number of Fourier terms.
}
  \item{sigsq}{
	innovation variance.
}
  \item{ars}{
	a vector of AR coefficients, where the time series is parameterized as X[t] = ars[1] * X[t-1] + ars[2]*X[t-2] + ... + e[t]+ mas[1]*e[t-1] + ...
}
  \item{mas}{
	a vector of MA coefficients, where the time series is parameterized as X[t] = ars[1] * X[t-1] + ars[2]*X[t-2] + ... + e[t]+ mas[1]*e[t-1] + ...
}
  \item{length.out}{
	number of equally spaced time points at which each gene is obsered.  
}
  \item{fourier.coefs}{
	A J by 2K+1 matrix with fourier coefficients.  The rows are clusters, and the columns are constant term, cosine term, sine term, cosine term, sine term,... where lower frequency terms come first.
}
  \item{tau}{
	A vector of length J specifying the frequencies for the different clusters.  
}
}
\details{
}
\value{
Returns a list with two components
\item{Y}{The simulated data}
\item{tm}{The corresponding time points}
}
\references{
	Ning Li, et al.  Functional clustering of periodic transcriptional profiles through ARMA(p,q)
}
\author{
Timothy McMurry and Arthur Berg
}
\note{
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
set.seed(100)
Data <- geneARMAsim(400, ars=c(.5, .1))
f1 <- geneARMAfit(Data$Y, Data$tm, 2, 2, 2, 0, eps.conv = .001, max.iter = 15, tau.init=c(.25, .45))
plot(f1, y=NULL, "all.means")
plot(f1, y=NULL, "single.cluster", j=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

