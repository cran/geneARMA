\name{geneARMAfit}
\Rdversion{1.1}
\alias{geneARMAfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a periodic model to time-course gene expression data.
}
\description{
Fits a Fourier series model to time-course gene expression data.  Data are classified into one of J groups, each with diffneret mean parameters for the gene expression profiles.  The signal mean is modeled by an order K Fourier series.  The residuals are modeled as a Gaussian ARMA(p,q) process.  Fitting is done by a modified EM algorithm, where some updates are estimated by numerical approximation.
}
\usage{
geneARMAfit(Y, times, J = 4, K = 2, p = 2, q = 0, arma.skip = 10, eps.conv = 0.01, max.iter = Inf, print.updates = TRUE, omega.init, fourier.init, tau.init, sigsq.init)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
	A matrix with the normalized gene expression data.  The rows should be genes and the columns observation times.
}
  \item{times}{
	A vector containing the times at which gene expression was measured.  The order must correspon to the order of the data in Y.
}
  \item{J}{
	The number of clusters to fit.
}
  \item{K}{
	The numer of terms in the Fourier sum.
}
  \item{p}{
	The autoregressive order.
}
  \item{q}{
	The moving-average order.	
}
  \item{arma.skip}{
	The number of iterations the algorithm should run before attempting to fit the time series model.  It is often helpful to let the algorithm do approximate clustering first, and then refine the time series fit.
}
  \item{eps.conv}{
	The convergence criterion.  If the change in all estimated terms is less than eps.conv, the program stops.  
}
  \item{max.iter}{
	The maximum number of iterations to run before stopping.
}
  \item{print.updates}{
  	Turn on/off output produced as the program runs.  
}
  \item{omega.init}{
	Initial cluster probabilities.  Uniform if left unspecified. 
}
  \item{fourier.init}{
	A J by 2K+1 matrix with initial values for the fourier coefficients.  The columns are constant term, cosine term, sine term, cosine term, sine term,... where lower frequency terms come first.  Filled with iid N(0,1) if left unspecified.
}
  \item{tau.init}{
	Initial values for the period of the signal.  Defaults to a range of values covering some fraction of the overall length of time measured.  
}
  \item{sigsq.init}{
	Initial value for the time series innovation variance.  The initial value should be very large, probably hundreds of times the actual, in order for the algorithm to find good clusters.  
}
}
\details{
	An EM algorithm is used to fit the model and cluster the data.  The results are somewhat sensitive to initial conditions.  While the program attempts to assign reasonable values, much better results can often be obtained by setting the starting values by hand.  
	
	Some of the parameters do not have closed form maximum likelihood estiamates.  On each iteration these parameters are updated by a 1-step Newton-Ralphson estimate.  This may cause slightly unsteady estimates near convergence.  
}
\value{
	Returns a list with two components.  
	\item{Theta}{is a list of all the model parameters estimated by the EM algorithm.}
	\item{Data}{is a list containing the original data and observation times.}
}
\references{
	Ning Li, et al.  Functional clustering of periodic transcriptional profiles through ARMA(p,q)
}
\author{
Timothy McMurry and Arthur Berg
}
\note{
%%  ~~further notes~~
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

