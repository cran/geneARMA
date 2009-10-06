\name{plot.geneARMA}
\Rdversion{1.1}
\alias{plot.geneARMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Plot method for "geneARMA" objects.  Displays fitted models and gene expression by group and time.  
}
\description{
	Produces two types of plots depending on "type" variable.  The first, "all.means," shows all periodic signals fitted to each of the clusters, with the last cluster, which is typically used for the group of genes with no abnormal expression, colored black.  The second type of plot shows the gene experssion values by time with the fitted mean overlayed.  
}
\usage{ 
\method{plot}{geneARMA}(x, y=NULL, type = c("all.means", "single.cluster"), j, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	A fitted model object of class "geneARMA" produced by geneARMAfit(...).
}
  \item{y}{
	Unused.  Supplied for compatability with plot method.   
}
  \item{type}{
	all.means produces a plot of all cluster means.  single.cluster plots the data for a single cluster with the mean overlayed.  
}
  \item{j}{
	is used to indicate which cluster to plot for single.cluster plots only.
}
  \item{\dots}{
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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

