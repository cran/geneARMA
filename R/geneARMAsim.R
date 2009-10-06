geneARMAsim <-
function(n = 100, J = 2, K = 2, sigsq = .01,  ars=numeric(0), mas=numeric(0), length.out = 21, fourier.coefs, tau)
{
	if(missing(fourier.coefs)){
		if(J==2 && K==2) { 
			fourier.coefs <- matrix(c(0,3,0,0,-2,1,0,1,0,.5), nrow=2, byrow=TRUE)
			tau <- c(.27, .4)
		}
		else stop("Fourier coefficients unspecified with no defaults")
	}
	if(missing(tau)) stop("Frequencies unspecified with no defaults")

	tm <- seq(0, 1, length.out=length.out)
	
	sim.errors <- function(ars, mas, sigsq) {
		if(length(ars) == 0 && length(mas) == 0) rnorm(length.out, sd=sqrt(sigsq))
		else arima.sim(model= list(ar = ars, ma = mas), length.out, sd=sqrt(sigsq))
	}
	if(!(J == nrow(fourier.coefs) &&
		J == length(tau) &&
		K == (ncol(fourier.coefs) - 1)/2 )) stop("Error: Inputs have incompatible dimensions")

	mu2 <- function(j, fourier.coefs){
		coefs <- fourier.coefs[j,]
		onetime <- function(s) {
			coefs[1] + sum( sapply(1:K, 
				function(k) (coefs[2*k] * cos(2*pi*k*s/tau[j]) + 
					coefs[2*k+1] * sin(2*pi*k*s/tau[j])) ) )
		}
		sapply(tm, onetime)
	}

	
	Y <- matrix(,nrow=n, ncol=length.out)
	
	MU <- sapply(1:J, mu2, fourier.coefs = fourier.coefs)
	for(i in 1L:n) {
		Y[i,] <- MU[,i %% J + 1] + sim.errors(ars, mas, sigsq)
	}
	
	return(list(Y=Y, tm=tm))
}

