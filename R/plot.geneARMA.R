plot.geneARMA <- function(x, y=NULL, type=c("all.means", "single.cluster"),j,...)
{
	fitted.model <- x
	Theta <- fitted.model$Theta
	Data <- fitted.model$Data
	J <- length(Theta$omega)
	
	mu2 <- function(j, Theta, tm){
		coefs <- Theta$Mu$cmat[j,]
		K <- (length(coefs)-1)/2
		onetime <- function(s) {
			coefs[1] + sum( sapply(1:K, 
				function(k) (coefs[2*k] * cos(2*pi*k*s/Theta$Mu$tau[j]) + 
					coefs[2*k+1] * sin(2*pi*k*s/Theta$Mu$tau[j])) ) )
		}
		sapply(tm, onetime)
	}

	means.plot <- function(Theta, Data)
	{
		MU <- sapply(1:J, mu2, Theta=Theta,tm=Data$tm)
		clrs <- c(rainbow(J-1, alpha=Theta$omega[-J]^.27), "#000000FF")
		plot(Data$tm, MU[,1], type="n", ylim=c(min(MU), max(MU)), 			xlab = "Time", 
			ylab="Normalized Expression", ...)	
		for(i in 1:J) points(Data$tm, MU[,i], type="l", lwd=3,
			col=clrs[i], ...)
	}
	
	cluster.plot <- function(Theta, Data, j)
	{
		ind <- which(Theta$P[,j] > .9)
		cldat <- Data$Y[ind,]
		MU <- mu2(j, Theta, Data$tm)
		plot(Data$tm, cldat[1,], ylim=c(min(cldat), 
			max(cldat)), type="n", 
			xlab = "Time",
			ylab = "Normalized Expression",
			main = paste("Cluster", j), ...)
		for(i in 1L:nrow(cldat)) {
			points(Data$tm, cldat[i,], pch = 19, 
				col="#0000FF30", , ...)
		}
		lines(Data$tm, MU, col="black", ...)
	}
				
	type <- match.arg(type)
	switch(type, 
		all.means=means.plot(Theta, Data),
		single.cluster=cluster.plot(Theta, Data, j)
	)
}

