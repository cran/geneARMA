geneARMAfit <-
function(Y, times, J = 4, K = 2, p = 2, q = 0, arma.skip = 10, eps.conv = .01, max.iter = Inf, print.updates = TRUE, omega.init, fourier.init, tau.init, sigsq.init)
{
	Coef <- function(rs)
	{
		p <- length(rs)
		if(p==0) return(numeric(0))
		
		one.coef <- function(k, zr)		
		{
			if(k == 0) return(1)
			tmp <- combn(1:p, k)
			(-1)^k*sum(apply(tmp, 2, function(ind) prod(zr[ind])))
		}		
		return(sapply(1:p, one.coef, zr=rs))
	}
	
	buildSig <- function(ars, mas, sigsq, dim, inverse=TRUE)
	{
		acfun <- rARMAacf(ar = ars, ma = mas, lag.max = dim)
		R <- toeplitz(acfun)
		Sig <- sigsq * R
		if(inverse)
			return(solve(Sig))
		else return(Sig)
	}
	
	
	checkCausal <- function(ars)
	{
		p <- length(ars)
		if(p == 0) return(ars)
		else{
			rts <- polyroot(c(1,-ars))
			if(all(Mod(rts) > 1)) return(ars)
			else {
				rtsinv <- 1/rts
				toobig <- which(Mod(rtsinv) > 1)
				rtsinv[toobig] <- .5*rtsinv[toobig]/Mod(rtsinv[toobig])
				return(-Re(Coef(rtsinv)))
			}
		}
	}
		
	checkInvertible <- function(mas)
	{
		q <- length(mas)
		if(q == 0) return(mas)
		else{
			rts <- polyroot(c(1,mas))
			if(all(Mod(rts) > 1)) return(mas)
			else {
				rtsinv <- 1/rts
				toobig <- which(Mod(rtsinv) > 1)
				rtsinv[toobig] <- .5*rtsinv[toobig]/Mod(rtsinv[toobig])
				return(Re(Coef(rtsinv)))
			}
		}
	}
	
	
	update.omega <- function(Theta) apply(Theta$P, 2, mean)
	
	update.sigsq <- function(Theta, Data, MP)
	{
		J <- MP$J
		ngene <- MP$ngene
	
		P <- Theta$P
		Y <- Data$Y
	
		acfun <- rARMAacf(ar = Theta$Var$ars, ma = Theta$Var$mas, lag.max = ncol(Data$Y)-1 )
		R <- toeplitz(acfun)
		RInv <- solve(R)	
		
		MU <- sapply(1:J, mu, Theta = Theta, Data = Data, MP = MP)
	
		oneterm <- function(i,j) {
			f <- function(i,j) P[i,j]*(Y[i,] - MU[,j]) %*% RInv %*% (Y[i,] - MU[,j])
			mapply(f, i, j)
		}
		
		sum(outer(1:ngene, 1:J, oneterm))/(ngene*ncol(Data$Y))
	}
	
	mu <- function(j, Theta, Data, MP)
	{
		K <- MP$K
		coefs <- Theta$Mu$cmat[j,]
		onetime <- function(s) 
		{
			coefs[1] + sum( sapply(1:K, 
				function(k) (coefs[2*k] * cos(2*pi*k*s/Theta$Mu$tau[j]) + 
					coefs[2*k+1] * sin(2*pi*k*s/Theta$Mu$tau[j])) ) )
		}
		
		sapply(Data$tm, onetime)
	}
	
		
	update.cj <- function(Theta, Data, MP)
	{
		if(max(abs(diff(Data$tm)[-1] - 
			diff(Data$tm)[-(length(Data$tm)-1)])) < 1e-13)
		{
			if( any(abs(Theta$Mu$tau/diff(Data$tm)[1] - 
				round(Theta$Mu$tau/diff(Data$tm)[1])) < 1e-13)  )
				warning("Tau is a multiple of tm[i+1]-tm[i].  If matrix inverses cause a crash, this may be the reason.")
		}
		
		J <- MP$J
		K <- MP$K
		ngene <- MP$ngene
	
		tau <- Theta$Mu$tau
		tm <- Data$tm
		P <- Theta$P
		Y <- Data$Y
	
		Dmat <- function(j) {	
			scpair <- function(k)
			{
				a <- 2*pi*tm*k/tau[j]
				cbind(cos(a), sin(a))
			}
			tmp <- matrix(,nrow=length(tm), ncol=2*K+1)
			tmp[,1] <- 1
			for(k in 1L:K) tmp[,(2*k):(2*k+1)] <- scpair(k)
			tmp
		}
		
		acfun <- rARMAacf(ar = Theta$Var$ars, ma = Theta$Var$mas, lag.max = ncol(Data$Y)-1)
		R <- toeplitz(acfun)
		Sig <- Theta$Var$sigsq * R
		SigInv <- solve(Sig)
	
		chat <- function(j)
		{
			D <- Dmat(j)
			oneterm.left <- function(i) P[i,j] * t(D) %*% SigInv %*% D
			
			left.bigarray <- array(,dim=c(2*K+1, 2*K+1, ngene))
			for(i in 1L:ngene) left.bigarray[,,i] <- oneterm.left(i)
			DSinvD <- apply(left.bigarray, c(1,2), sum)
	
			oneterm.right <- function(i) P[i,j] * Y[i,] %*% SigInv %*% D
			tmp <- sapply(1:ngene,oneterm.right)	
			right <- apply(tmp, 1, sum)
	
			solve(DSinvD, right)
		}
			
		t(sapply(1:J, chat))
	}
		
	update.tau <- function(Theta, Data, MP)
	{
		J <- MP$J
		K <- MP$K
		ngene <- MP$ngene
	
		tau <- Theta$Mu$tau
		P <- Theta$P
		Y <- Data$Y
		tm <- Data$tm
		
		acfun <- rARMAacf(ar = Theta$Var$ars, ma = Theta$Var$mas, lag.max = ncol(Data$Y)-1)
		R <- toeplitz(acfun)
		Sig <- Theta$Var$sigsq * R
		SigInv <- solve(Sig)
		
		MU <- sapply(1:J, mu, Theta = Theta, Data = Data, MP = MP)
		
		mu.d2 <- function(j, Theta, Data)
		{
			coefs <- Theta$Mu$cmat[j,]
			tau <- Theta$Mu$tau[j]
		
			tm <- Data$tm
		
			oneterm <- function(l, m) {
				f <- function(k,s) {
					-coefs[2*k]*(cos(2*pi*k*s/tau)*4*pi^2*s^2/tau^4 + 
					sin(2*pi*k*s/tau)*4*pi*k*s/tau^3 ) + 
					coefs[2*k+1]*(cos(2*pi*k*s/tau)*4*pi*k*s/tau^3 -
					sin(2*pi*k*s/tau)*4*pi^2*k^2*s^2/tau^4 )
				}
				mapply(f, l, m)
			}
			
		apply(outer(1:K, tm, oneterm), 2, sum)
		}
		
		muD2.mat <- sapply(1:J, mu.d2, Theta, Data)
		
		delta.j <- function(j)
		{
			fc <- Theta$Mu$cmat[j,]
			tauj <- tau[j]
			
			delta.ijl <- function(l)   	#NOTE: doesn't actually depend on i
			{ 
				oneterm <- function(k) {
					a <- 2*pi*k*tm[l]/tauj
					fc[2*k]*sin(a)*a/tauj - fc[2*k+1]*cos(a)*a/tauj
				}
				sum(sapply(1:K, oneterm))
			}
	
			sapply(1:ncol(Y), delta.ijl)
		}
		
		delta <- sapply(1:J, delta.j)
		
		d1 <- function(j) 
		{
			oneterm <- function(i) P[i,j] * as.vector(Y[i,] - MU[,j]) %*% 			SigInv %*% delta[,j]
			sum(sapply(1:ngene, oneterm))
		}
		
		
		d2 <- function(j) 
		{
			oneterm <- function(i) -P[i,j]*delta[,j] %*% SigInv %*% delta[,j] + 
				P[i,j]*(Y[i,] - MU[,j]) %*% SigInv %*% muD2.mat[,j]
			sum(sapply(1:ngene, oneterm))
		}
		
		
		onetau <- function(j) {
			fd <- d1(j)
			sd <- d2(j)
	#		cat("(Analytic)  First:",fd, "  Second:", sd, "\n")
			tau[j] - fd/sd
		}
		
		sapply(1:J, onetau)
	}
		
	update.arma <- function(Theta, Data, MP) 
	{
		J <- MP$J
		K <- MP$K
		ngene <- MP$ngene
	
		P <- Theta$P
		omega <- Theta$omega
		Y <- Data$Y
		ars <- Theta$Var$ars
		mas <- Theta$Var$mas
	
		hn <- 1/max(ngene, 1e6)
	
		MU <- sapply(1:J, mu, Theta = Theta, Data = Data, MP = MP)
	
		onecoef <- function(k, which = c("ar", "ma"))
		{
			if(which == "ar")
			{
				arsP <- arsM <- ars;
				masP <- masM <- mas;
				arsP[k] <- arsP[k] + hn
				arsM[k] <- arsM[k] - hn
			} else if(which == "ma") {
				arsP <- arsM <- ars;
				masP <- masM <- mas;
				masP[k] <- masP[k] + hn
				masM[k] <- masM[k] - hn
			} else stop("Error: Must specify either 'ar' or 'ma'.")
	
			detSig <- det(buildSig(ars, mas, Theta$Var$sigsq, 
				ncol(Data$Y)-1, inverse=FALSE))
			detSigP <- det(buildSig(arsP, masP, Theta$Var$sigsq, 
				ncol(Data$Y)-1, inverse=FALSE))
			detSigM <- det(buildSig(arsM, masM, Theta$Var$sigsq, 
				ncol(Data$Y)-1, inverse=FALSE))
			
			D1detSig <- (detSigP - detSig)/hn
			D2detSig <- (detSigP - 2*detSig + detSigM)/hn^2
	
			SigInv <- buildSig(ars, mas, 
				Theta$Var$sigsq, ncol(Data$Y)-1, inverse=TRUE)
			SigInvP <- buildSig(arsP, masP, 
				Theta$Var$sigsq, ncol(Data$Y)-1, inverse=TRUE)
			SigInvM <- buildSig(arsM, masM, 
				Theta$Var$sigsq, ncol(Data$Y)-1, inverse=TRUE)
			
			
			D1SigInv <- (SigInvP - SigInv)/hn
			D2SigInv <- (SigInvP - 2*SigInv + SigInvM)/hn^2
			
			d1ot <- function(i, j) {
				f <- function(i,j) (-1/2) * P[i,j] * (D1detSig/detSig +
					(Y[i,]-MU[,j]) %*% D1SigInv %*% (Y[i,]-MU[,j]))
				mapply(f, i, j)
			}
	
			d2ot <- function(i, j) {
				f <- function(i,j) (-1/2) * P[i,j] * ( -D1detSig^2/detSig^2 + 
					D2detSig/detSig + (Y[i,]-MU[,j]) %*% D2SigInv %*% (Y[i,]-MU[,j]))
				mapply(f, i, j)
			}
			
			d1tot <- sum(outer(1:ngene, 1:J, d1ot))
			d2tot <- sum(outer(1:ngene, 1:J, d2ot))

			if(which == "ar")
				return(ars[k] - d1tot/d2tot)
			else
				return(mas[k] - d1tot/d2tot)
		}
		
		arsNEW <- sapply(seq_along(ars), onecoef, which="ar")
		arsNEW <- checkCausal(arsNEW)
		masNEW <- sapply(seq_along(mas), onecoef, which="ma")
		masNEW <- checkInvertible(masNEW)
			
		return(list(ars = arsNEW, mas = masNEW))
	}
	
			
	BiglogL <- function(Theta, Data, MP) 
	{
		P <- Theta$P	
		omega <- Theta$omega
		J <- MP$J
		K <- MP$K
		Y <- Data$Y
		ngene <- MP$ngene
		
		MU <- sapply(1:J, mu, Theta = Theta, Data = Data, MP = MP)
	
		acfun <- rARMAacf(ar = Theta$Var$ars, ma = Theta$Var$mas, lag.max = ncol(Data$Y)-1)
		R <- toeplitz(acfun)
		Sig <- Theta$Var$sigsq * R
		SigInv <- solve(Sig)
	
		oneterm <- function(i, j) {
			f <- function(i,j) 
				P[i,j]*(log(omega[j]) - ncol(Y)*log(2*pi)/2 - log(det(Sig))/2
					- (Y[i,] - MU[,j]) %*% SigInv %*% (Y[i,] - MU[,j])/2)
			mapply(f, i, j)
			}
		tmp <- outer(1:ngene, 1:J, oneterm)
		tmp[!is.finite(tmp)] <- min(tmp[is.finite(tmp)])
		sum(tmp)
	}

	
	update.p <- function(Theta, Data, MP, Sig)
	{
		J <- MP$J
		ngene <- MP$ngene
		POld <- Theta$P
		ars <- Theta$Var$ars
		mas <- Theta$Var$mas
		sigsq <- Theta$Var$sigsq
		
		MU <- sapply(1:J, mu, Theta = Theta, Data = Data, MP = MP)
		omega <- Theta$omega
	
		Sig <- buildSig(ars, mas, sigsq, ncol(Data$Y)-1, inverse=FALSE)
	
		oneterm <- function(i,j)
		{
			f <- function(i,j) {
				y <- Data$Y[i,]
				omega[j] * dmvnorm(y, MU[,j], Sig)
			}
			mapply(f, i, j)
		}
		
		PRaw <- outer(1:ngene, 1:J, oneterm)
		denom <- apply(PRaw, 1, sum)
		PNew <- PRaw/denom
		PNew[is.nan(PNew)] <- POld[is.nan(PNew)]
		PNew
	}
	
	
	continue.em <- function(Theta.old, Theta, epsConv)
	{
		delta <- max(abs(
			unlist(Theta.old[-match("ll", names(Theta.old))]) - 
			unlist(Theta[-match("ll", names(Theta))])))
		cat("delta =", delta, "   epsConv =", epsConv, "\n")
		
	
		if ( delta < epsConv ) return(FALSE)
		else return(TRUE)
	}
	
	
	rARMAacf <- function(ar = numeric(0), ma = numeric(0), lag.max, pacf = FALSE)
	{
		if(length(ar) == 0 && length(ma) == 0)
		{
			if(missing(lag.max)) "Stop: LagMax must be specified."
			else return(c(1, rep(0, lag.max)))
		} else ARMAacf(ar, ma, lag.max, pacf=FALSE)
	}
	
	
	runEM <- function(Theta, Data, MP, CP)
	{	
		max.iter <- CP$max.iter
		eps.conv <- CP$eps.conv
		
		CONT <- TRUE
		iter <- 0
		
		while(CONT)
		{
			iter <- iter + 1
	
			Theta.old <- Theta
			Theta <- ee(Theta, Data, MP)
			Theta <- em(Theta, Data, MP, CP, iter)
	
			Theta$ll <- BiglogL(Theta, Data, MP)
	
			if(CP$print.updates) {
				cat("Theta:\n")
				print(Theta$omega)
				print(Theta$Mu)
				print(Theta$Var)
			}
	
			CONT <- continue.em(Theta.old, Theta, eps.conv)
			if(iter >= max.iter) CONT <- FALSE
	
			if(CP$print.updates) {
				cat("Log Likelihood:", Theta$ll, "\n")
				cat("iteration = ", iter, "\n\n\n")
			}
		}
		
		tmp <- list(Theta = Theta, Data = Data)
		class(tmp) <- "geneARMA"
		return(tmp)
	}
			
			
	ee <- function(Theta, Data, MP)
	{
		Theta$P <- update.p(Theta, Data, MP)
		return(Theta)
	}
	
	
	em <- function(Theta, Data, MP, CP, iter) 
	{
		Theta$omega <- update.omega(Theta)	
		Theta$Mu$cmat <- update.cj(Theta, Data, MP)  
		Theta$Var$sigsq <- update.sigsq(Theta,Data, MP)
		Theta$Mu$tau <- update.tau(Theta, Data, MP)
	
		if( iter > CP$arma.skip )
		{
			new.arma <- update.arma(Theta, Data, MP)
				Theta$Var$ars <- new.arma$ars
				Theta$Var$mas <- new.arma$mas
		}
		return(Theta)
	}


	ngene <- nrow(Y)
	
	if(missing(tau.init)) {
		mxt <- max(times)
		mnt <- min(times)
		tau.init <- seq((mxt - mnt)/4, 3*(mxt-mnt)/4, length.out=J)
	}
	if(missing(omega.init)) omega.init <- rep(1/J, J)
	if(missing(fourier.init)) {
		fourier.init <- rnorm(J*(2*K+1))
		dim(fourier.init) <- c(J, 2*K+1)
	}
	if(missing(sigsq.init)) sigsq.init <- max(Y) - min(Y)
	
	P.init <- runif(ngene*J)
	dim(P.init) <- c(ngene, J)
	for(i in 1:ngene) P.init[i,] <- P.init[i,]/sum(P.init[i,])

	ars <- numeric(p)
	mas <- numeric(q)
	
	Data <- list(Y = Y, tm = times)
	Theta <- list(
		omega = omega.init, 
		P = P.init, 
		Mu = list(	cmat = fourier.init, 
					tau = tau.init),
		Var = list(	ars = ars, 
					mas = mas, 
					sigsq = sigsq.init),
		ll = -Inf)

	MP <- list(J=J, K=K, ngene=nrow(Data$Y), p=p, q=q)
	CP <- list(arma.skip = arma.skip,
			eps.conv = eps.conv, 
			max.iter=max.iter, 
			print.updates=print.updates)
	
	return(runEM(Theta, Data, MP, CP))

}

