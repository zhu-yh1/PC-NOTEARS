#' Notears implementation with edge constraints
#'
#' This function performs causal discovery with continuos-optimization based method NOTEARS.
#'
#' @param X a numeric data matrix with samples as rows and variables as columns.
#' @param lambda1 lambda for l1 penalty, default 0.01.
#' @param losses "logistic" for binary variables and "l2" for continuous variables.
#' @param no_parent constraint that certain variables do not have parent nodes.
#' @param no_edge a 0,1 matrix where 1 species no_edge.
#' @param max.iter default 100.
#' @param h.tol default 1e-6.
#' @param rho.max default 1e+6.
#' @param w.threshold default 0.01.
#' 
#' @return A list containing adjacency matrix for causal graph and intercept.
#'
#'
#' @export
notearsInterceptMultiLoss = function(X, lambda1=0.01, losses=NULL, no_parents=NULL, no_edge=NULL,
									 max.iter=100, h.tol=1e-6, rho.max=1e+6, w.threshold=0.01, 
									 seed=1) {
	
	#variables are in columns of X  
	if(is.null(losses)){
		losses=character(ncol(X))
		for(i in 1:ncol(X)){
			if(all(X[,i] %in% c(0,1))){
				losses[i]="logistic"
			}
			else{
				losses[i]="l2"
			}
		}
		message("The losses are")
		show(table(losses))
	}
	else{
		stopifnot(all(losses %in% c("l2", "logistic", "poisson")))
	}
	losses=c("l2", losses) #first one is intercept
	
	n = nrow(X)
	d = ncol(X)
	nd=n*d
	d1=d+1
	
	G.loss=matrix(0,d1,d1)
	
	G.h=matrix(0, d1,d1)
	#add intercept 
	X=cbind(rep(1, nrow(X)),X)
	#thie frist row of W is not penalizs
	isPenalized=matrix(1, ncol=d+1, nrow=d+1)
	isPenalized[1,]=0
	isPenalizedV=as.vector(isPenalized)
	isPnealizedFrac=mean(isPenalizedV)
	
	#constraints
	upperLimit=matrix(Inf, ncol=d+1, nrow=d+1)
	upperLimit[,1]=0 #the intercept has no parents
	if(!is.null(no_parents)){
		upperLimit[, no_parents+1]=0
	}
	upperLimit[2:(d+1),2:(d+1)][no_edge==1]=0
	
	upperLimitVec=as.vector(upperLimit)
	upperLimitVec=rep(upperLimitVec,2)
	isZeroFrac=mean(upperLimitVec==0)
	# show(isZeroFrac)
	
	loss.func = function(W) { # ?
		M = X %*% W #predict X
		R = X - M #assuming we mostly have L2 loss
		loss=0
		for(i in 1:ncol(X)){
			loss.type=losses[i]
			if (loss.type=='l2') {
				
				loss =loss+ 0.5 / n * sum(R[,i] ** 2)
			} else if (loss.type=='logistic') {
				loss =loss+ 1.0 / n * sum(log(sum(exp(M[,i])+1)) - X[,i] * M[,i])
			} else if (loss.type=='poisson') {
				Si=exp(M[,i])
				loss=loss+ 1.0 / n * sum(Si - X[,i] * M[,i])
			}
		}
		return(loss)
	}
	
	G.loss.func = function(W) { # ?
		
		G.loss[]=0
		M = X %*% W
		R = X - M
		for(i in 1:ncol(X)){
			loss.type=losses[i]
			if (loss.type=='l2') {
				
				G.loss[,i] = -1.0 / n * t(X) %*% R[,i]
			} else if (loss.type=='logistic') {
				G.loss[,i] = 1.0 / n * t(X) %*% (1.0 / (1 + exp(-1 * M[,i])) - X[,i])
			} else if (loss.type=='poisson') {
				Si = exp(M[,i])
				G.loss[,i] = 1.0 / n * t(X) %*% (Si - X[,i])
			}
		}
		return(G.loss)
	}
	
	
	
	
	h.func = function(W) { # ?
		WnoInt=W[-1,-1]
		M = diag(1, d, d) + WnoInt * WnoInt / d
		E = matrixcalc::matrix.power(M, d-1)
		h.new = sum(t(E) * M) - d
		return(h.new)
	}
	
	G.h.func = function(W) { # ?
		WnoInt=W[-1,-1]
		M = diag(1, d, d) + WnoInt * WnoInt / d
		E = matrixcalc::matrix.power(M, d-1)
		G.hsub = t(E) * WnoInt * 2
		#need to reshape back to d+1
		G.h[-1,-1]=G.hsub
		return(G.h)
	}
	
	adj = function(w) { # ?
		w = as.matrix(w)
		w.pos = w[1:(length(w)/2),]
		w.neg = w[((length(w)/2)+1):(length(w)),]
		dim(w.pos) = c(d+1,d+1)
		dim(w.neg) = c(d+1,d+1)
		W = w.pos - w.neg
		return(W)
	}
	
	fn = function(w) { # ?
		W = adj(w)
		loss = loss.func(W)
		h.new = h.func(W)
		
		obj = loss + 0.5 * rho * h.new * h.new + alpha * h.new + lambda1 * sum(w*isPenalizedV)
		return(obj)
	}
	
	gr = function(w) { # ?
		W = adj(w)
		G.loss = G.loss.func(W)
		h.new = h.func(W)
		G.h = G.h.func(W)
		
		G.smooth = G.loss + (rho * h.new + alpha) * G.h
		
		G.obj = c(array(G.smooth + lambda1*isPenalizedV), array(-1 * G.smooth + lambda1*isPenalizedV))
		
		return(G.obj)
	}
	
	
	nodes = colnames(X)
	w.est = replicate(2*(d+1)*(d+1), 0)
	rho = 1
	alpha = 0
	h = Inf
	
	for (i in 1:max.iter) {
		w.new = NULL
		h.new = NULL
		while (rho < rho.max) {
			w.new = optim(w.est, fn, gr=gr, method='L-BFGS-B', lower=0, upper=upperLimitVec)$par # ?
			
			
			h.new = h.func(adj(w.new))
			if (h.new > (0.25*h)) {
				rho = rho*10
			} else {
				break
			}
		}
		w.est = w.new
		
		h = h.new
		alpha = alpha + (rho*h)
		if (h <= h.tol | rho >= rho.max) {
			break
		}
	}
	W.est = adj(w.est)
	
	
	
	W.est[abs(W.est) < w.threshold] = 0
	colnames(W.est) = rownames(W.est) = nodes
	rownames(W.est) = colnames(W.est) = colnames(X)
	return(list(graph=W.est[-1,-1], intercept=W.est[1,-1]))
}