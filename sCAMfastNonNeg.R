sCAMfastNonNeg <- function(X, A, Sest, iter=1, r = 1, lambda = 0, eta = 100, 
                           silent = FALSE, iteradmm=1000, eps = 1e-3, warm.start = NULL){
  require('nnls')
  require('corpcor')
  K <- nrow(Sest)
  L <- ncol(Sest)
  M <- nrow(X)
  c <- 2*(lambda + r)
  tau <- eta / r
  
  Si <- array(Sest, c(K, L, M))
  S0 <- matrix(c(t(Sest)), nrow=K * L, ncol=M)
  C3 <- rbind(S0, matrix(0, L * K, M))
  
  nuNorm <- c()
  epPriAll <- c()
  epDualAll <- c()
  iter0 <- 0
  while(iter0 < iter){
    j <- 0
    
    if (!is.null(warm.start)) {
        W <- rbind(warm.start, S0)
    } else {
        W <- rbind(matrix(0, L * K, M), S0)
    } 
    
    Z <- matrix(0, 2 * L * K, M)
    
    V <- matrix(0, L * K, M)
    for(i in 1:M){
      V[,i] = c(t(A[i,] * X[rep(i,K),] + 2 * lambda * Sest))
    }
    
    epPri <- Inf
    epDual <- Inf
    while((epPri > eps || epDual > eps) && j < iteradmm){
      H <- W + C3 + Z
      H <- H[1:(L*K),] + H[1:(L*K) + L*K,]
      
      U <- matrix(0, L * K, M)
      for(i in 1:M){
        U[,i] <- c(matrix(V[,i] + r * H[,i], ncol = K) %*% solve(A[i,]%o%A[i,] + c * diag(K)))
      }
      
      
      W1 <- U - S0 - Z[1:(L*K),]
      W2 <- U - Z[1:(L*K) + L*K,]
      
      nuNorm0 <- c()
      for(i in 1:K) {
        rsvd <- fast.svd(W1[1:L + (i-1)*L,])
        rd <- rsvd$d - tau
        rd[rd < 0] <- 0
        W1[1:L + (i-1)*L,] <- rsvd$u %*% diag(rd) %*% t(rsvd$v)
        nuNorm0 <- c(nuNorm0, sum(rsvd$d))
      }
      nuNorm <- rbind(nuNorm, nuNorm0)
      
      W2[W2 < 0] <- 0
      
      W0 <- W
      W <- rbind(W1, W2)
      
      Resi <- rbind(U, U) - W - C3
      Z <- Z - Resi
      
      epPri <- sqrt(sum(Resi^2))
      epDual <- sqrt(sum((W[1:(L*K),] - W0[1:(L*K),] + W[1:(L*K) + L*K,] - W0[1:(L*K) + L*K,])^2))
      epPriAll <- c(epPriAll, epPri)
      epDualAll <- c(epDualAll, epDual)
      if(!silent) cat("iter: ",j, ', pri: ',epPri, ', dual: ',epDual,"\n")
      j <- j + 1
    }
    iterrun <- j
    cat('ADMM pri: ',epPri,'dual: ',epDual, ", iter: ",j,"\n")
    
    for(i in 1:M){
      Si[,,i] <- matrix(W[1:(L*K) + L*K,i], nrow = K, byrow = TRUE)
    }
    
    
    for(i in 1:M){
      Phii2 <- Si[,,i] 
      Phii2[1,] <- Phii2[1,]/sqrt(sum(Phii2[1,]^2))
      for(j in 2:K){
        for(jj in 1:(j-1)){
          Phii2[j,] <- Phii2[j,] - Phii2[jj,] * (Phii2[j,] %*% Phii2[jj,]) 
        }
        Phii2[j,] <- Phii2[j,]/sqrt(sum(Phii2[j,]^2))
      }
      xp <- Phii2 %*% X[i,]
      Sp <- Phii2 %*% t(Si[,,i])
      A[i,] <- coef(nnls(Sp, xp))
    }
    
    iter0 <- iter0 + 1
  }
  
  return(list(A=A, S=Si, W=W, Resi=Resi, epPri=epPriAll, epDual=epDualAll, iterrun = iterrun, nuNorm=nuNorm))
}

