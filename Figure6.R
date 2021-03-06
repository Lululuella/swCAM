## Please download ".Rdata" and open "simu.Rproj" for fully replicating Figure 6.

rm(list=ls())
library(gtools)
library(nnls)
library(MASS)
library(truncnorm)
library(ggplot2)
library(scales)
library(gplots)


##########################################
# Simulation data
###########################################

mat <- as.matrix(read.table("GSE19380/plier-mm.matrix.txt",header=T,row.names=1))
Sraw <- mat[,c(1:4, 1:4+8, 1:4+12)]
label <- rep(1:3, each=4)
Smean0 <- sapply(1:3, function(x) rowMeans(Sraw[,label == x]))
L <- n <- 300

set.seed(111)
sampleIdx <- sample(1:nrow(Sraw),n)
Smean <- Smean0[sampleIdx,]

K <- 3
M <- m <- 50 #sample size
r<- 4 # 4 modules in each of 3 types  
q<-3
p<- c(0.06,0.06, 0.04,0.04)*n # number of genes involved in 4 modules

SnoiseR <- matrix(rnorm(r*q*m,0,1), r*q, m)
SnoiseR <- SnoiseR - rowMeans(SnoiseR) #center the noise signal
SnoiseL <- matrix(0, K*n, r*q)

ggco <- sample(1:n,q*sum(p)) #co-expressed gene index
A <- rdirichlet(m,c(1,1,1))
start <- 0
#1st cell type
for(i in 1:r){
    coeffi <- runif(p[i],0.15,0.3) 
    sign <- rep(c(1,-1), p[i]/2)
    SnoiseL[ggco[start + 1:p[i]], i] <- sign*coeffi
    start <- start + p[i]
}
#2nd cell type
for(i in 1:r){
    coeffi <- runif(p[i],0.15,0.3) 
    sign <- rep(c(1,-1), p[i]/2)
    SnoiseL[ggco[start + 1:p[i]] + n, i+r] <- sign*coeffi
    start <- start + p[i]
}
#3rd cell type
for(i in 1:r){
    coeffi <- runif(p[i],0.15,0.3) 
    sign <- rep(c(1,-1), p[i]/2)
    SnoiseL[ggco[start + 1:p[i]]  + 2 * n, i+2*r] <- sign*coeffi
    start <- start + p[i]
}

SnoiseVec <-  SnoiseL %*% SnoiseR
Snoise <- t(matrix(c(SnoiseVec),nrow=n))
S.sample <- t(Smean)[rep(1:K,m),] * (1+ Snoise)

X<- c()
for(i in 1:m){
    X <- rbind(X, A[i,,drop=F]%*%S.sample[K*(i-1) + 1:K,])
}
dX <- X - A%*%t(Smean)[1:K,]

Xsd <- runif(m*n, 0.02, 0.05)
Xnoise <- matrix(rnorm(m*n,0,1), m, n)
Xn <- X*(1+Xnoise*Xsd) #end of data simulation

# Assume ture A is known. Otherwize, A needs to be estimated
Sest <- apply(Xn,2, function(x) coef(nnls(A,x)))


################################################################################
# Figure 6a: swCAM with lambda = 5
################################################################################

Sest1 <- Sest
Aest1 <- A

source('sCAMfastNonNeg.R')

rsCAMdtrain <- sCAMfastNonNeg(Xn, Aest1, Sest1, iter = 1,  r = 1, lambda = 5, 
                              iteradmm=10000, silent = T, eps = 1e-10)
breaks = seq(-1, 1, length.out=129)

ggcoReorder <- ggco
ggcoReorder <- c(ggcoReorder, (1:n)[-ggco])

optvar0 <- Smean[,rep(1:K, each = M)]

optvar <- cbind(rsCAMdtrain$W[1:n,],rsCAMdtrain$W[n+1:n,],rsCAMdtrain$W[2*n+1:n,])
heatmap.2((t((optvar/optvar0)[c(ggcoReorder),])), trace="none",keysize = 1,
        Rowv = F,Colv = F,#col=redblue(64),
        add.expr = abline(v=c(0+0.5, L+0.5), h= c(0+0.5, M+0.5, 2*M +0.5, 3*M+0.5),lwd=1,col='black'),
        labRow = FALSE, labCol = FALSE,margins = c(2, 2),
        col=bluered(128), key=F, breaks=breaks, symkey=FALSE, density.info="none", 
        lmat=rbind(4:3,2:1), lhei=c(0.1,4), lwid=c(0.1,4))
    

################################################################################
# Figure 6b: swCAM with lambda = 5, delta = 1
################################################################################

source('sCAMfastNonNeg_NN_L21.R')

rsCAMdtrain <- sCAMfastNonNeg_NN_L21(Xn, Aest1, Sest1, iter = 1,  r = 1, lambda = 5, delta = 1,
                                     iteradmm=10000, silent = T, eps = 1e-10)
optvar0 <- Smean[,rep(1:K, each = M)]

optvar <- cbind(rsCAMdtrain$W[1:n,], rsCAMdtrain$W[n+1:n,], rsCAMdtrain$W[2*n+1:n,])
heatmap.2((t((optvar/optvar0)[c(ggcoReorder),])), trace="none",keysize = 1,
        Rowv = F,Colv = F,#col=redblue(64),
        add.expr = abline(v=c(0+0.5, L+0.5), h= c(0+0.5, M+0.5, 2*M +0.5, 3*M+0.5),lwd=1,col='black'),
        labRow = FALSE, labCol = FALSE,margins = c(2, 2),
        col=bluered(128), key=F, breaks=breaks, symkey=FALSE, density.info="none", 
        lmat=rbind(4:3,2:1), lhei=c(0.1,4), lwid=c(0.1,4))


################################################################################
# Figure 6c: 10-fold cross-validation to determine lambda
################################################################################

Kfold <- 10
set.seed(111)

nall <- ncol(Xn) * nrow(Xn)
tmpseq <- rep(1:Kfold, nall %/% Kfold)
if (nall %% Kfold > 0) {
    tmpseq <- c(tmpseq, 1:(nall %% Kfold))
}
sample_exp <- matrix(sample(tmpseq), nrow(Xn))

source('sCAMfastNonNegNA.R')

Sest1 <- Sest
Aest1 <- A

library(doSNOW)
cl <- makeCluster(Kfold, type = "SOCK")
registerDoSNOW(cl)

lambdaAll <- c(10000, 5000, 2000, 1500, 1000, 800, 500, 400, 300, 200, 150, 100, 80, 50, 40, 30, 20, 10 , 8, 5 , 2, 1, 0.5, 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001)

ptm <- proc.time()
res <- foreach (ifold = 1:Kfold, .combine=rbind) %dopar% {
    warmstart <- matrix(0, L * K, M)
    Xtrain <- Xn
    Xtrain[sample_exp == ifold] <- NA
    errlambda <- c()
    
    for (lambda in lambdaAll) {
        rsCAMdtrain <- sCAMfastNonNegNA(Xtrain, Aest1, Sest1, iter = 1,  r = 1, lambda = lambda, 
                                        iteradmm=10000, silent = T, eps = 1e-10, warm.start=warmstart)
        
        warmstart <- rsCAMdtrain$W[1:(L*K),]
        
        Xest <- c()
        for(i in 1:M){
            tmps1 <- Sest1 + t(matrix(rsCAMdtrain$W[1:(K*L),i], nrow = ncol(Xn)))
            tmps2 <- t(matrix(rsCAMdtrain$W[1:(K*L) + K*L,i], nrow = ncol(Xn)))
            tmps <- (tmps1 + tmps2) /2
            Xest <- rbind(Xest, Aest1[i,,drop=F]%*%tmps)
        }
        errlambda <- c(errlambda, sum((Xest[sample_exp == ifold] - Xn[sample_exp == ifold])^2),
                       rsCAMdtrain$epPri, rsCAMdtrain$epDual, rsCAMdtrain$iterrun)
    }
    errlambda
}
proc.time() - ptm

stopCluster(cl)

rmse <- res[,seq(1,ncol(res),4)]
rmse <- sqrt(rmse/ (M  * L /Kfold))

self_fun <- function(x) {sqrt(sum(x^2 *(M  * L /Kfold) ) / (M*L))}

df.cv <- data.frame(x=factor(rep(lambdaAll, each=nrow(rmse))),y=c(rmse))
ggplot(subset(df.cv, x %in% c(10000, 5000, 2000, 1500, 1000, 800, 500, 400, 300, 
                              200, 150, 100, 80, 50, 20, 10 , 8, 5 , 2, 1, 0.5, 
                              0.2, 0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001)), 
       aes(x=x, y=y)) + geom_boxplot() + theme_bw(base_size = 16) +
       stat_summary(fun.y=self_fun, geom="line", aes(group=1),  colour="blue") +
       theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=10))+ 
       xlab("lambda") + ylab("RMSE")
