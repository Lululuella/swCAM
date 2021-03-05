# K-fold cross validation to select lambda parameter

#input: 
#Xn: observation mixture matrix (sample * gene)
#Aest: estimated proporiton matrix from CAM (sample * subtype)
#Sest: estimated subtype expression matrix from CAM (subtype * gene)
#lambdaAll: candidate lambda parameter
#iteradmm: max iteration number; larger value will cause longer running time 
Xn <-
Aest <-
Sest <-
lambdaAll <- c(10000, 5000, 2000, 1500, 1000, 800, 500, 400, 300, 200, 150, 100, 80, 50, 40, 30, 20, 10 , 8, 5 , 2, 1, 0.5, 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001)
iteradmm <- 1000
Kfold <- 10 

L <- ncol(Xn)
K <- ncol(A)
M <- nrow(Xn)


source('sCAMfastNonNegNA.R')

#shuffle data
set.seed(111)
nall <- ncol(Xn) * nrow(Xn)
tmpseq <- rep(1:Kfold, nall %/% Kfold)
if (nall %% Kfold > 0) {
    tmpseq <- c(tmpseq, 1:(nall %% Kfold))
}
sample_exp <- matrix(sample(tmpseq), nrow(Xn))


library(doSNOW)
cl <- makeCluster(Kfold, type = "SOCK")
registerDoSNOW(cl)

ptm <- proc.time()
res <- foreach (ifold = 1:Kfold, .combine=rbind) %dopar% {
    warmstart <- matrix(0, L * K, M)
    Xtrain <- Xn
    Xtrain[sample_exp == ifold] <- NA
    errlambda <- c()
    
    for (lambda in lambdaAll) {
        rsCAMdtrain <- sCAMfastNonNegNA(Xtrain, Aest, Sest, iter = 1,  r = 1, lambda = lambda, 
                                        iteradmm=iteradmm, silent = T, eps = 1e-3, warm.start=warmstart)
        
        warmstart <- rsCAMdtrain$W[1:(L*K),]
        
        Xest <- c()
        for(i in 1:M){
            tmps1 <- Sest + t(matrix(rsCAMdtrain$W[1:(K*L),i], nrow = ncol(Xn)))
            tmps2 <- t(matrix(rsCAMdtrain$W[1:(K*L) + K*L,i], nrow = ncol(Xn)))
            tmps <- (tmps1 + tmps2) /2
            Xest <- rbind(Xest, Aest[i,,drop=F]%*%tmps)
        }
        errlambda <- c(errlambda, sum((Xest[sample_exp == ifold] - Xn[sample_exp == ifold])^2),
                    rsCAMdtrain$epPri, rsCAMdtrain$epDual, rsCAMdtrain$iterrun)
    }
    errlambda
}
proc.time() - ptm


#compute RMSE and plot the curve
rmse <- res[,seq(1,ncol(res),4)]
rmse <- sqrt(rmse/ (M  * L /Kfold))

self_fun <- function(x) {sqrt(sum(x^2 *(M  * L /Kfold) ) / (M*L))} 
df.cv <- data.frame(x=factor(rep(lambdaAll, each=nrow(rmse))),y=c(rmse))


library(ggplot2)
library(scales)
pict <- ggplot(subset(df.cv, x %in% lambdaAll), aes(x=x, y=y)) + 
    geom_boxplot() + 
    theme_bw(base_size = 24) +
    #stat_summary(fun.y=mean, geom="line", aes(group=1),  colour="blue") +
    stat_summary(fun.y=self_fun, geom="line", aes(group=1),  colour="blue") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=24))+ 
    xlab("lambda") + ylab("RMSE")
ggsave("cross-validation.png",plot = pict, width = 12, height = 6)
