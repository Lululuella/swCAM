# swCAM

#input: 
#Xn: observation mixture matrix (sample * gene)
#Aest: estimated proportion matrix from CAM (sample * subtype)
#Sest: estimated subtype expression matrix from CAM (subtype * gene)
#lambda: lambda parameter
#iteradmm: max iteration number in ADMM; larger value will cause longer running time 
Xn <-
Aest <-
Sest <-
lambda <- 
iteradmm <- 1000




source('sCAMfastNonNeg.R')


rsCAM <- sCAMfastNonNeg(Xn, Aest, Sest, lambda = lambda, iteradmm=iteradmm, silent = T)

# lambda parameter decides the penalty term of nuclear norm minimization.
# It can be decided by cross-validation scheme.
# From experience, large sample size need larger lambda parameter.


Swest <- rsCAM$S #sample-wise expression
sample_names <- rownames(Xn)
gene_names <- colnames(Xn)
save(Swest, sample_names, gene_names, file = "sample-wise S.RData")