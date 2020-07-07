# swCAM

#input: 
#Xn: observation mixture matrix (sample * gene)
#Aest: estimated proporiton matrix from CAM (sample * subtype)
#Sest: estimated subtype expression matrix from CAM (subtype * gene)
#eta: eta parameter
#iteradmm: max iteration number in ADMM; larger value will cause longer running time 
Xn <-
Aest <-
Sest <-
eta <- 
iteradmm <- 1000




source('sCAMfastNonNeg.R')


rsCAM <- sCAMfastNonNeg(Xn, Aest, Sest, eta = eta, iteradmm=iteradmm, silent = T)

# eta parameter decides the penalty term of nuclear norm minimization.
# It can be decided by cross-validation scheme.
# From experience, large sample size need larger eta parameter.


Swest <- rsCAM$S #sample-wise expression
sample_names <- rownames(Xn)
gene_names <- colnames(Xn)
save(Swest, sample_names, gene_names, file = "sample-wise S.RData")