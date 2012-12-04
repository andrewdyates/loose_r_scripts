# R CMD BATCH $HOME/qubic.R qubic.Rout
#

library(rqubic)
load("DCOR.int.RData")
typeof(DCOR.int)
dim(DCOR.int)
colnames(DCOR.int)[1]
rownames(DCOR.int)[1]
hist(DCOR.int, breaks=c(0,1,2,3,4,5,6,7,8))
seed <- generateSeeds(DCOR.int)
save(seed, file="DCOR.int.seed.RData")
DCOR.ES <- new("ExpressionSet", exprs=DCOR.int)
bic <- quBicluster(seed, DCOR.ES)
save(bic, file="DCOR.int.bic.RData")
save.image()
