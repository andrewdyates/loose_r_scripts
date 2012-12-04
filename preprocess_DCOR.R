# R CMD BATCH $HOME/qubic.R qubic.Rout
#

load("DCOR.R")
rmax<-apply(DCOR,1,max)
cmax<-apply(DCOR,2,max)
DCOR.sig<-DCOR[rmax>0.08,cmax>0.08]
print(dim(DCOR.sig))
#[1] 13076  9501
DCOR.sig <- (DCOR.sig-0.08)*10

# Set negative values to zero
DCOR.sig<-(DCOR.sig+abs(DCOR.sig))/2
save(DCOR.sig, file="DCOR.sig.RData")

DCOR.int<-as.matrix(round(DCOR.sig))
storage.mode(DCOR.int) <- "integer"
save(DCOR.int, file="DCOR.int.RData")

pdf("DCOR.int_hist.pdf", width=8, height=6")
h <- hist(DCOR.int, breaks=c(1,2,3,4,5,6,7))
dev.off()
