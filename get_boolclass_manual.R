#CLASSES_I = {
#  'UNL': 0,
#  'HIH': 1,
#  'PC': 2,
#  'LIL': 3,
#  'HIL': 4,
#  'NC': 5,
#  'LIH': 6,
#}

#install.packages("modeest")
library("RColorBrewer")
library(modeest)
library("gplots")

setwd("/nfs/01/osu6683")
source("density_merge_bicliques/bicliques_t.45_o.55_f.16_areamerge.R")
BC0<-DCOR[bcbest.0.row,bcbest.0.col]
source("density_merge_bicliques/revised_dec3/BCBig_t.65_o.80_f.16_areamerge.R")
BCBig<-DCOR[bcbig.0.row,bcbig.0.col]
#mRNA.expr <- read.table("gse15745_nov2012_experiments/gse15745_aligned_matrices_nov2/mRNA_correct_aligned.tab", header=TRUE, sep="\t", row.names=1);
#meth.expr <- read.table("gse15745_nov2012_experiments/gse15745_aligned_matrices_nov2/Methyl_correct_aligned.tab", header=TRUE, sep="\t", row.names=1);


#Classes
CLASSES <- read.table("Methyl_correct_aligned.tab.mRNA_correct_aligned.tab.stepminer.classes.tab", sep="\t");
rownames(CLASSES)<-rownames(meth.expr)
colnames(CLASSES)<-rownames(mRNA.expr)

BC0.cls<-CLASSES[bcbest.0.row,bcbest.0.col]
BCBig.cls<-CLASSES[bcbig.0.row,bcbig.0.col]

# calc mode
# mlv(as.matrix(BC0.cls))

int.hist <- function(x,ylab="Frequency",...) {
  barplot(table(factor(x,levels=min(x):max(x))),space=0,xaxt="n",ylab=ylab,...);axis(1)
}

r1<-c(0,1,1,1,1,1,1,4)
r2<-c(1,0,1,2,2,2,2,4)
r3<-c(1,1,0,1,2,3,2,4)
r4<-c(1,2,1,0,2,2,2,4)
r5<-c(1,2,2,2,0,1,2,4)
r6<-c(1,2,3,2,1,0,1,4)
r7<-c(1,2,2,2,2,1,0,4)
r8<-c(4,4,4,4,4,4,4,0)
glyph.dist<-matrix(c(r1,r2,r3,r4,r5,r6,r7,r8), 8,8)

setwd("dec18_biclq_class")
report.cls <- function(BC, i=0) {

  heatmap_cols <- brewer.pal(7,"Accent")
  sink(paste0("bc",i,"_cls_report.txt"))
  print(paste0("Biclique CLS ", i))
  print(dim(BC))

  modes.row<-lapply(apply(as.matrix(BC), 1, mlv),function(x){x$M})
  modes.col<-lapply(apply(as.matrix(BC), 2, mlv),function(x){x$M})

  title<-paste0("CLS Hist of all items, BC", i)
  png(paste0("hist_bc_cls",i,"_all.png"), width=10*72, height=8*72)
  h.all <- hist(unlist(BC), main=title, xlab="CLS", breaks=0:7-0.5)
  dev.off(); print(title); print(h.all)

  title<-paste0("CLS Mode Hist of rows (CpG), BC", i)
  png(paste0("hist_bc_cls",i,"_row.png"), width=10*72, height=8*72)
  h.row <- hist(unlist(modes.row), breaks=0:7-0.5)
  dev.off(); print(title); print(h.row)

  title<-paste0("CLS Mode Hist of cols (mRNA), BC", i)
  png(paste0("hist_bc_cls",i,"_col.png"), width=10*72, height=8*72)
  h.col <- hist(unlist(modes.col), breaks=0:7-0.5)
  dev.off(); print(title); print(h.col)
  
  write.table(BC, file=paste0("bc",i,".cls.tab"), quote=FALSE)
  png(paste0("clsmap_bc",i,".png"), width=dim(BC)[1]*2+1000, height=dim(BC)[1]*2+1000)
  my.hclust<-function(x, method="average", ...)
    hclust(x, method=method, ...)
  f.glyph.dist<- function(x, ...)
    glyph.dist[x]
  heatmap.2(as.matrix(BC), main=paste0("Class BC",i),
    col=heatmap_cols, ylab="CpG", xlab="mRNA", symm=FALSE, breaks=0:7-0.5,
    key=TRUE, symkey=FALSE, trace="none",
    distfun=my.dist,
    hclustfun=my.hclust);
  dev.off()
  png(paste0("img_bc_cls",i,".png"), width=dim(BC)[2]*2+600, height=dim(BC)[1]*2+600)
  image(as.matrix(BC), main=paste0("Class BC",i), col=heatmap_cols);
  dev.off()

  # Gene ontology.
  sink()
}


report.cls(BC0.cls,0)
report.cls(BCBig.cls,"Big")
