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


# LOCAL COPY LOAD
# ---------------
setwd("/Users/z/Dropbox/biostat/dec18_biclique")
BC0.cls <- read.table("bc0.cls.tab", header=TRUE, sep=" ", row.names=1);
BC0 <- read.table("bc0.dcor.tab", header=TRUE, sep=" ", row.names=1);
BCBig.cls <- read.table("bcBig.cls.tab", header=TRUE, sep=" ", row.names=1);
BCBig <- read.table("bcbig.dcor.tab", header=TRUE, sep=" ", row.names=1);
# ---------------

# calc mode
# mlv(as.matrix(BC0.cls))

int.hist <- function(x,ylab="Frequency",...) {
  barplot(table(factor(x,levels=min(x):max(x))),space=0,xaxt="n",ylab=ylab,...);axis(1)
}


# COMPUTE DISTANCE MATRIX
# ====================
# define glyph pairwise distance
r1<-c(0,1,1,1,1,1,1,4)
r2<-c(1,0,1,2,2,2,2,4)
r3<-c(1,1,0,1,2,3,2,4)
r4<-c(1,2,1,0,2,2,2,4)
r5<-c(1,2,2,2,0,1,2,4)
r6<-c(1,2,3,2,1,0,1,4)
r7<-c(1,2,2,2,2,1,0,4)
r8<-c(4,4,4,4,4,4,4,0)
dist.glyph <- matrix(c(r1,r2,r3,r4,r5,r6,r7,r8), 8,8)


# Compute distance between two vectors of glyphs
glyph.f <- function(A,B) {
  # glyphs are enumerated from zero already
  sqrt(sum(dist.glyph[as.numeric(A*8+(B+1))]**2))
}

# Compute distance matrix (this is rather slow, but faster than other options
gen.glyph.dist.m <- function(BC) {
  D <- matrix(0,nrow(BC),nrow(BC))
  for (i in 1:nrow(BC)) {
    for (j in (i+1):nrow(BC)) {
      if (j <= nrow(BC)) {
        D[j,i] <- glyph.f(BC[i,], BC[j,])
      }
    }
  }
  as.dist(D)
}

# Re-enumerate glyph symbols
#   1,2,3 to 0,1,2
#   0 to 3

renumerate <- function(BC) {
  BC[BC==0] <- -1
  BC[BC>=1 & BC<=3] <- BC[BC>=1 & BC<=3] -1
  BC[BC==-1] <- 3
  BC
}

# 

# Get indices from hierarchical clustering
D.DCOR.r <- as.dist(1-cor(t(BC0), method="pearson")) + dist(BC0)
D.DCOR.c <- as.dist(1-cor(BC0, method="pearson")) + dist(t(BC0))
D.cls.r <- dist(renumerate(BC0.cls))
D.cls.c <- dist(t(renumerate(BC0.cls)))
Rowv <- rowMeans(BC0, na.rm = TRUE)
Colv <- colMeans(BC0, na.rm = TRUE)

# how to combine D.DCOR and D.cls? norm between 0 and 1?
Rhclust <- as.dendrogram(hclust(D.DCOR.r*3+D.cls.r, method="average"))
Rhclust <- reorder(Rhclust, Rowv)
Chclust <- as.dendrogram(hclust(D.DCOR.c*3+D.cls.c, method="average"))
Chclust <- reorder(Chclust, Colv)

# To get final row and column row orders:
# order.dendrogram(Rhclust)
# order.dendrogram(Chclust)

# clusters glyphs by dCOR
heatmap.2(as.matrix(renumerate(BC0.cls)),
  col=heatmap_cols, ylab="CpG", xlab="mRNA", symm=FALSE, breaks=0:7-0.5,
  key=TRUE, symkey=FALSE, trace="none",
  Rowv=Rhclust, Colv=Chclust
);

# from 0.08 to 0.80 in 0.01 steps
heatmap_breaks_dcor <- seq(0.08,0.8,0.01)
heatmap_cols_dcor <- rev(colorRampPalette(brewer.pal(8,"RdYlBu"))(length(heatmap_breaks_dcor)-1))
# cluster dCOR by dCOR
heatmap.2(as.matrix(BC0), 
  col=heatmap_cols_dcor, ylab="CpG", xlab="mRNA", symm=FALSE, breaks=heatmap_breaks_dcor,
  key=TRUE, symkey=FALSE, trace="none",
  Rowv=Rhclust, Colv=Chclust
);

# There is some special ordering in the dendrogram here...
  my.hclust<-function(x, method="average", ...)
    hclust(x, method=method, ...)
  my.dist<- function(x, ...)
    as.dist(1-cor(t(x), method="pearson")) + dist(x)
  heatmap.2(as.matrix(BC0), 
    col=heatmap_cols_dcor, ylab="CpG", xlab="mRNA", symm=FALSE, breaks=heatmap_breaks_dcor,
    key=TRUE, symkey=FALSE, trace="none",
    distfun=my.dist,
    hclustfun=my.hclust);


hr<-hclust(dist(BC0.cls.r), method="average")
hr$order


setwd("dec18_biclq_class")
report.cls <- function(BC, i=0) {

  #heatmap_cols <- brewer.pal(7,"Accent")
  # divergent colors
  # heatmap_cols <- c("#000000", "#a00d42", "#d7424c", "#eb6532", "#40a185", "#2688bf", "#5b51a5")
  heatmap_cols <- c("#a00d42", "#d7424c", "#eb6532", "#000000", "#40a185", "#2688bf", "#5b51a5")
  # categorical, equal intensity
  #heatmap_cols <- c("#3f62f6", "#009426", "#a5a300", "#000000", "#cd4f00", "#dc1a00", "#b12ac1")

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
  
  write.table(BC, file=paste0("bc",i,".cls.tab"), quote=FALSE, sep="\t")
  png(paste0("clsmap_bc",i,".png"), width=dim(BC)[1]*2+1000, height=dim(BC)[1]*2+1000)

  my.hclust<-function(x, method="average", ...)
    hclust(x, method=method, ...)
  my.dist<- function(x, ...)
    dist(x)
  heatmap.2(as.matrix(BC),
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


# excerpt from custom distance matrix heatmap
  # Precompute dendrograms and row/column ordering
  BC.D <- gen.glyph.dist.m(BC)
  BC.Dt <- gen.glyph.dist.m(t(BC))
  Rhclust<-hclust(BC.D, method="average")
  Chclust<-hclust(BC.Dt, method="average")

  heatmap.2(as.matrix(BC), #main=paste0("Class BC",i),
    col=heatmap_cols, ylab="CpG", xlab="mRNA", symm=FALSE, breaks=0:7-0.5,
    key=TRUE, symkey=FALSE, trace="none",
    Rowv=as.dendrogram(Rhclust), Colv=as.dendrogram(Chclust)
    );
