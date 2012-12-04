library("RColorBrewer")
# TODO: load PCC, Spearman
library(illuminaHumanv2.db)
illuminaHumanv2()
library(IlluminaHumanMethylation27k.db)
IlluminaHumanMethylation27k()

load("DCOR.R")
#mRNA.expr <- read.table("gse15745_nov2012_experiments/gse15745_aligned_matrices_nov2/mRNA_correct_aligned.tab", header=TRUE, sep="\t", row.names=1);
#meth.expr <- read.table("gse15745_nov2012_experiments/gse15745_aligned_matrices_nov2/Methyl_correct_aligned.tab", header=TRUE, sep="\t", row.names=1);
#save(mRNA.expr, file="mRNA.expr.R")
#save(meth.expr, file="meth.expr.R")
load("mRNA.expr.R")
load("meth.expr.R")

source("density_merge_bicliques/bicliques_t.45_o.55_f.16_areamerge.R")
BC0<-DCOR[bcbest.0.row,bcbest.0.col]
BC1<-DCOR[bcbest.1.row,bcbest.1.col]
BC2<-DCOR[bcbest.2.row,bcbest.2.col]
BC3<-DCOR[bcbest.3.row,bcbest.3.col]

all.row<-c(bcbest.0.row,bcbest.1.row,bcbest.2.row,bcbest.3.row)
all.col<-c(bcbest.0.col,bcbest.1.col,bcbest.2.col,bcbest.3.col)
BC.all <- DCOR[all.row, all.col]

# TODO: Compare nearest mRNA to DNA methyl
# TODO: Gene ontology report

report <- function(BC, i=0) {

  heatmap_cols <- rev(colorRampPalette(brewer.pal(8,"RdYlBu"))(50))
  sink(paste0("bc",i,"_report.txt"))
  print(paste0("Biclique ", i))
  print(dim(BC))
  print(mean(unlist(BC)))

  title<-paste0("Hist of all items, BC", i)
  pdf(paste0("hist_bc",i,"_all.pdf"), width=10, height=8)
  h.all <- hist(unlist(BC), main=title, xlab="dCOR"); dev.off(); print(title); print(h.all)
  
  title<-paste0("Hist of column means, BC", i)
  pdf(paste0("hist_bc",i,"_col.pdf"), width=10, height=8)
  h.col <- hist(colMeans(BC), main=title, xlab="mean dCOR"); dev.off(); print(title); print(h.col)
  
  title<-paste0("Hist of row items, BC", i)
  pdf(paste0("hist_bc",i,"_row.pdf"), width=10, height=8)
  h.row <- hist(rowMeans(BC), main=title, xlab="mean dCOR"); dev.off(); print(title); print(h.row)

  write.table(BC, file=paste0("bc",i,".dcor.tab"), quote=FALSE)
  #pdf(paste0("heatmap_bc",i,".pdf"), width=dim(BC)[2]/10, height=dim(BC)[1]/10)
  png(paste0("img_bc",i,".png"), width=dim(BC)[1]*10+600, height=dim(BC)[1]*10+600)
  heatmap(as.matrix(BC), main=paste0("BC",i), col=heatmap_cols, ylab="CpG", xlab="mRNA");
  dev.off()
  png(paste0("img_bc",i,".png"), width=dim(BC)[2]*2+600, height=dim(BC)[1]*2+600)
  image(as.matrix(BC), main=paste0("BC",i), col=heatmap_cols);
  dev.off()

  bic.mrna.genes <- unlist(mget(colnames(BC), illuminaHumanv2SYMBOL, ifnotfound = NA))
  bic.meth.genes <- unlist(mget(rownames(BC), IlluminaHumanMethylation27kSYMBOL, ifnotfound = NA))
  write.table(bic.mrna.genes, file=paste0("bc",i,"_mRNA_syms.tab"), quote=FALSE)
  write.table(bic.meth.genes, file=paste0("bc",i,"_meth_syms.tab"), quote=FALSE)

  sink()
}

report(BC0,0)
report(BC1,1)
report(BC2,2)
report(BC3,3)
report(BC.all,"all")
