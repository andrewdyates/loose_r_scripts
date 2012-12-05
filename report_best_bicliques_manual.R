library("RColorBrewer")
# TODO: load PCC, Spearman
library(illuminaHumanv2.db)
illuminaHumanv2()
library(IlluminaHumanMethylation27k.db)
IlluminaHumanMethylation27k()
library(lumi)
library(annotate)
library(lumiHumanAll.db)
library(GOstats)
library("gplots")
setwd("/nfs/01/osu6683")


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
source("density_merge_bicliques/revised_dec3/BCBig_t.65_o.80_f.16_areamerge.R")
BCBig<-DCOR[bcbig.0.row,bcbig.0.col]

all.row<-unique(c(bcbest.0.row,bcbest.1.row,bcbest.2.row,bcbest.3.row))
all.col<-unique(c(bcbest.0.col,bcbest.1.col,bcbest.2.col,bcbest.3.col))
BC.all <- DCOR[all.row, all.col]

overlaps <- function(A, B) {
  r<-sum(!is.na(match(rownames(A),rownames(A))))
  c<-sum(!is.na(match(colnames(B),colnames(B))))
  r*c
  # how to best represent overlaps?
  }

# TODO: Compare nearest mRNA to DNA methyl
# TODO: Bicliques Overlaps report

report <- function(BC, i=0) {

  # from 0.08 to 0.80 in 0.01 steps
  heatmap_breaks <- seq(0.08,0.8,0.01)
  heatmap_cols <- rev(colorRampPalette(brewer.pal(8,"RdYlBu"))(length(heatmap_breaks)-1))
  sink(paste0("bc",i,"_report.txt"))
  print(paste0("Biclique ", i))
  print(dim(BC))
  print(mean(unlist(BC)))

  title<-paste0("dCOR Hist of all items, BC", i)
  png(paste0("hist_bc",i,"_all.png"), width=10*72, height=8*72)
  h.all <- hist(unlist(BC), main=title, xlab="dCOR"); dev.off(); print(title); print(h.all)
  
  title<-paste0("dCOR Hist of column means, BC", i)
  png(paste0("hist_bc",i,"_col.png"), width=10*72, height=8*72)
  h.col <- hist(colMeans(BC), main=title, xlab="mean dCOR"); dev.off(); print(title); print(h.col)
  
  title<-paste0("dCOR Hist of row items, BC", i)
  png(paste0("hist_bc",i,"_row.png"), width=10*72, height=8*72)
  h.row <- hist(rowMeans(BC), main=title, xlab="mean dCOR"); dev.off(); print(title); print(h.row)

  write.table(BC, file=paste0("bc",i,".dcor.tab"), quote=FALSE)
  png(paste0("heatmap_bc",i,".png"), width=dim(BC)[1]*2+1000, height=dim(BC)[1]*2+1000)
  my.hclust<-function(x, method="average", ...)
    hclust(x, method=method, ...)
  my.dist<- function(x, ...)
    as.dist(1-cor(t(x), method="pearson")) + dist(x)
  heatmap.2(as.matrix(BC), main=paste0("BC",i),
    col=heatmap_cols, ylab="CpG", xlab="mRNA", symm=FALSE, breaks=heatmap_breaks,
    key=TRUE, symkey=FALSE, trace="none",
    distfun=my.dist,
    hclustfun=my.hclust);
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
report(BCBig,"big")
report(BC.all,"all")


# ENRICHMENT
enrich <- function(probeIDs, type="mRNA", pv=0.01, fdr=0.2, ont="CC") {
  # ont in ('CC', 'BP', 'MF')
  if (type == "mRNA") {
    sigLL <- unlist(mget(probeIDs, illuminaHumanv2ENTREZID, ifnotfound = NA))
  } else { # DNA methylation
    sigLL <- unlist(mget(probeIDs, IlluminaHumanMethylation27kENTREZID, ifnotfound = NA))
  }
  sigLL <- unique(as.character(sigLL[!is.na(sigLL)]))
  params <- new("GOHyperGParams",
    geneIds= sigLL,
    annotation="lumiHumanAll.db",
    ontology=ont,
    pvalueCutoff=pv,
    conditional=FALSE,
    testDirection="over")
  
  hgOver <- hyperGTest(params)
  print(hgOver)
  ## Get the p-values of the test
  gGhyp.pv <- pvalues(hgOver)
  ## Adjust p-values for multiple test (FDR)
  gGhyp.fdr <- p.adjust(gGhyp.pv, 'fdr')
  sigGO.ID <- names(gGhyp.fdr[gGhyp.fdr < fdr])
  sigGO.Term <- getGOTerm(sigGO.ID)[[ont]]

  # save as table of:
  # sigGO.ID, sigGO.Term, pv, fdr
}