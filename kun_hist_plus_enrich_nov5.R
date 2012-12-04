# DCOR <- read.table("Methyl-correct-aligned.pkl_mRNA-correct-aligned.pkl.DCOR.values.tab", header=TRUE, sep="\t", row.names=1); save(DCOR, file="Methyl-correct-aligned.pkl_mRNA-correct-aligned.pkl.DCOR.values.RData")
# PCC <- read.table("Methyl-correct-aligned.pkl_mRNA-correct-aligned.pkl.PEARSON.values.tab", header=TRUE, sep="\t", row.names=1); save(PCC, file="Methyl-correct-aligned.pkl_mRNA-correct-aligned.pkl.PEARSON.values.RData")
# SP <- read.table("Methyl-correct-aligned.pkl_mRNA-correct-aligned.pkl.SPEARMAN.values.tab", header=TRUE, sep="\t", row.names=1); save(SP, file="Methyl-correct-aligned.pkl_mRNA-correct-aligned.pkl.SPEARMAN.values.RData")

library(illuminaHumanv2.db)
illuminaHumanv2()
library(IlluminaHumanMethylation27k.db)
IlluminaHumanMethylation27k()

load("Methyl-correct-aligned.pkl_mRNA-correct-aligned.pkl.DCOR.values.RData")
load("Methyl-correct-aligned.pkl_mRNA-correct-aligned.pkl.PEARSON.values.RData")
load("Methyl-correct-aligned.pkl_mRNA-correct-aligned.pkl.SPEARMAN.values.RData")
# HACK: manually load row and column names of bicliques generated manually
source("bicliques.R")

bic.dcor <- DCOR[meth,mrna]
bic.pcc <- PCC[meth,mrna]
bic.sp <- SP[meth,mrna]
write.table(bic.dcor, file="bic1.dcor.tab", quote=FALSE)
write.table(bic.pcc, file="bic1.pcc.tab", quote=FALSE)
write.table(bic.sp, file="bic1.sp.tab", quote=FALSE)

bic.mrna.genes <- unlist(mget(colnames(bic.dcor), illuminaHumanv2SYMBOL, ifnotfound = NA))
bic.meth.genes <- unlist(mget(rownames(bic.dcor), IlluminaHumanMethylation27kSYMBOL, ifnotfound = NA))

# Histograms
pdf("hist_dcor_bic1.pdf", width=10, height=8)
hist(as.matrix(bic.dcor))
dev.off()
pdf("hist_pcc_bic1.pdf", width=10, height=8)
hist(as.matrix(bic.pcc))
dev.off()
pdf("hist_sp_bic1.pdf", width=10, height=8)
hist(as.matrix(bic.sp))
dev.off()
pdf("hist_dcor_all.pdf", width=10, height=8)
hist(as.matrix(DCOR))
dev.off()
pdf("hist_pcc_all.pdf", width=10, height=8)
hist(as.matrix(PCC))
dev.off()
pdf("hist_sp_all.pdf", width=10, height=8)
hist(as.matrix(SP))
dev.off()

sum(!is.na(match(bic.meth.genes, bic.mrna.genes)))

meth.syms <- unlist(mget(rownames(DCOR),
     IlluminaHumanMethylation27kSYMBOL,
     ifnotfound = NA))
mrna.syms <- unlist(mget(colnames(DCOR),
     illuminaHumanv2SYMBOL,
     ifnotfound = NA))
name.match <- match(meth.syms, mrna.syms)

!is.na(name.match)

make.getD <- function(D, q) {
  function(i) {
    if (is.na(q[i])) {
      return(NA)
    } else {
      return(D[i, q[i]])
    }
  }
}
dcor.meth2gene<-sapply(1:length(name.match), make.getD(DCOR,name.match))
summary(dcor.meth2gene)
pcc.meth2gene<-sapply(1:length(name.match), make.getD(PCC,name.match))
summary(pcc.meth2gene)
sp.meth2gene<-sapply(1:length(name.match), make.getD(SP,name.match))
summary(sp.meth2gene)

pdf("hist_meth2gene_dcor_all.pdf", width=10, height=8)
hist(dcor.meth2gene)
dev.off()
pdf("hist_meth2gene_pcc_all.pdf", width=10, height=8)
hist(pcc.meth2gene)
dev.off()
pdf("hist_meth2gene_sp_all.pdf", width=10, height=8)
hist(sp.meth2gene)
dev.off()

# histograms of meth to nearest gene in biclique 1
q<-dcor.meth2gene[meth]
q[is.na(q)] <- 0
pdf("hist_meth2gene_dcor_bc1.pdf", width=10, height=8)
hist(q)
dev.off()

q<-pcc.meth2gene[meth]
q[is.na(q)] <- 0
pdf("hist_meth2gene_pcc_bc1.pdf", width=10, height=8)
hist(q)
dev.off()

q<-sp.meth2gene[meth]
q[is.na(q)] <- 0
pdf("hist_meth2gene_sp_bc1.pdf", width=10, height=8)
hist(q)
dev.off()

# CpG to random expressed mRNA probe


dcor.meth2randgene<-sapply(1:length(name.match), make.getD(DCOR,sample(name.match)))
summary(dcor.meth2randgene)
dcor.meth2gene<-sapply(1:length(name.match), make.getD(DCOR,name.match))
summary(dcor.meth2gene)

pcc.meth2randgene<-sapply(1:length(name.match), make.getD(PCC,sample(name.match)))
summary(pcc.meth2randgene)
pcc.meth2gene<-sapply(1:length(name.match), make.getD(PCC,name.match))
summary(pcc.meth2gene)


# ==============================
# 1. Get mRNA probe with maximum mean expression per unique gene symbol
# 1.1: load mRNA and methylation values
mRNA.expr <- read.table("gse15745_aligned_matrices/mRNA_correct_aligned.tab", header=TRUE, sep="\t", row.names=1);
save(mRNA.expr, file="mRNA_correct_aligned.RData")
meth.expr <- read.table("gse15745_aligned_matrices/Methyl_correct_aligned.tab", header=TRUE, sep="\t", row.names=1);
save(meth.expr, file="Methyl_correct_aligned.RData")

# 1.2: get unique mRNA names 
mrna.uniq <- unique(mrna.syms)
mrna.uniq <- mrna.uniq[!is.na(mrna.uniq)]

topProbe <- function(syms, expr) {
  u<-rowMeans(expr)
  function(s) {
    if (is.na(s)) { return(NA) }
    q<-syms==s
    q[is.na(q)] <- FALSE
    if (sum(q) == 0) { return(NA) }
    max.probe.name<-names(which.max(u[q]))[1]
    return(max.probe.name)
  }
}
topProbemRNA <- topProbe(mrna.syms, mRNA.expr)
topProbeMapmRNA<-sapply(mrna.syms, topProbemRNA)

# ------------------------------
# 2. map DNA methylation probes to most proximal, best expressed mRNA probes
# random gene test should be from set of all gene symbols covered by DNA methylation platform

# list of all DNA methylation probes to gene symbols
meth2genesym <- unlist(mget(rownames(DCOR), IlluminaHumanMethylation27kSYMBOL, ifnotfound = NA))
# no CpGs with missing genes!
# sum(is.na(meth2genesym))
# [1] 0
# note: maps probe to mRNA gene symbols expressed
topProbeMapMeth <- sapply(meth2genesym, topProbemRNA)

meth2mrna.colidx<-match(topProbeMapMeth, colnames(DCOR))
# TODO: MAKE RANDOM GENE MAP

dcor.meth2gene<-sapply(1:length(meth2mrna.colidx), make.getD(DCOR, meth2mrna.colidx))
summary(dcor.meth2gene)
# this is not quite correct
dcor.meth2rand<-sapply(1:length(meth2mrna.colidx), make.getD(DCOR, sample(meth2mrna.colidx)))
summary(dcor.meth2rand)

pcc.meth2gene<-sapply(1:length(meth2mrna.colidx), make.getD(PCC, meth2mrna.colidx))
summary(pcc.meth2gene)
# this is not quite correct
pcc.meth2rand<-sapply(1:length(meth2mrna.colidx), make.getD(PCC, sample(meth2mrna.colidx)))
summary(pcc.meth2rand)

pcc.meth2gene.0 <- pcc.meth2gene
pcc.meth2gene.0[is.na(pcc.meth2gene)] <- 0
pcc.meth2rand.0 <- pcc.meth2rand
pcc.meth2rand.0[is.na(pcc.meth2rand)] <- 0

pdf("hist_pcc_all_zeroed.pdf")
hist(pcc.meth2gene.0, breaks=200)
hist(pcc.meth2rand.0, breaks=200)
dev.off()
pdf("hist_pcc_all.pdf")
hist(pcc.meth2gene, breaks=200)
hist(pcc.meth2rand, breaks=200)
dev.off()

# Dr. Kun's improved histogram
# 10,415 of 24,334 (42.8%) CpG probes mapped to undetected mRNA probes excluded.
x <- hist(pcc.meth2gene, nclass=100)
y <- hist(pcc.meth2rand, nclass=100)

pdf("hist_all_meth2gene_meth2rand.pdf", width=10, height=10)
xrange=range(-1,1)
yrange=range(0,max(x$count, y$count))
plot(xrange, yrange,
  type="n", xlab="Pearson's Correlation",
  ylab="Count")

lines(x$mid, x$count, col="blue", type='l')
lines(y$mid, y$count, col="red", type='l')

# add a title and subtitle 
title("CpG vs Gene Expression", "10,415 of 24,334 (42.8%) CpG probes mapped to undetected mRNA probes excluded.")

# add a legend
D = 0.0736, p-value legend(xrange[1], yrange[2], c('nearest', 'random'), cex=0.8, col=c('blue','red'), lty=1)
dev.off()

# --------------------
# top 50 most variable cpgs
CpG_Probe_STD<-apply(meth.expr, 1, sd)
pdf("hist_std_cpg.pdf")
hist(CpG_Probe_STD)
dev.off()

#> sum(CpG_Probe_STD > 0.5031)
#[1] 6084
#> length(CpG_Probe_STD)
#[1] 24334
#> 6084/24334
#[1] 0.2500205

# HISTOGRAM OF TOP 25% BY VARIANCE OF CPGS
pcc.meth2gene.top <- pcc.meth2gene[CpG_Probe_STD > 0.5031]
pcc.meth2rand.top <- pcc.meth2rand[CpG_Probe_STD > 0.5031]
pcc.meth2gene.top[is.na(pcc.meth2gene.top)] <- 0
pcc.meth2rand.top[is.na(pcc.meth2rand.top)] <- 0
x <- hist(pcc.meth2gene.top, nclass=30)
y <- hist(pcc.meth2rand.top, nclass=30)

pdf("hist_top25pct_meth2gene_meth2rand_na0.pdf", width=10, height=10)
xrange=range(-1,1)
yrange=range(0,max(x$count, y$count))
plot(xrange, yrange,
  type="n", xlab="Pearson's Correlation",
  ylab="Count")

lines(x$mid, x$count, col="blue", type='l')
lines(y$mid, y$count, col="red", type='l')

# add a title and subtitle 
title("CpG vs Gene Expression", "Top 25% of CpGs (6,084 of 24,334) by variance. NA replaced by 0.")

# add a legend
legend(xrange[1], yrange[2], c('nearest', 'random'), cex=0.8, col=c('blue','red'), lty=1)
dev.off()



# HISTOGRAM OF TOP 50% BY VARIANCE OF CPGS
pcc.meth2gene.top <- pcc.meth2gene[CpG_Probe_STD > median(CpG_Probe_STD)]
pcc.meth2rand.top <- pcc.meth2rand[CpG_Probe_STD > median(CpG_Probe_STD)]
pcc.meth2gene.top[is.na(pcc.meth2gene.top)] <- 0
pcc.meth2rand.top[is.na(pcc.meth2rand.top)] <- 0
x <- hist(pcc.meth2gene.top, nclass=50)
y <- hist(pcc.meth2rand.top, nclass=50)

#pdf("hist_top50pct_meth2gene_meth2rand_nona.pdf", width=10, height=10)
pdf("hist_top50pct_meth2gene_meth2rand_na0.pdf", width=10, height=10)
xrange=range(-1,1)
yrange=range(0,max(x$count, y$count))
plot(xrange, yrange,
  type="n", xlab="Pearson's Correlation",
  ylab="Count")

lines(x$mid, x$count, col="blue", type='l')
lines(y$mid, y$count, col="red", type='l')

# add a title and subtitle 
#title("CpG vs Gene Expression", "Top 50% of CpGs (12,167 of 24,334) by variance. NA removed.")
title("CpG vs Gene Expression", "Top 50% of CpGs (12,167 of 24,334) by variance. NA replaced by 0.")

# add a legend
legend(xrange[1], yrange[2], c('nearest', 'random'), cex=0.8, col=c('blue','red'), lty=1)
dev.off()


# double check correlations: OK
cor(t(meth.expr[2,]), t(mRNA.expr[8349,]))

# ------------------------------
# NOVEMBER 8th SESSION
# ----------------------------------------

meth<-rows
mrna<-cols
bic.dcor <- DCOR[meth,mrna]
bic.pcc <- PCC[meth,mrna]
bic.sp <- SP[meth,mrna]
write.table(bic.dcor, file="bic1.dcor.tab", quote=FALSE)
write.table(bic.pcc, file="bic1.pcc.tab", quote=FALSE)
write.table(bic.sp, file="bic1.sp.tab", quote=FALSE)

# Histograms
pdf("hist_dcor_bic1.pdf", width=10, height=8)
hist(as.matrix(bic.dcor))
dev.off()
pdf("hist_pcc_bic1.pdf", width=10, height=8)
hist(as.matrix(bic.pcc))
dev.off()
pdf("hist_sp_bic1.pdf", width=10, height=8)
hist(as.matrix(bic.sp))
dev.off()
pdf("hist_dcor_all.pdf", width=10, height=8)
hist(as.matrix(DCOR))
dev.off()
pdf("hist_pcc_all.pdf", width=10, height=8)
hist(as.matrix(PCC))
dev.off()
pdf("hist_sp_all.pdf", width=10, height=8)
hist(as.matrix(SP))
dev.off()

# ==============================
# GENE ONTOLOGY
# ==============================
library(lumi)
library(annotate)
library(lumiHumanAll.db)
library(GOstats)
library(illuminaHumanv2.db)
illuminaHumanv2()
library(IlluminaHumanMethylation27k.db)
IlluminaHumanMethylation27k()


# ENTREZID
#> length(sigLL)
#[1] 123 -> 120 unique
#> length(colnames(bic.dcor))
#[1] 127
sigLL <- unlist(mget(colnames(bic.dcor), illuminaHumanv2ENTREZID, ifnotfound = NA))
sigLL <- unique(as.character(sigLL[!is.na(sigLL)]))

# BP
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
params <- new("GOHyperGParams",
  geneIds= sigLL,
  annotation="lumiHumanAll.db",
  ontology="CC",
  pvalueCutoff= 0.01,
  conditional=FALSE,
  testDirection="over")

hgOver <- hyperGTest(params)
## Get the p-values of the test
CC.gGhyp.pv <- pvalues(hgOver)
## Adjust p-values for multiple test (FDR)
CC.gGhyp.fdr <- p.adjust(CC.gGhyp.pv, 'fdr')
CC.sigGO.ID <- names(CC.gGhyp.fdr[gGhyp.fdr < 0.2])
CC.sigGO.Term <- getGOTerm(CC.sigGO.ID)[["CC"]]

# BP
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
params <- new("GOHyperGParams",
  geneIds= sigLL,
  annotation="lumiHumanAll.db",
  ontology="BP",
  pvalueCutoff= 0.01,
  conditional=FALSE,
  testDirection="over")

hgOver <- hyperGTest(params)
## Get the p-values of the test
BP.gGhyp.pv <- pvalues(hgOver)
## Adjust p-values for multiple test (FDR)
BP.gGhyp.fdr <- p.adjust(BP.gGhyp.pv, 'fdr')
BP.sigGO.ID <- names(BP.gGhyp.fdr[gGhyp.fdr < 0.2])
BP.sigGO.Term <- getGOTerm(BP.sigGO.ID)[["BP"]]

# MF
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
params <- new("GOHyperGParams",
  geneIds= sigLL,
  annotation="lumiHumanAll.db",
  ontology="MF",
  pvalueCutoff= 0.01,
  conditional=FALSE,
  testDirection="over")

hgOver <- hyperGTest(params)
## Get the p-values of the test
MF.gGhyp.pv <- pvalues(hgOver)
## Adjust p-values for multiple test (FDR)
MF.gGhyp.fdr <- p.adjust(MF.gGhyp.pv, 'fdr')
MF.sigGO.ID <- names(MF.gGhyp.fdr[gGhyp.fdr < 0.2])
MF.sigGO.Term <- getGOTerm(MF.sigGO.ID)[["MF"]]

geneList <- unlist(mget(colnames(bic.dcor), illuminaHumanv2SYMBOL, ifnotfound = NA))


write.table(CC.sigGO.Term, file="bic1.GO.CC.txt", quote=FALSE)
write.table(BP.sigGO.Term, file="bic1.GO.BP.txt", quote=FALSE)
write.table(MF.sigGO.Term, file="bic1.GO.MF.txt", quote=FALSE)
write.table(geneList, file="bic1.genesyms.txt", quote=FALSE)
