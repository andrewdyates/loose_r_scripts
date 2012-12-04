library(lumi)
library(IlluminaHumanMethylation27k.db)
IlluminaHumanMethylation27k()

# Use mappings based on hg19
IlluminaHumanMethylation27kCPG37 <- IlluminaHumanMethylation27kCHRLOC
IlluminaHumanMethylation27kCHR37 <- IlluminaHumanMethylation27kCHR

# Visualize data function
vis <- function(M, name, do_mds=TRUE) {
  # colors?
  pdf(paste0("plotSampleRelations_cluster_", name, ".pdf"), width=60, height=30)
  plotSampleRelation(M, method='cluster', cv.Th=0)
  dev.off()
  if (do_mds) {
    pdf(paste0("plotSampleRelations_mds_", name, ".pdf"), width=30, height=30)
    plotSampleRelation(M, method='mds', cv.Th=0)
    dev.off()
  }
  pdf(paste0("density_", name, ".pdf"), width=30, height=30)
  density(M, xlab="M-value")
  dev.off()
  pdf(paste0("plotColorBias1D_", name, ".pdf"), width=20, height=20)
  plotColorBias1D(M, channel='sum')
  dev.off()
  pdf(paste0("boxplotColorBias_", name, ".pdf"), width=80, height=12)
  boxplotColorBias(M, channel='sum')
  dev.off()
  pdf(paste0("plotColorBias2D_", name, ".pdf"), width=20, height=20)
  plotColorBias2D(M, selSample=1)
  dev.off()
  ## Because the distribution of M-value has two modes, we use a boxplot different from regular ones
  ## Figure 6: Multi-mode boxplot of M-value methylation levels before normalization
  pdf(paste0("boxplot_multimode_", name, ".pdf"), width=80, height=10)
  boxplot(M)
  dev.off()
  # intensity (methy + unmethy probe intensity) should have similar distribution over all arrays
  pdf(paste0("density_estIntens_", name, ".pdf"), width=30, height=70)
  density(estimateIntensity(M), xlab="log2(CpG-site Intensity)")
  dev.off()
  pdf(paste0("boxplot_estIntens_", name, ".pdf"), width=80, height=12)
  boxplot(estimateIntensity(M))
  dev.off()
  # TODO: Add heatmap
}


# LOAD FROM TEXT
samps <- read.table("GSE15745_GPL8490.samples.methylumi.tab",
  sep="\t",header=TRUE)

fileName <- "gse15745_gpl8490_methylumi2.tab"
M <- lumiMethyR(fileName,
  lib="IlluminaHumanMethylation27k.db", sep="\t",
  sampleDescriptions=samps)

# LOAD FROM BINARY
#load("M.methylumi.RData")

vis(M, "raw1")

# FILTER ARRAYS BY INTENSITY
intensity <- log2(exprs(estimateIntensity(M)))
pdf("hist_read_estIntensity_log_all.pdf", width=8, height=8)
hist(intensity); dev.off()
pdf("hist_array_estIntensity_log_all.pdf", width=8, height=8)
hist(colMeans(intensity)); dev.off()
x<-sd(colMeans(intensity))
y<-mean(colMeans(intensity))
array_filt <- y-colMeans(intensity) > x
print(sum(array_filt))
print(sort(colMeans(intensity)[array_filt]))
print(summary(colMeans(intensity)))
print(summary(colMeans(intensity)[!array_filt]))
pdf("hist_array_estIntensity_log_filt.pdf", width=8, height=8)
hist(colMeans(intensity)[!array_filt]); dev.off()
M.arrayfilt <- M[,!array_filt]

vis(M.arrayfilt, "arrayfilt2")


# FILTER PROBES
# by intensity
probe_intensity <- rowMeans(log2(exprs(estimateIntensity(M.arrayfilt))))
print(summary(probe_intensity))
pdf("hist_probe_estIntensity_log_all.pdf", width=8, height=8)
hist(probe_intensity); dev.off()
s<-sd(probe_intensity)
u<-mean(probe_intensity)
int_filt <- (u - probe_intensity) > s*2.5
print(sort(probe_intensity[int_filt])[1:30])
print(summary(probe_intensity[!int_filt]))
pdf("hist_probe_estIntensity_log_filt_int.pdf", width=8, height=8)
hist(probe_intensity[!int_filt]); dev.off()

# by missing chromosome
no_chr <- is.na(featureData(M.arrayfilt)$CHROMOSOME)
print(summary(probe_intensity[!no_chr]))

# by XY vs autosome (after filtering by chromosome and low intensity)
isX <- featureData(M.arrayfilt)$CHROMOSOME == "X"
isY <- featureData(M.arrayfilt)$CHROMOSOME == "Y"
isX[is.na(isX)] <- FALSE
isY[is.na(isY)] <- FALSE
print(summary(probe_intensity[!(isX|isY)]))

# combine filters
print(sum(int_filt))
print(sum(no_chr))
print(sum(isX))
print(sum(isY))
probe_filt <- int_filt | no_chr | isX | isY
print(sum(probe_filt))
print(summary(probe_intensity[!probe_filt]))
pdf("hist_probe_estIntensity_log_filt_all.pdf", width=8, height=8)
hist(probe_intensity[!probe_filt]); dev.off()

M.filt <- M.arrayfilt[!probe_filt,]
vis(M.filt, "filt3")


# COLOR BALANCE
# Filter arrays that cause unexplained problems (singular value error?) in color balance adjustment
# unfortunately, for now, this step must be executed manually
good <- !seq(1,dim(M.filt)[2]) %in% c(462)
M.adj <- lumiMethyC(M.filt[,good])
vis(M.adj, "colorAdj4")

# Outlier arrays by M-value distribution
sds <- apply(exprs(M.adj), 2, sd)
print(summary(sds))
pdf("hist_array_mvalue_std_adj.pdf", width=8, height=8)
hist(sds); dev.off()
u <- mean(sds)
s <- sd(sds)
filt_sd <- u-sds > s*2
print(sum(filt_sd))
pdf("hist_array_mvalue_std_adj_filt.pdf", width=8, height=8)
hist(sds[!filt_sd]); dev.off()
M.adj_sd <- M.adj[,!filt_sd]
vis(M.adj_sd, "stdfilt5")

 
## Background adjustment is skipped because the SSN normalization includes background adj
# NORMALIZATION
M.ssn <- lumiMethyN(M.adj_sd, method='ssn')
vis(M.ssn, "ssnNorm6")

# REMOVE BATCH EFFECTS
library(limma)
library(sva)

edata <- exprs(M.ssn)
pheno_select <- c(F,F,F,F,F,T,T,T,F,F,T,T,T,T)
pheno <- pData(M.ssn)[,pheno_select]
# names(pheno[,pheno_select])
# [1] "gender"           "age"              "prep_hyb_batch"   "pmi..hr."
# [5] "last_update_date" "source_name_ch1"  "tissuebank"

batch <- pheno$prep_hyb_batch
bank <- pheno$tissuebank
date <- pheno$last_update_date
pmi <- pheno$pmi..hr. # note: missing values, how to handle?
# impute missing values for pmi (mean)
pmi[is.na(pmi)] <- mean(pmi[!is.na(pmi)])

tissue_f <- as.factor(pheno$source_name_ch1)
sex_f <- as.factor(pheno$gender)
age_f <- pheno$age

# HOW TO MODEL COVARIATES??????
mod <- model.matrix(~tissue_f+sex_f+age_f, data=pheno)
combat_edata1 <- ComBat(dat=edata, batch=batch, mod=mod, numCovs=c(6), par.prior=TRUE, prior.plots=TRUE)
combat_edata2 <- ComBat(dat=combat_edata1, batch=bank, mod=mod, numCovs=c(6), par.prior=TRUE, prior.plots=TRUE)

# Cannot remove date.
# combat_edata3 <- ComBat(dat=combat_edata2, batch=date, mod=mod, numCovs=c(6), par.prior=TRUE, prior.plots=TRUE)
## Error in solve.default(t(design) %*% design) : 
##   system is computationally singular: reciprocal condition number = 3.52356e-21


M.correct <- removeBatchEffect(combat_edata2, covariates=pmi, design=mod)
vis(M.correct, "batchCorrect7")

tissue_color <- rep("black", length(tissue_f))
tissue_color[tissue_f == "Human Brain Tissue: TCTX"] <- "green"
tissue_color[tissue_f == "Human Brain Tissue: FCTX"] <- "red"
tissue_color[tissue_f == "Human Brain Tissue: PONS"] <- "blue"

name<-"batchCorrect7"
pdf(paste0("plotSampleRelations_cluster_", name, ".pdf"), width=60, height=30)
plotSampleRelation(M.correct, method='cluster', cv.Th=0, color=tissue_color)
dev.off()

pdf(paste0("plotSampleRelations_mds_", name, ".pdf"), width=60, height=30)
plotSampleRelation(M.correct, method='mds', cv.Th=0, color=tissue_color)
dev.off()

crtx <- tissue_f == "Human Brain Tissue: TCTX" | tissue_f == "Human Brain Tissue: FCTX"

# plot only cortex
pdf(paste0("plotSampleRelations_cortex_cluster_", name, ".pdf"), width=60, height=30)
plotSampleRelation(M.correct[,crtx], method='cluster', cv.Th=0, color=tissue_color[crtx])
dev.off()
pdf(paste0("plotSampleRelations_cortex_mds_", name, ".pdf"), width=60, height=30)
plotSampleRelation(M.correct[,crtx], method='mds', cv.Th=0, color=tissue_color[crtx])
dev.off()


# SAVE TO FILE
# --------
save(M.correct, file="Methylation_M.correct.RData")
save(M.ssn, file="Methylation_M.ssn.RData")

write.table(M.correct, file="Methylation_M.correct.tab", sep="\t", quote=FALSE)
write.table(t(pData(M.ssn)), file="Methylation_M.correct.samples.tab", sep="\t", quote=FALSE)
