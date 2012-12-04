library(lumi)
library(limma)
library(illuminaHumanv2.db)
illuminaHumanv2()


# 1: Load and import data.
# ====================
# 1.1: Load raw data.
rawsum.lm <- lumiR("GSE15745_GPL6104.raw.tab", sep="\t", columnNameGrepPattern=
list(exprs="AVG_Signal", se.exprs="BEAD_STDERR", beadNum=
"Avg_NBEADS", detection="Detection Pval", verbose=TRUE))
# 1.2: Load sample attributes.
# NOTE: This should be added to `phenoData` parameter of LumiBatch object
AttrT <- read.table("gse15745.gpl6104.samples.raw_aligned.tab", sep="\t",
header = TRUE, row.names = 1)

# 2: Quality control by detection call
# --------------------
arrayDetections<-detectionCall(rawsum.lm, type='sample', Th=0.01)
print(summary(arrayDetections))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   3228   12680   13340   12980   13800   14710
outliers<-(mean(arrayDetections) - arrayDetections) >= sd(arrayDetections)*2
print(sum(outliers))

# 3: Filter arrays, get attributes, and visualize batch effects
# --------------------
M <-rawsum.lm[,!outliers]
write.table(AttrT[,!outliers], "good_sample_attrs.tab",
  sep="\t", quote=FALSE)

# 4: Normalize array using neqc
rawsum.neqc <- new('LumiBatch', neqc(M, detection=detection(M)), se.exprs=se.exprs(M),
  beadNum=beadNum(M), detection=detection(M))


# 5: Filter Probes
# --------------------
# 5.1 Probe Quality
library(illuminaHumanv2.db)
illuminaHumanv2()

ids <- as.character(featureNames(M))
qual <- unlist(mget(ids, illuminaHumanv2PROBEQUALITY, ifnotfound = NA))
print(table(qual))
print("Na: ")
print(length(qual[is.na(qual)]))
AveSignal <- rowMeans(exprs(rawsum.neqc))
pdf("boxplot_probe_qual.pdf", width=12, height=8)
boxplot(log2(AveSignal) ~ qual)
dev.off()
rem <- qual == "No match" | qual == "Bad" | is.na(qual)


# 5.2 Probe Calls
# --------------------
tissue <- as.vector(as.matrix(AttrT)[7,!outliers])
cReads <- detectionCall(M[,tissue == "cerebellum"], type="probe", Th=0.01)
fcReads <- detectionCall(M[,tissue == "frontal cortex"], type="probe", Th=0.01)
tcReads <- detectionCall(M[,tissue == "temporal cortex"], type="probe", Th=0.01)
pReads <- detectionCall(M[,tissue == "pons"], type="probe", Th=0.01)

callP <- 0.95
cProbes <- cReads > round(callP*sum(tissue=="cerebellum"))
fcProbes <- fcReads > round(callP*sum(tissue=="frontal cortex"))
tcProbes <- tcReads > round(callP*sum(tissue=="temporal cortex"))
pProbes <- pReads > round(callP*sum(tissue=="pons"))

rowsIntersect <- cProbes & fcProbes & tcProbes & pProbes & !rem
rowsUnion <- (cProbes | fcProbes | tcProbes | pProbes) & !rem
sum(rowsIntersect)
sum(rowsUnion)
# 8205
# 10277

Mu <- rawsum.neqc[rowsUnion,]

# 6: REMOVE BATCH EFFECTS
# --------------------
tissue <- as.vector(as.matrix(AttrT)[7,!outliers])
prep_hyb_batch <- as.vector(as.matrix(AttrT)[2,!outliers])
tissuebank <- as.vector(as.matrix(AttrT)[9,!outliers])
pmi <- as.numeric(as.vector(as.matrix(AttrT[3,!outliers])))
age <- as.numeric(as.vector(as.matrix(AttrT[6,!outliers])))
gender <- as.vector(as.matrix(AttrT[10,!outliers]))
labels <- as.vector(as.matrix(AttrT[11,!outliers]))
titles <- as.vector(as.matrix(AttrT[4,!outliers]))
sample_nums <- as.vector(sapply(titles, function (x) unlist(strsplit(x, "-"))[2]))

mod <- model.matrix(~as.factor(tissue)+as.factor(gender)+age, data=AttrT[,!outliers])

batch <- as.factor(prep_hyb_batch)
bank <- as.factor(tissuebank)
edata <- exprs(Mu)
combat_edata1 <- ComBat(dat=edata, batch=batch, mod=mod, numCovs=c(6), par.prior=TRUE, prior.plots=TRUE)
combat_edata2 <- ComBat(dat=combat_edata1, batch=bank, mod=mod, numCovs=c(6), par.prior=TRUE, prior.plots=TRUE)
M.correct <- removeBatchEffect(combat_edata2, covariates=pmi, design=mod)


# 7: SAVE DATA
save(M.correct, file="mRNA_M.correct.RData")
write.table(M.correct, file="mRNA_M.correct.tab", sep="\t", quote=FALSE)
write.table(AttrT[,!outliers], file="mRNA_M.correct.samples.tab", sep="\t", quote=FALSE)


tissue_color <- rep("black", length(tissue))
tissue_color[tissue == "temporal cortex"] <- "green"
tissue_color[tissue == "frontal cortex"] <- "red"
tissue_color[tissue == "pons"] <- "blue"

name<-"mRNA_Correct"
pdf(paste0("plotSampleRelations_cluster_", name, ".pdf"), width=60, height=30)
plotSampleRelation(M.correct, method='cluster', cv.Th=0)
dev.off()

pdf(paste0("plotSampleRelations_mds_", name, ".pdf"), width=60, height=30)
plotSampleRelation(M.correct, method='mds', cv.Th=0, color=tissue_color)
dev.off()

library("RColorBrewer")
heatmap_cols <- rev(colorRampPalette(brewer.pal(8,"RdYlBu"))(50))

iqrange <- apply(M.correct, 1, IQR, na.rm = TRUE)
topVar <- order(iqrange, decreasing = TRUE)[1:500]
ids <- as.character(rownames(M.correct[topVar, ]))
syms <- unlist(mget(ids, illuminaHumanv2SYMBOL, ifnotfound = NA))
qual <- unlist(mget(ids, illuminaHumanv2PROBEQUALITY, ifnotfound = NA))
rowLabels <- as.list(paste(syms, qual))
rowColors <- rep("white", length(syms))
rowColors[is.na(syms)]<- "orange"
rowColors[qual == "Not Found" | is.na(qual) | qual == "Bad"] <- "red"

d <- dist(t(M.correct[topVar, ]))
hcluster <- hclust(d, method="average")

name<-"mRNA_M.correct.top500"
pdf(paste0("heatmap_", name, ".pdf"), width=60, height=60)
heatmap(M.correct[topVar, ], Colv=as.dendrogram(hcluster), labCol=colnames(M.correct), 
    col=heatmap_cols, ColSideColors=tissue_color, labRow=rowLabels,
    RowSideColors=rowColors)
  dev.off()
