library(lumi)
library(IlluminaHumanMethylation27k.db)
library(FDb.InfiniumMethylation.hg19)
library(IlluminaHumanMethylation450k.db)
IlluminaHumanMethylation450k()
IlluminaHumanMethylation27k()

# this is build 36!
#IlluminaHumanMethylation27kCPG37 <- IlluminaHumanMethylation27kCPGCOORDINATE
# based on hg19
IlluminaHumanMethylation27kCPG37 <- IlluminaHumanMethylation27kCHRLOC
IlluminaHumanMethylation27kCHR37 <- IlluminaHumanMethylation27kCHR


fileName <- "test4.txt"
example.lumiMethy <- lumiMethyR(fileName,
  lib="IlluminaHumanMethylation27k.db", sep="\t")

# Update build 36 positions to build 37 positions
#   something like: featureData(example.lumiMethy)$POSITION[1:2] <- c(4,5)
# ====
# get vector of chr, cpg for build 37 per probe id in original order
# use FDb.InfiniumMethylation.hg19
# hm27 <- get27k()


featureData(example.lumiMethy)$POSITION[1]


samps <- read.table(system.file("extdata/samples.txt",
  package = "methylumi"),sep="\t",header=TRUE)
## Perform the actual data reading
## This is an example of reading data from a
## Sentrix Array format file (actually two files,
## one for data and one for QC probes)

mldat <- methylumiR(system.file('extdata/exampledata.samples.txt',
   package='methylumi'),
   qcfile=system.file('extdata/exampledata.controls.txt',
   package="methylumi"),
   sampleDescriptions=samps)

## sampleDescriptions: A data.frame that contains at least one column,
##   SampleID (case insensitive).  This column MUST match the part
##   of the column headers before the ‘.Avg_Beta’, etc.  Also, if
MethyLumiSet



## experiment: get CPGs
cpg27k <- IlluminaHumanMethylation450k_get27k()
x <- as.vector(unlist(cpg27k))
