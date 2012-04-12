library(Rsamtools)
library(Biostrings)
library(GenomicRanges)
library(GenomicFeatures)
library(DESeq)
library(yaml)
source("diffExpFunctions.R")
p = yaml.load_file("params.yaml")

loadConfig(aligner=p$aligner,reference=p$reference,alignmentStringency=p$alignmentStringency,alignmentStrategy=p$alignmentStrategy)
GR<-getGRFromRefGene()
cList<-getCounts(GR,bamPaths, bamSamples, samples,minCount=p$minCount)

#fits: parameteric or local
#methods: pooled per-condition
#sharingMode="fit-only" gene-est-only maximum 
res<-doStats(cList$cds,fits='local',method='per-condition',sharingMode='maximum')
save(list=c("res","cList","p"),file="diffExpGenomic.RData",compress=TRUE)
