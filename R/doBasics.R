library(Rsamtools)
library(Biostrings)
library(GenomicRanges)
library(GenomicFeatures)
library(DESeq)
library(yaml)
source("diffExp/diffExpFunctions.R")
source("dirConfigTest.R")
p = yaml.load_file("diffExp/params.yaml")

loadConfig(aligner='novo',reference='hg19.ambig',alignmentStringency='tight',alignmentStrategy='random')
#GRgenome<-getGRFromGenome()
#GRrefgene<-getGRFromRefGene()
#GRmirna<-getGRFromGFF(concat(refsDir,"/hsa.chr.gff"))

load("savedGenomicRanges/refgene.gr.RData")
load("savedGenomicRanges/genome.gr.RData")
load("savedGenomicRanges/mirna.gr.RData")

grl <- GRangesList(mirna=GRmirna,refgene=GRrefgene,genome=GRgenome)
cList<-getCounts2(grl,bamPaths[19], bamSamples[19,], samples[19],minCount=0,mode=Raincatcher)

GR<-getGRFromRefGene()
cList<-getCounts(GR,bamPaths, bamSamples, samples,minCount=0)

#fits: parameteric or local
#methods: pooled per-condition
#sharingMode="fit-only" gene-est-only maximum 
res<-doStats(cList$cds,fits='local',method='per-condition',sharingMode='maximum')
save(list=c("res","cList","p"),file="diffExpGenomic.RData",compress=TRUE)
