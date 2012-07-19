library(rtracklayer)
library(stringr)

loadConfig <- function (aligner, alignmentStringency, reference, alignmentStrategy) {

  #bamDirectory<-concat("/bam",aligner,alignmentStringency,reference,alignmentStrategy,"/",sep="/")
  assign("bamDirectory", concat("/bam",aligner,alignmentStringency,reference,alignmentStrategy,"/",sep="/"), envir=.GlobalEnv)
  configFile<-'/nas/is1/leipzig/src/R/dirConfig.R'
  source(configFile)
}

seqlengths2gr <- function(x, strand="*"){
  ## convert 'seqlengths' of BSgenome to sGRanges
  GRanges(names(x), IRanges(1, x), strand=strand)
}
getGRFromGenome<-function(){
  session <- browserSession()
  genome(session) <- "hg19"
  GR <- seqlengths2gr(seqlengths(session))
  seqlevels(GR,force=TRUE)<-c(concat("chr",1:22),"chrX","chrY","chrM")
  elementMetadata(GR)<-NULL
  elementMetadata(GR)$id<-1:25
  elementMetadata(GR)$name<-c(concat("chr",1:22),"chrX","chrY","chrM")
  GR
}

getGRFromRefGene<-function(){
  knownGenes<-makeTranscriptDbFromUCSC(genome="hg19",
                                       tablename="refGene",
                                       transcript_ids=NULL,
                                       circ_seqs=DEFAULT_CIRC_SEQS,
                                       url="http://genome.ucsc.edu/cgi-bin/",
                                       goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath"
                                       )                        
  GR <- transcripts(knownGenes)
  names(elementMetadata(GR))<-c("id","name")
  seqlevels(GR,force=TRUE)<-c(concat("chr",1:22),"chrX","chrY","chrM")
  GR
}
getGRFromGFF<-function(gffFile){
  #concat(refsDir,"/hsa.chr.gff")
  gffRangedData<-import.gff(gffFile)
  #GRanges with 1523 ranges and 5 elementMetadata values:
  #  seqnames                 ranges strand   |     type   source    phase                                 group          tx_id
  #<Rle>              <IRanges>  <Rle>   | <factor> <factor> <factor>                           <character>    <character>
  #  NA     chr1     [  30366,   30503]      +   |    miRNA     <NA>     <NA> ACC="MI0006363"; ID="hsa-mir-1302-2"; hsa-mir-1302-2
  GR <- as(gffRangedData, "GRanges")
  accessions<-as.integer(str_match(as.character(elementMetadata(GR)$group),"MI([0-9]+)")[,2])
  mirnaNames<-as.character(str_match(as.character(elementMetadata(GR)$group),"hsa[a-zA-Z0-9-]+"))

  elementMetadata(GR)<-NULL
  elementMetadata(GR)$id<-accessions
  elementMetadata(GR)$name<-mirnaNames
  GR
}
getCounts <- function (GR,bamPaths, bamSamples, samples,minCount) {
  bamView<-BamViews(bamPaths=bamPaths,bamSamples=bamSamples,bamRanges=GR)
  bamcounts<-countBam(bamView)
  bfl <- BamFileList(bamPaths)
  olap <- summarizeOverlaps(GR, bfl)
  counts<-as.data.frame(assays(olap)$counts)
  colnames(counts)<-samples
  #technical replicates
  counts$WERI<-counts[,'31s8_WERI']+counts[,'36s2_WERI']+counts[,'WERI01_PGM']+counts[,'WERI02_PGM']+counts[,'WERI03_PGM']
  counts$Y79<-counts[,'36s1_Y79']+counts[,'36s2_Y79']
  counts$RB525T<-counts[,'RB525T']+counts[,'RB525T01_PGM']+counts[,'RB525T02_PGM']+counts[,'RB525T03_PGM']
  counts<-counts[,!str_detect(colnames(counts),'^3')]
  counts<-counts[,!str_detect(colnames(counts),'WERI.*PGM')]
  counts<-counts[,!str_detect(colnames(counts),'RB525T.*PGM')]
  #this messes up the id.var so matrix doesn't work
  
  #42s1_nrml  RB494N RB494T RB495N  RB495T RB498N RB498T RB517T  RB525T    WERI    Y79 (11 groups)
  conds<-c('N','N','T','N','T','N','T','T','T','T','T')
  individual<-c(1,2,2,3,3,4,4,5,6,7,8)
  replicate<-c(1,1,2,1,2,1,2,2,2,2,2)
  type<-c('IIfx','HS','HS','HS','HS','HS','HS','HS','HS','IIfx','IIfx')
  paired<-c(FALSE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE)
  pdata<-data.frame(condition=conds,replicate=replicate,type=type,individual=individual,paired=paired)
  rownames(pdata)<-c('42s1_nrml','RB494N','RB494T','RB495N','RB495T','RB498N','RB498T','RB517T','RB525T','WERI','Y79')
  
 
  rownames(counts)<-elementMetadata(rowData(olap))$tx_id
  
  #at least minCount from the paired lanes
  countsAboveThreshold<-subset(counts,rowSums(counts[,row.names(pdata[which(pdata$paired==TRUE),])])>=minCount)
  #subset(counts,rowSums(counts[-1])>minCount)
  
  countsMatrix<-as.matrix(countsAboveThreshold)
  
  
  colnames(countsMatrix)<-colnames(countsAboveThreshold)
  cds <- newCountDataSet( countsMatrix, pdata$condition )
  list(bamcounts=bamcounts,counts=counts,cds=cds,countsMatrix=countsMatrix)
}

getCounts2 <- function (GR,bamPaths, bamSamples, samples,minCount,mode) {
  bfl <- BamFileList(bamPaths)
  olap <- summarizeOverlaps(GR, bfl,mode=mode)
  precounts<-as.data.frame(assays(olap)$counts)
  colnames(counts)<-samples
  #technical replicates
  #compute from metadatalayerCake<-data.frame(class=c("exon","tx","intergenic","unmapped"),fraction=c(.02,.25,.50,.23))

#   counts$WERI<-counts[,'31s8_WERI']+counts[,'36s2_WERI']+counts[,'WERI01_PGM']+counts[,'WERI02_PGM']+counts[,'WERI03_PGM']
#   counts$Y79<-counts[,'36s1_Y79']+counts[,'36s2_Y79']
#   counts$RB525T<-counts[,'RB525T']+counts[,'RB525T01_PGM']+counts[,'RB525T02_PGM']+counts[,'RB525T03_PGM']
#   counts<-counts[,!str_detect(colnames(counts),'^3')]
#   counts<-counts[,!str_detect(colnames(counts),'WERI.*PGM')]
#   counts<-counts[,!str_detect(colnames(counts),'RB525T.*PGM')]
#   #this messes up the id.var so matrix doesn't work
#   
#   #42s1_nrml  RB494N RB494T RB495N  RB495T RB498N RB498T RB517T  RB525T    WERI    Y79 (11 groups)
#   conds<-c('N','N','T','N','T','N','T','T','T','T','T')
#   individual<-c(1,2,2,3,3,4,4,5,6,7,8)
#   replicate<-c(1,1,2,1,2,1,2,2,2,2,2)
#   type<-c('IIfx','HS','HS','HS','HS','HS','HS','HS','HS','IIfx','IIfx')
#   paired<-c(FALSE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE)
#   pdata<-data.frame(condition=conds,replicate=replicate,type=type,individual=individual,paired=paired)
#   rownames(pdata)<-c('42s1_nrml','RB494N','RB494T','RB495N','RB495T','RB498N','RB498T','RB517T','RB525T','WERI','Y79')
  
  
  rownames(counts)<-elementMetadata(rowData(olap))$tx_id
  
  #at least minCount from the paired lanes
  countsAboveThreshold<-subset(counts,rowSums(counts[,row.names(pdata[which(pdata$paired==TRUE),])])>=minCount)
  #subset(counts,rowSums(counts[-1])>minCount)
  
  countsMatrix<-as.matrix(countsAboveThreshold)
  
  
  colnames(countsMatrix)<-colnames(countsAboveThreshold)
  cds <- newCountDataSet( countsMatrix, pdata$condition )
  list(counts=counts,cds=cds,countsMatrix=countsMatrix)
}

doStats <- function (cds, fits, method, sharingMode) {
  cds <- estimateSizeFactors( cds )
  cds <- estimateDispersions( cds , fitType=fits, method=method, sharingMode=sharingMode)

  res <- nbinomTest( cds, "N", "T")
  
  #pds <- newCountDataSet( countsAboveThreshold, pdata )
  #pds <- estimateDispersions( cds , fitType="local", method='pooled', sharingMode='fit-only')
  #fit1 <- fitNbinomGLMs( pds, count ~  type + condition )
  #fit0 <- fitNbinomGLMs( pds, count ~  type )
  #pvalsGLM <- nbinomGLMTest( fit1, fit0 )
  #padjGLM <- p.adjust( pvalsGLM, method="BH" )
  #res$pvalsGLM<-pvalsGLM
  #res$padjGLM<-padjGLM
  write.csv(res,file="deseqResults.csv")
  res
}

