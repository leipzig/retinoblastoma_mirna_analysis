
loadConfig <- function (aligner, alignmentStringency, reference, alignmentStrategy) {
  alignmentParams<-list()
  alignmentParams[['loose']]<-'-l 17 -h 60 -t 65'
  alignmentParams[['tight']]<-'-l 17 -h 0 -t 0'
  
  bamDirectory<-concat("/bam",aligner,alignmentStringency,reference,alignmentStrategy,"/",sep="/")
  configFile<-'/nas/is1/leipzig/src/R/dirConfig.R'
  source(configFile)
}



getCounts <- function (bamPaths, bamSamples, samples) {
  knownGenes<-makeTranscriptDbFromUCSC(genome="hg19",
                                       tablename="refGene",
                                       transcript_ids=NULL,
                                       circ_seqs=DEFAULT_CIRC_SEQS,
                                       url="http://genome.ucsc.edu/cgi-bin/",
                                       goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath"
                                       )                        
  GR <- transcripts(knownGenes)
  seqlevels(GR,force=TRUE)<-c(concat("chr",1:22),"chrX","chrY","chrM")
  
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
  counts
}

getCds <- function (counts,minCount) {
  #42s1_nrml  RB494N RB494T RB495N  RB495T RB498N RB498T RB517T  RB525T    WERI    Y79 (11 groups)
  conds<-c('N','N','T','N','T','N','T','T','T','T','T')
  individual<-c(1,2,2,3,3,4,4,5,6,7,8)
  replicate<-c(1,1,2,1,2,1,2,2,2,2,2)
  type<-c('IIfx','HS','HS','HS','HS','HS','HS','HS','HS','IIfx','IIfx')
  paired<-c(FALSE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE)
  pdata<-data.frame(condition=conds,replicate=replicate,type=type,individual=individual,paired=paired)
  rownames(pdata)<-c('42s1_nrml','RB494N','RB494T','RB495N','RB495T','RB498N','RB498T','RB517T','RB525T','WERI','Y79')
  
  #at least minCount from the paired lanes
  
  countsAboveThreshold<-subset(counts,rowSums(counts[,row.names(pdata[which(pdata$paired==TRUE),])])>=minCount)
  #subset(counts,rowSums(counts[-1])>minCount)
  
  countsMatrix<-as.matrix(countsAboveThreshold)
  rownames(countsMatrix)<-elementMetadata(rowData(olap))$tx_name
  
  colnames(countsMatrix)<-colnames(countsAboveThreshold)
  cds <- newCountDataSet( countsMatrix, pdata$condition )
  cds
}

doStats <- function (cds, fits, method, sharingMode) {
  cds <- estimateDispersions( cds , fitType=fits, method=method, sharingMode=sharingMode)
  cds <- estimateSizeFactors( cds )
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