library(Rsamtools)
library(reshape)
library(mirbase.db)
library(stringr)
library(IRanges)
library(plyr)
library(doMC)
library(gtools)

#library(rredis)
#redisConnect()
registerDoMC()

##setup
#' This is data to be included in my package
#'
#' @name Ganguly Retinoblastoma HiSeq Runs
#' @docType data
#' @author Jeremy Leipzig \email{leipzigj@email.chop.edu}
#' @keywords bam
NULL

hiseqDataDir<-"/nas/is1/leipzig/Ganguly/gangulyRBhi/data/"
solexaDataDir<-"/nas/is1/leipzig/Ganguly/gangulySolexa/RB_miRNA_raw_data/reprocess/"
bamDirectory<-"bam/"
outputDir="/nas/is1/leipzig/Ganguly/gangulyRBhi/results/readlevelView/"
hiseqconds<-c("N","T","N","T","N","T","T","T")
hiseqsamples<-c("RB494N","RB494T","RB495N","RB495T","RB498N","RB498T","RB517T","RB525T")
solexaconds<-c("T","T","T","T","N")
solexasamples<-c("FGC0031_s_8.WERI","FGC0036_s_1.Y79","FGC0036_s_2.WERI","FGC0036_s_2.Y79","FGC0042_s_1.normal")
solexaprettynames<-c("31s8_WERI","36s1_Y79","36s2_WERI","36s2_Y79","42s1_nrml")
hiseqbamPaths<-concat(hiseqDataDir,bamDirectory,hiseqsamples,'.hairpin.sam.sorted.bam')
solexabamPaths<-concat(solexaDataDir,bamDirectory,solexasamples,'.bam')
bamPaths<-c(hiseqbamPaths,solexabamPaths)
samples<-c(hiseqsamples,solexaprettynames)
conds<-c(hiseqconds,solexaconds)
bamSamples<-DataFrame(conds=conds,row.names=samples)

#' Render a hairpin sequence properly placed mature sequences
#'
#' @param miRNA miRNA name
#' @keywords manip
#' @export
#' @examples
#' matureOnHairpin('hairpin-let-7a-1')
#' @return character[]
matureOnHairpin<-function(miRNA){
  seq<-hairpinLookup(miRNA)
  spacechar<-' '
  newline<-"<br>"
  myrange<-IRanges(start=matureFrom(mirnaMature[miRNA][[1]]),end=matureTo(mirnaMature[miRNA][[1]]),names=matureName(mirnaMature[miRNA][[1]]))
  dashes<-rep(spacechar,max(str_length(as.character(seq)),0))
  dashes[as.integer(myrange[str_sub(names(myrange),-1)=='*'])]<-'*'
  dashes[as.integer(myrange[str_sub(names(myrange),-1)!='*'])]<-'='
  dashstr<-concat(dashes,collapse="")
  concat(seq,dashstr,sep=newline)
}

#' Concatenate SQL-style
#' This is like paste but returns NULL anytime one of the given strings is NULL
#'
#' @param character
#' @export
#' @examples
#' concat
#' @return character
concat<-function(...,sep="",collapse=NULL){
  strings<-list(...)
  #NULL, NA
  if(
    all(unlist(llply(strings,length))>0)
    &&
    all(!is.na(unlist(strings)))
    ){
    do.call("paste", c(strings, list(sep = sep, collapse = collapse)))
  }else{
    NULL
  }
}

#' Placeholder for possible caching scenario
#' If we want to implement caching the edits we should do it here
fetchEdits<-function(targetseq,queryseq,...){
  #fetch<-redisGet(storedEdits[[targetseq]][[queryseq]])
  #if(is.null(fetch)){
    fetch<-calcEdits(targetseq,queryseq,...)
    #redisSet(storedEdits[[targetseq]][[queryseq]],fetch)
  #}
  fetch
}


calcEdits<-function(targetseq,queryseq,N_is_reference=TRUE){
  pwa<-pairwiseAlignment(targetseq,queryseq,type="local-global")

  #used for target positions
  target_poss_sub<-unlist(pwa@pattern@mismatch)
  #used for identifying the edits
  seq_poss_sub<-unlist(pwa@subject@mismatch)
  
  #pattern indels is given relative to the subject
  #pattern: [15] TCCCTGAGACCCTT-TAACCTGT 
  #subject:  [1] TCCCTGAGACCCTTTTAACCTGT 
  #"30i"
  target_poss_indel_ins<-concat(as.character(start(unlist(pwa@pattern@indel))+pwa@pattern@range@start),'i',sep='')
  seq_poss_indel_ins<-start(unlist(pwa@pattern@indel))
  
  
  #we won't see deletions on the target side
  #1234567890123456789012345678901234567
  #TGCCAGTCTCTAGGTCCCTGAGACCCTTTAACCTGTGAGGACATCCAGGGTCACAGGTGA
  #              TCCCTGAGACCCTT-AACCTGTG
  #pattern: [15] TCCCTGAGACCCTTTAACCTGTG 
  #subject:  [1] TCCCTGAGACCCTT-AACCTGTG 
  target_poss_indel_del<-concat(as.character(start(unlist(pwa@subject@indel))+pwa@pattern@range@start-1),'d',sep='')
  seq_poss_indel_del<-start(unlist(pwa@subject@indel))
  
  target_poss<-c(target_poss_sub,target_poss_indel_ins,target_poss_indel_del)
  seq_poss<-c(seq_poss_sub,seq_poss_indel_del,seq_poss_indel_ins)
  
  nucs<-strsplit(queryseq,'')[[1]][seq_poss]
  edits<-data.frame(nucs=nucs,poss=target_poss)
  editList<-list()
  if(N_is_reference){
    editList<-as.list(edits[edits$nucs!='N',])
  }else{
    editList<-as.list(edits)
  }
  #remove N's
  maskedSeq<-maskedEdits(as.character(subject(pwa)),as.character(pattern(pwa)))
  spacedSeq<-str_pad(maskedSeq,pwa@pattern@range@start+pwa@pattern@range@width-1)
  
  #formatEdits(nucs,target_poss)
  
  list(editList=editList,spacedSeq=spacedSeq)
}

maskedEdits<-function(subject,pattern){
  concat(apply(X=rbind(strsplit(subject,'')[[1]],strsplit(pattern,'')[[1]]),MARGIN=2,FUN=function(x){if(x[1]=='N'){x[2]}else{if(x[1]=='-'){''}else{x[1]}}}),collapse="")
}

formatEdits<-function(editList){
  edits<-vector()
  for(j in 1:length(editList$poss)){
   edits<-c(edits,concat(editList$poss[j],editList$nucs[j],sep=':'))
  }
  if(length(edits)==0){return("None")}else{return(concat(edits,collapse=" "))}
}

readlevel<-function(x){
   x$seq<-as.character(x$seq)
   x$qual<-as.character(x$qual)
   if(length(x$seq)>0){
       df<-as.data.frame(x)
       udf<-ddply(df,.(rname,seq),.fun=function(x){
            tmp <- x[1,]
            edits<-fetchEdits(hairpinLookup(as.character(tmp$rname)),as.character(tmp$seq))
            data.frame(spacedSeq=edits$spacedSeq,edits=formatEdits(edits$editList),nrow=nrow(x))
       })
  }
}
                   


##main
getBams<-function(){

  hairpinGR<-GRanges(names(hairpins),IRanges(1,width(hairpins)),strand="*")
  
  bamView<-BamViews(bamPaths=bamPaths,
  bamSamples=bamSamples,
  bamRanges=hairpinGR)
  
  what <- c("rname", "strand", "seq", "qual")
  param<-ScanBamParam(what=what,which=hairpinGR)
  #bams is a [[list of hiseqsamples[[list of genes]]]
  bams <- scanBam(bamView, param=param)
  unlist(bams)
}

#if(!exists("mirnaSEQ")){
  #x <- mirbaseSEQUENCE
  #mirnaSEQ <- mget(mappedkeys(x), x)
  #x <- mirbaseMATURE
  #mirnaMature<-mget(mappedkeys(x),x)


  hairpins <- read.DNAStringSet(paste(hiseqDataDir,"refs/hsa_hairpin.dna.fa",sep=""),use.names=TRUE)
  names(hairpins)<-str_match(names(hairpins),"^\\S+") 
  hairpinLookup<-function(x){as.character(hairpins)[x]}
  
  
  if(exists("hairpinSubset")){
    hairpins<-subset(hairpins,names(hairpins) %in% hairpinSubset)
  }   
  
  if(length(hairpins)==0){stop("hairpins cannot be of length 0")}
  
  if(!exists("skipBams") || skipBams!=TRUE){
    ubams<-getBams()
  }
  

  #ubams is a list of sample.genes i.e. $ RB498T.hsa-let-7i:1-60
  ureadlevel<-ldply(ubams,.fun=readlevel,.parallel=TRUE,.progress="none")
  if(nrow(ureadlevel)==0){stop("empty result set")}
  
  ureadlevel$sample<-str_split_fixed(ureadlevel$.id,'\\.',2)[,1]
  names(ureadlevel)<-c("id","rname","seq","spacedSeq","edits","value","variable")
  recasted<-cast(ureadlevel[,-1],fun.aggregate=sum)

  #for rstudio display
  res<-as.data.frame(recasted)
  d_ply(recasted,.(rname),.fun=function(x){
        save(x,file=concat("readlevels/",x$rname[1],".RData"),compress=TRUE)
        names(x)<-str_replace_all(names(x),'RB','')
        excelTable<-x[mixedorder(x$edits),]
        names(excelTable)[3]<-hairpinLookup(x$rname[1])
        write.csv(excelTable,file=concat("csv/",x$rname[1],"csv",sep="."))
  })
