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


outputDir="/nas/is1/leipzig/Ganguly/gangulyRBhi/results/readlevelView/"

aligner<-"novo"
reference<-"hairpin"
alignmentStringency<-"loose"
alignmentStrategy<-"none"
alignmentParams<-list()
alignmentParams[['loose']]<-'-l 17 -h 60 -t 65'
alignmentParams[['tight']]<-'-l 17 -h 0 -t 0'

bamDirectory<-concat("/bam",aligner,alignmentStringency,reference,alignmentStrategy,"/",sep="/")
configFile<-'/nas/is1/leipzig/src/R/dirConfig.R'
source(configFile)


if(!exists("mirnaSEQ")){
  x <- mirbaseSEQUENCE
  mirnaSEQ <- mget(mappedkeys(x), x)
  x <- mirbaseMATURE
  mirnaMature<-mget(mappedkeys(x),x)
}
#' Render a hairpin sequence properly placed mature sequences
#'
#' @param miRNA miRNA name
#' @keywords manip
#' @export
#' @examples
#' matureOnHairpin('hsa-let-7a-1')
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

getMature<-function(miRNA){
  mat<-NULL
  mat<-mirnaMature[miRNA][[1]]
  #some of these are listed without the -1
  if(is.null(mat)){mat<-mirnaMature[str_replace(miRNA,pattern="-[0-9]+$",replacement='')][[1]]}
  #some of these are listed without the a/b f
  if(is.null(mat)){mat<-mirnaMature[str_replace(miRNA,pattern="[ab]$",replacement='')][[1]]}
  #if(is.null(mat)){stop(concat("can't find the miRNA:",miRNA," in the mature lookup"))}
  mat
}
getRangesOfMature<-function(miRNA){
  if(is.null(getMature(miRNA))){
      myrange<-IRanges()
  }else{
      myrange<-IRanges(start=matureFrom(getMature(miRNA)),end=matureTo(getMature(miRNA)),names=matureName(getMature(miRNA)))
  }
  myrange
}

getReport<-function(hairpinName,querySeq){
  hairpin<-hairpinLookup(hairpinName)
  matRanges<-getRangesOfMature(hairpinName)
  #mutation at position 10
  editObj<-calcEdits(hairpin,querySeq)
  editsRelative(editObj,matRanges)
}

editsRelative<-function(editObj,matRanges){
  edits<-tas<-ntas<-isos<-status<-NULL
  if(length(matRanges)==0){
    tas<-edits<-ntas<-isos<-"N/A"
    status<-"No mature data"
  }else{
    editRanges<-IRanges(editObj$editList$poss,editObj$editList$poss)
    editOlap<-findOverlaps(editRanges,matRanges)
    if(length(editOlap)>0){
        edits<-formatInternalEdits(editOlap,editObj,matRanges)
        status<-"Mut"
    }else{
      edits<-"None"
    }
    
    #do we catch up to 3 past the end of the mature?
    flanks<-flank(matRanges,3,start=FALSE,both=FALSE)
    addOlap<-findOverlaps(editObj$refRanges,flanks)
    ntaOlap<-findOverlaps(editRanges,flanks)
    if(length(addOlap)>0){
      if(length(ntaOlap)>0){
        #if we have any nt_additions, the whole thing is nt, who are we kidding?
        tas<-"None"
        ntas<-formatNTAs(ntaOlap,editObj,matRanges)
      }else{
        ntas<-"None"
        tas<-formatTAs(addOlap,editObj,matRanges)
      }
      status<-"Mut"
    }else{
      tas<-ntas<-"None"
    }
    if(any(editObj$refRanges == matRanges)){
      status="Exact mature"
      isos<-"None"
    }else{
      #we might consider only exacts to be isomirs
      #if so add length(editOlap)==0 to these conditions
      if(length(addOlap)==0){
        #isomir or outer mature hairpin
        isoOlap<-findOverlaps(editObj$refRanges,matRanges)
        if(length(isoOlap)>0){
          isos<-formatIsomirs(isoOlap,editObj,matRanges)
          status<-"Mut"
        }else{
          isos<-"None"
          status<-"No mature olap"
        }
      }else{
        isos<-"None"
        stopifnot(!is.null(status))
      }
    }
  }#mature data available
  list(editObj=editObj,edits=edits,tas=tas,ntas=ntas,isos=isos,status=status)
}

formatIsomirs<-function(isomirs,editObj,matRanges){
  edits<-vector()
  for(j in seq(length=length(isomirs))){
    queryHits<-queryHits(isomirs)[j]
    subjectHits<-subjectHits(isomirs)[j]
    edits<-c(edits,concat('iso:',start(matRanges[subjectHits])-start(editObj$refRanges[1]),'/',end(editObj$refRanges[1])-end(matRanges[subjectHits])))
  }
  if(length(edits)==0){return("None")}else{return(concat(edits,collapse=" "))}
}

formatTAs<-function(additions,editObj,matRanges){
  edits<-vector()
  for(j in seq(length=length(additions))){
    queryHits<-queryHits(additions)[j]
    subjectHits<-subjectHits(additions)[j]
    edits<-c(edits,concat('TA:',end(editObj$refRanges[1])-end(matRanges[subjectHits])))
  }
  if(length(edits)==0){return("None")}else{return(concat(edits,collapse=" "))}
}
  
formatNTAs<-function(nt_additions,editObj,matRanges){
  edits<-vector()
  for(j in seq(length=length(nt_additions))){
    queryHits<-queryHits(nt_additions)[j]
    subjectHits<-subjectHits(nt_additions)[j]
    edits<-c(edits,concat(editObj$editList$poss[queryHits]-start(matRanges[subjectHits])+1,'NTA:',editObj$editList$nucs[queryHits],sep=''))
  }
  if(length(edits)==0){return("None")}else{return(concat(edits,collapse=" "))}
}

  
formatInternalEdits<-function(internalEdits,editObj,matRanges){
  edits<-vector()
  for(j in seq(length=length(internalEdits))){
    queryHits<-queryHits(internalEdits)[j]
    subjectHits<-subjectHits(internalEdits)[j]
   edits<-c(edits,concat(editObj$editList$poss[queryHits]-start(matRanges[subjectHits])+1,editObj$editList$muts[queryHits],':',editObj$editList$refs[queryHits],'>',editObj$editList$nucs[queryHits],sep=''))
  }
  if(length(edits)==0){return("None")}else{return(concat(edits,collapse=" "))}
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
  #target_poss_indel_ins<-concat(as.character(start(unlist(pwa@pattern@indel))+pwa@pattern@range@start),'i',sep='')
  target_poss_indel_ins<-start(unlist(pwa@pattern@indel))+pwa@pattern@range@start-2
  seq_poss_indel_ins<-start(unlist(pwa@pattern@indel))
  
  
  #we won't see deletions on the target side
  #1234567890123456789012345678901234567
  #TGCCAGTCTCTAGGTCCCTGAGACCCTTTAACCTGTGAGGACATCCAGGGTCACAGGTGA
  #              TCCCTGAGACCCTT-AACCTGTG
  #pattern: [15] TCCCTGAGACCCTTTAACCTGTG 
  #subject:  [1] TCCCTGAGACCCTT-AACCTGTG 
  #target_poss_indel_del<-concat(as.character(start(unlist(pwa@subject@indel))+pwa@pattern@range@start-1),'d',sep='')
  target_poss_indel_del<-start(unlist(pwa@subject@indel))+pwa@pattern@range@start-1
  seq_poss_indel_del<-start(unlist(pwa@subject@indel))
  
  
  target_poss<-c(target_poss_sub,target_poss_indel_ins,target_poss_indel_del)
  seq_poss<-c(seq_poss_sub,seq_poss_indel_del,seq_poss_indel_ins)
  
  mut_types<-c(rep('s',length(seq_poss_sub)),rep('i',length(seq_poss_indel_ins)),rep('d',length(seq_poss_indel_del)))

  
  nucs<-strsplit(queryseq,'')[[1]][seq_poss]
  refs<-strsplit(targetseq,'')[[1]][target_poss]
  edits<-data.frame(nucs=nucs,refs=refs,poss=target_poss,muts=mut_types)
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
  
  list(editList=editList,refRanges=pwa@pattern@range,spacedSeq=spacedSeq)
}

maskedEdits<-function(subject,pattern){
  concat(apply(X=rbind(strsplit(subject,'')[[1]],strsplit(pattern,'')[[1]]),MARGIN=2,FUN=function(x){if(x[1]=='N'){x[2]}else{if(x[1]=='-'){''}else{x[1]}}}),collapse="")
}

formatEdits<-function(editList){
  edits<-vector()
  for(j in 1:length(editList$poss)){
   edits<-c(edits,concat(editList$poss[j],':',editList$refs[j],'>',editList$nucs[j],sep=''))
  }
  if(length(edits)==0){return("None")}else{return(concat(edits,collapse=" "))}
}

genEditRep<-function(x){
            if(nrow(x)>0){
              tmp <- x[1,]
              #cat(paste(as.character(tmp$rname),as.character(tmp$seq),"\n"))
             #editReport<-getReport("hsa-let-7a-1","TGGAGGTAGTAGGTTGTATAGTTTTA")
              editReport<-getReport(as.character(tmp$rname),as.character(tmp$seq))
              return(data.frame(spacedSeq=editReport$editObj$spacedSeq,edits=editReport$edits,tas=editReport$tas,ntas=editReport$ntas,isos=editReport$isos,status=editReport$status,nrow=nrow(x)))
            }else{
              return(NULL)
            }
}
readlevel<-function(abam){
   abam$seq<-as.character(abam$seq)
   abam$qual<-as.character(abam$qual)
   udf<-NULL
   if(length(abam$seq)>0){
       df<-as.data.frame(abam)
       udf<-ddply(df,.(rname,seq),genEditRep)
  }
  return(udf)
}
                   
testRD<-function(i){
  cat(i)
  abam<-ubams[[i]]
  abam$seq<-as.character(abam$seq)
  abam$qual<-as.character(abam$qual)
  df<-as.data.frame(abam)
  if(nrow(df)>0){
    genEditRep(df[1,])
  }
}
#for(i in 1:length(ubams)){testRD(i)}

testReadLevel<-function(i){
  cat(i)
  abam<-ubams[[i]]
  udf<-readlevel(abam)
}
#for(i in 1:length(ubams)){testReadLevel(i)}  
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




  hairpins <- read.DNAStringSet(paste(refsDir,"/hairpin.fa",sep=""),use.names=TRUE)
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
  
  #value","variable")
  ureadlevel$sample<-str_split_fixed(ureadlevel$.id,'\\.',2)[,1]
  names(ureadlevel)<-c("id","rname","seq","spacedSeq","edits","tas","ntas","isos","status","value","variable")
  #recasted<-cast(ureadlevel[,-1],fun.aggregate=sum)
  #let's use spacedSeq
  recasted<-cast(ureadlevel[,c(-1,-3)],fun.aggregate=sum)
  notmirna<-DNAStringSet(scan("/nas/is1/leipzig/notmirna.txt",what=character()))
  tRNAs<-DNAStringSet(scan("/nas/is1/leipzig/tRNAs.txt",what=character()))

  recasted$query<-str_replace_all(recasted$spacedSeq,' ','')
  recasted$mismap<-FALSE
recasted$tRNA<-FALSE
  recasted$mismap[recasted$query %in% notmirna]<-TRUE
  recasted$tRNA[recasted$query %in% tRNAs]<-TRUE

  #find possible cross-mappings between miRNAs
  autoMismap<-recasted[,c('query','rname')]
#   possiblyAmbiguous<-autoMismap[autoMismap$query %in% which(table(autoMismap$query)>1),]
#   possiblyAmbiguous<-autoMismap[autoMismap$query %in%  autoMismap$query[which(table(autoMismap$query)>1)]
#   possiblyAmbiguous<-autoMismap[autoMismap$query %in%  autoMismap$query[which(table(autoMismap$query)>1)],]
#   possiblyAmbiguous<-autoMismap[autoMismap$query %in% which(table(autoMismap$query)>1),]
#   possiblyAmbiguous<-autoMismap[autoMismap$query %in%  autoMismap$query[which(table(autoMismap$query)>1)]
possiblyAmbiguous<-autoMismap[autoMismap$query %in% names(which((table(autoMismap$query)>1))),]

  possiblyAmbiguous$value<-1
  colnames(possiblyAmbiguous)<-c("query","variable","value")

  possiblyCast<-cast(possiblyAmbiguous,fun.aggregate=length)
  crossmirreport<-adply(possiblyCast,1,function(x){data.frame(query=x$query,hits=paste(colnames(possiblyCast)[which(x==1)],collapse=","))})[,c("query","hits")]
  recasted<-merge(recasted,crossmirreport,by="query",all.x=TRUE)



  save(recasted,file=concat(outputDir,"/recasted.RData"),compress=TRUE)
  #for rstudio display
  res<-as.data.frame(recasted)

  nomismaps<-subset(recasted,mismap==FALSE)
  nomismaps$rowsum<-rowSums(nomismaps[,samples])
  nomismaps<-subset(nomismaps,rowSums(nomismaps[,hiseqsamples]>0)>3)

  disp_df <- within(nomismaps, rm(query,status))
  disp_df$hits[is.na(disp_df$hits)]<-''
  #3 of the 6 matched samples must have this
  d_ply(disp_df,.(rname),.fun=function(x){
        if(sum(colSums(x[,hiseqsamples])>100)>3){
          save(x,file=concat(outputDir,"/readlevels/",x$rname[1],".RData"),compress=TRUE)
          names(x)<-str_replace_all(names(x),'RB','')
          excelTable<-x[mixedorder(x$edits),]
          names(excelTable)[2]<-hairpinLookup(x$rname[1])
          write.csv(excelTable,file=concat(outputDir,"csv/",x$rname[1],".csv",sep=""))
        }
  })
