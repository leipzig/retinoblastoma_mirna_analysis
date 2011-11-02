library("Rsamtools")
library(reshape)
library(mirbase.db)
library(stringr)
library(IRanges)
library(plyr)

library("doMC")
registerDoMC()

x <- mirbaseSEQUENCE
mirnaSEQ <- mget(mappedkeys(x), x)

x <- mirbaseMATURE
mirnaMature<-mget(mappedkeys(x),x)
#res$seq<-as.character(mirnaSEQ[as.character(res$miRNA)])

dataDir<-"/nas/is1/leipzig/Ganguly/gangulyRBhi/data/"
bamDirectory<-"bam/"
outputDir="/nas/is1/leipzig/Ganguly/gangulyRBhi/results/schugView/"

conds<-c("N","T","N","T","N","T")
samples<-c("RB494N","RB494T","RB495N","RB495T","RB498N","RB498T")
#conds<-c("N","T")
#samples<-c("RB494N","RB494T")
bamPaths<-paste(dataDir,bamDirectory,samples,'.hairpin.sam.sorted.bam',sep="")
bamSamples<-DataFrame(conds=conds,row.names=samples)



s <- read.DNAStringSet(paste(dataDir,"refs/hsa_hairpin.dna.fa",sep=""),use.names=TRUE)
hsa<-s[str_detect(names(s),"^hsa")]

names(hsa)<-str_match(names(hsa),"^\\S+")
hsaLookup<-as.character(hsa)
hsaGR<-GRanges(names(hsa),IRanges(1,width(hsa)),strand="*")


bamView<-BamViews(bamPaths=bamPaths,
bamSamples=bamSamples,
bamRanges=hsaGR)


what <- c("rname", "strand", "pos", "qwidth", "seq")

param<-ScanBamParam(what=what,which=hsaGR)

#bams is a [[list of samples[[list of genes]]]
bams <- scanBam(bamView, param=param)
#bams <- readBamGappedAlignments(bamView)


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

editColumn<-function(targetseq,queryseq){

  #submir<-str_sub(mirna,posOnMir,posOnMir+str_length(seq)-1)
  #pwa<-pairwiseAlignment(submir,seq,type="local-global")
  pwa<-pairwiseAlignment(targetseq,queryseq,type="local-global")
  
  #nucs<-strsplit(seq,'')[[1]][unlist(pwa@pattern@mismatch)]
  #poss<-unlist(pwa@subject@mismatch)+posOnMir-1
  
  #used for target positions
  target_poss_sub<-unlist(pwa@pattern@mismatch)
  target_poss_indel<-concat(as.character(start(unlist(pwa@pattern@indel))),'i',sep='')
  target_poss<-c(target_poss_sub,target_poss_indel)
  
  #used for identifying the edits
  seq_poss_sub<-unlist(pwa@subject@mismatch)
  seq_poss_indel<-start(unlist(pwa@subject@indel))
  seq_poss<-c(seq_poss_sub,seq_poss_indel)
  nucs<-strsplit(queryseq,'')[[1]][seq_poss]

  #formatEdits(nucs,target_poss)
  list(nucs=nucs,poss=target_poss)
}
formatEdits<-function(editList){
  edits<-vector()
  for(j in 1:length(editList$poss)){
   edits<-c(edits,concat(editList$poss[j],editList$nucs[j],sep=':'))
  }
  paste(edits,collapse=" ")
}

schug<-function(x){
   x$seq<-as.character(x$seq)
   if(length(x$seq)>0){
       df<-as.data.frame(x)
       cntdf<-ddply(df,.(rname,strand,seq,pos,qwidth),c("nrow"))
       udf<-ddply(cntdf,.(rname,strand,seq,pos,qwidth,nrow)
                  ,transform
                  ,spacedSeq=paste(str_pad(seq,pos+qwidth-1))
                  ,mirna=hsaLookup[as.character(rname)]
                ,edits=formatEdits(editColumn(hsaLookup[as.character(as.character(rname))],as.character(seq)))
                  )
        #udf$sample<-str_split_fixed(udf$.id,'\\.',2)[,1]
      udf[order(udf$edits,-udf$nrow),c("rname","spacedSeq","nrow","edits")]
  }
}
ubams<-unlist(bams)
#ubamsjr<-ubams[1:2]
uschug<-ldply(ubams,.fun=schug,.parallel=TRUE)
uschug$sample<-str_split_fixed(uschug$.id,'\\.',2)[,1]

recasted<-recast(uschug[,-1],rname+spacedSeq+edits~sample,sum)

#get sample name, i.e RB494N
#
#allbams<-ldply(ubams,.fun=function(x){ldply(x[1],.fun=schug,.parallel=TRUE)})
 
#allbams$mirna<-hsaLookup[as.character(allbams$rname)]

#allbams$edits<-apply(allbams,1,function(x){stringDist(c(x$V1,str_sub(x$mirna,x$pos,x$pos+x$qwidth)))})

#allbams$submir<-str_sub(allbams$mirna,allbams$pos,allbams$pos+allbams$qwidth-1)
#editDist<-function(query,submir){as.numeric(stringDist(c(query,submir)))}
 
#ddply(allbams[1:10,],.(rname,strand,seq,pos,qwidth,nrow,mirna),transform,edits=editColumn(as.character(seq),as.character(mirna),pos))