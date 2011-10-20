library("Rsamtools")
library(reshape)
library(mirbase.db)
library(stringr)
library(IRanges)

x <- mirbaseSEQUENCE
mirnaSEQ <- mget(mappedkeys(x), x)

x <- mirbaseMATURE
mirnaMature<-mget(mappedkeys(x),x)
res$seq<-as.character(mirnaSEQ[as.character(res$miRNA)])


s <- read.DNAStringSet("hsa_hairpin.dna.fa",use.names=TRUE)

names(hsa)<-str_match(names(hsa),"^\\S+")

hsaGR<-GRanges(names(hsa),IRanges(1,width(hsa)),strand="*")

hsa<-s[str_detect(names(s),"^hsa")]

what <- c("rname", "strand", "pos", "qwidth", "seq")

param<-ScanBamParam(what=what,which=hsaGR)
bams <- scanBam(bamFile, param=param)
  
schug<-function(x){
 x$seq<-as.character(x$seq)
 if(length(x$seq)>0){
 df<-as.data.frame(x)
 cntdf<-ddply(df,.(rname,strand,seq,pos,qwidth),c("nrow"))
 udf<-ddply(cntdf,.(rname,strand,seq,pos,qwidth,nrow),.fun=function(x){paste(str_pad(x$seq,x$pos+x$qwidth-1))})
 udf[order(udf$pos,-udf$nrow),]
}
}

allbams<-ldply(bams,.fun=schug)
