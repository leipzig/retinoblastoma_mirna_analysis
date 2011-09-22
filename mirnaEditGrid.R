library(reshape)
library(mirbase.db)
library(stringr)
library(IRanges)

x <- mirbaseSEQUENCE
mirnaSEQ <- mget(mappedkeys(x), x)

x <- mirbaseMATURE
mirnaMature<-mget(mappedkeys(x),x)


setwd("~/Documents/gangulyRBhi")

#1) forward ref alleles; 2) reverse ref; 3) forward non-ref; 4) reverse non-ref alleles
normal<-read.table("normal.hairpin.col",header=FALSE,col.names=c("miRNA","pos","ref","var","refnorm","revrefnorm","varnorm","revvarnorm"))
tumor<-read.table("tumor.hairpin.col",header=FALSE,col.names=c("miRNA","pos","ref","var","reftumor","revreftumor","vartumor","revvartumor"))


#questionable
normal$var<-gsub(",.","",normal$var)
tumor$var<-gsub(",.","",tumor$var)

all<-merge(normal,tumor,by=c("miRNA","pos","ref","var"),all=TRUE)
all[is.na(all)]<-0

all<-all[order(all$miRNA,all$pos),-grep('rev',names(all))]

getPval<-function(x){
  mirnaMat<-matrix(as.numeric(c(x$refnorm,x$varnorm,x$reftumor,x$vartumor)),nrow=2,dimnames=list(c("ref","var"),c("norm","tumor")))
  fisher.test(mirnaMat,alternative="two.sided")$p.val
}

res<-ddply(.data=all,.(miRNA,pos,ref,var,refnorm,varnorm,reftumor,vartumor),transform,fisherTwoSidedpval=fisher.test(matrix(c(refnorm,varnorm,reftumor,vartumor),nrow=2,dimnames=list(c("ref","var"),c("norm","tumor"))))$p.val)
res$seq<-as.character(mirnaSEQ[as.character(res$miRNA)])

matOnHair<-function(miRNA,var,pos,seq){
  myrange<-IRanges(start=matureFrom(mirnaMature[miRNA][[1]]),end=matureTo(mirnaMature[miRNA][[1]]),names=matureName(mirnaMature[miRNA][[1]]))
  dashes<-rep(" ",max(str_length(as.character(seq)),0))
  dashes[as.integer(myrange[str_sub(names(myrange),-1)=='*'])]<-'*'
  dashes[as.integer(myrange[str_sub(names(myrange),-1)!='*'])]<-'='
  dashstr<-paste(dashes,sep="",collapse="")
  #mutspace<-rep(" ",pos)
  #mutstr<-paste(mutspace,var,sep="",collapse="")
  mutstr<-str_pad(var,pos+str_length(var)-1)
  paste(dashstr,seq,mutstr,sep="\n")
}
res<-ddply(.data=res,c("miRNA","pos","ref","var"),transform,hairpin=matOnHair(as.character(miRNA),var,pos,seq))
res<-res[,-which(names(res)=='seq')]
write.csv(res,file="hairpin.edits.fisher.csv")


