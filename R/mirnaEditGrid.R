library(reshape)
library(mirbase.db)
library(stringr)
library(IRanges)

x <- mirbaseSEQUENCE
mirnaSEQ <- mget(mappedkeys(x), x)

x <- mirbaseMATURE
mirnaMature<-mget(mappedkeys(x),x)

workDir<-"/nas/is1/leipzig/Ganguly/gangulyRBhi/"
setwd("/nas/is1/leipzig/Ganguly/gangulyRBhi/")
variationDir<-"results/variation/"

#1) forward ref alleles; 2) reverse ref; 3) forward non-ref; 4) reverse non-ref alleles
normal<-read.table(paste(workDir,variationDir,"normal.hairpin.col",sep=""),header=FALSE,col.names=c("miRNA","pos","ref","var","refnorm","revrefnorm","varnorm","revvarnorm"))
tumor<-read.table(paste(workDir,variationDir,"tumor.hairpin.col",sep=""),header=FALSE,col.names=c("miRNA","pos","ref","var","reftumor","revreftumor","vartumor","revvartumor"))

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





matOnHair<-function(miRNA,var,pos,seq){
  spacechar<-' '
  newline<-"<br>"
  myrange<-IRanges(start=matureFrom(mirnaMature[miRNA][[1]]),end=matureTo(mirnaMature[miRNA][[1]]),names=matureName(mirnaMature[miRNA][[1]]))
  dashes<-rep(spacechar,max(str_length(as.character(seq)),0))
  dashes[as.integer(myrange[str_sub(names(myrange),-1)=='*'])]<-'*'
  dashes[as.integer(myrange[str_sub(names(myrange),-1)!='*'])]<-'='
  dashstr<-paste(dashes,sep="",collapse="")
  #mutspace<-rep(" ",pos)
  #mutstr<-paste(mutspace,var,sep="",collapse="")
  mutstr<-str_pad(var,pos+str_length(var)-1,pad=spacechar)
  paste("<pre>",dashstr,seq,mutstr,"</pre>",sep=newline)
}
res<-ddply(.data=all,.(miRNA,pos,ref,var,refnorm,varnorm,reftumor,vartumor),transform,fisherTwoSidedpval=fisher.test(matrix(c(refnorm,varnorm,reftumor,vartumor),nrow=2,dimnames=list(c("ref","var"),c("norm","tumor"))))$p.val)
res$seq<-as.character(mirnaSEQ[as.character(res$miRNA)])

rest<-ddply(.data=res,c("miRNA","pos","ref","var"),transform,hairpin=matOnHair(as.character(miRNA),var,pos,seq))

rest<-rest[,-which(names(rest)=='seq')]
sortable.html.table(rest,"res.html","/home/leipzig/public_html","miRNA edits")



write.csv(rest,file="hairpin.edits.fisher.csv")


