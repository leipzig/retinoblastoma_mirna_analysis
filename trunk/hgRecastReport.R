library(reshape)
setwd("/storage/Ganguly/gangulyRBhi")

#1) forward ref alleles; 2) reverse ref; 3) forward non-ref; 4) reverse non-ref alleles
normal<-read.table("normal.hg19.col",header=FALSE,col.names=c("miRNA","chr","pos","ref","var","refnorm","revrefnorm","varnorm","revvarnorm"))
#normal$cell<-"normal"
tumor<-read.table("tumor.hg19.col",header=FALSE,col.names=c("miRNA","chr","pos","ref","var","reftumor","revreftumor","vartumor","revvartumor"))
#tumor$cell<-"tumor"
#all<-rbind(normal,tumor)

#questionable
normal$var<-gsub(",.","",normal$var)
tumor$var<-gsub(",.","",tumor$var)

all<-merge(normal,tumor,by=c("miRNA","chr","pos","ref","var"),all=TRUE)
all[is.na(all)]<-0

all<-all[order(all$miRNA,all$pos),-grep('rev',names(all))]

#melted<-melt(all,measure.vars="cell")
#recasted<-cast(melted, miRNA+pos+ref+var+refnorm+varnorm ~ value)
#write.csv(recasted,file="hsa-let.edits.csv")

getPval<-function(x){
  #mirnaMat<-matrix(as.numeric(x[c("refnorm","varnorm","reftumor","vartumor")]),nrow=2,dimnames=list(c("ref","var"),c("norm","tumor")))
  mirnaMat<-matrix(as.numeric(c(x$refnorm,x$varnorm,x$reftumor,x$vartumor)),nrow=2,dimnames=list(c("ref","var"),c("norm","tumor")))
  fisher.test(mirnaMat,alternative="two.sided")$p.val
}

#res<-#ddply(.data=all,.(miRNA,pos,ref,var,refnorm,varnorm,reftumor,vartumor),summarize,foo=getPval)

res<-ddply(.data=all,.(miRNA,chr,pos,ref,var,refnorm,varnorm,reftumor,vartumor),transform,fisherTwoSidedpval=fisher.test(matrix(c(refnorm,varnorm,reftumor,vartumor),nrow=2,dimnames=list(c("ref","var"),c("norm","tumor"))))$p.val)
#ddply(.data=all,.(miRNA,pos,ref,var,refnorm,varnorm,reftumor,vartumor),transform,fisherTwoSidedpval=fisher.test(matrix(c(1,2,3,4),nrow=2,dimnames=list(c("ref","var"),c("norm","tumor")))),alternative="two.sided")$p.val)

print(res)
write.csv(res,file="hg19.edits.fisher.csv")