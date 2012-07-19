#bamDirectory should be set elsewhere
fastqDirectory<-"fastq/"

refsDir<-"/nas/is1/leipzig/refs/"
hiseqDataDir<-"/nas/is1/leipzig/doHiseq/"
solexaDataDir<-"/nas/is1/leipzig/doSolexa/"
ionDataDir<-"/nas/is1/leipzig/doPGM/"

outputDir="/nas/is1/leipzig/results/"
hiseqconds<-c("N","T","N","T","N","T","T","T")
hiseqsamples<-c("RB494N","RB494T","RB495N","RB495T","RB498N","RB498T","RB517T","RB525T")
solexaconds<-c("T","T","T","T","N")
solexasamples<-c("FGC0031_s_8.WERI","FGC0036_s_1.Y79","FGC0036_s_2.WERI","FGC0036_s_2.Y79","FGC0042_s_1.normal")
solexaprettynames<-c("31s8_WERI","36s1_Y79","36s2_WERI","36s2_Y79","42s1_nrml")
ionconds<-c("T","T","T","T","T","T")
ionsamples<-c("RB525T01","RB525T02","RB525T03","WERI01","WERI02","WERI03")
ionprettynames<-c("RB525T01_PGM","RB525T02_PGM","RB525T03_PGM","WERI01_PGM","WERI02_PGM","WERI03_PGM")

#replicate informations
metadata<-list()
nonreps<-c("42s1_nrml","RB494N","RB494T","RB495N","RB495T","RB498N","RB498T","RB517T")
metadata$replicates<-as.list(nonreps)
names(metadata$replicates)<-nonreps
metadata$replicates$RB525T<-c('RB525T','RB525T01_PGM','RB525T02_PGM','RB525T03_PGM')
metadata$replicates$WERI<-c('31s8_WERI','36s2_WERI','WERI01_PGM','WERI02_PGM','WERI03_PGM')
metadata$replicates$Y79<-c('36s1_Y79','36s2_Y79')
#42s1_nrml  RB494N RB494T RB495N  RB495T RB498N RB498T RB517T  RB525T    WERI    Y79 (11 groups)
metadata$conds<-c('N','T','N','T','N','T','T','T')
metadata$individual<-c(1,1,2,2,3,3,4,5)
metadata$replicate<-c(1,2,1,2,1,2,2,2)
metadata$type<-c('HS','HS','HS','HS','HS','HS','HS','HS')
metadata$paired<-c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,FALSE)
metadata$pdata<-data.frame(condition=metadata$conds,replicate=metadata$replicate,type=metadata$type,individual=metadata$individual,paired=metadata$paired)
rownames(metadata$pdata)<-c('RB494N','RB494T','RB495N','RB495T','RB498N','RB498T','RB517T','RB525T')

bamPaths<-c(concat(hiseqDataDir,bamDirectory,hiseqsamples,'.bam'))

readCountPaths<-c(concat(hiseqDataDir,fastqDirectory,hiseqsamples,'.cnt'))
readCounts<-sapply(readCountPaths,scan)

samples<-c(hiseqsamples)
bamconds<-c(hiseqconds)


names(bamPaths)<-samples
bamSamples<-DataFrame(conds=bamconds,row.names=samples,reads=readCounts)