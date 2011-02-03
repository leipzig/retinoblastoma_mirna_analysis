switch(Sys.info()[['sysname']],
Linux  = {gangulyRoot<-"/storage/Leipzig/ganguly"},
Darwin = {gangulyRoot<-"/Users/leipzig/Documents/ganguly"},
NULL = {stop("I can't figure out which machine you are on. Edit config.R accordingly")})
