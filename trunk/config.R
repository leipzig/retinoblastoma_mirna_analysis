mysys<-system2("uname",stdout=TRUE)
if(mysys=='Darwin'){
gangulyRoot<-"/Users/leipzig/Documents/ganguly"
}else{
    if(mysys=="Linux")
    {
    gangulyRoot<-"/storage/Leipzig/ganguly"
}else{
    stop("I can't figure out which machine you are on. Edit config.R accordingly")
}
}
