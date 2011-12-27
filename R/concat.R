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
