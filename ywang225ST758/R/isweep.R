#' Inverse Sweep
#'
#' @param A The input symmetric matrix
#' @param k Diagonal index entry set on which to sweep 
#' @export

isweep<-function(A,k=NULL){
       isweep_wrap=A;
    for(i in k){
      if(!is.character(isweep_wrap))
       isweep_wrap=isweep_k(isweep_wrap,i)
      else
       return("ERROR")
    }
    return(isweep_wrap);
}


  
  
