#' Sweep
#'
#' @param A The input symmetric matrix
#' @param k Diagonal index entry set on which to sweep 
#' @export

sweep<-function(A,k=NULL){
      sweep_wrap=A;
  for(i in k){
     if(!is.character(sweep_wrap))
       {sweep_wrap=sweep_k(sweep_wrap,i)
       #print(sweep_wrap);
       }
     else 
       return("ERROR")
  }
  return(sweep_wrap);
}

#A=matrix(c(1,2,3,2,1,4,3,4,5),nrow=3,ncol=3,byrow = TRUE)
#A1=sweep_k(A,1)
#A2=sweep_k(A1,2)
#A3=sweep_k(A2,3)
