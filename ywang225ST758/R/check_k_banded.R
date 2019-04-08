#' Check Matrix is k-banded
#' 
#' \code{check_k_banded} returns a Boolean variable indicating whether or
#' not a matrix is k-banded.
#' 
#' @param A The input matrix
#' @param k bandwidth parameter
#' @export
#' 
check_k_banded <- function(A, k) {
  # judge whether error
  if((as.integer(k)!=k)|k<0) {stop("ERROR")};
  # judge whether k band
  
  n=nrow(A);
  if(k>n-1) stop("Out of Range!")
  

  #check out of k bands
  for(i in 1:n)
  { for(j in 1:n)
      { if(abs(i-j)>k&&!(A[i,j]==0))
        return(FALSE) }
  }
  # check the k th band
  for(i in 1:(n-k))
    {if(A[i,i+k]!=0|A[i+k,i]!=0) 
      return(TRUE)}
  return(FALSE)
}







