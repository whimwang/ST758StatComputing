#' Inverse Sweep k
#'
#' \code{isweep_k} applies the inverse sweep operator to a symmetric matrix  on the kth diagonal entry if it is possible.
#'
#' @param A The input symmetric matrix
#' @param k Diagonal index on which to sweep
#' @export


isweep_k<-function(A,k){
  if(all.equal(A,t(A))!=TRUE|A[k,k]==0)
    return("ERROR");
  n=ncol(A);
  A_isweep_k=A;
  # calculate upper triangle
  for(i in 1:n)
    for(j in i:n)
    {
      if(i!=k&j!=k)  {A_isweep_k[i,j]=A[i,j]-A[i,k]*A[k,j]/A[k,k];#print(A_isweep_k);
                       next;}
      if(i==k&j!=k)  {A_isweep_k[i,j]=-A[i,j]/A[k,k];#print(A_isweep_k);
                       next;}
      if(i!=k&j==k)  {A_isweep_k[i,j]=-A[i,j]/A[k,k];#print(A_isweep_k);
                       next;}
      A_isweep_k[k,k]=-1/A[k,k];
    }
  
  # get the lower triangle based on symmetry
  for(i in 2:n)
    for(j in 1:(i-1))
  A_isweep_k[i,j]=A_isweep_k[j,i];
    
  return(A_isweep_k);
}


#A=matrix(c(1,2,3,2,1,4,3,4,5),nrow=3,ncol=3,byrow = TRUE)
#A1=sweep_k(A,1)
#A2=sweep_k(A1,2)
#A3=sweep_k(A2,3)

#B=matrix(c(1,3,3,1),nrow=2,ncol=2,byrow=TRUE)
#B1=sweep_k(B,1)
#B2=sweep_k(B1,2)
