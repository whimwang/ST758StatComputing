#' Sweep k
#'
#' \code{sweep_k} applies the sweep operator to a symmetric matrix  on the kth diagonal entry if it is possible.
#'
#' @param A The input symmetric matrix
#' @param k Diagonal index on which to sweep
#' @export
#' 
sweep_k<-function(A,k){
    if(all.equal(A,t(A))!=TRUE|A[k,k]==0)
        return("ERROR");
    n=ncol(A);
    A_sweep_k=A;
    # calculate upper triangle
   for(i in 1:n)
     for(j in i:n)
     {
       if(i!=k&j!=k)  {A_sweep_k[i,j]=A[i,j]-A[i,k]*A[k,j]/A[k,k];#print(A_sweep_k);
                         next;}
       if(i==k&j!=k)  {A_sweep_k[i,j]=A[i,j]/A[k,k];#print(A_sweep_k);
                         next;}
       if(i!=k&j==k)  {A_sweep_k[i,j]=A[i,j]/A[k,k];#print(A_sweep_k);
                         next;}
       A_sweep_k[k,k]=-1/A[k,k];
     }
    
    # get lower triangle based on symmetry
    for(i in 2:n)
      for(j in 1:(i-1))
          A_sweep_k[i,j]=A_sweep_k[j,i];
          
    return(A_sweep_k);
}
 
#A=matrix(c(1,2,3,2,1,4,3,4,5),nrow=3,ncol=3,byrow = TRUE)
#A1=sweep_k(A,1)
#A2=sweep_k(A1,2)
#A3=sweep_k(A2,3)





