#' Banded Cholesky
#'
#' \code{chol_banded} computes the Cholesky decomposition of a banded
#' positive definite matrix.
#'
#' @param A The input matrix
#' @param k bandwidth parameter
#' @export
#'
chol_banded <- function(A, k) {
  if(!isTRUE(all.equal(A,Matrix::t(A)))) stop("not symmetric matrix"); # check whether symmetric
  if(!check_k_banded(A,k)) stop("not k-banded");

  n=nrow(A);
  if(k==0) return(diag(sqrt(diag(A)),nrow=n,ncol=n));

  R=matrix(numeric(n^2),nrow=n,ncol=n);# upper triangle
  A_current=A

for(i in 1:(n-k)){
R[i,i]=sqrt(A_current[1,1])
R[i,(i+1):(i+k)]=A_current[1,2:(k+1)]/R[i,i]#i+k<=n
RR_star=A_current[-1,-1] #nrow of A_current>=2
if(k>=2) { RR_star[1:k,1:k]=A_current[2:(k+1),2:(k+1)]-R[i,(i+1):(i+k)]%*%t(R[i,(i+1):(i+k)])} #k>=1 i+k<=n
if(k==1) { RR_star=A_current[2:(k+1),2:(k+1)]-R[i,(i+1):(i+k)]%*%t(R[i,(i+1):(i+k)])}
A_current=RR_star
}

for(i in (n-k+1):(n-1)){
  R[i,i]=sqrt(A_current[1,1])
  R[i,(i+1):n]=A_current[1,2:(n-i+1)]/R[i,i]
  RR_star=A_current[-1,-1]-R[i,(i+1):n]%*%t(R[i,(i+1):n])#nrow of A_current>=2
  A_current=RR_star
}
i=i+1
if(i==n)# last step
  R[i,i]=sqrt(RR_star)
return(R)
}









