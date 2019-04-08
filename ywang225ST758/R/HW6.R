#' Softthreshold operator
#' 
#' @param x input vector
#' @param lambda threshold
#' @export
#' @examples
#' softthreshold(seq(-10, 10, length.out=10), 1)
softthreshold <- function(x, lambda) {
  temp=abs(x)-lambda;
  return(ifelse(temp>0,temp,-temp)*sign(x));
}
#' Proximal mapping of the scaled nuclear norm
#' 
#' @param X input matrix
#' @param gamma scale parameter
#' @export
#' @examples
#' set.seed(12345)
#' n <- 10; p <- 20
#' X <- matrix(rnorm(n*p), n, p)
#' Y <- prox_nuc(X, 1)
prox_nuc <- function(X, gamma) {
  decomposition=svd(X);
  new_d=softthreshold(decomposition$d,gamma);
  X_map=decomposition$u%*%diag(new_d)%*%t(decomposition$v)
  return(X_map);
}
#' Project matrix onto observed index set
#' 
#' @param X input matrix
#' @param Omega index set
#' @export
#' @examples
#' set.seed(12345)
#' n <- 10; p <- 20
#' X <- matrix(rnorm(n*p), n, p)
#' numel <- n*p
#' Omega <- sample(1:numel, floor(0.5*numel), replace=FALSE)
#' PX <- project_omega(X, Omega)
project_omega <- function(X, Omega) {
  p=ncol(X);n=nrow(X);
  res=Omega%%n;int=Omega%/%n;
  row=ifelse(res==0,n,res);
  col=ifelse(res==0,int,int+1);
  Result=matrix(0,nrow(X),ncol(X));
  for(i in 1:length(Omega))
    Result[row[i],col[i]]=X[row[i],col[i]];
  return(Result);
}
#' Proximal-Gradient Step
#' 
#' @param X input matrix
#' @param U current guess
#' @param Omega index set
#' @param lambda regularization parameter
#' @param t step size
#' @export
#' @examples
#' set.seed(12345)
#' n <- 10; p <- 20
#' X <- matrix(rnorm(n*p), n, p)
#' numel <- n*p
#' Omega <- sample(1:numel, floor(0.5*numel), replace=FALSE)
#' U <- matrix(rnorm(n*p), n, p)
#' U[Omega] <- X[Omega]
#' lambda <- 1
#' t <- 1
#' U1 <- prox_step(X, U, Omega, lambda, t)
prox_step <- function(X, U, Omega, lambda, t) {
  #X_half=X-(project_omega(X,Omega)-project_omega(U,Omega))*t;
  X_half=U-t*project_omega(U-X,Omega);
  out=prox_nuc(X_half,lambda)
  return(out);
}
#' Matrix Completion Loss
#' 
#' @param X input matrix
#' @param U current guess
#' @param Omega index set
#' @param lambda regularization parameter
#' @export
#' @examples
#' set.seed(12345)
#' n <- 10; p <- 20
#' X <- matrix(rnorm(n*p), n, p)
#' numel <- n*p
#' Omega <- sample(1:numel, floor(0.5*numel), replace=FALSE)
#' U <- matrix(rnorm(n*p), n, p)
#' U[Omega] <- X[Omega]
#' lambda <- 1
#' matrix_completion_loss(X, U, Omega, lambda)
matrix_completion_loss <- function(X, U, Omega, lambda) {
  diff_matrix=project_omega(X,Omega)-project_omega(U,Omega);
  nuclear_norm=sum(svd(U)$d);
  return((norm(diff_matrix,type="F"))^2*0.5+lambda*nuclear_norm)
}
#' Proximal Gradient Descent
#' 
#' @param X input matrix
#' @param U current guess
#' @param Omega index set
#' @param lambda regularization parameter
#' @param t step-size parameter
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @export
#' @examples
#' set.seed(12345)
#' n <- 10; p <- 20
#' X <- matrix(rnorm(n*p), n, p)
#' numel <- n*p
#' Omega <- sample(1:numel, floor(0.5*numel), replace=FALSE)
#' U <- matrix(rnorm(n*p), n, p)
#' U[Omega] <- X[Omega]
#' lambda <- 1
#' pg_sol <- prox_grad(X, U, Omega, lambda,1)
prox_grad <- function(X, U, Omega, lambda, t, max_iter=1e2, tol=1e-3) {
  times=0;
  U_new=prox_step(X,U,Omega,lambda,t);
  obj_new=matrix_completion_loss(X,U_new,Omega,lambda);
  obj=obj_new;
  norm_X_F=norm(X,type="F");
  relative_change_function_value=tol+1;
  while(times<max_iter&(relative_change_function_value>tol)){#|(obj_new<obj)
    times=times+1;
    obj=obj_new;
    U_new=prox_step(X,U_new,Omega,lambda,t);
    obj_new=matrix_completion_loss(X,U_new,Omega,lambda);
    #print(paste0(times," loss: ",obj_new,"change",relative_change_function_value));#
    relative_change_function_value=obj/obj_new-1;
    relative_change_iterate_value=norm(U_new-X,type="F")^2/norm_X_F^2;
  }
  return(list(iterate_value=U_new,
              objective_function_value=matrix_completion_loss(X,U_new,Omega,lambda),
              relative_change_function_value=relative_change_function_value,
              relative_change_iterate_value=relative_change_iterate_value));
  
}
#' Soft-Impute Algorithm
#' 
#' @param X input matrix
#' @param U current guess
#' @param Omega index set
#' @param lambda regularization parameter
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @export
#' @examples
#' set.seed(12345)
#' n <- 10; p <- 20
#' X <- matrix(rnorm(n*p), n, p)
#' numel <- n*p
#' Omega <- sample(1:numel, floor(0.5*numel), replace=FALSE)
#' U <- matrix(rnorm(n*p), n, p)
#' U[Omega] <- X[Omega]
#' lambda <- 1
#' si_sol <- soft_impute(X, U, Omega, lambda)
soft_impute <- function(X, U, Omega, lambda, max_iter=1e2, tol=1e-3) {
  times=0;n=nrow(X);p=ncol(X);total_seq=1:(n*p);
  U_new=U;
  Omega_compliment=total_seq[-Omega];
  obj_new=matrix_completion_loss(X,U_new,Omega,lambda);
  norm_X_F=norm(X,"F");
  relative_change_function_value=1;
  relative_change_iterate_value=1;
  obj=obj_new;
  
  while(times<max_iter&(relative_change_function_value>tol|obj_new>obj)){
    times=times+1; 
    obj=obj_new;#renew
    Y_matrix=project_omega(X,Omega)+project_omega(U_new,Omega_compliment);
    decomposition=svd(Y_matrix);
    new_d=softthreshold(decomposition$d,lambda);
    U_new=decomposition$u%*%diag(new_d)%*%t(decomposition$v);
    obj_new=matrix_completion_loss(X,U_new,Omega,lambda);
    #loss=matrix_completion_loss(X,U_new,Omega,lambda);
    relative_change_function_value=obj/obj_new-1;
    relative_change_iterate_value=norm(U_new-X,type="F")^2/norm_X_F^2;
    #print(paste0(times," function ",relative_change_function_value," change_function ",relative_change_function_value));
  }
  return(list(iterate_value=U_new,
              objective_function_value=matrix_completion_loss(X,U_new,Omega,lambda),
              relative_change_function_value=relative_change_function_value,
              relative__change_iterate_value=relative_change_iterate_value));
}