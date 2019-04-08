#' Gradient Descent (Backtracking Step-Size)
#' 
#' @param y binary response
#' @param X design matrix
#' @param b0 initial parameter estimate
#' @param delta shape parameter
#' @param lambda regularization parameter
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @export
#' 
gradient_descent_backtrack <- function(y, X, b0, delta, lambda=1e-4, max_iter=1e2, tol=1e-3) { 
beta=0.9;alpha=1;r=0.5;
times=0;

x_0=beta0;
x_1=beta0;
alpha_current=alpha;
while(times<max_iter){
  times=times+1;
  print(times);
  alpha_current=alpha;
  x_1_grad=gradf_hsvm(y,X,x_1,delta,lambda)
  while(fx_hsvm(y,X,x_1-alpha_current*x_1_grad,delta,lambda)>=fx_hsvm(y,X,x_1,delta,lambda)-r*alpha_current*sum(x_1_grad^2))
  {alpha_current=alpha_current*beta;
    print(paste("alpha_current",alpha_current,times));}
  
  x_0=x_1;
  x_1=x_1-alpha_current*x_1_grad;
  
  if(abs(fx_hsvm(y,X,x_1,delta,lambda)-fx_hsvm(y,X,x_0,delta,lambda))<tol)
    return( list(iterate_value=x_1,
                 objective_function_value=fx_hsvm(y, X, beta=x_1, delta, lambda),
                 norm2_gradient_value=sqrt(sum(x_1^2)),
                 relative_change_in_function=abs((fx_hsvm(y, X, beta=x_1, delta, lambda)-fx_hsvm(y, X, beta=x_0, delta, lambda))/fx_hsvm(y,X,beta=x_0,delta,lambda)),
                 relative_change_in_iterate=sqrt(sum((x_1-x_0)^2))/sqrt(sum(x_0^2)) ))  
  
}
return( list(iterate_value=x_1,
             objective_function_value=fx_hsvm(y, X, beta=x_1, delta, lambda),
             norm2_gradient_value=sqrt(sum(x_1^2)),
             relative_change_in_function=abs((fx_hsvm(y, X, beta=x_1, delta, lambda)-fx_hsvm(y, X, beta=x_0, delta, lambda))/fx_hsvm(y,X,beta=x_0,delta,lambda)),
             relative_change_in_iterate=sqrt(sum((x_1-x_0)^2))/sqrt(sum(x_0^2)) ))  
}









#===================================================
gradient_descent_backtrack <- function(y, X, b0, delta, lambda=1e-4, max_iter=1e2, tol=1e-3) {
  print("backtrack");
  beta=0.8;alpha=1;r=0.5;#chose
  times=0;
  alpha_current=alpha;
  x_0=b0;
  x_1=x_0-alpha_current*gradf_hsvm(y, X, x_0, delta);
  while(times<max_iter){
    times=times+1;
    x_1_grad=gradf_hsvm(y, X, x_1, delta)
    x_1_fx=fx_hsvm(y, X, x_1, delta)
    x_0_fx=fx_hsvm(y, X, x_0, delta)
    x_1_grad_norm2=norm(x_1_grad,type="2")
    if(x_1_grad_norm2<tol)# solution
      return(list(iterate_value=x_1,
                  objective_function_value=x_1_fx,
                  norm2_gradient_value=norm(x_1,type="2"),
                  relative_change_in_function=abs((x_1_fx-x_0_fx)/fx_hsvm(x_0)),
                  relative_change_in_iterate=norm(x_1-x_0,type="2")/norm(x_0,type="2") ))
    else{#not solution
    while( (fx_hsvm(y, X, x_1-alpha_current*x_1_grad, delta)>=
          fx_hsvm(y, X, x_1, delta)-1/2*alpha_current*x_1_grad_norm2^2))
        alpha_current=alpha_current*beta; # get the right step
    
    x_0=x_1;
    x_1=x_1-alpha_current*x_1_grad;
    alpha_current=alpha;
    
    }}
    if(times==max_iter)
      return("not converge")
}







