
# test-ridge

#library("testthat")
#Part I: check errors
#================
lambda1=c(1,1,3,4);
lambda2=c(-1,1,3,4);
y1=c(1,2,3)
y2=c(1,2,3,4)
X=matrix(data=c(1:6),nrow=3,ncol=2)

test_that("whether error",
          { expect_error(ridge_regression(y2,X,lambda1),"not conformable");
            expect_error(ridge_regression(y1,X,lambda2),"negative tuning parameter exists") }
          )

#Part II:check correctness of the estimated regression coefficients.
#======================

beta_est=ridge_regression(y1,X,lambda1)
test_that("whether equal",
          {for(i in 1:length(lambda1))
            expect_equal((t(X)%*%X+lambda1[i]*diag(1,ncol(X)))%*%beta_est[,i],t(X)%*%y1)
          }
          )
