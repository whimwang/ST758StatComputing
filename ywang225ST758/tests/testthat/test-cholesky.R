
# test check_k_banded
#=======================================


A=matrix(data=c(1,4,0,0,9,2,5,0,0,7,3,6,0,0,8,4),nrow=4,ncol=4,byrow=TRUE)
B=matrix(data=c(1,2,3,0,0,2,4,5,6,0,3,5,7,8,9,0,6,8,10,11,0,0,9,11,12),nrow=5,ncol = 5,byrow=TRUE)
library("testthat")
test_that("whether error",
          { expect_equal(check_k_banded(A,0),FALSE) # not "FALSE" the former is logical the latter is character
            expect_error(check_k_banded(A,0.5),"ERROR")# not Error:ERROR
            expect_error(check_k_banded(A,-1),"ERROR")
           expect_equal(check_k_banded(B,0),FALSE)
            expect_equal(check_k_banded(B,0),FALSE)
            expect_equal(check_k_banded(B,2),TRUE)
         })


# test chol_banded
#==========================================

temp_R=matrix(data=c(2,2,3,rep(0,2),
                     0,1:3,0,
                     0,0,1:3,
                     rep(0,3),1:2,
                     rep(0,4),1),
              nrow=5,ncol=5,byrow=TRUE)
temp_A=t(temp_R)%*%temp_R+2*diag(rep(1,5),nrow=5,ncol =5)




temp_R=matrix(data=c(1:4,0,
                     0,rep(1,4),
                     rep(0,2),rep(1,3),
                     rep(0,3),rep(1,2),
                     rep(0,4),1),nrow = 5,ncol=5,byrow = TRUE)

temp_B=t(temp_R)%*%temp_R+diag(rep(1,5),nrow=5,ncol=5)


temp_R=matrix(data=c(rep(1,5),rep(0,2),
                     0,rep(1,5),0,
                     rep(0,2),rep(1,5),
                     rep(0,3),rep(1,4),
                     rep(0,4),rep(1,3),
                     rep(0,5),rep(1,2),
                     rep(0,6),1),nrow=7,ncol=7,byrow = TRUE)
temp_C=t(temp_R)%*%temp_R+3*diag(rep(1,7),nrow=7,ncol=7)


test_that("whether chol_banded right",
          { expect_error(chol_banded(temp_A,1),"not k-banded")
            expect_equal(t(chol_banded(temp_A,2))%*%(chol_banded(temp_A,2)),temp_A)
            expect_equal(chol_banded(temp_A,2),chol(temp_A))

            expect_error(chol_banded(temp_B,2),"not k-banded")
            expect_equal(t(chol_banded(temp_B,3))%*%(chol_banded(temp_B,3)),temp_B)
            expect_equal(chol_banded(temp_B,3),chol(temp_B))

            expect_error(chol_banded(temp_R,4),"not symmetric matrix")
            expect_equal(t(chol_banded(temp_C,4))%*%(chol_banded(temp_C,4)),temp_C)
            expect_equal(chol_banded(temp_C,4),chol(temp_C))
          })



