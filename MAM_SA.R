###############################################################################################################
# Example of Multi-Race Model - Mating System
# Developers:[
# J. Augusto II, 
# Cherlynn Silva <cdnsprado@gmail.com https://github.com/CDsilvaA>
# ]
# Source: Mauricio Elzo
#
# Documentation
# Data provided by the user:
# X: incidence matrix
# r: a matrix of variances and residual covariances 
# y: vector of the response variable
# P: pedigree matrix (who is whose father)
# da: matrix of variances and additive genetic covariances
# dn: matrix of variances and non-additive genetic covariances
# KA: matrix of coefficients of additive genetic predications
# KN: matrix of coefficients of non-additive genetic predications
# KT: matrix of coefficients of total genetic prediction
#
#
# Data calculated by the script:
# Ga: inverse matrix of variances and additive genetic covariances
# Gn: inverse matrix of variances and nonadditive genetic covariances
# lhs: matrix left hand
# rhs: right hand matrix
# y_hat: vector of the solutions
# se_y_hat: standard solutions mistakes
# AMBV: additive genetic prediction solutions (deviations of breed B from breed A)
# SEP_AMBV: standard errors of the matrix of coefficients of additive genetic predications
# NMBV: non-additive genetic prediction solutions (assuming males are paired with females 1/2A1/2B and vice versa)
# SEP_NMBV: standard errors of the matrix of coefficients of nonadditive genetic predications
# TMBV: total genetic prediction solutions (sum of AMBV and NMBV)
# SEP_TMBV: standard errors of the matrix of coefficients of total genetic predications
#
#
# Necessary packages (calculation of the generalized inverse)

library(MASS)
#
#
###############################################################################################################
X<-matrix(c(1,	1.0,	0.0,	0,	1,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
            1,	0.0,	1.0,	0,	0,	1,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
            1,	0.5,	0.5,	1,	0,	1,	0,	0,	1,	0,	0,	0,	0,	1,	0,	0,	0,	0,
            1,	0.5,	0.5,	1,	0,	1,	0,	0,	0,	1,	0,	0,	1,	0,	0,	0,	0,	0,
            1,	0.5,	0.5,	1,	1,	0,	0,	0,	0,	0,	1,	0,	1,	1,	0,	0,	0,	0,
            1,	0.75,	0.25,	0.5,	1,	0,	0,	0,	0,	0,	0,	1,	0.5,	0,	0.5,	0,	0,	0), byrow=T, nrow=6)
colnames(X)<-c("mu","A","B","h","M","F","a_1","a_2","a_3","a_4","a_5","a_6",
               "na_1","na_2","na_3","na_4","na_5","na_6")
X

y<-c(289,245,256,261,292,286)
r<-diag(c(49,16,32.5,32.5,32.5,47))
da<-diag(c(36,9,20.5,13.5,11.25,15.625))
dn<-diag(c(16,16,12,12,8,8))

P<-matrix(c(0,0,0,0,0,0,
            0,0,0,0,0,0,
            0,1,0,0,0,0,
            1,0,0,0,0,0,
            1,1,0,0,0,0,
            1,0,1,0,0,0), byrow=T, ncol=6)
P

Ga<-solve(da)-0.5*(solve(da)%*%P)-0.5*(t(P)%*%solve(da))+0.25*(t(P)%*%solve(da)%*%P)
Ga

Gn<-solve(dn)-0.5*(solve(dn)%*%P)-0.5*(t(P)%*%solve(dn))+0.25*(t(P)%*%solve(dn)%*%P)
Gn


lhs<-(t(X)%*%solve(r)%*%X)
lhs[7:12,7:12]<-(t(X)%*%solve(r)%*%X)[7:12,7:12]+(Ga)
lhs[13:18,13:18]<-(t(X)%*%solve(r)%*%X)[13:18,13:18]+(Gn)
round(lhs,3)


rhs<-t(X)%*%solve(r)%*%y
round(rhs,2)


y_hat<-ginv(lhs)%*%rhs
se_y_hat<-sqrt(diag(ginv(lhs)))
solucoes<-cbind(round(y_hat,2), round(se_y_hat,2))
rownames(solucoes)<-c("mu","A","B","h","M","F","a_1","a_2","a_3","a_4","a_5","a_6",
                   "na_1","na_2","na_3","na_4","na_5","na_6")
colnames(solucoes)<-c("SOL","SESOL")
solucoes



KA<-matrix(c(0,0,0,0,0,0,
             1,0,0.5,0.5,0.5,0.75,
             -1,0,-0.5,-0.5,-0.5,-0.75,
             0,0,0,0,0,0,
             0,0,0,0,0,0,
             0,0,0,0,0,0,
             1,0,0,0,0,0,
             0,1,0,0,0,0,
             0,0,1,0,0,0,
             0,0,0,1,0,0,
             0,0,0,0,1,0,
             0,0,0,0,0,1,
             0,0,0,0,0,0,
             0,0,0,0,0,0,
             0,0,0,0,0,0,
             0,0,0,0,0,0,
             0,0,0,0,0,0,
             0,0,0,0,0,0),nrow=6)
KA
AMBV<-KA%*%y_hat
SEP_AMBV<-sqrt(diag(KA%*%ginv(lhs)%*%t(KA)))



KN<-matrix(c(0,0,0,0,0,0,
             0,0,0,0,0,0,
             0,0,0,0,0,0,
             0.5,0.5,0.5,0.5,0.5,0.5,
             0,0,0,0,0,0,
             0,0,0,0,0,0,
             0,0,0,0,0,0,
             0,0,0,0,0,0,
             0,0,0,0,0,0,
             0,0,0,0,0,0,
             0,0,0,0,0,0,
             0,0,0,0,0,0,
             0.5,0,0,0,0,0,
             0,0.5,0,0,0,0,
             0,0,0.5,0,0,0,
             0,0,0,0.5,0,0,
             0,0,0,0,0.5,0,
             0,0,0,0,0,0.5),nrow=6)
KN
NMBV<-KN%*%y_hat
SEP_NMBV<-sqrt(diag(KN%*%ginv(lhs)%*%t(KN)))



KT<-matrix(c(0,0,0,0,0,0,
             1,0,0.5,0.5,0.5,0.75,
             -1,0,-0.5,-0.5,-0.5,-0.75,
             0.5,0.5,0.5,0.5,0.5,0.5,
             0,0,0,0,0,0,
             0,0,0,0,0,0,
             1,0,0,0,0,0,
             0,1,0,0,0,0,
             0,0,1,0,0,0,
             0,0,0,1,0,0,
             0,0,0,0,1,0,
             0,0,0,0,0,1,
             0.5,0,0,0,0,0,
             0,0.5,0,0,0,0,
             0,0,0.5,0,0,0,
             0,0,0,0.5,0,0,
             0,0,0,0,0.5,0,
             0,0,0,0,0,0.5),nrow=6)
KT
TMBV<-KT%*%y_hat
SEP_TMBV<-sqrt(diag(KT%*%ginv(lhs)%*%t(KT)))



sol_k<-cbind(round(AMBV,2), round(SEP_AMBV,2), round(NMBV,2), round(SEP_NMBV,2), round(TMBV,2), round(SEP_TMBV,2))
rownames(sol_k)<-c("1","2","3","4","5","6")
colnames(sol_k)<-c("AMBV","SEP_AMBV","NMBV","SEP_NMBV","TMBV","SEP_TMBV")
sol_k
