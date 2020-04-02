#start_time <- Sys.time()
#source("LSVCE.R")
#end_time <- Sys.time()
#end_time-start_time
#!/usr/bin/Rscript --slave

#start_time <- Sys.time()
#library("readr")
#dat1a<-read.table("../RA.out",nrow=21996)
#dat2a<-read.table("../RA.out",skip=21996)
#dat1 <- dat1a[sample(1:nrow(dat1a), 2000, replace=FALSE),]
#dat2 <- dat2a[sample(1:nrow(dat2a), 2000, replace=FALSE),]

dat1<-read.table("../RA.out",nrow=1500)
date1=dat1[,c(1)]
deriv1=as.matrix(dat1[,c(4:15)])
oc1 <- dat1[,c(3)]

dat2=read.table("../RA.out",skip=21996,nrow=600)
date2=dat2[,c(1)]
deriv2=as.matrix(dat2[,c(4:15)])
oc2 <- dat2[,c(3)]

ua=1.5e11

#print(deriv2)

#print("ALKHQP==============")
#result2 <- lsfit(deriv1,oc1, intercept = FALSE)
#sol  <- coef(result2)
#print(sol)
#sigma <- ls.diag(result2)$std.err
#print(sigma)

A=rbind(deriv1,deriv2)
n1=length(oc1)
n2=length(oc2)

# simulation pour valider
#oc1=rnorm(n1,sd=10)
oc1=oc1/ua

#oc2=rnorm(n2,sd=7)
oc2=oc2/ua
yobs=c(oc1,oc2)
# LSVCE ===================

#library(matlib)
library("matrixcalc")
# inv(A) = inverse
# t(A) = transpose

c1=1
c2=1

Q11=diag(n1)
Q10=matrix(0L, nrow = n2, ncol = n1)
Q111=rbind(Q11,Q10)
Q01=matrix(0L, nrow = n2+n1, ncol = n2)
Q1=cbind(Q111,Q01)

Q22=diag(n2)
Q20=matrix(0L, nrow = n1, ncol = n2)
Q111=rbind(Q20,Q22)
Q01=matrix(0L, nrow = n2+n1, ncol = n1)
Q2=cbind(Q01,Q111)

Qy=c1*Q1+c2*Q2
invQy=(1./c1)*Q1+(1./c2)*Q2
#print(invQy)
#print(t(A))
#print("avant inv ==========")
#print(t(A)%*%invQy%*%A)

library("Rfast")

Nqy=matrix.inverse(t(A)%*%invQy%*%A)%*%t(A)%*%invQy
#print("NQY ==========")
#print(matrix.inverse(t(A)%*%invQy%*%A))

paortho=diag(n1+n2)-A%*%Nqy
e=paortho%*%(yobs)
matr=invQy%*%paortho
#print("paortho")
#print(paortho)
#print("e")
#print(e)
#print("matr")
#print(matr)
#N11=0.5*matrix.trace(Q1%*%matr%*%(Q1)%*%(matr))
#N22=0.5*matrix.trace(Q2%*%matr%*%(Q2)%*%(matr))
#N21=0.5*matrix.trace(Q2%*%matr%*%(Q1)%*%(matr))
aa=matr%*%(Q1%*%matr)
ab=matr%*%(Q2%*%matr)
#print("aa =======")
#print(Q1%*%aa)
#print("ab =======")
#print(Q1%*%ab)

##p1=Q1%*%aa
#p2=Q1%*%ab
#p3=Q2%*%aa
#p4=Q2%*%ab

p1=mat.mult(Q1,aa)
p2=mat.mult(Q1,ab)
p3=mat.mult(Q2,aa)
p4=mat.mult(Q2,ab)

N11=0.5*matrix.trace(p1)
N12=0.5*matrix.trace(p2)
N22=0.5*matrix.trace(p4) #matr%*%(Q2)%*%(matr))
N21=0.5*matrix.trace(p3) #matr%*%(Q1)%*%(matr))

N <- matrix( c(N11,N12,
	N12,N22),nrow=2, byrow=TRUE)

l1=0.5*(t(e)%*%(invQy)%*%(Q1)%*%invQy)%*%(e)
l2=0.5*(t(e)%*%(invQy)%*%(Q2)%*%invQy)%*%(e)

l=c(l1,l2)
#print(l)

newsigv=sqrt(matrix.inverse(N)%*%(l))*ua
print(newsigv*ua)
##end_time <- Sys.time()
##print(end_time-start_time)

