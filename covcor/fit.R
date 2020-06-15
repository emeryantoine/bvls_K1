system("paste Residuals.txt Jacobian.txt > RA.out")
dat=as.matrix(read.table("../../transfert/cas_complet/040520/RAW.out.412"))
deriv=dat[,c(4:415)]

oc2 <- dat[,c(1)]*1.0
result2 <- lsfit(deriv,oc2, intercept = FALSE)
sol  <- coef(result2)
print(sol)
sigma <- ls.diag(result2)$std.err
print(sigma)

matrcorr <- ls.diag(result2)$corr
print(matrcorr)

#print(matrcorr)
#corrplot(matrcorr,method="pie")

#array([2.24125218e+01, 2.33252756e+02, 1.23250823e+03, 2.55331392e-07,
#     1.25508987e+05])
       
#MTM = np.dot(M.transpose(),M)
##MTMINV = np.linalg.inv(MTM)
#MTY = np.dot(M.transpose(),Y)
#P = np.dot(MTMINV,MTY)

#Residuals = Y - np.dot(M,P)
#ChiSq = np.dot(Residuals.transpose(),Residuals)

#C = MTMINV * ChiSq/(nobs-npar)

library("matrixcalc")
#library("Matrix")

mat=t(deriv)%*%deriv
invmat=matrix.inverse(mat)
MTY=t(deriv)%*%oc2
P=invmat%*%MTY
Residuals = oc2 - deriv%*%P
ChiSq = t(Residuals)%*%Residuals
fact=ChiSq/(nrow(deriv)-ncol(deriv))
C=as.vector(fact)*diag(invmat)
sqrt(as.vector(fact)*diag(invmat))


# avec ponderation on a quand meme besoin de S/(nobs-nparam)
# pour cov

deriv=dat[,c(2:6)]
oc2 <- dat[,c(1)]
W = diag(length(oc2))
sigma=10e-3
W=((1.0/sigma**2))*W

# calcul de delta
mat=t(deriv)%*%W%*%deriv
invmat=matrix.inverse(mat)
MTY=t(deriv)%*%W%*%oc2
P=invmat%*%MTY

# calcul des residus
Residuals = oc2 - deriv%*%P

# calcul covariance
ChiSq = t(Residuals)%*%W%*%Residuals
fact=ChiSq/(nrow(deriv)-ncol(deriv))
C=as.vector(fact)*diag(invmat)
sqrt(as.vector(fact)*diag(invmat))
