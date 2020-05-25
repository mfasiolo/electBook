#knots function
default.knots=function(x,num.knots)
{  if (missing(num.knots)) num.knots=max(5,min(floor(length(unique(x))/4),35))
   return(quantile(unique(x),seq(0,1,length=(num.knots+2))[-c(1,(num.knots+2))]))}

#variance smoother
newtonr=function(lambda2, eta, W2, R, tol=10^(-2), maxiter=25){
dif=1; k=0
while(dif>tol & k<maxiter){
k=k+1
sigma=exp(W2%*%eta)
nR=R/sqrt(sigma)
H1=0
for(i in 1:n){ H1=H1+t(W2[ID==i,]*nR[ID==i])%*%(nR[ID==i]*W2[ID==i,])}
u=drop(apply(W2, 2, sum)-t(nR^2)%*%W2)+drop(2*lambda2*D2%*%eta)
H=H1+2*lambda2*D2
eta1=eta-solve(H)%*%u
lambda2=1/(mean((diag(solve(H)))[(p2):(p2-1+K2)])+mean((eta1[(p2):(p2-1+K2)])^2))
dif=sum((eta1-eta)^2)
eta=eta1}
print(lambda2)
return(list(eta=eta1,lambda=lambda2))}

#update smoothing parameter for variance smoother
eta_f=function(lambda2, eta, W2, R, tol=10^(-2), maxiter=2){
k=0; dif=1
while(dif>tol && k < maxiter){
k=k+1
temp=newtonr(lambda2,eta, W2, R)
nlambda2=temp$lambda
dif=abs(nlambda2-lambda2)
lambda2=nlambda2
}
return(list(eta=temp$eta, lambda2=lambda2))}

 
#EM algorithm
EM_f=function(beta.hat, sigma, lambda1, lambda3, lambdanu, D,tol=10^(-2), maxiter=25){
Dnu=solve(solve(D)+lambdanu*D3)
dif=1; k=0
while(dif>tol & k<maxiter){
k=k+1
a=d=0
for(i in 1:n){
varm=diag(sigma[ID==i])
invcovm=solve(varm+ Znu[ID==i,]%*%Dnu%*%t(Znu[ID==i,]))
a=a+t(M[ID==i,])%*%invcovm%*%M[ID==i,]
d=d+t(M[ID==i,])%*%invcovm%*%Y[ID==i]}
nbeta.hat=solve(a+ c(rep(lambda1, p+K1), rep(lambda3,p+K1))*D12)%*%d

#b=matrix(rep(0,n*(pnu+Knu+1)),nrow=n)
b=matrix(rep(0,n*(pnu+Knu)),nrow=n)
R=Y-cbind(W1, G)%*%nbeta.hat
for(i in 1:n){
varm=diag(sigma[ID==i])
invcovm=solve(varm+Znu[ID==i,]%*%Dnu%*%t(Znu[ID==i,]))
b[i,]=Dnu%*%t(Znu[ID==i,])%*%invcovm%*%(R[ID==i])}

D0=0
for(i in 1:n){
varm=diag(sigma[ID==i])
invcovm=solve(varm+Znu[ID==i,]%*%Dnu%*%t(Znu[ID==i,]))
nQ=invcovm-invcovm%*%M[ID==i,]%*%solve(a+lambda1*D12)%*%t(M[ID==i,])%*%invcovm
D0=D0+b[i,]%*%t(b[i,])+Dnu-Dnu%*%t(Znu[ID==i,])%*%nQ%*%Znu[ID==i,]%*%Dnu}

dif=mean((beta.hat-nbeta.hat)^2)
beta.hat=nbeta.hat}

return(list(D=D0/n, nbeta.hat=nbeta.hat, b.hat=b))}

#choice of smoothing parameter  for mean function
lambda1f=function(lambda, lambda3, lambdanu, sigma, beta, maxiter=3){
r=Y-M1%*%c(beta[1:p], beta[(p+K1+1):(2*p+K1)])
n.row=0
V1=matrix(rep(0,N1^2),nrow=N1)
Dnu=solve(solve(D)+lambdanu*D3)
for(i in 1:n){
temp1=length(Y[ID==i])
varm=diag(sigma[ID==i])
V1[(n.row+1):(n.row+temp1),(n.row+1):(n.row+temp1)]=varm+Znu[ID==i,]%*%Dnu%*%t(Znu[ID==i,])
n.row=n.row+temp1}

dif=1; k=0
while(dif>.01 & k<maxiter){
V=V1+Z1%*%t(Z1)*c(lambda)+GZ1%*%t(GZ1)/lambda3 
inV=solve(V)
dV=Z1%*%t(Z1)  
dinV=-inV%*%dV%*%inV
tempm2=solve(t(M1)%*%inV%*%M1)
u=sum(diag(inV%*%dV))+sum(diag(tempm2%*%t(M1)%*%dinV%*%M1))+t(r)%*%dinV%*%r
H=sum(diag(dinV%*%dV))-2*t(r)%*%dinV%*%dV%*%inV%*%r-sum(diag(tempm2%*%t(M1)%*%dinV%*%M1%*%tempm2%*%t(M1)%*%dinV%*%M1))-2*sum(diag(tempm2%*%t(M1)%*%dinV%*%dV%*%inV%*%M1))
nlambda=lambda-u/H;
dif=abs(nlambda-lambda)
lambda=nlambda
k=k+1}

return(c(1/nlambda))}

#choice of smoothing parameter  for varying coefficient
lambda3f=function(nlambda=25, lambda1, lambdanu, sigma, beta){
REML=rep(0, nlambda) 
r=Y-M1%*%c(beta[1:p], beta[(p+K1+1):(2*p+K1)]) 
V1=matrix(rep(0,N1^2),nrow=N1) 
Dnu=solve(solve(D)+lambdanu*D3)
for(j in 1:nlambda){ lambda3=1.5^(-1/2*nlambda+j) 
n.row=0 
 for(i in 1:n){temp1=length(Y[ID==i]) 
 varm=diag(sigma[ID==i])
 V1[(n.row+1):(n.row+temp1),(n.row+1):(n.row+temp1)]=varm+Znu[ID==i,]%*%Dnu%*%t(Znu[ID==i,]) 
 n.row=n.row+temp1}
V=V1+Z1%*%t(Z1)/lambda1+GZ1%*%t(GZ1)/lambda3 
inV=solve(V) 
a=det(V);  if (a<10^(-300)) loga=-300*log(10) else loga=log(a)
REML[j]=loga+t(r)%*%inV%*%r+log(det(t(M1)%*%inV%*%M1))
if(is.na(REML[j]) | REML[j]>10^5) REML[j]=10^5 }
vec=cbind(REML, (1:nlambda))[order(REML),]
return(1.5^(-1/2*nlambda+vec[1,2]))}


#choice of smoothing parameter for subject specific curves
lambdanuf=function(nlambda=25, sigma, lambda1){
REML=rep(0,nlambda)
for(j in 1: nlambda){
lambdanu=1.5^(-1/2*nlambda+j) 
r=Y-M%*%nbeta.hat
Dnu=solve(solve(D)+lambdanu*D3) 

temp1=0
 for(i in 1:n){
 varm=diag(sigma[ID==i])
 covm=(varm+Znu[ID==i,]%*%Dnu%*%t(Znu[ID==i,]))  
 invcovm=solve(covm)
 REML[j]=REML[j]+log(det(covm))+t(r[ID==i])%*%invcovm%*%r[ID==i]
 temp1=temp1+t(M[ID==i,])%*%invcovm%*%M[ID==i,]}

REML[j]=REML[j]+log(det(temp1))
if(is.na(REML[j])) REML[j]=10^5 }
vec=cbind(REML, (1:nlambda))[order(REML),]
return(1.5^(-1/2*nlambda+vec[1,2]))}


