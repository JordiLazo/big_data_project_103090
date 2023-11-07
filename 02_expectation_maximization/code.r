data <- read.table("Batman.data",header=TRUE)


X<-data$X

mean(X)
var(X)
n<-length(X)
X_image<-matrix(X,nrow=sqrt(n),ncol=sqrt(n),byrow=F)

image(X_image)

hist(X_image)


Mu1<-4.0
Mu2<-2.0
V1<-2.15
V2<-2.15
pi1<-0.75
pi2<-1-pi1

values<-1 
break_s<-0 
pi0<-10.0
n_simul<-1000
mx1<-matrix(rep(1,n),nrow=n, ncol=1)
epsi<-0.0001

plot(0,pi1,ylim=c(0,1),xlim=c(0,200))

for(j in 1:n_simul){
if(break_s==0){
  if(values==0){

   Mu1<-round((t(Z1) %*% X)/sum(Z1),4)
   Mu2<-round((t(Z2) %*% X)/sum(Z2),4)
   B1<-X-mx1 %*% Mu1
   B2<-X-mx1 %*% Mu2
   V1<- round((t(B1) %*% diag(array(Z1)) %*% B1)/sum(Z1),4)
   V2<- round((t(B2) %*% diag(array(Z2)) %*% B2)/sum(Z2),4)
   pi1<-round(sum(Z1)/n,4)
   pi2<-round(sum(Z2)/n,4)
  }

values<-0 
Exp1<-c()
for(i in 1:n){
f1<-dnorm(X[i],mean=Mu1,sd=sqrt(V1))
f2<-dnorm(X[i],mean=Mu2,sd=sqrt(V2))
Exp1[i]<-(pi1*f1)/(pi1*f1+(1.0-pi1)*f2)
}

Exp2<-c()
for(i in 1:n){
f1<-dnorm(X[i],mean=Mu1,sd=sqrt(V1))
f2<-dnorm(X[i],mean=Mu2,sd=sqrt(V2))
Exp2[i]<-(pi2*f2)/(pi2*f2+(1.0-pi2)*f1)
}

Z1<-round(Exp1,4)
Z2<-round(Exp2,4)
print(pi1)

 if(abs(pi1-pi0)<epsi){break_s<-1}

 } else if(break_s==1){
           break }
pi0<-pi1 
print(pi1)
points(j,pi1,pch=19)

}

Mu1; Mu2; V1; V2; pi1

x11()
X_image2<-matrix(t(Z1),nrow=sqrt(n),ncol=sqrt(n),byrow=F)
image(X_image2)

Z11<-array(rep(0,n))
for(i in 1:n){
if(Z1[i]>0.5){Z11[i]=1}
}
Z11_image<-matrix(Z11,nrow=sqrt(n),ncol=sqrt(n),byrow=F)
x11()
image(Z11_image)

postscript(file="Image3.ps")
pdf(file="Image3.pdf")
jpeg()
layout(matrix(c(1,2),1,2,byrow=TRUE),respect=FALSE)
image(X_image)
image(Z11_image)
dev.off()
dev.off()
dev.off()