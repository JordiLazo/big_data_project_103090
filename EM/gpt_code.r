
data <- read.table("Batman.data",header=TRUE)


X<-data$X

mean(X)
var(X)
n<-length(X) #n1+n2

### Transform X variable in a matrix, you can rotate image by changing byrow=T o byrow=F
X_image<-matrix(X,nrow=sqrt(n),ncol=sqrt(n),byrow=F)

image(X_image)

###observe histogram of the values ##
hist(X_image)


### Start with some initial values for Mu1, Mu2, V1, V2, pi1, pi2
## We could have also considered starting with new values for the Z1, Z2 new random variables, but I prefer to give values to Mu1, Mu2, etc.
Mu1<-4.0  #180.0
Mu2<-2.0   #170.0
V1<-2.15    #272.0
V2<-2.15    #272.0
pi1<-0.75 #0.65
pi2<-1-pi1
####

##Working variables
values<-1 
break_s<-0 
pi0<-10.0
n_simul<-1000 ###numbMu1<-round((t(Z1) %*% X)/sum(Z1),4)
mx1<-matrix(rep(1,n),nrow=n, ncol=1)
epsi<-0.0001 #estimators difference
######


#Plot to show the evolution of estimates of the parameter pi1
plot(0,pi1,ylim=c(0,1),xlim=c(0,200))
###########################LOOP############################
for(j in 1:n_simul){ #First For loop
if(break_s==0){
  if(values==0){
##MLE estimators for a Gaussian mixture model
   Mu1<-round((t(Z1) %*% X)/sum(Z1),4)
   Mu2<-round((t(Z2) %*% X)/sum(Z2),4)
   B1<-X-mx1 %*% Mu1
   B2<-X-mx1 %*% Mu2
   V1<- round((t(B1) %*% diag(array(Z1)) %*% B1)/sum(Z1),4)
   V2<- round((t(B2) %*% diag(array(Z2)) %*% B2)/sum(Z2),4)
   pi1<-round(sum(Z1)/n,4)
   pi2<-round(sum(Z2)/n,4)
  }
##We always start here for the first iteration
values<-0 ##parameter values controls that for the first iteration we start here the code
 
### calcul new Z1, values of new variables ###
Exp1<-c() ##Here we set vector Exp1, expectation of the Z1 values
for(i in 1:n){
f1<-dnorm(X[i],mean=Mu1,sd=sqrt(V1)) #density function for a Normal variable
f2<-dnorm(X[i],mean=Mu2,sd=sqrt(V2))
Exp1[i]<-(pi1*f1)/(pi1*f1+(1.0-pi1)*f2) #Expectation of Z1 values
}
###we can show the Exp1 values in an histogram
#hist(p1)


### calcul new Z2, values of new variables ###
##Note that Exp2 can be just computed as Exp2=1-Exp1
Exp2<-c()
for(i in 1:n){
f1<-dnorm(X[i],mean=Mu1,sd=sqrt(V1))
f2<-dnorm(X[i],mean=Mu2,sd=sqrt(V2))
Exp2[i]<-(pi2*f2)/(pi2*f2+(1.0-pi2)*f1) #Expectation of Z2 values
}

####now rewrite Z1 and Z2, these are the values that the new used for new iteration, to obtain new mu1, mu2, V1, v2 and pi1, pi2###
Z1<-round(Exp1,4)
Z2<-round(Exp2,4)
print(pi1)

##loop to check if we can stop the iteration procedure
 if(abs(pi1-pi0)<epsi){break_s<-1}

 } else if(break_s==1){
           break }
pi0<-pi1 
print(pi1)
points(j,pi1,pch=19)

}#close First For loop



### Print estimated parameters ###
Mu1; Mu2; V1; V2; pi1

####### Genereta new image ####
### We can display Z1 or Z2 ####
x11()
X_image2<-matrix(t(Z1),nrow=sqrt(n),ncol=sqrt(n),byrow=F)
image(X_image2)

### assign 0 or 1 in terms of Z1 probabilities ###
Z11<-array(rep(0,n))
for(i in 1:n){
if(Z1[i]>0.5){Z11[i]=1}
}
Z11_image<-matrix(Z11,nrow=sqrt(n),ncol=sqrt(n),byrow=F)
x11()
image(Z11_image)


### To generate a graphic ###
postscript(file="Image3.ps")
pdf(file="Image3.pdf")
jpeg()
layout(matrix(c(1,2),1,2,byrow=TRUE),respect=FALSE)
image(X_image)
image(Z11_image)
dev.off()
dev.off()
dev.off()