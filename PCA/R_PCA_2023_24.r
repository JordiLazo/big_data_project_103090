#                 ####Data Heredity of Headform in Man####

# X <- matrix(c(191, 155,
#               195, 149,
#               181, 148, 
#               183, 153, 
#               176, 144, 
#               208, 157, 
#               189, 150, 
#               197, 159, 
#               188, 152, 
#               192, 150, 
#               179, 158, 
#               183, 147, 
#               174, 150, 
#               190, 159, 
#               188, 151, 
#               163, 137, 
#               195, 155, 
#               186, 153, 
#               181, 145, 
#               175, 140, 
#               192, 154, 
#               174, 143, 
#               176, 139, 
#               197, 167, 
#               190, 163),
#             nrow=25,     ## number of rows 
#             ncol=2,      ## number of columns 
#             byrow = T)   ## fill matrix by rows
# n <- nrow(X); m <- ncol(X)    # number of rows and columns
# name_v<-c("X1","X2") ### name of variables 
# name_vl<-c(seq(1,n,1))
# print(X)







                   ### If we call an external file ###

# data <- read.table("School_qualifications.data",header=TRUE)
# data <- read.table("School_qualifications2.csv",sep=",",header=TRUE)
# data <- read.table("forestfires.csv",sep=",",header=TRUE)
# data <- read.table("Car.csv",sep=",",header=TRUE)
data <- read.table("Fire.csv",sep=",",header=TRUE)
# data <- read.table("Survey.csv",sep=",",header=TRUE)
# data <- read.table("wineDataset.csv",sep=",",header=TRUE)
# data <- read.table("BodyFat.csv",sep=",",header=TRUE)

n<-length(data[,1]) #number of subjects
m<-length(data[1,]) #number of variables
name_v<-c() #column names, variables
name_vl<-c() #row names, individuals
for(i in 1:m){name_v[i]<-colnames(data)[i]}
for(i in 1:n){name_vl[i]<-rownames(data)[i]}
i1<-0
k<-c()
for(i in 1:n){
 for(j in 1:m){
   i1<-i1+1
  k[i1]<-data[i,j]
 }
}

X <- matrix(k,
            nrow=n,     ## number of rows 
            ncol=m,      ## number of columns 
            byrow = T)   ## fill matrix by rows





               ########## Assuming that we have the X matrix ############### 


## obtain the 1 matrix####
mx1<-matrix(rep(1,n),nrow=n, ncol=1)
#### 1 transpose matrix ####
mx1T<-t(mx1)
#### Identity matrix ####
Id <- diag(n)
#### Weight matrix D ###
D <- Id/(n) 

#### Compute the  B matrix (centrered X matrix) ####

B<-(Id-mx1 %*% mx1T %*% D) %*% X

##Obtain the average value of each column
X[,1] #first column
mean(X[,1]) #mean first column
mean(X[,2]) #mean second column
X[,1]-mean(X[,1]) #corrected column
B[,1] #compare
X[,2]-mean(X[,2]) #corrected column
B[,2] #compare


### Compute the S matrix (variance-covariance matrix of X) ###

S<-(t(B) %*% B)/(n)

##Remenber the variance-covariance matrix are the poblational parameters (/n), try sum(B[,1]^2)/n=S[1,1], sum(B[,2]^2)/n=S[2,2], or
##sum(B[,1]*B[,2])/n=S[1,2]=S[2,1]

### Obtain a diagonal matrix D_1s of order m with the inver of standard deviation (using S)

D_1s<- diag(sqrt(1/diag(S)), nrow=m, ncol=m)

### Now D_1s is a diagonal matrix of order m with the inver of standard deviation

### Obtain the Y matrix (centrered and reduced X matrix) ###

Y<-B %*% D_1s

##Population variance var[
### obtain correlation matrix Directly with R ###
R<-cor(X,method="pearson")

###Note that variance-covariance matrix of Y is equal to R
###Try sum(Y[,1]^2)/n=R[1,1], sum(Y[,2]^2)/n=R[2,2], or sum(Y[,1]*Y[,2])/n=R[1,2]=R[2,1]

##if you write var(Y[,1]) you will obtain the sample variance (we divide by n-1), if we want the population variance, we need to obtain 
## var(Y[,1])*(n-1)/n

#plot(Y[,1],Y[,2])

### Obtain Eigenvalues and Eigenvector from matrix R ###

Eig<-eigen(R)

### Eigenvalues ###

Eig$values

### Eigenvectors ### Each column is a eigenvector

Eig$vectors

### First component ###
FC<-matrix(Eig$vector[,1],
                  nrow=m,
                  ncol=1,
                  byrow=T)
Z1<-Y %*% (1*FC)
print("JORDIIII")
print(FC)

### We can change the signe of the components, as long we change all the values
### If we want to compare Z1 and Z2 in a graph we better change the sign of all the components
### Second component ###
SC<-matrix(Eig$vector[,2],
                  nrow=m,
                  ncol=1,
                  byrow=T)
Z2<-Y %*% (1*SC)

TC<-matrix(Eig$vector[,3],
                  nrow=m,
                  ncol=1,
                  byrow=T)
Z3<-Y %*% (1*TC)


### Explained variability for each component ###
Z_var<-round(Eig$values/m*100,2)

### Scree Plot ###
pdf(file="Scree.pdf")
plot(seq(1,m,1), Eig$values, asp=1, type = "b", xlab = "Components", ylab = "Variance",
     main="Scree plot",cex.main=2.0,cex.axis=1.50,cex.lab=1.5,pch=19)
dev.off()

### Biplot ###
pdf(file="Biplot.pdf")
plot(Z1, Z2, asp=1,,cex.axis=1.30,cex.lab=1.30,main="Biplot",cex.main="1.5",pch=19,
     xlab=paste("Factor 1 (", Z_var[1],"% )"),
     ylab=paste("Factor 2 (", Z_var[2],"% )"))
abline(h=0,col="blue");abline(v=0,col="blue")
for(i in 1:n){
#text(Z1[i]+runif(1,0,max(Z1)*0.1),Z2[i]+runif(1,0,max(Z2)*0.1),name_vl[i],cex=1.50)
text(Z1[i]+max(Z1)*0.05,Z2[i]+max(Z2)*0.05,name_vl[i],cex=1.0)
}
dev.off()


### To obtain the matrix of the correlation between components and variables
### Here we consider only 3 components
Rcx <- diag(0, nrow=2, ncol=m)
for (j in 1:2){
  for (k in 1:m) {
    Rcx[j,k] <- sqrt(Eig$values[j])*Eig$vector[k,j]
  }
}
#rows are components; columns are variables
#For instace, for the case of head measurements, the correlation between the first component and variable Y1 is 0.9312
#The correlation between the second component and variable Y1 is negative and with the second variable Y2 is positive

### Correlation circle ###
pdf(file="Corr_plot.pdf")
t <- seq(0,1,0.01)
x <- cos(2*pi*t)
y <- sin(2*pi*t)
plot(x, y, asp=1, type = "l", xlab = "C1", ylab = "C2",
     main="Correlacions between components and variables",axes=FALSE,xlim=c(-1,1),ylim=c(-1,1),cex.main=1.50,
cex.lab=1.30)
axis(1,seq(-1,1,0.5),crt=360,cex.axis=1.30,pos=c(-1,-1))
axis(2,seq(-1,1,0.5),cex.axis=1.30,pos=c(-1,-1),las=2,hadj=1.0)
lines(seq(-1,1,0.5),rep(0,5),col="blue")
lines(rep(0,5),seq(-1,1,0.5),col="blue")
points(Rcx[1,],Rcx[2,], pch=19)
for(i in 1:m){
text(Rcx[1,i]+0.05,Rcx[2,i]+0.05,name_v[i],cex=1.50)
}
rect(-1,-1,1,1)
dev.off()


#####Using the prcomp package of R#######


X.pca<-prcomp(X, center = TRUE,scale.=TRUE)
print(X.pca)
X.pca$sdev
sum((X.pca$sdev)^2)
X.pca$rotation
#Stand deviation are the sqrt(lambda)
#rotation matrix is the matrix of all the Principal components, a square matrix with dimension the number of variables

#The vector with the first PC is 
X.pca$rotation[,1]

#Scree plot
plot(X.pca,type="l")

##The output for the Z new variables are (this a matrix with row the number of individuals and columns de number of variables
X.pca$x




##########Football player data set######
X <- matrix(c(15.5, 60.0, 21.1, 10.3, 13.4, 12.4,
              15.4, 59.7, 20.0, 12.8, 14.5, 11.3,
              15.1, 59.7, 20.2, 11.4, 14.1, 12.1,
              14.3, 56.9, 18.9, 11.0, 13.4, 11.0,
              14.8, 58.0, 20.1,  9.6, 11.1, 11.7,
              15.2, 57.5, 18.5,  9.9, 12.8, 11.4,
              15.4, 58.0, 20.8, 10.2, 12.8, 11.9, 
              16.3, 58.0, 20.1,  8.8, 13.0, 12.9,
              15.5, 57.0, 19.6, 10.5, 13.9, 11.8,
              15.0, 56.0, 19.6, 10.4, 14.5, 12.0,  
              15.5, 57.2, 20.0, 11.2, 13.4, 12.4, 
              15.5, 56.5, 19.8,  9.2, 12.8, 12.2,
              15.7, 57.5, 19.8, 11.8, 12.6, 12.5,
              14.4, 57.0, 20.4, 10.2, 12.7, 12.3,
              14.9, 54.8, 18.5, 11.2, 13.8, 11.3,
              16.5, 59.8, 20.2,  9.4, 14.3, 12.2, 
              15.5, 56.1, 18.8,  9.8, 13.8, 12.6,
              15.3, 55.0, 19.0, 10.1, 14.2, 11.6, 
              14.5, 55.6, 19.3, 12.0, 12.6, 11.6,
              15.5, 56.5, 20.0,  9.9, 13.4, 11.5, 
              15.2, 55.0, 19.3,  9.9, 14.4, 11.9,
              15.3, 56.5, 19.3,  9.1, 12.8, 11.7,
              15.3, 56.8, 20.2,  8.6, 14.2, 11.5,
              15.8, 55.5, 19.2,  8.2, 13.0, 12.6,  
              14.8, 57.0, 20.2,  9.8, 13.8, 10.5,
              15.2, 56.9, 19.1,  9.6, 13.0, 11.2,
              15.9, 58.8, 21.0,  8.6, 13.5, 11.8,
              15.5, 57.3, 20.1,  9.6, 14.1, 12.3,
              16.5, 58.0, 19.5,  9.0, 13.9, 13.3,
              17.3, 62.6, 21.5, 10.3, 13.8, 12.8,
              14.9, 56.5, 20.4,  7.4, 13.0, 12.0,
              15.4, 57.5, 19.5, 10.5, 13.8, 11.5,
              15.3, 55.4, 19.2,  9.7, 13.3, 11.5,
              14.6, 56.0, 19.8,  8.5, 12.0, 11.5,
              16.2, 56.5, 19.5, 11.5, 14.5, 11.8,
              14.6, 58.0, 19.9, 13.0, 13.4, 11.5,
              15.9, 56.7, 18.7, 10.8, 12.8, 12.6,
              14.7, 55.8, 18.7, 11.1, 13.9, 11.2,
              15.5, 58.5, 19.4, 11.5, 13.4, 11.9,
              16.1, 60.0, 20.3, 10.6, 13.7, 12.2,
              15.2, 57.8, 19.9, 10.4, 13.5, 11.4,
              15.1, 56.0, 19.4, 10.0, 13.1, 10.9,
              15.9, 59.8, 20.5, 12.0, 13.6, 11.5,
              16.1, 57.7, 19.7, 10.2, 13.6, 11.5,
              15.7, 58.7, 20.7, 11.3, 13.6, 11.3,
              15.3, 56.9, 19.6, 10.5, 13.5, 12.1,
              15.3, 56.9, 19.5,  9.9,  14.0, 12.1,
              15.2, 58.0, 20.6, 11.0, 15.1, 11.7,
              16.6, 59.3, 19.9, 12.1, 14.6, 12.1,
              15.5, 58.2, 19.7, 11.7, 13.8, 12.1,
              15.8, 57.5, 18.9, 11.8, 14.7, 11.8,
              16.0, 57.2, 19.8, 10.8, 13.9, 12.0,
              15.4, 57.0, 19.8, 11.3, 14.0, 11.4,
              16.0, 59.2, 20.8, 10.4, 13.8, 12.2,
              15.4, 57.6, 19.6, 10.2, 13.9, 11.7,
              15.8, 60.3, 20.8, 12.4, 13.4, 12.1,
              15.4, 55.0, 18.8, 10.7, 14.2, 10.8,
              15.5, 58.4, 19.8, 13.1, 14.5, 11.7,
              15.7, 59.0, 20.4, 12.1, 13.0, 12.7,
              17.3, 61.7, 20.7, 11.9, 13.3, 13.3),
            nrow=60,     ## number of rows 
            ncol=6,      ## number of columns 
            byrow = T)   ## fill matrix by rows
n <- nrow(X); m <- ncol(X)    # number of rows and columns
name_v<-c("A","B","C","D","E","F") ### name of variables 
name_vl<-c(seq(1,n,1))






### Representation of data and first component ###

postscript(file="Ortho_fig.ps")
par(mar=2.50+c(1.50,8.80,1.50,8.80))
plot(Y[,1],Y[,2])
#lines(c(-21,0),c(-14.4,0), type = "l", col = "blue")
lines(c(-21,21)*0.1024398,c(-14.4,14.4)*0.1356314, type = "l", col = "blue")
lines(c(-11.72,0)*0.1024398,c(-1.12,0)*0.1356314,col="yellow",lwd=3)
lines(c(-11.72,-7.05)*0.1024398,c(-1.12,-4.85)*0.1356314,col="red",lwd=3)
lines(c(-7.05,0)*0.1024398,c(-4.85,0)*0.1356314,col="green",lwd=3)
text(-13.0*0.1024398,-2*0.1356314,"x1")
text(-10.0*0.1024398,-4*0.1356314,"r1")
text(-1.50*0.1024398,-2.5*0.1356314,"z1")
text(-5.0*0.1024398,1.0*0.1356314,"h1")
abline(h=0)
abline(v=0)
dev.off()



