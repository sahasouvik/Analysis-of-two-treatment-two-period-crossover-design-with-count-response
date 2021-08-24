#Here G is the number of categories
u=NULL
for(i in 1:50){
a=runif(1,0,1)
if(a<=1/3){
u[i]=0
}else if(a<=2/3){
u[i]=1
}else{
u[i]=2
}
}
#Generating the vactor containing the categories
u
a=0.9
lambda_A=5
lambda_B=4
lambda_AB=4.5
lambda_BA=3.5
beta=4
X_A=NULL
X_B=NULL
for(i in 1:50){
x=rpois(1,lambda=lambda_A*a^(2-u[i]))
r=rpois(1,lambda=lambda_B*a^(2-u[i]))
X_A=c(X_A,x)
X_B=c(X_B,r)
}
X_A
X_B
X_AB=NULL
X_BA=NULL
for(i in 1:50){
x=rpois(1,lambda=lambda_AB+beta*(lambda_AB-lambda_A)*a^(2-u[i])+beta*X_A[i])
r=rpois(1,lambda=lambda_BA+beta*(lambda_BA-lambda_B)*a^(2-u[i])+beta*X_B[i])
X_AB=c(X_AB,x)
X_BA=c(X_BA,r)
}
X_AB
X_BA
cor(X_A,X_AB)
cor(X_B,X_BA)