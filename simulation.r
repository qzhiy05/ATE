
m <- 1000


simipw <- function(x1,x2,x3){

n<-3000

x1 <-rnorm(n)
x2 <- rnorm(n)
x3 <- rbinom(n,1,0.3)


p1 <- exp(0.3+x1+x2+x3)/(exp(0.3+x1+x2+x3)+exp(0.8+2*x1+2*x2+2*x3)+exp(0.1+3*x1+3*x2+3*x3))
p2 <- exp(0.8+2*x1+2*x2+2*x3)/(exp(0.3+x1+x2+x3)+exp(0.8+2*x1+2*x2+2*x3)+exp(0.1+3*x1+3*x2+3*x3))
p3 <- exp(0.1+3*x1+3*x2+3*x3)/(exp(0.3+x1+x2+x3)+exp(0.8+2*x1+2*x2+2*x3)+exp(0.1+3*x1+3*x2+3*x3))

  
  A <- 1:n;
  for (i in 1:n){ 
    p<- c(p1[i],p2[i],p3[i])
    A[i]<- sample(1:3, 1, replace=T, prob=p)
  }


wheat <- data.frame(A,x1,x2,x3)


library(nnet)

mod.fitr<-multinom(A~ x1 + x2 + x3, data = wheat)
summary(mod.fitr)

library(broom)
g<-tidy(mod.fitr) 


a1<-g[[3]][[1]]
a2<-g[[3]][[2]]
a3<-g[[3]][[3]]
a4<-g[[3]][[4]]
a5<-g[[3]][[5]]
a6<-g[[3]][[6]]
a7<-g[[3]][[7]]
a8<-g[[3]][[8]]


XC1<-rep(0,times=n)
XC2<-rep(0,times=n)
phat1<-rep(0,times=n)
phat2<-rep(0,times=n)
phat3<-rep(0,times=n)

err <- rnorm(n)

y <- 1:n;
ave1 <- rep(0,times=n)
ave2 <- rep(0,times=n)
ave3 <- rep(0,times=n)
r1<-rep(0,times=n)
r2<-rep(0,times=n)
r3<-rep(0,times=n)

z0 <- rep(0,times=n)
z1 <- rep(0,times=n)
z2 <- rep(0,times=n)
z3 <- rep(0,times=n)

X <- model.matrix(~x1+x2+x3)
C1 <- matrix(c(a1,a2,a3,a4),ncol=1)
C2 <- matrix(c(a5,a6,a7,a8),ncol=1)

for (i in 1:n){


XC1[i] <- exp(X[i,]%*%C1)
XC2[i] <- exp(X[i,]%*%C2)


phat1[i] <- 1/(1+XC1[i]+XC2[i])
phat2[i] <- (XC1[i])*phat1[i]
phat3[i] <- (XC2[i])*phat1[i]



I<-rmultinom(1,1,c(phat1[i],phat2[i],phat3[i]))

alpha <- c(0,1,2)
beta1 <- c(1,0.5,1)
beta2 <- c(2,1,2)
beta3 <- c(3,2,1)

y[i] <- t(I)%*%alpha+t(I)%*%beta1*x1[i]+t(I)%*%beta2*x2[i]+t(I)%*%beta3*x3[i]+err[i]

z0[i] <- t(I)%*%alpha
z1[i] <- t(I)%*%beta1*x1[i]
z2[i] <- t(I)%*%beta2*x2[i]
z3[i] <- t(I)%*%beta3*x3[i]

if (I[1,]==1){
   r1[i] <- phat1[i]^(-1)
   ave1[i] <- (y[i]/phat1[i])

}
else if (I[2,]==1){
   r2[i] <- phat2[i]^(-1)
   ave2[i] <- (y[i]/phat2[i])

}
else if (I[3,]==1){
   r3[i] <- phat3[i]^(-1)
   ave3[i] <- (y[i]/phat3[i])

}


}

#####ANCOVA
data <- data.frame(y,z0,z1,z2,z3)
data$z0 <- factor(data$z0)

lm.out <- lm(y~z0+z1+z2+z3, data=data)
anova(lm.out)

summary(lm.out)


w<-tidy(lm.out)

b1<-w[[2]][[1]]
b2<-w[[2]][[2]]
b3<-w[[2]][[3]]
b4<-w[[2]][[4]]
b5<-w[[2]][[5]]
b6<-w[[2]][[6]]


mean(data$z1)
mean(data$z2)
mean(data$z3)

y1hat <- (b1)+(b4)*mean(data$z1)+(b5)*mean(data$z2)+(b6)*mean(data$z3)
y2hat <- (b1)+(b4)*mean(data$z1)+(b5)*mean(data$z2)+(b6)*mean(data$z3)+(b2)
y3hat <- (b1)+(b4)*mean(data$z1)+(b5)*mean(data$z2)+(b6)*mean(data$z3)+(b3)


delta21<-mean(y2hat)-mean(y1hat)
delta32<-mean(y3hat)-mean(y2hat)
delta31<-mean(y3hat)-mean(y1hat)


####Inverse weighting
est1 <- sum(ave1)/sum(r1)
est2 <- sum(ave2)/sum(r2)
est3 <- sum(ave3)/sum(r3)

ATE21 <- est2-est1


ATE32 <- est3-est2


ATE31 <- est3-est1

UU<-matrix(c(ATE21,ATE32,ATE31,delta21,delta32,delta31),nrow=1)


return(UU)
}


resultp <- matrix(NA, nrow=m, ncol=6) #creat a matrix to hold outcome
for(i in 1:m) resultp[i,]<-simipw(x1,x2,x3)


mean(resultp[,1]) ##ATE21--IPTW
var(resultp[,1]) ##ATE21--IPTW

mean(resultp[,2]) ##ATE32--IPTW
var(resultp[,2]) ##ATE32--IPTW

mean(resultp[,3]) ##ATE31--IPTW
var(resultp[,3]) ##ATE31--IPTW

mean(resultp[,4]) ##ATE21--ANCOVA
var(resultp[,4]) ##ATE21--ANCOVA

mean(resultp[,5]) ##ATE32--ANCOVA
var(resultp[,5]) ##ATE32--ANCOVA

mean(resultp[,6]) ##ATE31--ANCOVA
var(resultp[,6]) ##ATE31--ANCOVA



mse21<-mean((resultp[,1]-0.7)^2) ##IPTW
mse21

mse32<-mean((resultp[,2]-0.7)^2) ##IPTW
mse32

mse31<-mean((resultp[,3]-1.4)^2) ##IPTW
mse31

amse21<-mean((resultp[,4]-0.7)^2) ##ANCOVA
amse21

amse32<-mean((resultp[,5]-0.7)^2) ##ANCOVA
amse32

amse31<-mean((resultp[,6]-1.4)^2) ##ANCOVA
amse31

