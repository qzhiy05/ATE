library(nnet)
library(broom)
library(speff2trial)
library(tidyr)

data(ACTG175)

mydata <- ACTG175 %>% drop_na(cd496)
summary(mydata)

attach(mydata)

sd(arms)
sd(cd40)
sd(cd80)
sd(strat)
sd(offtrt)
sd(cd496)
sd(cd820)
sd(days)

mul <-multinom(arms~strat+cd40+cd80+offtrt+cd496+cd820,data=mydata)
g <-tidy(mul)
g

print(g,n=21)


in1<-g[[3]][[1]]
a11<-g[[3]][[2]]
a12<-g[[3]][[3]]
a13<-g[[3]][[4]]
a14<-g[[3]][[5]]
a15<-g[[3]][[6]]
a16<-g[[3]][[7]]
in2<-g[[3]][[8]]
a21<-g[[3]][[9]]
a22<-g[[3]][[10]]
a23<-g[[3]][[11]]
a24<-g[[3]][[12]]
a25<-g[[3]][[13]]
a26<-g[[3]][[14]]
in3<-g[[3]][[15]]
a31<-g[[3]][[16]]
a32<-g[[3]][[17]]
a33<-g[[3]][[18]]
a34<-g[[3]][[19]]
a35<-g[[3]][[20]]
a36<-g[[3]][[21]]

#####IPW
n<-1342
XC1<-rep(0,times=n)
XC2<-rep(0,times=n)
XC3<-rep(0,times=n)
phat1<-rep(0,times=n)
phat2<-rep(0,times=n)
phat3<-rep(0,times=n)
phat4<-rep(0,times=n)


ave1 <- rep(0,times=n)
ave2 <- rep(0,times=n)
ave3 <- rep(0,times=n)
ave4 <- rep(0,times=n)

r1<-rep(0,times=n)
r2<-rep(0,times=n)
r3<-rep(0,times=n)
r4<-rep(0,times=n)

chy<-rep(0,times=n)

T<-mydata[,27]


tr1 <- data.frame(dplyr::filter(mydata,arms==0))
coff1 <- lm(log(days)~strat+cd40+cd80+offtrt+cd496+cd820, data=tr1)
stratco1 <- tidy(coff1)
alpha0 <- stratco1[[2]][[1]]
beta10 <- stratco1[[2]][[2]]
beta20 <- stratco1[[2]][[3]]
beta30 <- stratco1[[2]][[4]]
beta40 <- stratco1[[2]][[5]]
beta50 <- stratco1[[2]][[6]]
beta60 <- stratco1[[2]][[7]]


tr2 <- data.frame(dplyr::filter(mydata,arms==1))
coff2 <- lm(log(days)~strat+cd40+cd80+offtrt+cd496+cd820, data=tr2)
stratco2 <- tidy(coff2)
alpha1 <- stratco2[[2]][[1]]
beta11 <- stratco2[[2]][[2]]
beta21 <- stratco2[[2]][[3]]
beta31 <- stratco2[[2]][[4]]
beta41 <- stratco2[[2]][[5]]
beta51 <- stratco2[[2]][[6]]
beta61 <- stratco2[[2]][[7]]


tr3 <- data.frame(dplyr::filter(mydata,arms==2))
coff3 <- lm(log(days)~strat+cd40+cd80+offtrt+cd496+cd820, data=tr3)
stratco3 <- tidy(coff3)
alpha2 <- stratco3[[2]][[1]]
beta12 <- stratco3[[2]][[2]]
beta22 <- stratco3[[2]][[3]]
beta32 <- stratco3[[2]][[4]]
beta42 <- stratco3[[2]][[5]]
beta52 <- stratco3[[2]][[6]]
beta62 <- stratco3[[2]][[7]]


tr4 <- data.frame(dplyr::filter(mydata,arms==3))
coff4 <- lm(log(days)~strat+cd40+cd80+offtrt+cd496+cd820, data=tr4)
stratco4 <- tidy(coff4)
alpha3 <- stratco4[[2]][[1]]
beta13 <- stratco4[[2]][[2]]
beta23 <- stratco4[[2]][[3]]
beta33 <- stratco4[[2]][[4]]
beta43 <- stratco4[[2]][[5]]
beta53 <- stratco4[[2]][[6]]
beta63 <- stratco4[[2]][[7]]



for (i in 1:n){

X <- model.matrix(~strat+cd40+cd80+offtrt+cd496+cd820,data=mydata)
C1 <- matrix(c(in1,a11,a12,a13,a14,a15,a16),ncol=1)
C2 <- matrix(c(in2,a21,a22,a23,a24,a25,a26),ncol=1)
C3 <- matrix(c(in3,a31,a32,a33,a34,a35,a36),ncol=1)


XC1[i] <- exp(X[i,]%*%C1)
XC2[i] <- exp(X[i,]%*%C2)
XC3[i] <- exp(X[i,]%*%C3)

phat1[i] <- 1/(1+XC1[i]+XC2[i]+XC3[i])
phat2[i] <- (XC1[i])*phat1[i]
phat3[i] <- (XC2[i])*phat1[i]
phat4[i] <- (XC3[i])*phat1[i]

coef0 <- c(alpha0,beta10,beta20,beta30,beta40,beta50,beta60)
coef1 <- c(alpha1,beta11,beta21,beta31,beta41,beta51,beta61)
coef2 <- c(alpha2,beta12,beta22,beta32,beta42,beta52,beta62)
coef3 <- c(alpha3,beta13,beta23,beta33,beta43,beta53,beta63)


if (T[i]==0){
   chy[i] <- X[i,]%*%coef0
   r1[i] <- phat1[i]^(-1)
   ave1[i] <- (chy[i]/phat1[i])
}
else if (T[i]==1){
   chy[i] <- X[i,]%*%coef1
   r2[i] <- phat2[i]^(-1)
   ave2[i] <- (chy[i]/phat2[i])
}
else if (T[i]==2){
   chy[i] <- X[i,]%*%coef2
   r3[i] <- phat3[i]^(-1)
   ave3[i] <- (chy[i]/phat3[i])
}
else if (T[i]==3){
   chy[i] <- X[i,]%*%coef3
   r4[i] <- phat4[i]^(-1)
   ave4[i] <- (chy[i]/phat4[i])
}


}


est1 <- sum(ave1)/sum(r1)
est2 <- sum(ave2)/sum(r2)
est3 <- sum(ave3)/sum(r3)
est4 <- sum(ave4)/sum(r4)

ATE21 <- est2-est1
ATE21

ATE31 <- est3-est1
ATE31

ATE41 <- est4-est1
ATE41

ATE32 <- est3-est2
ATE32

ATE42 <- est4-est2
ATE42

ATE43 <- est4-est3
ATE43


######ANCOVA
arms1 <- factor(mydata$arms)

lm.out <-lm(log(days)~arms1+strat+cd40+cd80+offtrt+cd496+cd820,data=mydata)
anova(lm.out)

summary(lm.out)


w<-tidy(lm.out)
w

b1<-w[[2]][[1]]
b2<-w[[2]][[2]]
b3<-w[[2]][[3]]
b4<-w[[2]][[4]]
b5<-w[[2]][[5]]
b6<-w[[2]][[6]]
b7<-w[[2]][[7]]
b8<-w[[2]][[8]]
b9<-w[[2]][[9]]
b10<-w[[2]][[10]]

y1hat <- (b1)+(b5)*mean(strat)+(b6)*mean(cd40)+(b7)*mean(cd80)+(b8)*mean(offtrt)+(b9)*mean(cd496)+(b10)*mean(cd820)
y2hat <- (b1)+(b5)*mean(strat)+(b6)*mean(cd40)+(b7)*mean(cd80)+(b8)*mean(offtrt)+(b9)*mean(cd496)+(b10)*mean(cd820)+(b2)
y3hat <- (b1)+(b5)*mean(strat)+(b6)*mean(cd40)+(b7)*mean(cd80)+(b8)*mean(offtrt)+(b9)*mean(cd496)+(b10)*mean(cd820)+(b3)
y4hat <- (b1)+(b5)*mean(strat)+(b6)*mean(cd40)+(b7)*mean(cd80)+(b8)*mean(offtrt)+(b9)*mean(cd496)+(b10)*mean(cd820)+(b4)

delta21<-mean(y2hat)-mean(y1hat)
delta31<-mean(y3hat)-mean(y1hat)
delta41<-mean(y4hat)-mean(y1hat)
delta32<-mean(y3hat)-mean(y2hat)
delta42<-mean(y4hat)-mean(y2hat)
delta43<-mean(y4hat)-mean(y3hat)


delta21
delta31
delta41
delta32
delta42
delta43

