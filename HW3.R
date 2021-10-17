#Install Package
library("survival")
library("survminer")
library(dplyr)

#K-M Function
S1 <- 1-3/24
S2 <- S1*(1-2/20)
S3 <- S2*(1-3/15)
S4 <- S3*(1-3/8)
S5 <- S4*(1-2/5)
S <- c(S1,S2,S3,S4,S5)

#Log-Rank Test Data
d1 <- c(2,1,2,2,1)
d2 <- c(1,1,1,1,1)
d <- d1+d2
n1 <- c(12,10,8,4,2)
n2 <- c(12,10,7,4,3)
n <- n1+n2

#Stamdard Log-Rank Test
T1 <- 0
T2 <- 0
for (i in 1:5) {
  T1 <- T1 + d1[i]-d[i]*n1[i]/n[i]
  T2 <- T2 + (n1[i]*n2[i]*d[i]*(n[i]-d[i]))/(n[i]*n[i]*(n[i]-1))
}
T <- (T1)^2/T2

#Peto Log-Rank Test
T1_P <- 0
T2_P <- 0
for (i in 1:5) {
  T1_P <- T1_P + (d1[i]-d[i]*n1[i]/n[i])*S[i]
  T2_P <- T2_P + ((n1[i]*n2[i]*d[i]*(n[i]-d[i]))/(n[i]*n[i]*(n[i]-1)))*S[i]^2
}
T_P <- (T1_P)^2/T2_P

#Fleming and Harrington Log-Rank Test
T1_FH <- 0
T2_FH <- 0
for (i in 1:5) {
  T1_FH <- T1_FH + (d1[i]-d[i]*n1[i]/n[i])*(1-S[i])
  T2_FH <- T2_FH + ((n1[i]*n2[i]*d[i]*(n[i]-d[i]))/(n[i]*n[i]*(n[i]-1)))*(1-S[i])^2
}
T_FH <- (T1_FH)^2/T2_FH

#Sample Size Calculation
beta <- log(3/2)
p1 <- 1/2
Z_1_Alpha <- 2.57
Z_1_beta <- 2.33
Sample_Size <- (Z_1_Alpha+Z_1_beta)^2/(beta^2*p1*(1-p1))

#Newton Method:
lamba1 <- log(2)/2
lamba2 <- log(2)/3
xn <- 6
for (i in 1:10) {
  xn = xn-(xn-(1-exp(-lamba1*xn))/(2*lamba1)-(1-exp(-lamba2*xn))/(2*lamba2)-73/15)/(1-exp(-lamba1*xn)/2-exp(-lamba2*xn)/2)
}

#Data In Text
setwd('/Users/Mac/Desktop')
DATA = read.table('/Users/Mac/Desktop/bmt.txt',
                  row.name=NULL, col.name = c("Group","t1","t2",
                                              "d1", "d2", "d3","ta","da","tc","dc","tp","dp",
                                              "z1","z2","z3","z4","z5","z6","z7","z8","z9","z10"))

#Survival Function
fit <- survfit(Surv(ta, da)~Group, data = DATA, conf.type = 'log-log')
print(fit)

#Log Rank Test
surv_diff <- survdiff(Surv(ta, da)~Group, data = DATA)

#Survival Function Plot
plot(fit, lty=1:3, col=1:3, xlab="Time", ylab="Kaplan-Meier Estimates",ylim = c(0.7, 1))

#Survival Function Plot
plot(fit, lty=1:3, col=1:3, xlab="Time", ylab="Kaplan-Meier Estimates",
     xlim = c(0, 100),ylim = c(0.7, 1))

#Survival Function of Relapse
fit1 <- survfit(Surv(t2, d2)~Group, data = DATA, conf.type = 'log-log')
print(fit1)

#Log Rank Test
surv_diff_2 <- survdiff(Surv(t2, d2)~Group, data = DATA)

#Survival Function Plot
plot(fit1, lty=1:3, col=1:3, xlab="Time", ylab="Kaplan-Meier Estimates",
     ylim = c(0.2, 1))

#Stratified Log-rank test
surv_diff3 <- survdiff(Surv(t2, d2) ~ Group +strata(da), data =  DATA)
