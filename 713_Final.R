# Down load R Packages Needed
library(dplyr)
library(PASWR)
library(survival)
library(survminer)
library(ggpubr)
library(ggplot2)

#Import Data
setwd('/Users/Mac/Desktop')
DATA = read.table('/Users/Mac/Desktop/pbc.txt',
      row.name=NULL, col.name = c("id","time","Status",
                                  "Treatment", "Age", "Sex","Ascites","Hepato","Spiders","Edema","Bili","Chol",
                                  "Albumin","Copper","Alk.phos","Ast","Trig","Platelet","Protime","Stage"))

#Data Cleaning
DATA <- DATA %>% filter(id!="NA", time!="NA", Status!="NA", Treatment!="NA", Age!="NA", Sex!="NA", Ascites!="NA",
                        Hepato!="NA", Spiders!="NA", Edema!="NA", Bili!="NA", Chol!="NA", Albumin!="NA", Copper!="NA",
                        Alk.phos!="NA", Ast!="NA", Trig!="NA", Platelet!="NA", Protime!="NA", Stage!="NA")

#Statistical Methods for Purpose 1
DATA_T <- DATA %>% filter(Treatment!="2")
DATA_NT <- DATA %>% filter(Treatment!="1")

x1 <- DATA_T$Age
y1 <- DATA_NT$Age
T_test_1 <- t.test(x1,y1)

x2 <- DATA_T$Bili
y2 <- DATA_NT$Bili
T_test_2 <- t.test(x2,y2)

x3 <- DATA_T$Chol
y3 <- DATA_NT$Chol
T_test_3 <- t.test(x3,y3)

x4 <- DATA_T$Albumin
y4 <- DATA_NT$Albumin
T_test_4 <- t.test(x4,y4)

x5 <- DATA_T$Copper
y5 <- DATA_NT$Copper
T_test_5 <- t.test(x5,y5)

x6 <- DATA_T$Alk.phos
y6 <- DATA_NT$Alk.phos
T_test_6 <- t.test(x6,y6)

x7 <- DATA_T$Ast
y7 <- DATA_NT$Ast
T_test_7 <- t.test(x7,y7)

x8 <- DATA_T$Trig
y8 <- DATA_NT$Trig
T_test_8 <- t.test(x8,y8)

x9 <- DATA_T$Platelet
y9 <- DATA_NT$Platelet
T_test_9 <- t.test(x9,y9)

x10 <- DATA_T$Protime
y10 <- DATA_NT$Protime
T_test_10 <- t.test(x10,y10)

Chi_Square_1 <- chisq.test(DATA$Treatment, DATA$Sex)
Chi_Square_2 <- chisq.test(DATA$Treatment, DATA$Ascites)
Chi_Square_3 <- chisq.test(DATA$Treatment, DATA$Hepato)
Chi_Square_4 <- chisq.test(DATA$Treatment, DATA$Spiders)
Chi_Square_5 <- chisq.test(DATA$Treatment, DATA$Edema)
Chi_Square_6 <- chisq.test(DATA$Treatment, DATA$Stage)

#Statistical Methods For Purpose 2
DATA2 <- DATA
DATA2$Status <- as.integer(DATA$Status == 2)
Survival_All <- survfit(Surv(time, Status)~1, data = DATA2)
plot(Survival_All, lty=1:3, col=1:3, xlab="Time", ylab="Kaplan-Meier Estimator")

#Nelson-Aalen estimate of the cumulative hazard function and its 95% CI
ggsurvplot(Survival_All,
           conf.int = TRUE,
           risk.table.col = "strata", 
           ggtheme = theme_bw(), 
           palette = c("#E7B800", "#2E9FDF"),
           fun = "cumhaz")

#Log-Rank test
Survival_Treatment <- survfit(Surv(time, Status)~Treatment, data = DATA2)
plot(Survival_Treatment, lty=1:3, col=1:3, xlab="Time", ylab="Kaplan-Meier Estimator")
survdiff(Surv(time, Status) ~ Treatment, data = DATA2)

Survival_Sex <- survfit(Surv(time, Status)~Sex, data = DATA2)
plot(Survival_Sex, lty=1:3, col=1:3, xlab="Time", ylab="Kaplan-Meier Estimator")
survdiff(Surv(time, Status) ~ Sex, data = DATA2)

#Cox Model
cox_model1 <- coxph(Surv(time, Status) ~ Treatment + Sex, data = DATA2)
cox_model2 <- coxph(Surv(time, Status) ~ Treatment + Ascites, data = DATA2)
cox_model3 <- coxph(Surv(time, Status) ~ Treatment + Hepato, data = DATA2)
cox_model4 <- coxph(Surv(time, Status) ~ Treatment + Spiders, data = DATA2)
cox_model5 <- coxph(Surv(time, Status) ~ Treatment + factor(Edema), data = DATA2)
cox_model6 <- coxph(Surv(time, Status) ~ Treatment + Bili, data = DATA2)
cox_model7 <- coxph(Surv(time, Status) ~ Treatment + Chol, data = DATA2)
cox_model8 <- coxph(Surv(time, Status) ~ Treatment + Albumin, data = DATA2)
cox_model9 <- coxph(Surv(time, Status) ~ Treatment + Copper, data = DATA2)
cox_model10 <- coxph(Surv(time, Status) ~ Treatment + Alk.phos, data = DATA2)
cox_model11 <- coxph(Surv(time, Status) ~ Treatment + Ast, data = DATA2)
cox_model12 <- coxph(Surv(time, Status) ~ Treatment + Trig, data = DATA2)
cox_model13 <- coxph(Surv(time, Status) ~ Treatment + Platelet, data = DATA2)
cox_model14 <- coxph(Surv(time, Status) ~ Treatment + Protime, data = DATA2)
cox_model15 <- coxph(Surv(time, Status) ~ Treatment + factor(Stage), data = DATA2)

#Proportional Hazard Assumption - g(t) = t
DATA2$V1 <- DATA2$time*DATA2$Ascites
cox_model_t1 <- coxph(Surv(time, Status) ~ Treatment + Ascites + V1, data = DATA2)
DATA2$V2 <- DATA2$time*DATA2$Hepato
cox_model_t2 <- coxph(Surv(time, Status) ~ Treatment + Hepato + V2, data = DATA2)
DATA2$V3 <- DATA2$time*DATA2$Spiders
cox_model_t3 <- coxph(Surv(time, Status) ~ Treatment + Spiders + V3, data = DATA2)
DATA2$V4 <- DATA2$time*DATA2$Edema
cox_model_t4 <- coxph(Surv(time, Status) ~ Treatment + factor(Edema) + V4, data = DATA2)
DATA2$V5 <- DATA2$time*DATA2$Bili
cox_model_t5 <- coxph(Surv(time, Status) ~ Treatment + Bili + V5, data = DATA2)
DATA2$V6 <- DATA2$time*DATA2$Chol
cox_model_t6 <- coxph(Surv(time, Status) ~ Treatment + Chol + V6, data = DATA2)
DATA2$V7 <- DATA2$time*DATA2$Albumin
cox_model_t7 <- coxph(Surv(time, Status) ~ Treatment + Albumin + V7, data = DATA2)
DATA2$V8 <- DATA2$time*DATA2$Copper
cox_model_t8 <- coxph(Surv(time, Status) ~ Treatment + Copper + V8, data = DATA2)
DATA2$V9 <- DATA2$time*DATA2$Alk.phos
cox_model_t9 <- coxph(Surv(time, Status) ~ Treatment + Alk.phos + V9, data = DATA2)
DATA2$V10 <- DATA2$time*DATA2$Ast
cox_model_t10 <- coxph(Surv(time, Status) ~ Treatment + Ast + V10, data = DATA2)
DATA2$V11 <- DATA2$time*DATA2$Trig
cox_model_t11 <- coxph(Surv(time, Status) ~ Treatment + Trig + V11, data = DATA2)
DATA2$V12 <- DATA2$time*DATA2$Platelet
cox_model_t12 <- coxph(Surv(time, Status) ~ Treatment + Platelet + V12, data = DATA2)
DATA2$V13 <- DATA2$time*DATA2$Protime
cox_model_t13 <- coxph(Surv(time, Status) ~ Treatment + Protime + V13, data = DATA2)
DATA2$V14 <- DATA2$time*DATA2$Stage
cox_model_t14 <- coxph(Surv(time, Status) ~ Treatment + factor(Stage) + V14, data = DATA2)

#Proportional Hazard Assumption - g(t) = ln(t)
DATA2$V1 <- log(DATA2$time)*DATA2$Ascites
cox_model_t1 <- coxph(Surv(time, Status) ~ Treatment + Ascites + V1, data = DATA2)
DATA2$V2 <- log(DATA2$time)*DATA2$Hepato
cox_model_t2 <- coxph(Surv(time, Status) ~ Treatment + Hepato + V2, data = DATA2)
DATA2$V3 <- log(DATA2$time)*DATA2$Spiders
cox_model_t3 <- coxph(Surv(time, Status) ~ Treatment + Spiders + V3, data = DATA2)
DATA2$V4 <- log(DATA2$time)*DATA2$Edema
cox_model_t4 <- coxph(Surv(time, Status) ~ Treatment + factor(Edema) + V4, data = DATA2)
DATA2$V5 <- log(DATA2$time)*DATA2$Bili
cox_model_t5 <- coxph(Surv(time, Status) ~ Treatment + Bili + V5, data = DATA2)
DATA2$V6 <- log(DATA2$time)*DATA2$Chol
cox_model_t6 <- coxph(Surv(time, Status) ~ Treatment + Chol + V6, data = DATA2)
DATA2$V7 <- log(DATA2$time)*DATA2$Albumin
cox_model_t7 <- coxph(Surv(time, Status) ~ Treatment + Albumin + V7, data = DATA2)
DATA2$V8 <- log(DATA2$time)*DATA2$Copper
cox_model_t8 <- coxph(Surv(time, Status) ~ Treatment + Copper + V8, data = DATA2)
DATA2$V9 <- log(DATA2$time)*DATA2$Alk.phos
cox_model_t9 <- coxph(Surv(time, Status) ~ Treatment + Alk.phos + V9, data = DATA2)
DATA2$V10 <- log(DATA2$time)*DATA2$Ast
cox_model_t10 <- coxph(Surv(time, Status) ~ Treatment + Ast + V10, data = DATA2)
DATA2$V11 <- log(DATA2$time)*DATA2$Trig
cox_model_t11 <- coxph(Surv(time, Status) ~ Treatment + Trig + V11, data = DATA2)
DATA2$V12 <- log(DATA2$time)*DATA2$Platelet
cox_model_t12 <- coxph(Surv(time, Status) ~ Treatment + Platelet + V12, data = DATA2)
DATA2$V13 <- log(DATA2$time)*DATA2$Protime
cox_model_t13 <- coxph(Surv(time, Status) ~ Treatment + Protime + V13, data = DATA2)
DATA2$V14 <- log(DATA2$time)*DATA2$Stage
cox_model_t14 <- coxph(Surv(time, Status) ~ Treatment + factor(Stage) + V14, data = DATA2)

#Cox Model For all Variables
cox_model_whole <- coxph(Surv(time, Status) ~ Treatment + Ascites + Sex + Age + Hepato + Spiders 
                         + factor(Edema) + Bili + Chol + Albumin + Copper + Alk.phos + Ast + Trig + Platelet
                         + Protime + factor(Stage), data = DATA2)

cox_model_whole_1 <- coxph(Surv(time, Status) ~ Treatment + Sex + Age + Hepato + Spiders 
                         + factor(Edema) + Bili + Chol + Albumin + Copper + Alk.phos + Ast + Trig + Platelet
                         + Protime + factor(Stage), data = DATA2)

cox_model_whole_2 <- coxph(Surv(time, Status) ~ Treatment + Sex + Age + Hepato + Spiders 
                           + factor(Edema) + Bili + Chol + Albumin + Copper + Ast + Trig + Platelet
                           + Protime + factor(Stage), data = DATA2)

cox_model_whole_3 <- coxph(Surv(time, Status) ~ Treatment + Sex + Age + Spiders 
                           + factor(Edema) + Bili + Chol + Albumin + Copper + Ast + Trig + Platelet
                           + Protime + factor(Stage), data = DATA2)

cox_model_whole_4 <- coxph(Surv(time, Status) ~ Treatment + Sex + Age 
                           + factor(Edema) + Bili + Chol + Albumin + Copper + Ast + Trig + Platelet
                           + Protime + factor(Stage), data = DATA2)

cox_model_whole_5 <- coxph(Surv(time, Status) ~ Treatment + Sex + Age 
                           + factor(Edema) + Bili + Chol + Albumin + Copper + Ast + Platelet
                           + Protime + factor(Stage), data = DATA2)

cox_model_whole_6 <- coxph(Surv(time, Status) ~ Treatment + Sex + Age 
                           + factor(Edema) + Bili + Chol + Albumin + Copper + Ast
                           + Protime + factor(Stage), data = DATA2)

cox_model_whole_7 <- coxph(Surv(time, Status) ~ Treatment + Age 
                           + factor(Edema) + Bili + Chol + Albumin + Copper + Ast
                           + Protime + factor(Stage), data = DATA2)

cox_model_whole_8 <- coxph(Surv(time, Status) ~ Treatment + Age 
                           + factor(Edema) + Bili + Albumin + Copper + Ast
                           + Protime + factor(Stage), data = DATA2)

#Assumption for Selected Variables
est_final = cox.zph(cox_model_whole_8, transform = "identity")