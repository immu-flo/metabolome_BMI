#### Import libraries ####
setwd("H:/R/Multicohort/")			     	# Enter location of data file to import
library(openxlsx)
library(ggplot2)
library(survival)
library(moments)
library(haven)
library(glmnet)
library(ggsci)
library(cowplot)
library(survminer)
library(tidyverse)
library(tableone)
library(pROC)
library(haven)
#### Import Data ####
Data <- read.xlsx("Metabolic_BMI_dataset.xlsx") 	# Specify data file
start <- which(names(Data)=="SUGAR_MOS")+1

#Change R-unfriendly metabolite names 
annotations <- read.csv("metabolite_names.csv")
colnames(Data) <- c(colnames(Data[1:(start-1)]),annotations$Name)

#### Data tidy-up ####
Data$pr_dm[is.na(Data$pr_dm)] <-0
Data <- subset(Data,pr_dm!=1)
Data <- subset(Data,Glucose<7.1)
Data <- subset(Data,BMI>18)

MOS <- subset(Data,Cohort=="MOS")
MDC <- subset(Data,Cohort=="MDC")
CIAO <- subset(Data,Cohort=="CIAO")

X_MOS <- data.matrix(scale(MOS[start:ncol(MOS)]))
Y_MOS <- scale(MOS$BMI)

X_MDC <- data.matrix(scale(MDC[start:ncol(MDC)]))
Y_MDC <- scale(MDC$BMI)

X_CIAO <- data.matrix(scale(CIAO[start:ncol(CIAO)]))
Y_CIAO <- scale(CIAO$BMI)

#### Create and export Table 1 ####
myVars <- c("age", "female", "BMI", "Waist", "Glucose", "HDL", "TG","current_smoker")
catVars <- c("female","current_smoker")
tab1 <- tableone::CreateTableOne(data = Data,strata = "Cohort",vars = myVars, factorVars = catVars)
tab1exp <- print(tab1, quote = FALSE, noSpaces = TRUE)
write.csv2(tab1exp,"table1.csv")

#### Ridge regression BMI MOS #####
set.seed(1000)
random_20_percent_of_samples<-rownames(MOS)[sample(1:length(rownames(MOS)),dim(MOS)[1]*0.2)]
MOS_val<-MOS[random_20_percent_of_samples,]
dim(MOS_val)

## Define training set
MOS_train<-MOS[-match(random_20_percent_of_samples,rownames(MOS)),]
dim(MOS_train)

X_MOS_train <- data.matrix(scale(MOS_train[start:ncol(MOS_train)]))
X_MOS_val <- data.matrix(scale(MOS_val[start:ncol(MOS_val)]))
Y_MOS_train <- scale(MOS_train$BMI)
Y_MOS_val <- scale(MOS_val$BMI)

lambdas <- 10^seq(3, -2, by = -.1)

fit <- glmnet(X_MOS_train, Y_MOS_train, alpha = 0, lambda = lambdas)

cv_fit <- cv.glmnet(X_MOS_train, Y_MOS_train, alpha = 0, lambda = lambdas)
opt_lambda <- cv_fit$lambda.min

# Export plot of lambda optimization
pdf("Figure S2.pdf")
plot(cv_fit)
dev.off()

y_predicted_mos_train <-predict(fit, s = opt_lambda, newx = X_MOS_train,type="response")
y_predicted_mos <- predict(fit, s = opt_lambda, newx = X_MOS_val,type="response")

rsq <- function (x, y) cor(x, y) ^ 2
rsq_MOS <- rsq(Y_MOS_val,y_predicted_mos)
rsq_MOS

MOS <- rbind(MOS_train,MOS_val)
MOS <- MOS %>% add_column(predicted_BMI=c(y_predicted_mos_train,y_predicted_mos),.before=10)
MOS <- MOS %>% add_column(delta=MOS$predicted_BMI-scale(MOS$BMI),.before=11)
MOS <- MOS %>% add_column(prediction_group=ifelse(MOS$delta>5/sd(MOS$BMI),"Overestimated",
                                           ifelse(MOS$delta< -5/sd(MOS$BMI),"Underestimated",
                                                  ifelse(scale(MOS$BMI)>(30-mean(MOS$BMI))/sd(MOS$BMI),"Predicted Obesity",  
                                                         ifelse(scale(MOS$BMI)>(25-mean(MOS$BMI))/sd(MOS$BMI),
                                                                "Predicted Overweight","Predicted Normalweight")))),
                          .before=11)

MOS <- add_column(.data=MOS,metabolic_BMI=as.numeric(as.character((MOS$BMI+MOS$delta*sd(MOS$BMI)))),.before=11)
MOS$prediction_group <- factor(MOS$prediction_group,levels=c("Overestimated","Predicted Normalweight",
                                                             "Predicted Overweight","Predicted Obesity",
                                                             "Underestimated"))
# Prediction of OB/OW in MOS validation set
#OB
MOS_val["metabolic_BMI"] <- y_predicted_mos
MOS_NWOB <- subset(MOS_val,BMI >30 | BMI <25)
MOS_NWOB["OB"] <- ifelse(MOS_NWOB$BMI>30,1,0)
roc(MOS_NWOB$OB,MOS_NWOB$metabolic_BMI)
#OW
MOS_val["metabolic_BMI"] <- y_predicted_mos
MOS_NWOW <- subset(MOS_val,BMI <30)
MOS_NWOW["OW"] <- ifelse(MOS_NWOW$BMI>25,1,0)
roc(MOS_NWOW$OW,MOS_NWOW$metabolic_BMI)

#### Identify prediction outliers MOS ####

p_ridge_mos <- ggplot(MOS,aes(x=metabolic_BMI,y=BMI,color=prediction_group))+
geom_point()+
theme_classic()+
scale_color_jama()+
xlab("Metabolic BMI")+
theme(legend.title = element_blank())+
xlim(c(18,40))+
ylim(c(18,50))+
annotate(geom="text", x=35, y=40, label=paste("r2=",round(rsq_MOS,2),sep=""))
  
#### Extract coefficients from Ridge ####
tmp_coeffs <- coef(cv_fit, s = "lambda.min")
ridge_coeff <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
ridge_coeff["direction"] <- ifelse(ridge_coeff$coefficient>0,"Pos","Neg")
ridge_coeff <- subset(ridge_coeff,name!="(Intercept)")
ridge_coeff <- ridge_coeff[order(abs(ridge_coeff$coefficient),decreasing = TRUE),]
ridge_coeff <- ridge_coeff[1:50,]

p_ridge_coeff <- ggplot(ridge_coeff,aes(x=reorder(name,abs(coefficient)),y=coefficient,fill=direction))+
  geom_bar(stat="identity")+
  theme_classic()+
  scale_fill_jama()+
  theme(axis.text.x = element_text(angle = 90),
        legend.title = element_blank())+
  xlab("")+
  coord_flip()

#### Ridge regression BMI validation in MDC ####

y_predicted_MDC <- predict(fit, s = opt_lambda, newx = X_MDC,type="response")

rsq_MDC <- rsq(Y_MDC,y_predicted_MDC)
rsq_MDC

MDC <- MDC %>% add_column(delta=as.numeric(as.character(y_predicted_MDC-scale(MDC$BMI))),.before=10)
MDC <- MDC %>% add_column(metabolic_BMI=MDC$BMI+MDC$delta*sd(MDC$BMI),.before=11)
MDC <- MDC %>% add_column(prediction_group=ifelse(MDC$delta>5/sd(MDC$BMI) ,"Overestimated",
                                                  ifelse(MDC$delta< -5/sd(MDC$BMI),"Underestimated",
                                                         ifelse(scale(MDC$BMI)>(30-mean(MDC$BMI))/sd(MDC$BMI),"Predicted Obesity",  
                                                                ifelse(scale(MDC$BMI)>(25-mean(MDC$BMI))/sd(MDC$BMI),
                                                                       "Predicted Overweight","Predicted Normalweight")))),
                          .before=11)

MDC$prediction_group <- factor(MDC$prediction_group,levels=c("Overestimated","Predicted Normalweight",
                                                             "Predicted Overweight","Predicted Obesity",
                                                             "Underestimated"))

p_ridge_mdc <- ggplot(MDC,aes(x=metabolic_BMI,y=BMI,color=prediction_group))+
  geom_point()+
  theme_classic()+
  scale_color_jama()+
  theme(legend.title = element_blank())+
  xlab("Metabolic BMI")+
  xlim(c(18,40))+
  ylim(c(18,50))+
  annotate(geom="text", x=35, y=40, label=paste("r2=",round(rsq_MDC,2),sep=""))

#### Ridge regression BMI validation in CIAO ####

y_predicted_CIAO <- predict(fit, s = opt_lambda, newx = X_CIAO,type="response")

rsq_CIAO <- rsq(Y_CIAO,y_predicted_CIAO)
rsq_CIAO

CIAO <- CIAO %>% add_column(delta=as.numeric(as.character(y_predicted_CIAO-scale(CIAO$BMI))),.before=10)
CIAO <- CIAO %>% add_column(metabolic_BMI=CIAO$BMI+CIAO$delta*sd(CIAO$BMI),.before=11)
CIAO <- CIAO %>% add_column(prediction_group=ifelse(CIAO$delta>5/sd(CIAO$BMI) ,"Overestimated",
                                                    ifelse(CIAO$delta< -5/sd(CIAO$BMI),"Underestimated",
                                                           ifelse(scale(CIAO$BMI)>(30-mean(CIAO$BMI))/sd(CIAO$BMI),"Predicted Obesity",  
                                                                  ifelse(scale(CIAO$BMI)>(25-mean(CIAO$BMI))/sd(CIAO$BMI),
                                                                         "Predicted Overweight","Predicted Normalweight")))),
                            .before=11)

CIAO$prediction_group <- factor(CIAO$prediction_group,levels=c("Overestimated","Predicted Normalweight",
                                                             "Predicted Overweight","Predicted Obesity",
                                                             "Underestimated"))

p_ridge_ciao <- ggplot(CIAO,aes(x=metabolic_BMI,y=BMI,color=prediction_group))+
  geom_point()+
  theme_classic()+
  scale_color_jama()+
  theme(legend.title = element_blank())+
  xlab("Metabolic BMI")+
  xlim(c(18,40))+
  ylim(c(18,50))+
  annotate(geom="text", x=35, y=40, label=paste("r2=",round(rsq_CIAO,2),sep=""))

#### Print combined plots ####
ridge_plots <- plot_grid(plot_grid(p_ridge_mos,p_ridge_mdc,p_ridge_ciao,ncol = 1,nrow = 3,labels=c("A","B","C"))
          ,p_ridge_coeff,ncol = 2,labels = c("","D"))

ggsave("Figure1.pdf",ridge_plots, height = 9,width = 16)

#### Outlier analysis in MOS ####
## Outlier analysis
OE_MOS <- subset(MOS,prediction_group=="Overestimated")
UE_MOS <- subset(MOS,prediction_group=="Underestimated")
NW_MOS <- subset(MOS,prediction_group=="Predicted Normalweight")
OW_MOS <- subset(MOS,prediction_group=="Predicted Overweight")
OB_MOS <- subset(MOS,prediction_group=="Predicted Obesity")

#Proportion of correct predictions
correct_prediction_MOS <- 1-(nrow(OE_MOS)+nrow(UE_MOS))/nrow(MOS)
proportion_OE_MOS <- nrow(OE_MOS)/nrow(MOS)
proportion_UE_MOS <- nrow(UE_MOS)/nrow(MOS)

# ANOVA risk factors
aov_MOS_BMI <- TukeyHSD(aov(MOS$BMI~MOS$prediction_group))
aov_MOS_metabolicBMI <- TukeyHSD(aov(MOS$metabolic_BMI~MOS$prediction_group))
aov_MOS_HDL <- TukeyHSD(aov(MOS$HDL~MOS$prediction_group))
aov_MOS_Glucose <- TukeyHSD(aov(MOS$Glucose~MOS$prediction_group))
aov_MOS_TG <- TukeyHSD(aov(MOS$TG~MOS$prediction_group))
aov_MOS_Waist <- TukeyHSD(aov(MOS$Waist~MOS$prediction_group))

# Case-rate of categorical variables
case_rate <- function(k)(sum(k==1,na.rm=TRUE)/length(k))

# smoking
case_rate(OE_MOS$current_smoker)
case_rate(NW_MOS$current_smoker)
case_rate(OW_MOS$current_smoker)
case_rate(OW_MOS$current_smoker)
case_rate(OB_MOS$current_smoker)

# Sex
case_rate(OE_MOS$female)
case_rate(NW_MOS$female)
case_rate(UE_MOS$female)
case_rate(OW_MOS$female)
case_rate(OB_MOS$female)

chisq_smok_MOS <- chisq.test(MOS$current_smoker,MOS$prediction_group)
chisq_female_MOS <- chisq.test(MOS$female,MOS$prediction_group)

#Boxplots 
box_BMI_MOS <- ggplot(MOS,aes(x=prediction_group,y=BMI,color=prediction_group))+
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+
  geom_jitter(alpha=0.3,width=0.1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+
  scale_color_jama()+
  ylim(c(18,50))

my_legend <- get_legend(
  box_BMI_MOS +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

box_BMI_MOS <- ggplot(MOS,aes(x=prediction_group,y=BMI,color=prediction_group))+
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+
  geom_jitter(alpha=0.3,width=0.1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+
  scale_color_jama()+
  theme(legend.position = "none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ylim(c(18,50))

box_Glucose_MOS <- ggplot(MOS,aes(x=prediction_group,y=Glucose,color=prediction_group))+
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+
  geom_jitter(alpha=0.3,width=0.1)+
  theme(axis.text.x = element_blank())+
  xlab("")+
  scale_color_jama()+
  theme(legend.position = "none")+
  ylim(c(3.5,7))

box_HDL_MOS <- ggplot(MOS,aes(x=prediction_group,y=HDL,color=prediction_group))+
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+
  geom_jitter(alpha=0.3,width=0.1)+
  theme(axis.text.x = element_blank())+
  xlab("")+
  scale_color_jama()+
  theme(legend.position = "none")+
  ylim(c(0.4,4))

box_Waist_MOS <- ggplot(MOS,aes(x=prediction_group,y=Waist,color=prediction_group))+
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+
  geom_jitter(alpha=0.3,width=0.1)+
  theme(axis.text.x = element_blank())+
  xlab("")+
  scale_color_jama()+
  theme(legend.position = "none")+
  ylim(c(50,150))

box_age_MOS <- ggplot(MOS,aes(x=prediction_group,y=age,color=prediction_group))+
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+
  geom_jitter(alpha=0.3,width=0.1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+
  scale_color_jama()+
  theme(legend.position = "none")

box_TG_MOS <- ggplot(MOS,aes(x=prediction_group,y=TG,color=prediction_group))+
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+
  geom_jitter(alpha=0.3,width=0.1)+
  theme(axis.text.x = element_blank())+
  xlab("")+
  scale_color_jama()+
  theme(legend.position = "none")+
  ylim(c(0.3,5.5))

box_mBMI_MOS <- ggplot(MOS,aes(x=prediction_group,y=metabolic_BMI,color=prediction_group))+
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+
  geom_jitter(alpha=0.3,width=0.1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+
  scale_color_jama()+
  theme(legend.position = "none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ylim(c(18,50))

# Prediction of OW/OB
MOS_NWOB <-  subset(MOS_val,BMI>30 | BMI<25)
MOS_NWOB$OB <- ifelse(MOS_NWOB$BMI>30,1,0)
MOS_NWOW <- subset(MOS_val, BMI<30)
MOS_NWOW$OW <- ifelse(MOS_NWOW$BMI>25,1,0)

roc_MOS_OB <-  roc(MOS_NWOB$OB,MOS_NWOB$metabolic_BMI)
roc_MOS_OW <-  roc(MOS_NWOW$OW,MOS_NWOW$metabolic_BMI)

MOS_diet <- subset(MOS,!is.na(SFA_MOS))

#Diet analyses
aov_MOS_sfa <- TukeyHSD(aov(MOS_diet$SFA_MOS~MOS_diet$prediction_group))
aov_MOS_pufa <- TukeyHSD(aov(MOS_diet$PUFA_MOS~MOS_diet$prediction_group))
aov_MOS_wgra <- TukeyHSD(aov(MOS_diet$WGRA_MOS~MOS_diet$prediction_group))
aov_MOS_fruitveg <- TukeyHSD(aov(MOS_diet$FRVEG_MOS~MOS_diet$prediction_group))
aov_MOS_fish <- TukeyHSD(aov(MOS_diet$FISH_MOS~MOS_diet$prediction_group))
aov_MOS_sugar <- TukeyHSD(aov(MOS_diet$SUGAR_MOS~MOS_diet$prediction_group))
aov_MOS_meat <- TukeyHSD(aov(MOS_diet$MEAT_MOS~MOS_diet$prediction_group))

linreg_MOS_fruitveg <- summary(lm(FRVEG_MOS~age+female+current_smoker+prediction_group,data=MOS_diet))
 
# Physical activity MOS
aov_MOS_pa <- TukeyHSD(aov(MOS$pa_MOS~MOS$prediction_group))
MOS$prediction_group <- factor(MOS$prediction_group,ordered = FALSE)
MOS_diet$prediction_group <- factor(MOS_diet$prediction_group,ordered = FALSE)
summary(lm(pa_MOS~prediction_group+age+female+current_smoker,data=MOS))
summary(lm(pa_MOS~prediction_group+age+female+current_smoker,data=MOS_diet))
summary(lm(pa_MOS~prediction_group+age+female+current_smoker+FRVEG_MOS,data=MOS_diet))

write.xlsx(data.frame(aov_MOS_pa$`MOS$prediction_group`),"aov_MOS_pa.xlsx",overwrite = TRUE)

#### Outlier analysis in CIAO ####

#Proportion of correct predictions
OE_CIAO <- subset(CIAO,prediction_group=="Overestimated")
UE_CIAO <- subset(CIAO,prediction_group=="Underestimated")
NW_CIAO <- subset(CIAO,prediction_group=="Predicted Normalweight")
OW_CIAO <- subset(CIAO,prediction_group=="Predicted Overweight")
OB_CIAO <- subset(CIAO,prediction_group=="Predicted Obesity")
correct_prediction_CIAO <- 1-(nrow(OE_CIAO)+nrow(UE_CIAO))/nrow(CIAO)
proportion_OE_CIAO <- nrow(OE_CIAO)/nrow(CIAO)
proportion_UE_CIAO <- nrow(UE_CIAO)/nrow(CIAO)

#ANOVA risk factors CIAO
aov_CIAO_BMI <- TukeyHSD(aov(CIAO$BMI~CIAO$prediction_group))
aov_CIAO_metabolicBMI <- TukeyHSD(aov(CIAO$metabolic_BMI~CIAO$prediction_group))
aov_CIAO_Glucose <- TukeyHSD(aov(CIAO$Glucose~CIAO$prediction_group))
aov_CIAO_TG <- TukeyHSD(aov(CIAO$TG~CIAO$prediction_group))
aov_CIAO_Waist <- TukeyHSD(aov(CIAO$Waist~CIAO$prediction_group))
aov_CIAO_HDL <- TukeyHSD(aov(CIAO$HDL~CIAO$prediction_group))

# Case-rate of categorical variables
case_rate(OE_CIAO$current_smoker)
case_rate(NW_CIAO$current_smoker)
case_rate(UE_CIAO$current_smoker)
case_rate(OW_CIAO$current_smoker)
case_rate(OB_CIAO$current_smoker)

chisq_smok_CIAO <- chisq.test(CIAO$current_smoker,CIAO$prediction_group)
chisq_sex_CIAO <- chisq.test(CIAO$female,CIAO$prediction_group)

# Boxplots
box_BMI_CIAO <- ggplot(CIAO,aes(x=prediction_group,y=BMI,color=prediction_group))+
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+
  geom_jitter(alpha=0.3,width=0.1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+
  scale_color_jama()+
  theme(legend.position = "none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ylim(c(18,50))

box_Glucose_CIAO <- ggplot(CIAO,aes(x=prediction_group,y=Glucose,color=prediction_group))+
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+
  geom_jitter(alpha=0.3,width=0.1)+
  theme(axis.text.x = element_blank())+
  xlab("")+
  scale_color_jama()+
  theme(legend.position = "none")+
  ylim(c(3.5,7))

box_HDL_CIAO <- ggplot(CIAO,aes(x=prediction_group,y=HDL,color=prediction_group))+
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+
  geom_jitter(alpha=0.3,width=0.1)+
  theme(axis.text.x = element_blank())+
  xlab("")+
  scale_color_jama()+
  theme(legend.position = "none")+
  ylim(c(0.4,4))

box_Waist_CIAO <- ggplot(CIAO,aes(x=prediction_group,y=Waist,color=prediction_group))+
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+
  geom_jitter(alpha=0.3,width=0.1)+
  theme(axis.text.x = element_blank())+
  xlab("")+
  scale_color_jama()+
  theme(legend.position = "none")+
  ylim(c(50,150))

box_age_CIAO <- ggplot(CIAO,aes(x=prediction_group,y=age,color=prediction_group))+
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+
  geom_jitter(alpha=0.3,width=0.1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+
  scale_color_jama()+
  theme(legend.position = "none")

box_TG_CIAO <- ggplot(CIAO,aes(x=prediction_group,y=TG,color=prediction_group))+
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+
  geom_jitter(alpha=0.3,width=0.1)+
  theme(axis.text.x = element_blank())+
  xlab("")+
  scale_color_jama()+
  theme(legend.position = "none")+
  ylim(c(0.3,5.5))

box_mBMI_CIAO <- ggplot(CIAO,aes(x=prediction_group,y=metabolic_BMI,color=prediction_group))+
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+
  geom_jitter(alpha=0.3,width=0.1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+
  scale_color_jama()+
  theme(legend.position = "none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ylim(c(18,50))

# Prediction of OW/OB
CIAO_NWOB <-  subset(CIAO,BMI>30 | BMI<25)
CIAO_NWOB$OB <- ifelse(CIAO_NWOB$BMI>30,1,0)
CIAO_NWOW <- subset(CIAO, BMI<30)
CIAO_NWOW$OW <- ifelse(CIAO_NWOW$BMI>25,1,0)

roc_CIAO_OB <-  roc(CIAO_NWOB$OB,CIAO_NWOB$metabolic_BMI)
roc_CIAO_OW <-  roc(CIAO_NWOW$OW,CIAO_NWOW$metabolic_BMI)

#### Outlier analysis in MDC ####
#Proportions of correct predictions
OE_MDC <- subset(MDC,prediction_group=="Overestimated")
UE_MDC <- subset(MDC,prediction_group=="Underestimated")
NW_MDC <- subset(MDC,prediction_group=="Predicted Normalweight")
OW_MDC <- subset(MDC,prediction_group=="Predicted Overweight")
OB_MDC <- subset(MDC,prediction_group=="Predicted Obesity")

correct_prediction_MDC <- 1-(nrow(OE_MDC)+nrow(UE_MDC))/nrow(MDC)
proportion_OE_MDC <- nrow(OE_MDC)/nrow(MDC)
proportion_UE_MDC <- nrow(UE_MDC)/nrow(MDC)

# ANOVA risk factors MDC
aov_MDC_BMI <- TukeyHSD(aov(MDC$BMI~MDC$prediction_group))
aov_MDC_metabolicBMI <- TukeyHSD(aov(MDC$metabolic_BMI~MDC$prediction_group))
aov_MDC_Glucose <- TukeyHSD(aov(MDC$Glucose~MDC$prediction_group))
aov_MDC_TG <- TukeyHSD(aov(MDC$TG~MDC$prediction_group))
aov_MDC_HDL <- TukeyHSD(aov(MDC$HDL~MDC$prediction_group))
aov_MDC_Waist <- TukeyHSD(aov(MDC$Waist~MDC$prediction_group))

# Prediction cross-sectional Obesity
MDC_NWOB <- subset(MDC,BMI >30 | BMI <25)
MDC_NWOB["OB"] <- ifelse(MDC_NWOB$BMI>30,1,0)
roc_MDC_OB <-  roc(MDC_NWOB$OB,MDC_NWOB$metabolic_BMI)

# Prediction cross-sectional Overweight
MDC_NWOW <- subset(MDC,BMI <30)
MDC_NWOW["OW"] <- ifelse(MDC_NWOW$BMI>25,1,0)
roc_MDC_OW <-  roc(MDC_NWOW$OW,MDC_NWOW$metabolic_BMI)

# Case-rate of categorical variables
# Incident Diabetes
case_rate(OE_MDC$inc_dm_2015)
case_rate(NW_MDC$inc_dm_2015)
case_rate(OW_MDC$inc_dm_2015)
case_rate(OB_MDC$inc_dm_2015)
case_rate(UE_MDC$inc_dm_2015)

# smoking
case_rate(OE_MDC$current_smoker)
case_rate(NW_MDC$current_smoker)
case_rate(UE_MDC$current_smoker)
case_rate(OW_MDC$current_smoker)
case_rate(OB_MDC$current_smoker)

# Sex
case_rate(OE_MDC$female)
case_rate(NW_MDC$female)
case_rate(UE_MDC$female)
case_rate(OW_MDC$female)
case_rate(OB_MDC$female)
# Chisquare-tests of categorical variables

chisq_smok_MDC <- chisq.test(MDC$current_smoker,MDC$prediction_group)
chisq_sex_MDC <- chisq.test(MDC$female,MDC$prediction_group)

# Diet and life-style
MDC_diet <- subset(MDC,!is.na(MDC$SFA_MDC))
MDC_pa <- subset(MDC,!is.na(pa_MDC))

aov_MDC_sfa <- TukeyHSD(aov(MDC_diet$SFA_MDC~MDC_diet$prediction_group))
aov_MDC_pufa <- TukeyHSD(aov(MDC_diet$PUFA_MDC~MDC_diet$prediction_group))
aov_MDC_wgra <- TukeyHSD(aov(MDC_diet$WGRA_MDC~MDC_diet$prediction_group))
aov_MDC_fruitveg <- TukeyHSD(aov(MDC_diet$FRVEG_MDC~MDC_diet$prediction_group))
aov_MDC_fish <- TukeyHSD(aov(MDC_diet$FISH_MDC~MDC_diet$prediction_group))
aov_MDC_sugar <- TukeyHSD(aov(MDC_diet$SUGAR_MDC~MDC_diet$prediction_group))
aov_MDC_meat <- TukeyHSD(aov(MDC_diet$MEAT_MDC~MDC_diet$prediction_group))

aov_MDC_pa <- TukeyHSD(aov(MDC_pa$pa_MDC~MDC_pa$prediction_group))

linreg_MDC_fruitveg <- summary(lm(FRVEG_MDC~age+female+current_smoker+prediction_group,data=MDC_diet))

linreg_MDC_pa <- summary(lm(pa_MDC~prediction_group+age+female+current_smoker,data=MDC,
           subset=!is.na(pa_MDC)))


box_pa <- ggplot(MDC_pa,aes(x=prediction_group,y=pa_MDC,color=prediction_group))+
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+
  geom_jitter(alpha=0.3,width=0.1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+
  scale_color_jama()

# Boxplots

box_BMI_MDC <- ggplot(MDC,aes(x=prediction_group,y=BMI,color=prediction_group))+
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+
  geom_jitter(alpha=0.3,width=0.1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+
  scale_color_jama()+
  ylim(c(18,50))+
  theme(legend.position = "none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


box_Glucose_MDC <- ggplot(MDC,aes(x=prediction_group,y=Glucose,color=prediction_group))+
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+
  geom_jitter(alpha=0.3,width=0.1)+
  theme(axis.text.x = element_blank())+
  xlab("")+
  scale_color_jama()+
  theme(legend.position = "none")+
  ylim(c(3.5,7))

box_HDL_MDC <- ggplot(MDC,aes(x=prediction_group,y=HDL,color=prediction_group))+
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+
  geom_jitter(alpha=0.3,width=0.1)+
  theme(axis.text.x = element_blank())+
  xlab("")+
  scale_color_jama()+
  theme(legend.position = "none")+
  ylim(c(0.4,4))

box_Waist_MDC <- ggplot(MDC,aes(x=prediction_group,y=Waist,color=prediction_group))+
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+
  geom_jitter(alpha=0.3,width=0.1)+
  theme(axis.text.x = element_blank())+
  xlab("")+
  scale_color_jama()+
  theme(legend.position = "none")+
  ylim(c(50,150))

box_bodyfatp_MDC <- ggplot(MDC,aes(x=prediction_group,y=bodyfatp,color=prediction_group))+
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+
  geom_jitter(alpha=0.3,width=0.1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+
  scale_color_jama()+
  theme(legend.position = "none")

box_age_MDC <- ggplot(MDC,aes(x=prediction_group,y=age,color=prediction_group))+
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+
  geom_jitter(alpha=0.3,width=0.1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+
  scale_color_jama()+
  theme(legend.position = "none")

box_TG_MDC <- ggplot(MDC,aes(x=prediction_group,y=TG,color=prediction_group))+
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+
  geom_jitter(alpha=0.3,width=0.1)+
  theme(axis.text.x = element_blank())+
  xlab("")+
  scale_color_jama()+
  theme(legend.position = "none")+
  ylim(c(0.3,5.5))

box_mBMI_MDC <- ggplot(MDC,aes(x=prediction_group,y=metabolic_BMI,color=prediction_group))+
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+
  geom_jitter(alpha=0.3,width=0.1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+
  scale_color_jama()+
  theme(legend.position = "none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ylim(c(18,50))

# Cox regression models
MDC$inc_dm_2017 <- as.numeric(as.character(MDC$inc_dm_2017))
MDC_diet$inc_dm_2017 <- as.numeric(as.character(MDC_diet$inc_dm_2017))


cox_DM_mod1 <- coxph(Surv(fudm_2017,inc_dm_2017)~prediction_group+age+female,data=MDC)
cox_DM_mod2 <- coxph(Surv(fudm_2017,inc_dm_2017)~prediction_group+age+female+Glucose+TG+current_smoker+HDL,data=MDC)
cox_DM_mod3 <- coxph(Surv(fudm_2017,inc_dm_2017)~prediction_group+age+female+Glucose+TG+current_smoker+HDL+
                       FRVEG_MDC+WGRA_MDC+SFA_MDC+pa_MDC,data=MDC_diet)

cox_mort_mod1 <- coxph(Surv(fuend_2017,dead_2017)~prediction_group+age+female,data=MDC)
cox_mort_mod2 <- coxph(Surv(fuend_2017,dead_2017)~prediction_group+age+female+Glucose+TG+current_smoker+HDL,data=MDC)
cox_mort_mod3 <- coxph(Surv(fuend_2017,dead_2017)~prediction_group+age+female+Glucose+TG+current_smoker+HDL+
                         FRVEG_MDC+WGRA_MDC+SFA_MDC+pa_MDC,data=MDC_diet)

#UE as reference group
MDC$prediction_group = factor(MDC$prediction_group,
                              levels=c("Underestimated","Overestimated","Predicted Normalweight",
                                       "Predicted Overweight","Predicted Obesity"))

cox_DM_mod1_UEref <- coxph(Surv(fudm_2017,inc_dm_2017)~prediction_group+age+female,data=MDC)
cox_DM_mod2_UEref <- coxph(Surv(fudm_2017,inc_dm_2017)~prediction_group+age+female+Glucose+TG+current_smoker+HDL,data=MDC)
cox_mort_mod1_UEref <- coxph(Surv(fuend_2017,dead_2017)~prediction_group+age+female,data=MDC)
cox_mort_mod2_UEref <- coxph(Surv(fuend_2017,dead_2017)~prediction_group+age+female+Glucose+TG+current_smoker+HDL,data=MDC)

# Forest plots of incident T2D and all-cause mortality
plotdata_cox <- as.data.frame(rbind(cbind(exp(cox_DM_mod2$coefficients),exp(confint(cox_DM_mod2)),
                                          summary(cox_DM_mod2)$coefficients[,4],rep("T2D",10)),
                                    cbind(exp(cox_mort_mod2$coefficients),exp(confint(cox_mort_mod2)),
                                          summary(cox_mort_mod2)$coefficients[,4],rep("Mortality",10))))
plotdata_cox <- add_column(.data=plotdata_cox,prediction_group=rep(names(cox_DM_mod2$coefficients),2),.before=1)


colnames(plotdata_cox) <- c("prediction_group","HR","lowint","highint","p","Endpoint")
plotdata_cox$prediction_group <- gsub("prediction_group","",plotdata_cox$prediction_group)
plotdata_cox <- plotdata_cox[c(1:4,11:14),]

plotdata_cox <- add_row(.data=plotdata_cox,.before = 5)
plotdata_cox <- add_row(.data=plotdata_cox)
rownames(plotdata_cox) <- 1:nrow(plotdata_cox)
plotdata_cox$prediction_group <- c(plotdata_cox$prediction_group[1:4],"Overestimated",
                                   plotdata_cox$prediction_group[6:9],"Overestimated")

plotdata_cox$HR <- c(as.numeric(as.character(plotdata_cox$HR))[1:4],1,as.numeric(as.character(plotdata_cox$HR))[6:9],1)
plotdata_cox$lowint <- c(as.numeric(as.character(plotdata_cox$lowint))[1:4],1,as.numeric(as.character(plotdata_cox$lowint))[6:9],1)
plotdata_cox$highint <- c(as.numeric(as.character(plotdata_cox$highint))[1:4],1,as.numeric(as.character(plotdata_cox$highint))[6:9],1)
plotdata_cox$p <- c(as.numeric(as.character(plotdata_cox$p))[1:4],1,as.numeric(as.character(plotdata_cox$p))[6:9],1)
plotdata_cox$Endpoint <- c(rep("T2D",5),rep("Mortality",5))

plotdata_cox_T2D <- subset(plotdata_cox,Endpoint=="T2D")
plotdata_cox_Mort <- subset(plotdata_cox,Endpoint=="Mortality")

p_incdm <- ggplot(plotdata_cox_T2D,aes(x=prediction_group,y=HR,min=lowint,max=highint,color=prediction_group))+
  geom_pointrange()+
  theme_classic()+
  scale_color_jama()+
  xlab("")+
  ylab("Hazard ratio")+
  ylim((min(plotdata_cox$lowint)),max(plotdata_cox$highint))+
  geom_hline(yintercept = 1,linetype=2)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

p_mort <- ggplot(plotdata_cox_Mort,aes(x=prediction_group,y=HR,min=lowint,max=highint,color=prediction_group))+
  geom_pointrange()+
  theme_classic()+
  scale_color_jama()+
  xlab("")+
  ylab("Hazard ratio")+
  ylim((min(plotdata_cox$lowint)),max(plotdata_cox$highint))+
  geom_hline(yintercept = 1,linetype=2)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank())

Figure3 <- plot_grid(p_incdm,p_mort,rel_widths = c(1,1.75),labels=c("A","B"))

ggsave("Figure3.pdf",Figure3)

# BMI re-examination
MDC_re <- subset(MDC,!is.na(BMI_AUS))
MDC_re_OE <- subset(MDC_re,prediction_group=="Overestimated")
MDC_re_NW <- subset(MDC_re,prediction_group=="Predicted Normalweight")
MDC_re_OB <- subset(MDC_re,prediction_group=="Predicted Obesity")
MDC_re_OW <- subset(MDC_re,prediction_group=="Predicted Overweight")
MDC_re_UE <- subset(MDC_re,prediction_group=="Underestimated")

BMI_AUS_OE <- ggpaired(MDC_re_OE, cond1 = "BMI", cond2 = "BMI_AUS",fill=pal_jama()(5)[1],
        palette = "jama",line.color = "gray",ylim=c(18,50),ylab="BMI",xlab = "")

BMI_AUS_NW <- ggpaired(MDC_re_NW, cond1 = "BMI", cond2 = "BMI_AUS",fill=pal_jama()(5)[2],
                       palette = "jama",line.color = "gray",ylim=c(18,50),ylab="BMI",xlab = "")

BMI_AUS_OB <- ggpaired(MDC_re_OB, cond1 = "BMI", cond2 = "BMI_AUS",fill=pal_jama()(5)[4],
                       palette = "jama",line.color = "gray",ylim=c(18,50),ylab="BMI",xlab = "")

BMI_AUS_OW <- ggpaired(MDC_re_OW, cond1 = "BMI", cond2 = "BMI_AUS",fill=pal_jama()(5)[3],
                       palette = "jama",line.color = "gray",ylim=c(18,50),ylab="BMI",xlab = "")

BMI_AUS_UE <- ggpaired(MDC_re_UE, cond1 = "BMI", cond2 = "BMI_AUS",fill=pal_jama()(5)[5],
                       palette = "jama",line.color = "gray",ylim=c(18,50),ylab="BMI",xlab = "")


BMI_AUS_box <- plot_grid(BMI_AUS_OE,BMI_AUS_NW,BMI_AUS_OW,BMI_AUS_OB,BMI_AUS_UE,ncol=5,nrow=1,
          labels=c("Overestimated","Predicted Normalweight","Predicted Overweight",
                   "Predicted Obesity","Underestimated"),label_size = 10)
ggsave("FigureS6.pdf",BMI_AUS_box,width = 16,height = 6)


aov_MDC_deltaBMI <- TukeyHSD(aov(MDC_re$BMI_delta~MDC_re$prediction_group))
write.xlsx(cbind(rownames(aov_MDC_deltaBMI$`MDC_re$prediction_group`),
                 aov_MDC_deltaBMI$`MDC_re$prediction_group`),"aov_BMI_delta.xlsx",overwrite = TRUE)

#### Metabolite differences between outlier groups MOS####
anova_NWOE=anova_OWOE=anova_UEOE=anova_UEOB=anova_NWOE_p=anova_OWOE_p=anova_UEOE_p=anova_UEOB_p=NULL
a <- (which(names(MOS)=="SUGAR_MOS")+1):ncol(MOS)
for (i in 1:length(a)){
  anova <- TukeyHSD(aov(MOS[,a[i]]~MOS$prediction_group))
  anova_NWOE[i] <- anova$`MOS$prediction_group`[1]
  anova_OWOE[i] <- anova$`MOS$prediction_group`[3]
  anova_UEOE[i] <- anova$`MOS$prediction_group`[4]
  anova_UEOB[i] <- anova$`MOS$prediction_group`[9]
  anova_NWOE_p[i] <- anova$`MOS$prediction_group`[31]
  anova_OWOE_p[i] <- anova$`MOS$prediction_group`[33]
  anova_UEOE_p[i] <- anova$`MOS$prediction_group`[34]
  anova_UEOB_p[i] <- anova$`MOS$prediction_group`[39]
}
anova_out <- data.frame(names(MOS)[a],
                        anova_NWOE,anova_NWOE_p,ifelse(anova_NWOE_p<0.05/109,1,0),
                        anova_OWOE,anova_OWOE_p,ifelse(anova_OWOE_p<0.05/109,1,0),
                        anova_UEOE,anova_UEOE_p,ifelse(anova_UEOE_p<0.05/109,1,0),
                        anova_UEOB,anova_UEOB_p,ifelse(anova_UEOB_p<0.05/109,1,0))
anova_out <- anova_out %>% rename(Metabolites="names.MOS..a.",NWOE_sig="ifelse.anova_NWOE_p...0.05.109..1..0.",
                                  OWOE_sig="ifelse.anova_OWOE_p...0.05.109..1..0.",UEOE_sig="ifelse.anova_UEOE_p...0.05.109..1..0.",
                                  UEOB_sig="ifelse.anova_UEOB_p...0.05.109..1..0.")

write.xlsx(anova_out,"MOS_ANOVA_metabolites_prediction_groups.xlsx",overwrite = TRUE)

# barplot metabolite difference
selected_metabolites <- subset(anova_out,NWOE_sig==1 | OWOE_sig  ==1 | UEOE_sig==1 | UEOB_sig==1)

selected_MOS <- MOS %>% select(c(selected_metabolites$Metabolites,prediction_group))

#NWOE
anova_NWOE <- subset(anova_out,NWOE_sig==1)
anova_NWOE["direction"] <- ifelse(anova_NWOE$anova_NWOE>0,"Higher in NW","Higher in OE")

NWOE_plot <- ggplot(anova_NWOE,aes(x=anova_NWOE,y=reorder(Metabolites,-anova_NWOE_p),fill=direction))+
  geom_bar(stat="identity")+
  theme_classic()+
  theme(legend.title = element_blank())+
  scale_fill_jama()+
  ylab("Metabolites")+
  xlab("Difference NW-OE")+
  theme(legend.title = element_blank())

#UEOB
anova_UEOB <- subset(anova_out,UEOB_sig==1)
anova_UEOB["direction"] <- ifelse(anova_UEOB$anova_UEOB>0,"Higher in UE","Higher in OB")

UEOB_plot <- ggplot(anova_UEOB,aes(x=anova_UEOB,y=reorder(Metabolites,-anova_UEOB_p),fill=direction))+
  geom_bar(stat="identity")+
  theme_classic()+
  theme(legend.title = element_blank())+
  scale_fill_jama()+
  ylab("Metabolites")+
  xlab("Difference UE-OB")+
  theme(legend.title = element_blank())

#### Metabolite differences between outlier groups MDC####
anova_NWOE_MDC=anova_OWOE_MDC=anova_UEOE_MDC=anova_UEOB_MDC=anova_NWOE_p_MDC=anova_OWOE_p_MDC=anova_UEOE_p_MDC=anova_UEOB_p_MDC=NULL
a <- (which(names(MDC)=="SUGAR_MOS")+1):ncol(MDC)
for (i in 1:length(a)){
  anova_MDC <- TukeyHSD(aov(MDC[,a[i]]~MDC$prediction_group))
  anova_NWOE_MDC[i] <- anova_MDC$`MDC$prediction_group`[1]
  anova_OWOE_MDC[i] <- anova_MDC$`MDC$prediction_group`[3]
  anova_UEOE_MDC[i] <- anova_MDC$`MDC$prediction_group`[4]
  anova_UEOB_MDC[i] <- anova_MDC$`MDC$prediction_group`[9]
  anova_NWOE_p_MDC[i] <- anova_MDC$`MDC$prediction_group`[31]
  anova_OWOE_p_MDC[i] <- anova_MDC$`MDC$prediction_group`[33]
  anova_UEOE_p_MDC[i] <- anova_MDC$`MDC$prediction_group`[34]
  anova_UEOB_p_MDC[i] <- anova_MDC$`MDC$prediction_group`[39]
}
anova_out_MDC <- data.frame(names(MDC)[a],
                            anova_NWOE_MDC,anova_NWOE_p_MDC,ifelse(anova_NWOE_p_MDC<0.05/109,1,0),
                            anova_OWOE_MDC,anova_OWOE_p_MDC,ifelse(anova_OWOE_p_MDC<0.05/109,1,0),
                            anova_UEOE_MDC,anova_UEOE_p_MDC,ifelse(anova_UEOE_p_MDC<0.05/109,1,0),
                            anova_UEOB_MDC,anova_UEOB_p_MDC,ifelse(anova_UEOB_p_MDC<0.05/109,1,0))
anova_out_MDC <- anova_out_MDC %>% rename(Metabolites="names.MDC..a.",NWOE_sig="ifelse.anova_NWOE_p_MDC...0.05.109..1..0.",
                                          OWOE_sig="ifelse.anova_OWOE_p_MDC...0.05.109..1..0.",UEOE_sig="ifelse.anova_UEOE_p_MDC...0.05.109..1..0.",
                                          UEOB_sig="ifelse.anova_UEOB_p_MDC...0.05.109..1..0.")

write.xlsx(anova_out_MDC,"MDC_ANOVA_metabolites_prediction_groups.xlsx",overwrite = TRUE)

# barplot metabolite difference
selected_metabolites <- subset(anova_out_MDC,NWOE_sig==1 | OWOE_sig  ==1 | UEOE_sig==1 | UEOB_sig==1)

selected_MDC <- MDC %>% select(c(selected_metabolites$Metabolites,prediction_group))

#NWOE
anova_NWOE_MDC <- subset(anova_out_MDC,NWOE_sig==1)
anova_NWOE_MDC["direction"] <- ifelse(anova_NWOE_MDC$anova_NWOE>0,"Higher in NW","Higher in OE")

NWOE_plot_MDC <- ggplot(anova_NWOE_MDC,aes(x=anova_NWOE,y=reorder(Metabolites,-anova_NWOE_p_MDC),fill=direction))+
  geom_bar(stat="identity")+
  theme_classic()+
  theme(legend.title = element_blank())+
  scale_fill_jama()+
  ylab("Metabolites")+
  xlab("Difference NW-OE")+
  theme(legend.title = element_blank())

#UEOB
anova_UEOB <- subset(anova_out,UEOB_sig==1)
anova_UEOB["direction"] <- ifelse(anova_UEOB$anova_UEOB>0,"Higher in UE","Higher in OB")

UEOB_plot <- ggplot(anova_UEOB,aes(x=anova_UEOB,y=reorder(Metabolites,-anova_UEOB_p),fill=direction))+
  geom_bar(stat="identity")+
  theme_classic()+
  theme(legend.title = element_blank())+
  scale_fill_jama()+
  ylab("Metabolites")+
  xlab("Difference UE-OB")+
  theme(legend.title = element_blank())


#### Metabolite differences between outlier groups CIAO####
anova_NWOE_CIAO=anova_OWOE_CIAO=anova_UEOE_CIAO=anova_UEOB_CIAO=anova_NWOE_p_CIAO=anova_OWOE_p_CIAO=anova_UEOE_p_CIAO=anova_UEOB_p_CIAO=NULL
a <- (which(names(CIAO)=="SUGAR_MOS")+1):ncol(CIAO)
for (i in 1:length(a)){
  anova_CIAO <- TukeyHSD(aov(CIAO[,a[i]]~CIAO$prediction_group))
  anova_NWOE_CIAO[i] <- anova_CIAO$`CIAO$prediction_group`[1]
  anova_OWOE_CIAO[i] <- anova_CIAO$`CIAO$prediction_group`[3]
  anova_UEOE_CIAO[i] <- anova_CIAO$`CIAO$prediction_group`[4]
  anova_UEOB_CIAO[i] <- anova_CIAO$`CIAO$prediction_group`[9]
  anova_NWOE_p_CIAO[i] <- anova_CIAO$`CIAO$prediction_group`[31]
  anova_OWOE_p_CIAO[i] <- anova_CIAO$`CIAO$prediction_group`[33]
  anova_UEOE_p_CIAO[i] <- anova_CIAO$`CIAO$prediction_group`[34]
  anova_UEOB_p_CIAO[i] <- anova_CIAO$`CIAO$prediction_group`[39]
}
anova_out_CIAO <- data.frame(names(CIAO)[a],
                             anova_NWOE_CIAO,anova_NWOE_p_CIAO,ifelse(anova_NWOE_p_CIAO<0.05/109,1,0),
                             anova_OWOE_CIAO,anova_OWOE_p_CIAO,ifelse(anova_OWOE_p_CIAO<0.05/109,1,0),
                             anova_UEOE_CIAO,anova_UEOE_p_CIAO,ifelse(anova_UEOE_p_CIAO<0.05/109,1,0),
                             anova_UEOB_CIAO,anova_UEOB_p_CIAO,ifelse(anova_UEOB_p_CIAO<0.05/109,1,0))
anova_out_CIAO <- anova_out_CIAO %>% rename(Metabolites="names.CIAO..a.",NWOE_sig="ifelse.anova_NWOE_p_CIAO...0.05.109..1..0.",
                                            OWOE_sig="ifelse.anova_OWOE_p_CIAO...0.05.109..1..0.",UEOE_sig="ifelse.anova_UEOE_p_CIAO...0.05.109..1..0.",
                                            UEOB_sig="ifelse.anova_UEOB_p_CIAO...0.05.109..1..0.")

write.xlsx(anova_out_CIAO,"CIAO_ANOVA_metabolites_prediction_groups.xlsx",overwrite = TRUE)

plotdata_anova_all <- read.xlsx("H:/R/Multicohort/BMI/anova_OENW_allcohorts.xlsx")

OENW_plot_all <- ggplot(plotdata_anova_all,aes(x=anova_NWOE,y=as.factor(Order),fill=Cohort,label=Metabolites))+
  geom_bar(stat="identity")+
  theme(panel.background = element_rect(fill="white"),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        legend.title = element_blank())+
  scale_fill_jama()+
  ylab("Metabolites")+
  xlab("Difference NW-OE")+
  theme(legend.title = element_blank())+
  geom_hline(data=plotdata_anova_all,aes(yintercept=ifelse(Cohort=="CIAO",Order+0.5,0)),color="black",alpha=0.01)+
  geom_text(data=plotdata_anova_all, aes(y = Order, x= -1.5, 
                                         label = ifelse(Cohort=="MDC"
                                                        ,Metabolites,""),angle=0 ), 
            size = 3, show.legend = FALSE, colour = "black")+
  xlim(-2,1)

ggsave("FigureS5.pdf",OENW_plot_all)

#### Boxplots & ANOVA Diet ####
box_Fruitveg_MDC <- ggplot(MDC_diet,aes(x=prediction_group,y=FRVEG_MDC,color=prediction_group))+
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+
  geom_jitter(alpha=0.3,width=0.1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+
  ylab("Fruit and vegetable intake (g/MJ)")+
  ylim(c(0,150))+
  scale_color_jama()+
  theme(legend.position = "none")+
  annotate("text", x = 1.5, y = 100, label = "p<0.001")+
  annotate("segment", x = 1, xend = 2, y = 95, yend = 95, colour = "black")+
  annotate("text", x = 4.5, y = 130, label = "p=0.02")+
  annotate("segment", x = 4, xend = 5, y = 125, yend = 125, colour = "black")+
  theme(legend.position = "none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

box_Fruitveg_MOS <- ggplot(MOS_diet,aes(x=prediction_group,y=FRVEG_MOS,color=prediction_group))+
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+
  geom_jitter(alpha=0.3,width=0.1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+
  ylab("Fruit and vegetable intake (g/MJ)")+
  ylim(c(0,150))+
  scale_color_jama()+
  theme(legend.position = "none")+
  annotate("text", x = 1.5, y = 100, label = "p=0.003")+
  annotate("segment", x = 1, xend = 2, y = 95, yend = 95, colour = "black")+
  theme(legend.position = "none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p_combined_Diet<- plot_grid(plot_grid(box_Fruitveg_MDC,box_Fruitveg_MOS,
                         labels = c("A", "B"),
                         ncol = 2, nrow = 1),my_legend,
                         nrow=2,rel_heights = c(1,0.1))

ggsave("FigureS7.pdf",p_combined_Diet, height = 6,width = 10)

anova_diet_MDC <- rbind(aov_MDC_sfa$`MDC_diet$prediction_group`,aov_MDC_pufa$`MDC_diet$prediction_group`,
      aov_MDC_fruitveg$`MDC_diet$prediction_group`,aov_MDC_fish$`MDC_diet$prediction_group`,
      aov_MDC_wgra$`MDC_diet$prediction_group`,aov_MDC_meat$`MDC_diet$prediction_group`,
      aov_MDC_sugar$`MDC_diet$prediction_group`)

anova_diet_MOS <- rbind(aov_MOS_sfa$`MOS_diet$prediction_group`,aov_MOS_pufa$`MOS_diet$prediction_group`,
                        aov_MOS_fruitveg$`MOS_diet$prediction_group`,aov_MOS_fish$`MOS_diet$prediction_group`,
                        aov_MOS_wgra$`MOS_diet$prediction_group`,aov_MOS_meat$`MOS_diet$prediction_group`,
                        aov_MOS_sugar$`MOS_diet$prediction_group`)

write.csv2(anova_diet_MDC, file="Table_S8_MDC.csv")
write.csv2(anova_diet_MOS, file="Table_S8_MOS.csv")

#### Export ANOVA groupwise comparisons risk factors ####

# MOS
aov_MOS_export <- cbind(rownames(aov_MOS_BMI$`MOS$prediction_group`),aov_MOS_BMI$`MOS$prediction_group`,
                        aov_MOS_Waist$`MOS$prediction_group`,aov_MOS_metabolicBMI$`MOS$prediction_group`,
                        aov_MOS_Glucose$`MOS$prediction_group`,aov_MOS_TG$`MOS$prediction_group`,
                        aov_MOS_HDL$`MOS$prediction_group`)

# MDC
aov_MDC_export <- cbind(rownames(aov_MDC_BMI$`MDC$prediction`),aov_MDC_BMI$`MDC$prediction`,
                        aov_MDC_Waist$`MDC$prediction`,aov_MDC_metabolicBMI$`MDC$prediction`,
                        aov_MDC_Glucose$`MDC$prediction`,aov_MDC_TG$`MDC$prediction`,
                        aov_MDC_HDL$`MDC$prediction`)

# CIAO
aov_CIAO_export <- cbind(rownames(aov_CIAO_BMI$`CIAO$prediction_group`),aov_CIAO_BMI$`CIAO$prediction_group`,
                        aov_CIAO_Waist$`CIAO$prediction_group`,aov_CIAO_metabolicBMI$`CIAO$prediction_group`,
                        aov_CIAO_Glucose$`CIAO$prediction_group`,aov_CIAO_TG$`CIAO$prediction_group`,
                        aov_CIAO_HDL$`CIAO$prediction_group`)

write.xlsx(aov_MOS_export,"ANOVA_MOS_riskfactors.xlsx",overwrite = TRUE)
write.xlsx(aov_MDC_export,"ANOVA_MDC_riskfactors.xlsx",overwrite = TRUE)
write.xlsx(aov_CIAO_export,"ANOVA_CIAO_riskfactors.xlsx",overwrite = TRUE)

#### export plots ####
FigureS3 <- plot_grid(box_Waist_MOS,box_Waist_MDC,box_Waist_CIAO,
                        box_TG_MOS,box_TG_MDC,box_TG_CIAO,
                        nrow=2,ncol=3,labels=c("A","B","C","D","E","F"))
FigureS3 <- plot_grid(FigureS3,my_legend,ncol=1,rel_heights = c(1,0.1))

FigureS4 <- plot_grid(box_HDL_MOS,box_HDL_MDC,box_HDL_CIAO,
                        box_Glucose_MOS,box_Glucose_MDC,box_Glucose_CIAO,
                        nrow=2,ncol=3,labels=c("A","B","C","D","E","F"))
FigureS4 <- plot_grid(FigureS4,my_legend,ncol=1,rel_heights = c(1,0.1))

ggsave("FigureS3.pdf",FigureS3,width = 10,height = 6)
ggsave("FigureS4.pdf",FigureS4,width = 10,height = 6)

plotBMI <- plot_grid(box_BMI_MOS,box_BMI_MDC,box_BMI_CIAO,nrow=1,ncol = 3, labels = c("A","B","C"))
plotmBMI <- plot_grid(box_mBMI_MOS,box_mBMI_MDC,box_mBMI_CIAO,nrow=1,ncol = 3, labels = c("D","E","F"))

Figure2 <- plot_grid(plotBMI,plotmBMI,nrow=2,ncol=1)
Figure2 <- plot_grid(Figure2,my_legend,ncol = 1, rel_heights = c(1, .1))
ggsave("Figure2.pdf",Figure2,height = 6,width = 10)
