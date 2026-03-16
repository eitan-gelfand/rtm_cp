# upload libraries
library(Matrix)
library(lme4)
library("lm.beta")  
library(arm)
library(rsq)
library(sjPlot)
library(sjmisc)
library(emmeans)
center_scale <- function(x) {
  scale(x, scale = FALSE)
}

#call data
data <- read.csv("Final4 Clean for Updating.csv")

data$Range <- as.numeric(data$Range)
data$ACC <- as.numeric(data$ACC)
data$New_t1 <- as.numeric(data$New_t1)
data$New_t2 <- as.numeric(data$New_t2)
data$New_tinf <- as.numeric(data$New_tinf)

#CP subgroup data 
own_CP <- subset(data,(ExperimentName=="Caucasian")&(Range>10)& (Range<50)&(Group==1))
other_CP <- subset(data,(ExperimentName=="Asian")& (Range>10)& (Range<50)&(Group==1))

#TD subgroup data 
own_TD <- subset(data, (ExperimentName=="Caucasian")& (Range>10) & (Range<50) & (Group==2))
other_TD <- subset(data,(ExperimentName=="Asian")& (Range>10)& (Range<50) & (Group==2))

own_L_TD <- subset(own_TD,(CFMTgroup1=="L"))
other_L_TD <- subset(other_TD,(CFMTgroup1=="L"))
own_H_TD <- subset(own_TD,(CFMTgroup1=="H"))
other_H_TD <- subset(other_TD,(CFMTgroup1=="H"))

###############################################################################

#create model for own-race Low CFMT TD

ownGLM_LTD0 <- glmer(ACC ~ center_scale(New_tinf) + center_scale(Range) + (1|Subject),data=own_L_TD,family = binomial(link = "probit"),nAGQ = 0)
ownGLM_LTD1 <- glmer(ACC ~ center_scale(New_tinf) + center_scale(New_t1) + center_scale(Range) +(1|Subject),data=own_L_TD,family = binomial(link = "probit"),nAGQ = 0)


#compare models
anova(ownGLM_LTD0,ownGLM_LTD1,test="Chisq") #reject Null model

#show model summary
summary(ownGLM_LTD0)
summary(ownGLM_LTD1)

#variance explained
rsq(ownGLM_LTD1)

#show beta
own_LTD_z0<-standardize(ownGLM_LTD0)
summary(own_LTD_z0)
own_LTD_z1<-standardize(ownGLM_LTD1)
summary(own_LTD_z1)

#create model for other-race Low CFMT TD

otherGLM_LTD0 <- glmer(ACC ~ center_scale(New_tinf) + center_scale(Range) + (1|Subject),data=other_L_TD,family = binomial(link = "probit"),nAGQ = 0)
otherGLM_LTD1 <- glmer(ACC ~ center_scale(New_tinf) + center_scale(New_t1) + center_scale(Range) +(1|Subject),data=other_L_TD,family = binomial(link = "probit"),nAGQ = 0)

#compare models
anova(otherGLM_LTD0,otherGLM_LTD1,test="Chisq") #reject Null model

#show model summary
summary(otherGLM_LTD0)
summary(otherGLM_LTD1)

#variance explained
rsq(otherGLM_LTD1)

#show beta
other_LTD_z0<-standardize(otherGLM_LTD0)
summary(other_LTD_z0)
other_LTD_z1<-standardize(otherGLM_LTD1)
summary(other_LTD_z1)

#create model for own-race High  CFMT TD

ownGLM_HTD0 <- glmer(ACC ~ center_scale(New_tinf) + center_scale(Range) + (1|Subject),data=own_H_TD,family = binomial(link = "probit"),nAGQ = 0)
ownGLM_HTD1 <- glmer(ACC ~ center_scale(New_tinf) + center_scale(New_t1) + center_scale(Range) +(1|Subject),data=own_H_TD,family = binomial(link = "probit"),nAGQ = 0)


#compare models
anova(ownGLM_HTD0,ownGLM_HTD1,test="Chisq") #reject Null model

#show model summary
summary(ownGLM_HTD0)
summary(ownGLM_HTD1)

#variance explained
rsq(ownGLM_HTD1)

#show beta
own_HTD_z0<-standardize(ownGLM_HTD0)
summary(own_HTD_z0)
own_HTD_z1<-standardize(ownGLM_HTD1)
summary(own_HTD_z1)

#create model for other-race High CFMT TD

otherGLM_HTD0 <- glmer(ACC ~ center_scale(New_tinf) + center_scale(Range) + (1|Subject),data=other_H_TD,family = binomial(link = "probit"),nAGQ = 0)
otherGLM_HTD1 <- glmer(ACC ~ center_scale(New_tinf) + center_scale(New_t1) + center_scale(Range) +(1|Subject),data=other_H_TD,family = binomial(link = "probit"),nAGQ = 0)

#compare models
anova(otherGLM_HTD0,otherGLM_HTD1,test="Chisq") #reject Null model

#show model summary
summary(otherGLM_HTD0)
summary(otherGLM_HTD1)

#variance explained
rsq(otherGLM_HTD1)

#show beta
other_HTD_z0<-standardize(otherGLM_HTD0)
summary(other_HTD_z0)
other_HTD_z1<-standardize(otherGLM_HTD1)
summary(other_HTD_z1)

######################################################################################################

#create model for own CP

ownGLM_CP0 <- glmer(ACC ~ center_scale(New_tinf) + center_scale(Range) + (1|Subject),data=own_CP ,family = binomial(link = "probit"),nAGQ = 0)
ownGLM_CP1 <- glmer(ACC ~ center_scale(New_tinf) + center_scale(New_t1) + center_scale(Range) +(1|Subject),data=own_CP ,family = binomial(link = "probit"),nAGQ = 0)

#compare models
anova(ownGLM_CP0,ownGLM_CP1,test="Chisq") #reject Null model

#show model summary
summary(ownGLM_CP0)
summary(ownGLM_CP1)

#variance explained
rsq(ownGLM_CP1)

#show beta
own_CP_z0<-standardize(ownGLM_CP0)
summary(own_CP_z0)
own_CP_z1<-standardize(ownGLM_CP1)
summary(own_CP_z1)

#create model for other race CP

otherGLM_CP0 <- glmer(ACC ~ center_scale(New_tinf) + center_scale(Range) + (1|Subject),data=other_CP ,family = binomial(link = "probit"),nAGQ = 0)
otherGLM_CP1 <- glmer(ACC ~ center_scale(New_tinf) + center_scale(New_t1) + center_scale(Range) +(1|Subject),data=other_CP,family = binomial(link = "probit"),nAGQ = 0)

#compare models
anova(otherGLM_CP0,otherGLM_CP1,test="Chisq") #reject Null model

#show model summary
summary(otherGLM_CP0)
summary(otherGLM_CP1)

#variance explained
rsq(otherGLM_CP1)

#show beta
other_CP_z0<-standardize(otherGLM_CP0)
summary(other_CP_z0)
other_CP_z1<-standardize(otherGLM_CP1)
summary(other_CP_z1)

###############################################################################################################################################
