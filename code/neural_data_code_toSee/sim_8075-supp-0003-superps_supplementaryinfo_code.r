############################################################################################
################################# SUPPLEMENTARY MATERIALS FOR ##############################
########################### Should a Propensity Score Model be Super? ######################
################ The Utility of Machine Learning Procedures for Causal Adjustment ##################
############## Authors: Shomoita Alam, Erica E. M. Moodie, David A. Stephens ###############
############################################################################################

##############################################
###############Reading data###################
##############################################

rm(list=ls())

kc.dat <- read.table("KingCounty2001_data.txt",header=T,sep=",")
set.seed(198603)
summary(kc.dat)

head(kc.dat)

##############################################
###########Recoding the covariates############
##############################################

kc.dat$sex_coded[kc.dat$sex=="F"]<- 0
kc.dat$sex_coded[kc.dat$sex=="M"]<- 1

kc.dat$race_coded[kc.dat$race=="white"]<- 0
kc.dat$race_coded[kc.dat$race=="asian"]<- 1
kc.dat$race_coded[kc.dat$race=="hispanic" | kc.dat$race=="other"
                  | kc.dat$race=="black"]<- 3

kc.dat$smoker_coded[kc.dat$smoker=="N"]<- 0
kc.dat$smoker_coded[kc.dat$smoker=="Y"]<- 1

kc.dat$drinker_coded[kc.dat$drinker=="N"]<- 0
kc.dat$drinker_coded[kc.dat$drinker=="Y"]<- 1

kc.dat$parity_coded[kc.dat$parity==0]<- 0
kc.dat$parity_coded[kc.dat$parity==1]<- 1
kc.dat$parity_coded[kc.dat$parity>=2]<- 2

head(kc.dat)

cov.data<- data.frame(sex=as.factor(kc.dat$sex_coded), age=kc.dat$age,
                      race=as.factor(kc.dat$race_coded), parity=as.factor(kc.dat$parity_coded),
                      married=kc.dat$married, welfare=kc.dat$welfare,
                      smoker=as.factor(kc.dat$smoker_coded),
                      drinker=as.factor(kc.dat$drinker_coded), wpre=kc.dat$wpre, edu=kc.dat$edu,
                      firstep=kc.dat$firstep, bwt=kc.dat$bwt)

##############################################
########Generating simulated datasets#########
##############################################

data_final<- list()
ran.seed<- sample(2000:200000, 500, replace=F)
for (i in 1:500){
  set.seed(ran.seed[i])
  boot_ind<- sample(1:2500, 500, replace=T)
  data_final[[i]]<- cov.data[boot_ind,]
  rownames(data_final[[i]])<- 1:500
  
  n <- dim(data_final[[i]])[1]
  outcome.mod <- lm(bwt ~ firstep+sex+race+parity+smoker,data=data_final[[i]])
  summary(outcome.mod)
  data.gen.coefs <- round(coef(outcome.mod))
  data.gen.coefs[2] <- 0
  
  des.X <- model.matrix(outcome.mod)
  
  data_final[[i]]$newY <- des.X %*% data.gen.coefs
  + rnorm(n,0,summary(outcome.mod)$sigma)
  summary(data_final[[i]])
}

save.image("KC_sample.RData")

##############################################
#################End Here#####################
##############################################

##################################################################################
#############Codes for running Super Learner in Compute Canada Clusters###########
##################################################################################

rm(list=ls(all=TRUE))

args=commandArgs();
####use this to divide to different jobs to parallel run ####
start=as.numeric(args[4]);
end=as.numeric(args[5]);

library(twang)
library(SuperLearner)
library(Hmisc)
library(gam)
library(glmnet)
library(earth)
library(caret)
library(randomForest)
library(nnet)
library(rpart)
library(ipred)
library(gbm)

load("KC_sample.RData")

ps<- NULL
for (i in start:end){
  outc_SL<-data_final[[i]][,11] # outc_SL: treatment; data_final: dataset
  expl_SL<-(data_final[[i]][,c(1,3:4,7)]) # explanatory variables
  
  my.library <- c("SL.glm", "SL.glm.interaction", "SL.step", "SL.stepAIC",
                  "SL.step.interaction", "SL.glmnet", "SL.bayesglm", "SL.caret", "SL.earth",
                  "SL.nnet","SL.knn", "SL.randomForest", "SL.rpart", "SL.ipredbagg",
                  "SL.rpartPrune", "SL.gbm", "SL.gam")
  
  cv.fit <- CV.SuperLearner(Y = outc_SL, X = expl_SL, family="binomial",
                            SL.library=my.library, verbose=FALSE, cvControl=list(stratifyCV=TRUE,
                                                                                 shuffle=TRUE,V=10), method ="method.NNLS")
  
  predictions <- cbind(cv.fit$SL.predict,cv.fit$library.predict)
  ps<- rbind(ps, predictions)
}

devout=paste("SLCV","_",start,"_",end,".txt",sep="")
write.table(ps, devout, row.names = F, sep="\t")


########################################################################
#############Codes for running GBM in Compute Canada Clusters###########
########################################################################

rm(list=ls(all=TRUE))

args=commandArgs();
####use this to divide to different jobs to parallel run ####
start=as.numeric(args[4]);
end=as.numeric(args[5]);

library(twang)
library(gbm)

load("KC_sample.RData")

ps<- NULL
for (i in start:end){
  gbm.fit <- ps(firstep~sex+race+parity+smoker, data=data_final[[i]])
  ps.gbm <- gbm.fit$ps$ks.mean.ATE
  ps<- cbind(ps, ps.gbm)
}

devout=paste("GBM_PS","_",start,"_",end,".txt",sep="")
write.table(ps, devout, row.names = F, sep="\t")

####################################################################
################Null treatment with four confounders################
####################################################################

rm(list=ls())

load("KC_sample.RData")

library(Hmisc)
library(tableone)
library(survey)
library(Matching)
library(MatchIt)

#######################################
#####Getting propensity scores, PS#####
#######################################

auc1<- list()
auc2<- list()
auc3<- list()
auc4<- list()
pa1<- list()
pa2<- list()
pa3<- list()
pa4<- list()

#Balance measures: ASAM

smds<-list()

ps.lr.match<- list()
matched.samp<- list()
dropped.match<- list()
smds.Matched<- list()
ps.match.lr<- list()
ps.lr.match2<- list()
matched.samp2<- list()
dropped.match2<- list()
smds.Matched2<- list()
ps.match.lr2<- list()
ps.SL.match<- list()
matched.samp.SL<- list()
dropped.match.SL<- list()
smds.Matched.SL<- list()
ps.match.SL<- list()
ps.GBM.match<- list()
matched.samp.GBM<- list()
dropped.match.GBM<- list()
smds.Matched.GBM<- list()
ps.match.GBM<- list()

ps.IPW.lr<- list()
smds.IPW<- list()
ps.IPW.lr2<- list()
smds.IPW2<- list()
ps.IPW.SL<- list()
smds.IPW.SL<- list()
ps.IPW.GBM<- list()
smds.IPW.GBM<- list()

keep.out.estimates<- list()
keep.out.estimates.adjw<- list()
keep.out.estimates.adjps<- list()
keep.out.estimates.adjps2<- list()
keep.out.estimates.ps.sl<- list()
keep.out.estimates.ps.gbm<- list()
mean.est.match<- list()
se.est.match<- list()
mean.est.match2<- list()
se.est.match2<- list()
mean.est.match.SL<- list()
se.est.match.SL<- list()
mean.est.match.GBM<- list()
se.est.match.GBM<- list()
IPW.est<- list()
IPW.est2<- list()
IPW.est.SL<- list()
IPW.est.GBM<- list()

for (i in 1:500){
  ps.mod<- glm(data_final[[i]][,11] ~ . , data = data_final[[i]][,c(1,3:4,7)],
               family = binomial(link = "logit"))
  ps.mod2<- glm(data_final[[i]][,11] ~ .^2, data = data_final[[i]][,c(1,3:4,7)],
                family = binomial(link = "logit"))
  
  #PS.LR and PS.LR2
  data_final[[i]]$ps.lr <- predict(ps.mod,type="response")
  data_final[[i]]$ps.lr2 <- predict(ps.mod2,type="response")
  
  #PS.SL
  all_SL_pred<- read.table(paste0("SLCV_",i,"_",i,".txt"), header=T)
  data_final[[i]]$ps.sl<- all_SL_pred[,1]
  
  #PS.GBM
  
  all_GBM_pred<- read.table(paste0("GBM_PS_",i,"_",i,".txt"), header=T)
  data_final[[i]]$ps.gbm<- all_GBM_pred[,1]
  
  auc1[[i]]<- as.numeric(somers2(data_final[[i]]$ps.lr, data_final[[i]]$firstep)["C"])
  auc2[[i]]<- as.numeric(somers2(data_final[[i]]$ps.lr2, data_final[[i]]$firstep)["C"])
  auc3[[i]]<- as.numeric(somers2(data_final[[i]]$ps.sl, data_final[[i]]$firstep)["C"])
  auc4[[i]]<- as.numeric(somers2(data_final[[i]]$ps.gbm, data_final[[i]]$firstep)["C"])
  
  pa1[[i]]<- (length(which(data_final[[i]]$ps.lr>=0.5 & data_final[[i]]$firstep==1))
              + length(which(data_final[[i]]$ps.lr<0.5 & data_final[[i]]$firstep==0)))/n
  pa2[[i]]<- (length(which(data_final[[i]]$ps.lr2>=0.5 & data_final[[i]]$firstep==1))
              + length(which(data_final[[i]]$ps.lr2<0.5 & data_final[[i]]$firstep==0)))/n
  pa3[[i]]<- (length(which(data_final[[i]]$ps.sl>=0.5 & data_final[[i]]$firstep==1))
              + length(which(data_final[[i]]$ps.sl<0.5 & data_final[[i]]$firstep==0)))/n
  pa4[[i]]<- (length(which(data_final[[i]]$ps.gbm>=0.5 & data_final[[i]]$firstep==1))
              + length(which(data_final[[i]]$ps.gbm<0.5 & data_final[[i]]$firstep==0)))/n
  
  #Original
  vars <- names(data_final[[i]])[c(1,3:4,7)]
  tabUnmatched <- CreateTableOne(vars = vars, strata = "firstep",
                                 data = data_final[[i]], test = FALSE)
  smds[[i]]<- ExtractSmd(tabUnmatched)*100
  
  #PS.LR.Matching
  ps.lr.match[[i]] <- Match(Y=data_final[[i]]$newY, Tr=data_final[[i]]$firstep,
                            X=data_final[[i]]$ps.lr,estimand="ATE", caliper=0.2, M=1, ties=TRUE,
                            version="standard", replace=TRUE)
  matched.samp[[i]] <- data_final[[i]][c(ps.lr.match[[i]]$index.control,
                                         ps.lr.match[[i]]$index.treated),]
  matched.samp[[i]]$wts <- ps.lr.match[[i]]$weights[c(ps.lr.match[[i]]$index.control,
                                                      ps.lr.match[[i]]$index.treated)]
  
  dropped.match[[i]]<- ps.lr.match[[i]]$ndrops.matches
  
  dim(matched.samp[[i]])
  matched.samp[[i]]$wts[which(matched.samp[[i]]$wts %in% NA)]<- 1
  dim(matched.samp[[i]])
  
  ps.match.lr[[i]] <- svydesign(ids=~0, data=matched.samp[[i]],
                                weights=matched.samp[[i]]$wts)
  tabMatched <- svyCreateTableOne(vars = vars, strata = "firstep",
                                  data = ps.match.lr[[i]], test = FALSE)
  smds.Matched[[i]]<- ExtractSmd(tabMatched)*100
  
  #PS.LR2.Matching
  ps.lr.match2[[i]] <- Match(Y=data_final[[i]]$newY, Tr=data_final[[i]]$firstep,
                             X=data_final[[i]]$ps.lr2,estimand="ATE",
                             caliper=0.2, M=1, ties=TRUE, version="standard", replace=TRUE)
  matched.samp2[[i]] <- data_final[[i]][c(ps.lr.match2[[i]]$index.control,
                                          ps.lr.match2[[i]]$index.treated),]
  matched.samp2[[i]]$wts <- ps.lr.match2[[i]]$weights[c(ps.lr.match2[[i]]$index.control,
                                                        ps.lr.match2[[i]]$index.treated)]
  
  dropped.match2[[i]]<- ps.lr.match2[[i]]$ndrops.matches
  
  dim(matched.samp2[[i]])
  matched.samp2[[i]]$wts[which(matched.samp2[[i]]$wts %in% NA)]<- 1
  dim(matched.samp2[[i]])
  
  ps.match.lr2[[i]] <- svydesign(ids=~0, data=matched.samp2[[i]],
                                 weights=matched.samp2[[i]]$wts)
  tabMatched2 <- svyCreateTableOne(vars = vars, strata = "firstep",
                                   data = ps.match.lr2[[i]], test = FALSE)
  smds.Matched2[[i]]<- ExtractSmd(tabMatched2)*100
  
  #PS.SL.Matching
  ps.SL.match[[i]] <- Match(Y=data_final[[i]]$newY, Tr=data_final[[i]]$firstep,
                            X=data_final[[i]]$ps.sl,estimand="ATE",
                            caliper=0.2, M=1, ties=TRUE, version="standard", replace=TRUE)
  
  matched.samp.SL[[i]] <- data_final[[i]][c(ps.SL.match[[i]]$index.control,
                                            ps.SL.match[[i]]$index.treated),]
  matched.samp.SL[[i]]$wts <- ps.SL.match[[i]]$weights[c(ps.SL.match[[i]]$index.control,
                                                         ps.SL.match[[i]]$index.treated)]
  
  dropped.match.SL[[i]]<- ps.SL.match[[i]]$ndrops.matches
  
  dim(matched.samp.SL[[i]])
  matched.samp.SL[[i]]$wts[which(matched.samp.SL[[i]]$wts %in% NA)]<- 1
  dim(matched.samp.SL[[i]])
  
  ps.match.SL[[i]] <- svydesign(ids=~0, data=matched.samp.SL[[i]],
                                weights=matched.samp.SL[[i]]$wts)
  tabMatched.SL <- svyCreateTableOne(vars = vars, strata = "firstep",
                                     data = ps.match.SL[[i]], test = FALSE)
  smds.Matched.SL[[i]]<- ExtractSmd(tabMatched.SL)*100
  
  #PS.GBM.Matching
  ps.GBM.match[[i]] <- Match(Y=data_final[[i]]$newY, Tr=data_final[[i]]$firstep,
                             X=data_final[[i]]$ps.gbm,estimand="ATE",
                             caliper=0.2, M=1, ties=TRUE, version="standard", replace=TRUE)
  
  matched.samp.GBM[[i]] <- data_final[[i]][c(ps.GBM.match[[i]]$index.control,
                                             ps.GBM.match[[i]]$index.treated),]
  matched.samp.GBM[[i]]$wts<-
    ps.GBM.match[[i]]$weights[c(ps.GBM.match[[i]]$index.control,
                                ps.GBM.match[[i]]$index.treated)]
  
  dropped.match.GBM[[i]]<- ps.GBM.match[[i]]$ndrops.matches
  
  dim(matched.samp.GBM[[i]])
  matched.samp.GBM[[i]]$wts[which(matched.samp.GBM[[i]]$wts %in% NA)]<- 1
  dim(matched.samp.GBM[[i]])
  
  ps.match.GBM[[i]] <- svydesign(ids=~0, data=matched.samp.GBM[[i]],
                                 weights=matched.samp.GBM[[i]]$wts)
  tabMatched.GBM <- svyCreateTableOne(vars = vars, strata = "firstep",
                                      data = ps.match.GBM[[i]], test = FALSE)
  smds.Matched.GBM[[i]]<- ExtractSmd(tabMatched.GBM)*100
  
  #PS.LR.IPW
  data_final[[i]]$ps.lr.weight<- (data_final[[i]]$firstep/data_final[[i]]$ps.lr
                                  + (1-data_final[[i]]$firstep)/(1-data_final[[i]]$ps.lr))
  data_final[[i]]$ps.lr.weight[data_final[[i]]$ps.lr.weight==NaN
                               | data_final[[i]]$ps.lr.weight==Inf]<- NA
  data_final[[i]]$ps.lr.weight[which(is.na(data_final[[i]]$ps.lr.weight))]
  <- max(1/.05,max(data_final[[i]]$ps.lr.weight, na.rm = T))
  ps.IPW.lr[[i]] <- svydesign(ids=~0, data=data_final[[i]],
                              weights=data_final[[i]]$ps.lr.weight)
  IPW.mod<- lm(newY ~ firstep, weights=ps.lr.weight, data=data_final[[i]])
  tabIPW <- svyCreateTableOne(vars = vars, strata = "firstep",
                              data = ps.IPW.lr[[i]], test = FALSE)
  smds.IPW[[i]]<- ExtractSmd(tabIPW)*100
  
  #PS.LR2.IPW
  data_final[[i]]$ps.lr.weight2<- (data_final[[i]]$firstep/data_final[[i]]$ps.lr2
                                   + (1-data_final[[i]]$firstep)/(1-data_final[[i]]$ps.lr2))
  data_final[[i]]$ps.lr.weight2[data_final[[i]]$ps.lr.weight2==NaN
                                | data_final[[i]]$ps.lr.weight2==Inf]<- NA
  data_final[[i]]$ps.lr.weight2[which(is.na(data_final[[i]]$ps.lr.weight2))]
  <- max(1/.05,max(data_final[[i]]$ps.lr.weight2, na.rm = T))
  ps.IPW.lr2[[i]] <- svydesign(ids=~0, data=data_final[[i]],
                               weights=data_final[[i]]$ps.lr.weight2)
  IPW.mod2<- lm(newY ~ firstep, weights=ps.lr.weight2, data=data_final[[i]])
  tabIPW2 <- svyCreateTableOne(vars = vars, strata = "firstep",
                               data = ps.IPW.lr2[[i]], test = FALSE)
  smds.IPW2[[i]]<- ExtractSmd(tabIPW2)*100
  
  #PS.SL.IPW
  data_final[[i]]$ps.SL.weight<- (data_final[[i]]$firstep/data_final[[i]]$ps.sl
                                  + (1-data_final[[i]]$firstep)/(1-data_final[[i]]$ps.sl))
  data_final[[i]]$ps.SL.weight[data_final[[i]]$ps.SL.weight==NaN
                               | data_final[[i]]$ps.SL.weight==Inf]<- NA
  data_final[[i]]$ps.SL.weight[which(is.na(data_final[[i]]$ps.SL.weight))]
  <- max(1/.05,max(data_final[[i]]$ps.SL.weight, na.rm = T))
  ps.IPW.SL[[i]] <- svydesign(ids=~0, data=data_final[[i]],
                              weights=data_final[[i]]$ps.SL.weight)
  IPW.mod.SL<- lm(newY ~ firstep, weights=ps.SL.weight, data=data_final[[i]])
  tabIPW.SL <- svyCreateTableOne(vars = vars, strata = "firstep",
                                 data = ps.IPW.SL[[i]], test = FALSE)
  smds.IPW.SL[[i]]<- ExtractSmd(tabIPW.SL)*100
  
  #PS.GBM.IPW
  data_final[[i]]$ps.GBM.weight<- (data_final[[i]]$firstep/data_final[[i]]$ps.gbm
                                   + (1-data_final[[i]]$firstep)/(1-data_final[[i]]$ps.gbm))
  data_final[[i]]$ps.GBM.weight[data_final[[i]]$ps.GBM.weight==NaN
                                | data_final[[i]]$ps.GBM.weight==Inf]<- NA
  data_final[[i]]$ps.GBM.weight[which(is.na(data_final[[i]]$ps.GBM.weight))]
  <- max(1/.05,max(data_final[[i]]$ps.GBM.weight, na.rm = T))
  ps.IPW.GBM[[i]] <- svydesign(ids=~0, data=data_final[[i]],
                               weights=data_final[[i]]$ps.GBM.weight)
  IPW.mod.GBM<- lm(newY ~ firstep, weights=ps.GBM.weight, data=data_final[[i]])
  tabIPW.GBM <- svyCreateTableOne(vars = vars, strata = "firstep",
                                  data = ps.IPW.GBM[[i]], test = FALSE)
  smds.IPW.GBM[[i]]<- ExtractSmd(tabIPW.GBM)*100
  
  #ATE Estimation
  
  #Logistic
  #Naive
  naive.mod<- lm(data_final[[i]]$newY ~ data_final[[i]]$firstep)
  keep.out.estimates[[i]]<- summary(naive.mod)$coef[2, c(1,2)]
  
  #Adjusted with Ws
  adj.w.mod<- lm(data_final[[i]]$newY ~ data_final[[i]]$firstep + data_final[[i]][,1]
                 + data_final[[i]][,3] + data_final[[i]][,4] + data_final[[i]][,7])
  keep.out.estimates.adjw[[i]]<- summary(adj.w.mod)$coef[2, c(1,2)]
  
  #Adjusted with PS.LR
  adj.ps.lr.mod<- lm(data_final[[i]]$newY ~ data_final[[i]]$firstep +
                       data_final[[i]]$ps.lr)
  keep.out.estimates.adjps[[i]]<- summary(adj.ps.lr.mod)$coef[2, c(1,2)]
  
  #Adjusted with PS.LR2
  adj.ps.lr2.mod<- lm(data_final[[i]]$newY ~ data_final[[i]]$firstep +
                        data_final[[i]]$ps.lr2)
  keep.out.estimates.adjps2[[i]]<- summary(adj.ps.lr2.mod)$coef[2, c(1,2)]
  
  #Adjusted by PS.SL
  
  adjps.sl.mod<- lm(data_final[[i]]$newY ~ data_final[[i]]$firstep +
                      data_final[[i]]$ps.sl)
  keep.out.estimates.ps.sl[[i]]<- summary(adjps.sl.mod)$coef[2, c(1,2)]
  
  #Adjusted by PS.GBM
  
  adjps.gbm.mod<- lm(data_final[[i]]$newY ~ data_final[[i]]$firstep +
                       data_final[[i]]$ps.gbm)
  keep.out.estimates.ps.gbm[[i]]<- summary(adjps.gbm.mod)$coef[2, c(1,2)]
  
  #PS.LR.Matching Model
  
  mean.est.match[[i]]<- ps.lr.match[[i]]$est
  se.est.match[[i]]<- ps.lr.match[[i]]$se
  
  #PS.LR2.Matching Model
  
  mean.est.match2[[i]]<- ps.lr.match2[[i]]$est
  se.est.match2[[i]]<- ps.lr.match2[[i]]$se
  
  #PS.SL.Matching.Model
  
  mean.est.match.SL[[i]]<- ps.SL.match[[i]]$est
  se.est.match.SL[[i]]<- ps.SL.match[[i]]$se
  
  #PS.GBM.Matching.Model
  
  mean.est.match.GBM[[i]]<- ps.GBM.match[[i]]$est
  se.est.match.GBM[[i]]<- ps.GBM.match[[i]]$se
  
  #PS.LR.IPW.Model
  
  IPW.est[[i]]<- summary(IPW.mod)$coef[2, c(1:2)]
  
  #PS.LR2.IPW.Model
  
  IPW.est2[[i]]<- summary(IPW.mod2)$coef[2, c(1:2)]
  
  #PS.SL.IPW.Model
  
  IPW.est.SL[[i]]<- summary(IPW.mod.SL)$coef[2, c(1:2)]
  
  #PS.GBM.IPW.Model
  
  IPW.est.GBM[[i]]<- summary(IPW.mod.GBM)$coef[2, c(1:2)]
  
  print(paste("Simulation Complete: ",i))
}

save.image("Loop_ends_all.RData")


auc.ps.lr<- do.call(rbind, auc1); dim(auc.ps.lr)
auc.ps.lr2<- do.call(rbind, auc2); dim(auc.ps.lr2)
auc.ps.sl<- do.call(rbind, auc3); dim(auc.ps.sl)
auc.ps.gbm<- do.call(rbind, auc4); dim(auc.ps.gbm)

pa.ps.lr<- do.call(rbind, pa1); dim(pa.ps.lr)
pa.ps.lr2<- do.call(rbind, pa2); dim(pa.ps.lr2)
pa.ps.sl<- do.call(rbind, pa3); dim(pa.ps.sl)
pa.ps.gbm<- do.call(rbind, pa4); dim(pa.ps.gbm)

Average_AUC_PS1<- mean(auc.ps.lr); Average_AUC_PS1
Average_AUC_PS2<- mean(auc.ps.lr2); Average_AUC_PS2
Average_AUC_SL<- mean(auc.ps.sl); Average_AUC_SL
Average_AUC_GBM<- mean(auc.ps.gbm); Average_AUC_GBM
Average_PA_PS1<- mean(pa.ps.lr); Average_PA_PS1
Average_PA_PS2<- mean(pa.ps.lr2); Average_PA_PS2
Average_PA_SL<- mean(pa.ps.sl); Average_PA_SL
Average_PA_GBM<- mean(pa.ps.gbm); Average_PA_GBM

AUCs.A1<- cbind(Average_AUC_PS1, Average_AUC_PS2, Average_AUC_SL, Average_AUC_GBM)
AUCs.A1

PAs.A1<- cbind(Average_PA_PS1, Average_PA_PS2, Average_PA_SL, Average_PA_GBM)
PAs.A1

#Original.ASAM
smds.r<- as.data.frame(do.call(rbind, smds))
smds.r$Mean<- rowMeans(smds.r, na.rm = T)
head(smds.r)
mean(smds.r$Mean)

#PS.LR.ASAM
mean(do.call(rbind,dropped.match))

smds.r.match<- as.data.frame(do.call(rbind, smds.Matched))
smds.r.match$Mean<- rowMeans(smds.r.match, na.rm = T)
head(smds.r.match)
mean(smds.r.match$Mean)

#PS.LR2.ASAM
mean(do.call(rbind,dropped.match2))

smds.r.match2<- as.data.frame(do.call(rbind, smds.Matched2))
smds.r.match2$Mean<- rowMeans(smds.r.match2, na.rm = T)
head(smds.r.match2)
mean(smds.r.match2$Mean)

#SL.Match.ASAM
mean(do.call(rbind,dropped.match.SL))

smds.r.match.SL<- as.data.frame(do.call(rbind, smds.Matched.SL))
smds.r.match.SL$Mean<- rowMeans(smds.r.match.SL, na.rm = T)
head(smds.r.match.SL)
mean(smds.r.match.SL$Mean)

#GBM.Match.ASAM
mean(do.call(rbind,dropped.match.GBM))

smds.r.match.GBM<- as.data.frame(do.call(rbind, smds.Matched.GBM))
smds.r.match.GBM$Mean<- rowMeans(smds.r.match.GBM, na.rm = T)
head(smds.r.match.GBM)
mean(smds.r.match.GBM$Mean)

#PS.LR.IPTW

smds.r.IPW<- as.data.frame(do.call(rbind, smds.IPW))
smds.r.IPW$Mean<- rowMeans(smds.r.IPW, na.rm = T)
head(smds.r.IPW)
mean(smds.r.IPW$Mean)

#PS.LR2.IPTW
smds.r.IPW2<- as.data.frame(do.call(rbind, smds.IPW2))
smds.r.IPW2$Mean<- rowMeans(smds.r.IPW2, na.rm = T)
head(smds.r.IPW2)
mean(smds.r.IPW2$Mean)

#PS.SL.IPW
smds.r.IPW.SL<- as.data.frame(do.call(rbind, smds.IPW.SL))
smds.r.IPW.SL$Mean<- rowMeans(smds.r.IPW.SL, na.rm = T)
head(smds.r.IPW.SL)
mean(smds.r.IPW.SL$Mean)

#PS.GBM.IPW
smds.r.IPW.GBM<- as.data.frame(do.call(rbind, smds.IPW.GBM))
smds.r.IPW.GBM$Mean<- rowMeans(smds.r.IPW.GBM, na.rm = T)
head(smds.r.IPW.GBM)
mean(smds.r.IPW.GBM$Mean)

#ATE estimation

#Naive
mean.sd.A<- do.call(rbind, keep.out.estimates); dim(mean.sd.A)
mse.ps.est<- sqrt(var(mean.sd.A[,1])+ mean(mean.sd.A[,1]-0)^2)

#Adjusted W

mean.sd.A.adjw<- do.call(rbind, keep.out.estimates.adjw); dim(mean.sd.A.adjw)
mse.ps.est.adjw<-  sqrt(var(mean.sd.A.adjw[,1])+ mean(mean.sd.A.adjw[,1]-0)^2)

##Adjusted PS-LR

mean.sd.A.adjps<- do.call(rbind, keep.out.estimates.adjps); dim(mean.sd.A.adjps)
mse.ps.est.adjps<- sqrt(var(mean.sd.A.adjps[,1])+ mean(mean.sd.A.adjps[,1]-0)^2)

##Adjusted PS-LR2

mean.sd.A.adjps2<- do.call(rbind, keep.out.estimates.adjps2); dim(mean.sd.A.adjps2)
mse.ps.est.adjps2<- sqrt(var(mean.sd.A.adjps2[,1])+ mean(mean.sd.A.adjps2[,1]-0)^2)

#Adjusted by PS.SL

mean.sd.A.sl<- do.call(rbind, keep.out.estimates.ps.sl); dim(mean.sd.A.sl)
mse.ps.est.sl<- sqrt(var(mean.sd.A.sl[,1])+ mean(mean.sd.A.sl[,1]-0)^2)

#Adjusted by PS.GBM

mean.sd.A.gbm<- do.call(rbind, keep.out.estimates.ps.gbm); dim(mean.sd.A.gbm)
mse.ps.est.gbm<- sqrt(var(mean.sd.A.gbm[,1])+ mean(mean.sd.A.gbm[,1]-0)^2)

#PS.LR.Matching

mean.A.match<- do.call(rbind, mean.est.match); dim(mean.A.match)
sd.A.match<- do.call(rbind, se.est.match); dim(sd.A.match)
mse.est.match<- sqrt(var(mean.A.match[,1])+ mean(mean.A.match[,1]-0)^2)

#PS.LR2.Matching

mean.A.match2<- do.call(rbind, mean.est.match2); dim(mean.A.match2)
sd.A.match2<- do.call(rbind, se.est.match2); dim(sd.A.match2)
mse.est.match2<- sqrt(var(mean.A.match2[,1])+ mean(mean.A.match2[,1]-0)^2)

#PS.SL.Matching

mean.A.match.SL<- do.call(rbind, mean.est.match.SL); dim(mean.A.match.SL)
sd.A.match.SL<- do.call(rbind, se.est.match.SL); dim(sd.A.match.SL)
mse.est.match.SL<- sqrt(var(mean.A.match.SL[,1])+ mean(mean.A.match.SL[,1]-0)^2)

#PS.GBM.Matching

mean.A.match.GBM<- do.call(rbind, mean.est.match.GBM); dim(mean.A.match.GBM)
sd.A.match.GBM<- do.call(rbind, se.est.match.GBM); dim(sd.A.match.GBM)
mse.est.match.GBM<- sqrt(var(mean.A.match.GBM[,1])+ mean(mean.A.match.GBM[,1]-0)^2)

#PS.LR.IPTW

IPW.est.all<- do.call(rbind, IPW.est)
mse.IPW.all<- sqrt(var(IPW.est.all[,1])+ mean(IPW.est.all[,1]-0)^2)

#PS.LR2.IPTW

IPW.est.all2<- do.call(rbind, IPW.est2)
mse.IPW.all2<- sqrt(var(IPW.est.all2[,1])+ mean(IPW.est.all2[,1]-0)^2)

#PS.SL.IPW

IPW.est.all.SL<- do.call(rbind, IPW.est.SL)
mse.IPW.all.SL<- sqrt(var(IPW.est.all.SL[,1])+ mean(IPW.est.all.SL[,1]-0)^2)

#PS.GBM.IPW

IPW.est.all.GBM<- do.call(rbind, IPW.est.GBM)
mse.IPW.all.GBM<- sqrt(var(IPW.est.all.GBM[,1])+ mean(IPW.est.all.GBM[,1]-0)^2)

###############Combining results#############

#Naive

Estimate<- mean(mean.sd.A[,1])
Absolute_Bias<- abs(Estimate-(0))
Empirical_SE<- sd(mean.sd.A[,1])
rMSE<- mse.ps.est
ASAM<- mean(smds.r$Mean)
Discarded<- NA

Naive<- cbind(Estimate, Absolute_Bias, Empirical_SE,  rMSE, ASAM, Discarded)
Naive

#Adjusted by Ws

Estimate<- mean(mean.sd.A.adjw[,1])
Absolute_Bias<- abs(Estimate-(0))
Empirical_SE<- sd(mean.sd.A.adjw[,1])
rMSE<- mse.ps.est.adjw
ASAM<- mean(smds.r$Mean)
Discarded<- NA

Adjusted.W<- cbind(Estimate, Absolute_Bias, Empirical_SE,  rMSE, ASAM, Discarded)
Adjusted.W

#Adjusted by PS-LR

Estimate<- mean(mean.sd.A.adjps[,1])
Absolute_Bias<- abs(Estimate-(0))
Empirical_SE<- sd(mean.sd.A.adjps[,1])
rMSE<- mse.ps.est.adjps
ASAM<- mean(smds.r$Mean)
Discarded<- NA

Adjusted.PS<- cbind(Estimate, Absolute_Bias, Empirical_SE,  rMSE, ASAM, Discarded)
Adjusted.PS

#Adjusted by PS-LR2

Estimate<- mean(mean.sd.A.adjps2[,1])
Absolute_Bias<- abs(Estimate-(0))
Empirical_SE<- sd(mean.sd.A.adjps2[,1])
rMSE<- mse.ps.est.adjps2
ASAM<- mean(smds.r$Mean)
Discarded<- NA

Adjusted.PS2<- cbind(Estimate, Absolute_Bias, Empirical_SE,  rMSE, ASAM, Discarded)
Adjusted.PS2

#Adjusted by PS-SL

Estimate<- mean(mean.sd.A.sl[,1])
Absolute_Bias<- abs(Estimate-(0))
Empirical_SE<- sd(mean.sd.A.sl[,1])
rMSE<- mse.ps.est.sl
ASAM<- mean(smds.r$Mean)
Discarded<- NA

Adjusted.PS.SL<- cbind(Estimate, Absolute_Bias, Empirical_SE,  rMSE, ASAM, Discarded)
Adjusted.PS.SL


#Adjusted by PS-GBM

Estimate<- mean(mean.sd.A.gbm[,1])
Absolute_Bias<- abs(Estimate-(0))
Empirical_SE<- sd(mean.sd.A.gbm[,1])
rMSE<- mse.ps.est.gbm
ASAM<- mean(smds.r$Mean)
Discarded<- NA

Adjusted.PS.GBM<- cbind(Estimate, Absolute_Bias, Empirical_SE,  rMSE, ASAM, Discarded)
Adjusted.PS.GBM

#Logit Matching

Estimate<- mean(mean.A.match[,1])
Absolute_Bias<- abs(Estimate-(0))
Empirical_SE<- sd(mean.A.match[,1])
rMSE<- mse.est.match
ASAM<- mean(smds.r.match$Mean)
Discarded<- mean(do.call(rbind,dropped.match))

Logit_Matching<- cbind(Estimate, Absolute_Bias, Empirical_SE,  rMSE, ASAM, Discarded)
Logit_Matching

#Logit Matching2

Estimate<- mean(mean.A.match2[,1])
Absolute_Bias<- abs(Estimate-(0))
Empirical_SE<- sd(mean.A.match2[,1])
rMSE<- mse.est.match2
ASAM<- mean(smds.r.match2$Mean)
Discarded<- mean(do.call(rbind,dropped.match2))

Logit_Matching2<- cbind(Estimate, Absolute_Bias, Empirical_SE,  rMSE, ASAM, Discarded)
Logit_Matching2

#SL Matching

Estimate<- mean(mean.A.match.SL[,1])
Absolute_Bias<- abs(Estimate-(0))
Empirical_SE<- sd(mean.A.match.SL[,1])
rMSE<- mse.est.match.SL
ASAM<- mean(smds.r.match.SL$Mean)
Discarded<- mean(do.call(rbind,dropped.match.SL))

SL_Matching<- cbind(Estimate, Absolute_Bias, Empirical_SE,  rMSE, ASAM, Discarded)
SL_Matching

#GBM Matching

Estimate<- mean(mean.A.match.GBM[,1])
Absolute_Bias<- abs(Estimate-(0))
Empirical_SE<- sd(mean.A.match.GBM[,1])
rMSE<- mse.est.match.GBM
ASAM<- mean(smds.r.match.GBM$Mean)
Discarded<- mean(do.call(rbind,dropped.match.GBM))

GBM_Matching<- cbind(Estimate, Absolute_Bias, Empirical_SE,  rMSE, ASAM, Discarded)
GBM_Matching

#Logit IPTW

Estimate<- mean(IPW.est.all[,1])
Absolute_Bias<- abs(Estimate-(0))
Empirical_SE<- sd(IPW.est.all[,1])
rMSE<- mse.IPW.all
ASAM<- mean(smds.r.IPW$Mean)
Discarded<- NA

Logit_IPTW<- cbind(Estimate, Absolute_Bias, Empirical_SE,  rMSE, ASAM, Discarded)
Logit_IPTW

#Logit IPTW2

Estimate<- mean(IPW.est.all2[,1])
Absolute_Bias<- abs(Estimate-(0))
Empirical_SE<- sd(IPW.est.all2[,1])
rMSE<- mse.IPW.all2
ASAM<- mean(smds.r.IPW2$Mean)
Discarded<- NA

Logit_IPTW2<- cbind(Estimate, Absolute_Bias, Empirical_SE,  rMSE, ASAM, Discarded)
Logit_IPTW2

#SL IPTW

Estimate<- mean(IPW.est.all.SL[,1])
Absolute_Bias<- abs(Estimate-(0))
Empirical_SE<- sd(IPW.est.all.SL[,1])
rMSE<- mse.IPW.all.SL
ASAM<- mean(smds.r.IPW.SL$Mean)
Discarded<- NA

SL_IPTW<- cbind(Estimate, Absolute_Bias, Empirical_SE,  rMSE, ASAM, Discarded)
SL_IPTW

#GBM IPTW

Estimate<- mean(IPW.est.all.GBM[,1])
Absolute_Bias<- abs(Estimate-(0))
Empirical_SE<- sd(IPW.est.all.GBM[,1])
rMSE<- mse.IPW.all.GBM
ASAM<- mean(smds.r.IPW.GBM$Mean)
Discarded<- NA

GBM_IPTW<- cbind(Estimate, Absolute_Bias, Empirical_SE,  rMSE, ASAM, Discarded)
GBM_IPTW

ScenarioA1.out<- rbind(Naive, Adjusted.W, Adjusted.PS, Adjusted.PS2, Adjusted.PS.SL,
                       Adjusted.PS.GBM, Logit_Matching, Logit_Matching2, SL_Matching, GBM_Matching,
                       Logit_IPTW, Logit_IPTW2, SL_IPTW, GBM_IPTW)

rownames(ScenarioA1.out)<- c("Naive", "Adjusted by Ws", "Adjusted by PS-LR",
                             "Adjusted by PS-LR2",  "Adjusted by PS-SL", "Adjusted by PS-GBM",
                             "Logit PS-LR Matching", "Logit PS-LR2 Matching", "SL Matching", "GBM Matching",
                             "Logit PS-LR IPTW", "Logit PS-LR2 IPTW", "SL IPTW", "GBM IPTW")
ScenarioA1.out

Summ.l<- NULL
for (i in 1:500){
  Summ.l<- rbind(Summ.l, summary(data_final[[i]]$ps.lr.weight))
}
Summ.A1.lr<- colMeans(Summ.l)
Summ.A1.lr

Summ.l2<- NULL
for (i in 1:500){
  Summ.l2<- rbind(Summ.l2, summary(data_final[[i]]$ps.lr.weight2))
}
Summ.A1.lr2<- colMeans(Summ.l2)
Summ.A1.lr2

Summ.sl<- NULL
for (i in 1:500){
  Summ.sl<- rbind(Summ.sl, summary(data_final[[i]]$ps.SL.weight))
}
Summ.A1.sl<- colMeans(Summ.sl)
Summ.A1.sl

Summ.gbm<- NULL
for (i in 1:500){
  Summ.gbm<- rbind(Summ.gbm, summary(data_final[[i]]$ps.GBM.weight))
}
Summ.A1.gbm<- colMeans(Summ.gbm)
Summ.A1.gbm

fin.AUC.PA<- rbind(AUCs.A1, PAs.A1)
fin.weight<- rbind(Summ.A1.lr, Summ.A1.lr2, Summ.A1.sl, Summ.A1.gbm)
rownames(fin.weight)<- c("PS-LR", "PS-LR2", "PS-SL", "PS-GBM")

save.image("NonBootResults_all.RData")

##############################################
#################End Here#####################
##############################################
