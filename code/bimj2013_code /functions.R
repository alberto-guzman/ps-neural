
#########################################################################################

# In this file there are functions:
# funcov and funsim for generating a simulated dataset and estimating on the simulated data, respectively;
# funcov_cs and funsim_cs are the analogous version using real covariates;
# funcov_tp and funsim_opt: utilities for setting tuning parameters.

#########################################################################################

# Main references

# ref 1: Setoguchi et al. Pharmacoepidemiology and Drug Safety (discrete outcome + matching)
# ref 2: Lee, Lessler and Stuart Statist & Med 2010, 29 337--346 (outcome continuo + weighting)

#######   Steps

#1) generate covariate Wi,  i=1, ... 10

#2) generate treatment T and outcome Y

#3) estimate ps

#4) estimate ATT  using inverse probability weighting and ps matching 

#  funcov for steps 1 and 2; funsim for steps 3 and 4


###### install and load packages

list.of.packages = c("Matching","rpart","randomForest","gbm","twang","ipred","neuralnet",
"nnet","e1071","klaR","xtable","flexmix","AUC","Hmisc","Kendall","lattice") # replace xx and yy with package names
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0) {install.packages(new.packages)}
lapply(list.of.packages, require, character.only=T)

############## utiity function: F.sample.cor

#inputs: x, rho
#output: y, i.e., a vector correlated to x by rho 

F.sample.cor <- function(x, rho) {
y <- (rho * (x - mean(x)))/sqrt(var(x)) + sqrt(1 - rho^2) * rnorm(length(x))
#cat("Sample corr = ", cor(x, y), "\n")
return(y)
}

#########################################
############## funzione: funcov (size,scenario)
#########################################

#inputs: sample size n, scenario
#output: a data set of size n in the chosen scenario, containing:
    #covariates w1...w10
    #treatment T
    #outcome Y

funcov<-function(size,scenarioT,scenarioY){

#~~ generate individual level covariates
w1 <- rnorm(size, mean=0, sd=1)
w2 <- rnorm(size, mean=0, sd=1)
w3 <- rnorm(size, mean=0, sd=1)
w4 <- rnorm(size, mean=0, sd=1)
w5 <- F.sample.cor(w1, 0.2)#0.2,0.9
#w5b <- F.sample.cor(w1, 0.0001)
w6 <- F.sample.cor(w2, 0.9)#0.2,0.9
w7 <- rnorm(size, mean=0, sd=1)
w8 <- F.sample.cor(w3, 0.2)
w9 <- F.sample.cor(w4, 0.9)
#w9b <- F.sample.cor(w4, 0.0001)
w10 <- rnorm(size, mean=0, sd=1)

#~~ dichotomize variables (will attenuate correlations above)
w1 <- ifelse(w1 > mean(w1), 1, 0)
w3 <- ifelse(w3 > mean(w3), 1, 0)
w5 <- ifelse(w5 > mean(w5), 1, 0)
#w5b <- ifelse(w5b > mean(w5b), 1, 0)
w6 <- ifelse(w6 > mean(w6), 1, 0)
w8 <- ifelse(w8 > mean(w8), 1, 0)
w9 <- ifelse(w9 > mean(w9), 1, 0)
#w9b <- ifelse(w9 > mean(w9), 1, 0)


#~~ scenarios for data generation models
# A: model with additivity and linearity
# B: mild non-linearity
# C: moderate non-linearity
# D: mild non-additivity
# E: mild non-additivity and non-linearity
# F: moderate non-additivity
# G: moderate non-additivity and non-linearity
# binary exposure modeling


#~~~~~~~~~~~~~~Global Variables~~~~~~~~~~~~~~~~~~~~#

#~~ coefficients for the treatment and potential outcomes equations

#~~ coefficients for treatment equation (from Setoguchi (reference 1))
 b0 <-  0
 b1 <-  0.8
 b2 <- -0.25 
 b3 <-  0.6 
 b4 <- -0.4 
 b5 <- -0.8 
 b6 <- -0.5 
 b7 <-  0.7 
 
 #~~ coefficients for potential outcomes equation Y(0) and Y(1)
 #~~ coefficients for Y(0): from Setoguchi (reference 1)
 #~~ coefficients for Y(1): same as above plus delta. In particular 
 #~~ the intercept in Y(1) has been fixed to obtain Y(1)-Y(0)=-0.4 in each scenario
 a00 <- -3.85   
 a01 <-  0.3     
 a02 <- -0.36   
 a03 <- -0.73  
 a04 <- -0.2   
 a05 <-  0.71   
 a06 <- -0.19 
 a07 <-  0.26 

 #~~ scenarios for treatment assignment


if (scenarioT == "A") {
trueps <- (1 + exp( -(0 + b1*w1 + b2*w2 + b3*w3 + b4*w4 + b5*w5 + b6*w6 + b7*w7) ) ) ^ -1
} else
if (scenarioT == "B") {
trueps <- (1 + exp( -(0 + b1*w1 + b2*w2 + b3*w3 + b4*w4 + b5*w5 + b6*w6 + b7*w7
+ b2*w2*w2) ) )^-1
} else
if (scenarioT == "C") {
trueps <- (1 + exp( -(0 + b1*w1 + b2*w2 + b3*w3 + b4*w4 + b5*w5 + b6*w6 + b7*w7
+ b2*w2*w2 +b4*w4*w4 + b7*w7*w7) ) )^-1
} else 
if (scenarioT == "D") {
trueps <- (1 + exp( -(0 + b1*w1 + b2*w2 + b3*w3 + b4*w4 + b5*w5 + b6*w6 + b7*w7
+ b1*0.5*w1*w3 + b2*0.7*w2*w4 + b4*0.5*w4*w5 + b5*0.5*w5*w6) ) )^-1
} else
if (scenarioT == "E") {
trueps <- (1 + exp( -(0 + b1*w1 + b2*w2 + b3*w3 + b4*w4 + b5*w5 + b6*w6 + b7*w7
+ b2*w2*w2 + b1*0.5*w1*w3 + b2*0.7*w2*w4 + b4*0.5*w4*w5 + b5*0.5*w5*w6) ) )^-1
} else
if (scenarioT == "F") {
trueps <- (1 + exp( -(0 + b1*w1 + b2*w2 + b3*w3 + b4*w4 + b5*w5 + b6*w6 + b7*w7
+ b1*0.5*w1*w3 + b2*0.7*w2*w4 + b3*0.5*w3*w5 + b4*0.7*w4*w6 + b5*0.5*w5*w7
+ b1*0.5*w1*w6 + b2*0.7*w2*w3 + b3*0.5*w3*w4 + b4*0.5*w4*w5 + b5*0.5*w5*w6) ) )^-1
} else

{
# scenario G
trueps <- (1 + exp( -(0 + b1*w1 + b2*w2 + b3*w3 + b4*w4 + b5*w5 + b6*w6 + b7*w7
+ b2*w2*w2 + b4*w4*w4 + b7*w7*w7 + b1*0.5*w1*w3 + b2*0.7*w2*w4 +b3*0.5*w3*w5
+ b4*0.7*w4*w6 + b5*0.5*w5*w7 + b1*0.5*w1*w6 + b2*0.7*w2*w3 + b3*0.5*w3*w4
+ b4*0.5*w4*w5 + b5*0.5*w5*w6) ) )^-1
}

#~~ binary treatment T 
unif1 <- runif(size,0,1)
T     <- ifelse(trueps > unif1, 1, 0)	# there is a probability of unif1 that T=1

#~~ scenarios for outcome

if (scenarioY == "a") {
Y0    <-       a00 +       a01*w1 +       a02*w2 +       a03*w3 +       a04*w4       +
    a05*w8 +       a06*w9 +        a07*w10                
Y1    <- (a00-0.4) + (a01 + 0)*w1 + (a02 + 0)*w2 + (a03 + 0)*w3 + (a04 + 0)*w4 +      (a05+0)*w8 +  (a06+0)*w9 +     (a07+0)*w10
} else
if (scenarioY == "d") {# 1 inter
Y0    <-       a00 +       a01*w1 +       a02*w2 +       a03*w3 +       a04*w4       +
    a05*w8 +       a06*w9 +        a07*w10                
Y1    <- (a00-0.4) + (a01 + 0)*w1 + (a02 + 0.25*a02)*w2 + (a03 + 0)*w3 + (a04 + 0)*w4 +      (a05+0)*w8 +  (a06+0)*w9 +     (a07+0)*w10
} else
if (scenarioY == "f") {# 2 inter
Y0    <-       a00 +       a01*w1 +       a02*w2 +       a03*w3 +       a04*w4       +
    a05*w8 +       a06*w9 +        a07*w10                
Y1    <- (a00-0.4) + (a01 + 0)*w1 + (a02 + 0.25*a02)*w2 + (a03 + 0)*w3 + (a04 + 0.25*a04)*w4 +      (a05+0)*w8 +  (a06+0)*w9 +     (a07+0)*w10
} else 
if (scenarioY == "g") {  # (2 quad, 2 inter)
Y0    <-       a00 +       a01*w1 +       a02*w2 +       a03*w3 +       a04*w4       +
    a05*w8 +       a06*w9 +        a07*w10 
    + 0.1*a02*w2^2 + 0.1*a04*w4^2                
Y1    <- (a00-0.4) + (a01 + 0)*w1  + (a03 + 0)*w3   +   (a05+0)*w8 +  (a06+0)*w9 + (a07+0)*w10 
+ 0.25*0.1*a02*w2^2 + 0.25*0.1*a04*w4^2 
+ (a02 + 0.25*a02)*w2 + (a04 + 0.25*a04)*w4 
} else
if (scenarioY == "b") { #(1 quad terms)
Y0    <-       a00 +       a01*w1 +       a02*w2 +       a03*w3 +       a04*w4       +
    a05*w8 +       a06*w9 +        a07*w10  
    + 0.5*a02*w2^2              
Y1    <- (a00-0.2) + (a01 + 0)*w1 + (a02 + 0)*w2 + (a03 + 0)*w3 + (a04 + 0)*w4 +      (a05+0)*w8 +  (a06+0)*w9 +     (a07+0)*w10 
    + a02*w2^2
} else
if (scenarioY == "c") { #(2 quad terms)
Y0    <-       a00 +       a01*w1 +       a02*w2 +       a03*w3 +       a04*w4       +
    a05*w8 +       a06*w9 +        a07*w10  
    + 0.5*a02*w2^2 + 0.5*a04*w4^2                
Y1    <- (a00-0.1) + (a01 + 0)*w1 + (a02 + 0)*w2 + (a03 + 0)*w3 + (a04 + 0)*w4 +      (a05+0)*w8 +  (a06+0)*w9 +     (a07+0)*w10 
     + a02*w2^2 + a04*w4^2 
} else
if (scenarioY == "e") 
{
# scenario e (1 quad term + 1 interaction)
Y0    <-       a00 +       a01*w1 +       a02*w2 +       a03*w3 +       a04*w4       +
    a05*w8 +       a06*w9 +        a07*w10  
    + 0.5*a02*w2^2              
Y1    <- (a00-0.2) + (a01 + a01)*w1 + (a02 + a02)*w2 + (a03 + a03)*w3 + (a04 + a04)*w4 +      (a05+0)*w8 +  (a06+0)*w9 +     (a07+0)*w10 
+ a02*w2^2 + (a04 + 0.25*a04)*w4
}

# continuous outcome Y

   Y    <- T*Y1 + (1-T)*Y0
   
# true individual effect of T on Y

   indeff <- Y1-Y0

  
# create simulation dataset

 sim <- as.data.frame(cbind(w1, w2, w3 ,w4, w5, w6, w7, w8, w9, w10, T, Y, indeff ))
 return(sim)
} 


 ############## funzione: funsim (x, psmethod)
 
 ############## funsim uses psmethod to estimate ps; then use 1) weighting and 2) matching to estimate treatment effect. Finally, it calculates performance metrics.

#inputs: x = a simulated dataset created by funcov, a psmethod
#output: a list containing performance metrics (bias and var of estimate, covariate balance diagnostics)
 
 
 
 funsim <- function(x,psmethod,par="ATT"){

#~~ estimate ps           
       
    
if (psmethod =="truelogit")
{   	
	  mod = glm(T~ w1 + w2 + w3 + w4, data=x, family=binomial)
	  ps = mod$fitted
	  
} else if (psmethod == "logit"){
     
    mod = glm(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, data=x ,family=binomial)
	 ps = mod$fitted
        
} else if (psmethod == "tree"){

    mod = rpart(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10,method="class",data=x)
	 ps = predict(mod)[,2]     
        
}  else if (psmethod == "randomforest"){
	
      	mod = randomForest(factor(T)~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, data=x, ntree= 500)
 	   ps<-predict(mod , type="prob")[,2]
    
}  else if (psmethod == "gbm"){
              
     mod = gbm(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, data=x, distribution = "bernoulli", interaction.depth = 1, n.trees = 100)  
	  ps = predict.gbm(mod,data=x,n.trees=100, type="response") 
	 
}  else if (psmethod == "gbmtwang"){
	
   mod = ps(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, data=x, n.trees = 10000,interaction.depth = 3,verbose=FALSE,shrinkage = 0.0005)
   ps = as.vector(mod$ps[,1]) 
     
} else if (psmethod == "bag"){
	
   mod = bagging(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, data=x)
   ps =  predict(mod,newdata=x,type="prob")    
}

    else if (psmethod == "nnet")    {	
    	                         
   mod= nnet(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10,
   data=x, entropy=T,size=10,decay=0,maxit=2000,trace=F)
   ps = as.numeric(predict(mod, type='raw')) #nb: anche predizioni fuori da [0,1]
#   ps=exp(ps)/(1+exp(ps))          
                               }
    else if (psmethod == "nb")     {
    
	 mod = naiveBayes(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, data=x)
      ps = predict(mod,newdata=x,type="raw")[,2] 
}

#### fit measure for ps

auc<-auc(sensitivity(ps, factor(x$T), perc.rank = TRUE))

#### true (sample) att:

g <- mean(x$indeff[x$T==1]) 

### true (sample) ate

if(par=="ATE") {g <- mean(x$indeff) }

#### estimating ATT via propensity score weighting

weights     <- ifelse(x$T==1,1,ps/(1-ps))

#### estimating ATE via propensity score weighting

if(par=="ATE") { weights     <- ifelse(x$T==1,1/ps,1/(1-ps)) }

hatg          <- wtd.mean(x$Y[x$T==1],weights=weights[x$T==1])-wtd.mean(x$Y[x$T==0],weights=weights[x$T==0])

absrbias    <- abs((hatg -g)/g)*100 
varhatg      <- (hatg-g)^2 

modw     <- lm( Y ~  T, data=x, weights=weights)
hatgsew <- summary(modw)$coefficients[c("T"),c("Std. Error")]
covw <- ifelse(g > hatg-2*hatgsew  & g < hatg + 2*hatgsew , 1, 0)



# estimating ATT via matching with caliper

rr = Match(Y=x$Y, Tr=x$T, X=ps, caliper=0.25, M=1, replace=TRUE,ties=FALSE)

if(par=="ATE"){
rr = Match(Y=x$Y, Tr=x$T, X=ps, caliper=0.25,M=1,estimand=par,replace=TRUE,ties=FALSE)
}

# estimating ATE via matching with caliper

allmdata<-rbind(x[rr$index.treated,],x[rr$index.control,])

hatgm          <-rr$est
absrbiasm    <- abs((hatgm -g)/g)*100 
varhatgm      <- (hatgm-g)^2 

 hatgsem          <-rr$se.standard
 covm <- ifelse(g > hatgm-2*hatgsem  & g < hatgm + 2*hatgsem , 1, 0)
# = mean(mdata$Y[mdata$Tr==1])-mean(mdata$Y[mdata$Tr==0])

# ~  size of matched dataset

orig.nobs   <-rr$orig.nobs
orig.tnobs  <-rr$orig.treated.nobs
#match.nobs   <-rr$nobs # this is wrong! it is the nobs inthe unmatched dataset
match.nobs   <-length(rr$mdata$Tr)
match.tnobs  <-length(rr$mdata$Tr[rr$mdata$Tr==1])	      




# Balance BEFORE for w1 ... w10

bb = MatchBalance(T ~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10
+w1*w2+w1*w3+w1*w4+w1*w5+w1*w6+w1*w7+w1*w8+w1*w9+w1*w10
+w2*w3+w2*w4+w2*w5+w2*w6+w2*w7+w2*w8+w2*w9+w2*w10
+w3*w4+w3*w5+w3*w6+w3*w7+w3*w8+w3*w9+w3*w10
+w4*w5+w4*w6+w4*w7+w4*w8+w4*w9+w4*w10
+w5*w6+w5*w7+w5*w8+w5*w9+w5*w10
+w6*w7+w6*w8+w6*w9+w6*w10
+w7*w8+w7*w9+w7*w10
+w8*w9+w8*w10
+w9*w10, data=x, ks = FALSE, nboots=0, print.level=0) 
# e.g.  


# ASAM
asb_b    <-vector();for(i in 1:10){asb_b[[i]] <- bb$BeforeMatching[[i]]$sdiff}

# ASAM with interactions
interasbb<-vector();for(i in 1:(length(bb$BeforeMatching))){interasbb[[i]] <- bb$BeforeMatching[[i]]$sdiff}

# qq mean difference raw
vecqqmeanrawb    <-vector();for(i in 1:10){vecqqmeanrawb[[i]] <- bb$BeforeMatching[[i]]$qqsummary.raw$meandiff}

# qq mean diff
vecqqmeanb<-vector();for(i in 1:10){vecqqmeanb[[i]] <- bb$BeforeMatching[[i]]$qqsummary$meandiff}
interqqb<-vector();for(i in 1:length(bb$BeforeMatching)){interqqb[[i]] <- bb$BeforeMatching[[i]]$qqsummary$meandiff}

# qq max difference
vecqqmaxb<-vector();for(i in 1:10){vecqqmaxb[[i]] <- bb$BeforeMatching[[i]]$qqsummary$maxdiff}

# qq max with difference
interqqmaxb<-vector();for(i in 1:(length(bb$BeforeMatching))){interqqmaxb[[i]] <- bb$BeforeMatching[[i]]$qqsummary$maxdiff} 

# variance ratio 
# i.e. before log (var ratio) for ps; reference: Imbens-Rubin Ch 10, eq. 14.5 and tobacco litigation       
varratio_b =abs(log(var(ps[x$T==1])/var(ps[x$T==0])))  
 
#~~ performance metrics: balance AFTER weighting for w1 ...w10

ba = MatchBalance(T ~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10
+w1*w2+w1*w3+w1*w4+w1*w5+w1*w6+w1*w7+w1*w8+w1*w9+w1*w10
+w2*w3+w2*w4+w2*w5+w2*w6+w2*w7+w2*w8+w2*w9+w2*w10
+w3*w4+w3*w5+w3*w6+w3*w7+w3*w8+w3*w9+w3*w10
+w4*w5+w4*w6+w4*w7+w4*w8+w4*w9+w4*w10
+w5*w6+w5*w7+w5*w8+w5*w9+w5*w10
+w6*w7+w6*w8+w6*w9+w6*w10
+w7*w8+w7*w9+w7*w10
+w8*w9+w8*w10
+w9*w10, data=x,weights=weights, ks = FALSE, nboots=0, print.level=0)

# ASAM
asb_a    <-vector()
for(i in 1:10){asb_a[[i]] <- ba$BeforeMatching[[i]]$sdiff}

# ASAM with interactions 
interasba<-vector()
for(i in 1:(length(ba$BeforeMatching))){interasba[[i]] <- ba$BeforeMatching[[i]]$sdiff}


# QQ mean raw
# # mean difference of empirical quantiles for w1 with weights applied
# # (cannot use weighths with matchBalance forr qq measure)
covnames          <- colnames(x[,-which(names(x)%in%c("Y","T"))])
vecqqmeanrawaw <- vector() 
for(i in 1:length(covnames) ) {
	vecqqmeanrawaw[i] <- abs(
 mean(wtd.quantile(x[x$T==1,] [, covnames[i] ],probs=seq(0.01,1-0.01,0.01),weights=weights[x$T==1])-mean(wtd.quantile(x[x$T==0,] [, covnames[i] ],probs=seq(0.01,1-0.01,0.01),weights=weights[x$T==0]))))
	}

# # equivalent to:
# W1qqrawa=abs(mean(wtd.quantile(x[x$T==1,]$w1,probs=seq(0.01,1-0.01,0.01),weights=weights[x$T==1])-mean(wtd.quantile(x[x$T==0,]$w1,probs=seq(0.01,1-0.01,0.01),weights=weights[x$T==0]))))

#mean and max differences of ecdf after weighting (nb: non si possono usare pesi con matchbalance per qq)
covnames          <- colnames(x[,-which(names(x)%in%c("Y","T"))])
vecqqmeanaw <- vecqqmaxaw<-vector() 
for(i in 1:length(covnames) ) {
	vals <- sort(unique(x[, covnames[i] ]))
	wt <- stepfun(wtd.Ecdf(x[x$T==1,][, covnames[i] ],weights=weights[x$T==1])$x[-1],wtd.Ecdf(x[x$T==1,][, covnames[i] ],weights=weights[x$T==1])$ecdf)
swt <- wt(vals)
wc <- stepfun(wtd.Ecdf(x[x$T==0,][, covnames[i] ],weights=weights[x$T==0])$x[-1],wtd.Ecdf(x[x$T==0,][, covnames[i] ],weights=weights[x$T==0])$ecdf)
swc <- wc(vals)
	vecqqmeanaw[i]        <- mean(abs(swt - swc))
	vecqqmaxaw[i] <- max(abs(swt - swc))
	}


# variance ratio
varratio_a =abs(log(wtd.var(ps[x$T==1],weights=weights[x$T==1])/wtd.var(ps[x$T==0],weights=weights[x$T==0])))  # log (var ratio) after weighting; reference in Inbens-Rubin Ch 10, eq. 14.5 and tobacco litigation




#~~ performance metrics: balance AFTER matching for w1 ...w10

bam = MatchBalance(T ~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10
+w1*w2+w1*w3+w1*w4+w1*w5+w1*w6+w1*w7+w1*w8+w1*w9+w1*w10
+w2*w3+w2*w4+w2*w5+w2*w6+w2*w7+w2*w8+w2*w9+w2*w10
+w3*w4+w3*w5+w3*w6+w3*w7+w3*w8+w3*w9+w3*w10
+w4*w5+w4*w6+w4*w7+w4*w8+w4*w9+w4*w10
+w5*w6+w5*w7+w5*w8+w5*w9+w5*w10
+w6*w7+w6*w8+w6*w9+w6*w10
+w7*w8+w7*w9+w7*w10
+w8*w9+w8*w10
+w9*w10, data=allmdata, ks = FALSE, nboots=0, print.level=0) 

# ASAM after matching
asb_am    <-vector()
for(i in 1:10){asb_am[[i]] <- bam$BeforeMatching[[i]]$sdiff}

# ASAM with interactions after matching
interasbam<-vector()
for(i in 1:(length(bam$BeforeMatching))){interasbam[[i]] <- bam$BeforeMatching[[i]]$sdiff}

vecqqmeanrawam    <-vector()
for(i in 1:10){vecqqmeanrawam[[i]] <- bam$BeforeMatching[[i]]$qqsummary.raw$meandiff}

# mean diff in ecdf after matching
vecqqmeanam    <-vector()
for(i in 1:10){vecqqmeanam[[i]] <- bam$BeforeMatching[[i]]$qqsummary$meandiff}

interqqam<-vector()
for(i in 1:(length(bam$BeforeMatching))){interqqam[[i]] <- bam$BeforeMatching[[i]]$qqsummary$meandiff}

# standardized max difference of ecdfs after matching
vecqqmaxam    <-vector()
for(i in 1:10){vecqqmaxam[[i]] <- bam$BeforeMatching[[i]]$qqsummary$maxdiff}

interqqmaxam<-vector()
for(i in 1:(length(bam$BeforeMatching))){interqqmaxam[[i]] <- bam$BeforeMatching[[i]]$qqsummary$maxdiff} 

 # variance ratio
 varratio_am =abs(log(var(ps[rr$index.treated])/var(ps[rr$index.control])))  # abs log (var ratio) after matching

output<-list( 
"auc"=auc,
"hatgw"=hatg,"absrbiasw"=absrbias,"varhatgw"=varhatg,
"hatgm"=hatgm,"absrbiasm"=absrbiasm,"varhatgm"=varhatgm,
"hatgsem"=hatgsem,"hatgsew"=hatgsew,
"covw"=covw,"covm"=covm,
#  mean (over covariates) of balance summaries : 
# asam
"masb_b"=mean(abs(asb_b)),
"masb_aw"=mean(abs(asb_a)),
"masb_am"=mean(abs(asb_am)),
# mean asam>20
"over20masb_b"=sum(abs(asb_b[which(abs(asb_b)>20)])-20)/length(asb_b),
"over20masb_aw"=sum(abs(asb_a[which(abs(asb_a)>20)])-20)/length(asb_a),
"over20masb_am"=sum(abs(asb_am[which(abs(asb_am)>20)])-20)/length(asb_am),
# mean asam>20
"over10masb_b"=sum(abs(asb_b[which(abs(asb_b)>10)])-10)/length(asb_b),
"over10masb_aw"=sum(abs(asb_a[which(abs(asb_a)>10)])-10)/length(asb_a),
"over10masb_am"=sum(abs(asb_am[which(abs(asb_am)>10)])-10)/length(asb_am),
# asam with interactions
"masbinter_b"=mean(abs(interasbb)),
"masbinter_aw"=mean(abs(interasba)),
"masbinter_am"=mean(abs(interasbam)),
# 
"qqmeanraw_b"=mean(vecqqmeanrawb),
"qqmeanraw_aw"=mean(vecqqmeanrawaw),
"qqmeanraw_am"=mean(vecqqmeanam),
#
"qqmean_b"=mean(vecqqmeanb),
"qqmean_aw"=mean(vecqqmeanaw),
"qqmean_am"=mean(vecqqmeanam),
#"qqmeaninter_b"=mean(abs(interqqb)),"qqmeaninter_am"=mean(abs(interqqam)),
"qqmax_b"=mean(vecqqmaxb),
"qqmax_aw"=mean(vecqqmaxaw),
"qqmax_am"=mean(vecqqmaxam),
# variance ratio
"varratio_b"=varratio_b,
"varratio_aw"=varratio_a,
"varratio_am"=varratio_am,
#
"match.tnobs"=match.tnobs,"match.nobs"=match.nobs)
 return(output)
 }
 
 ##################
 # function funsim_opt: funsim with optimally tuned parameters
 ##################
 
 # scenarios<-c("Aa" ,"Ad", "Ae", "Ag", "Ca", "Cd", "Ce", "Cg","Fa", "Fd", "Fe", "Fg","Ga", "Gd", "Ge", "Gg")  
  # tunopt<-c(cp,  mtry,  sizem, decaym, sizew, decayw)

 
 funsim_opt <- function(x,psmethod,par="ATT",tunopt=NULL){

#~~ estimate ps           
       
    
if (psmethod =="truelogit")
{   	
	  mod = glm(T~ w1 + w2 + w3 + w4, data=x, family=binomial)
	  ps = mod$fitted
	  
} else if (psmethod == "logit"){
     
    mod = glm(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, data=x ,family=binomial)
	 ps = mod$fitted
        
} else if (psmethod == "tree"){

    mod = rpart(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10,
    cp=tunopt[2],
    method="class",data=x)
	 ps = predict(mod)[,2]     
        
}  else if (psmethod == "randomforest"){
	
      	mod = randomForest(factor(T)~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, 
      	mtry=tunopt[1], 
      	data=x, ntree= 500)
 	   ps<-predict(mod , type="prob")[,2]
    
}  else if (psmethod == "gbm"){
              
     mod = gbm(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, data=x, distribution = "bernoulli", interaction.depth = 1, n.trees = 100)  
	  ps = predict.gbm(mod,data=x,n.trees=100, type="response") 
	 
}  else if (psmethod == "gbmtwang"){
	
   mod = ps(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, data=x, n.trees = 10000,interaction.depth = 3,verbose=FALSE,shrinkage = 0.0005)
   ps = as.vector(mod$ps[,1]) 
     
} else if (psmethod == "bag"){
	
   mod = bagging(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, data=x)
   ps =  predict(mod,newdata=x,type="prob")    
}
    else if (psmethod == "nn")    {	
    # mod for psm	                         
   mod= nnet(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10,
   data=x, entropy=T,
   size=tunopt[3],
   decay=tunopt[4],maxit=2000,trace=F)
   ps = as.numeric(predict(mod, type='raw')) 
   #ps=exp(ps)/(1+exp(ps))
   # mod for psw
      mod= nnet(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10,
   data=x, entropy=T,
   size=tunopt[5],
   decay=tunopt[6],maxit=2000,trace=F)
   psw = as.numeric(predict(mod, type='raw')) 
   #psw=exp(psw)/(1+exp(psw))          
        
                               }
    else if (psmethod == "nb")     {
    
	 mod = naiveBayes(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, laplace=0,
	 data=x)
      ps = predict(mod,newdata=x,type="raw")[,2] 
}

#### fit measure for ps

auc<-auc(sensitivity(ps, factor(x$T), perc.rank = TRUE))

#### true (sample) att:

g <- mean(x$indeff[x$T==1]) 

### true (sample) ate

if(par=="ATE") {g <- mean(x$indeff) }

#### estimating ATT via propensity score weighting

weights     <- ifelse(x$T==1,1,ps/(1-ps))
if (psmethod == "nn")    {weights     <- ifelse(x$T==1,1,psw/(1-psw))}

#### estimating ATE via propensity score weighting

if(par=="ATE") { weights     <- ifelse(x$T==1,1/ps,1/(1-ps)) 
	if (psmethod == "nn")    {weights     <- ifelse(x$T==1,1/psw,1/(1-psw))}
	 }

hatg          <- wtd.mean(x$Y[x$T==1],weights=weights[x$T==1])-wtd.mean(x$Y[x$T==0],weights=weights[x$T==0])

absrbias    <- abs((hatg -g)/g)*100 
varhatg      <- (hatg-g)^2 

modw     <- lm( Y ~ T , data=x, weight=weights)
hatgsew <- summary(modw)$coefficients[c("T"),c("Std. Error")]
covw <- ifelse(g > hatg-2*hatgsew  & g < hatg + 2*hatgsew , 1, 0)



# estimating ATT via matching with caliper

rr = Match(Y=x$Y, Tr=x$T, X=ps, caliper=0.25, M=1, replace=TRUE,ties=FALSE)

if(par=="ATE"){
rr = Match(Y=x$Y, Tr=x$T, X=ps, caliper=0.25,M=1,estimand=par,replace=TRUE,ties=FALSE)
}

allmdata<-rbind(x[rr$index.treated,],x[rr$index.control,])

hatgm          <-rr$est
absrbiasm    <- abs((hatgm -g)/g)*100 
varhatgm      <- (hatgm-g)^2 

 hatgsem          <-rr$se.standard
 covm <- ifelse(g > hatgm-2*hatgsem  & g < hatgm + 2*hatgsem , 1, 0)



# ~  size of matched dataset

orig.nobs   <-rr$orig.nobs
orig.tnobs  <-rr$orig.treated.nobs
match.nobs   <-length(rr$mdata$Tr)
match.tnobs  <-length(rr$mdata$Tr[rr$mdata$Tr==1])	      




# Balance BEFORE for w1 ... w10

bb = MatchBalance(T ~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10
+w1*w2+w1*w3+w1*w4+w1*w5+w1*w6+w1*w7+w1*w8+w1*w9+w1*w10
+w2*w3+w2*w4+w2*w5+w2*w6+w2*w7+w2*w8+w2*w9+w2*w10
+w3*w4+w3*w5+w3*w6+w3*w7+w3*w8+w3*w9+w3*w10
+w4*w5+w4*w6+w4*w7+w4*w8+w4*w9+w4*w10
+w5*w6+w5*w7+w5*w8+w5*w9+w5*w10
+w6*w7+w6*w8+w6*w9+w6*w10
+w7*w8+w7*w9+w7*w10
+w8*w9+w8*w10
+w9*w10, data=x, ks = FALSE, nboots=0, print.level=0) 
# e.g.  
#ASAM Before for W1: abs(bb$BeforeMatching[[1]]$sdiff) 
#ASAM Before for W1: abs(mean(data$w1[data$T==1])-mean(data$w1[data$T==0]))/sd(data$w1[data$T==1])*100

# ASAM
asb_b    <-vector();for(i in 1:10){asb_b[[i]] <- bb$BeforeMatching[[i]]$sdiff}

# ASAM with interactions
interasbb<-vector();for(i in 1:(length(bb$BeforeMatching))){interasbb[[i]] <- bb$BeforeMatching[[i]]$sdiff}

# qq mean difference raw
vecqqmeanrawb    <-vector();for(i in 1:10){vecqqmeanrawb[[i]] <- bb$BeforeMatching[[i]]$qqsummary.raw$meandiff}

# qq mean diff
vecqqmeanb<-vector();for(i in 1:10){vecqqmeanb[[i]] <- bb$BeforeMatching[[i]]$qqsummary$meandiff}
interqqb<-vector();for(i in 1:length(bb$BeforeMatching)){interqqb[[i]] <- bb$BeforeMatching[[i]]$qqsummary$meandiff}

# qq max difference
vecqqmaxb<-vector();for(i in 1:10){vecqqmaxb[[i]] <- bb$BeforeMatching[[i]]$qqsummary$maxdiff}

# qq max with difference
interqqmaxb<-vector();for(i in 1:(length(bb$BeforeMatching))){interqqmaxb[[i]] <- bb$BeforeMatching[[i]]$qqsummary$maxdiff} 

# variance ratio 
# i.e. before log (var ratio) for ps; reference: Imbens-Rubin Ch 10, eq. 14.5 and tobacco litigation       
varratio_b =abs(log(var(ps[x$T==1])/var(ps[x$T==0])))  
  
    
#~~ performance metrics: balance AFTER weighting for w1 ...w10

ba = MatchBalance(T ~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10
+w1*w2+w1*w3+w1*w4+w1*w5+w1*w6+w1*w7+w1*w8+w1*w9+w1*w10
+w2*w3+w2*w4+w2*w5+w2*w6+w2*w7+w2*w8+w2*w9+w2*w10
+w3*w4+w3*w5+w3*w6+w3*w7+w3*w8+w3*w9+w3*w10
+w4*w5+w4*w6+w4*w7+w4*w8+w4*w9+w4*w10
+w5*w6+w5*w7+w5*w8+w5*w9+w5*w10
+w6*w7+w6*w8+w6*w9+w6*w10
+w7*w8+w7*w9+w7*w10
+w8*w9+w8*w10
+w9*w10, data=x,weights=weights, ks = FALSE, nboots=0, print.level=0)

# ASAM
asb_a    <-vector()
for(i in 1:10){asb_a[[i]] <- ba$BeforeMatching[[i]]$sdiff}

# ASAM with interactions 
interasba<-vector()
for(i in 1:(length(ba$BeforeMatching))){interasba[[i]] <- ba$BeforeMatching[[i]]$sdiff}


# QQ mean raw
# # mean difference of empirical quantiles for w1 with weights applied
# # (cannot use weighths with matchBalance forr qq measure)
covnames          <- colnames(x[,-which(names(x)%in%c("Y","T"))])
vecqqmeanrawaw <- vector() 
for(i in 1:length(covnames) ) {
	vecqqmeanrawaw[i] <- abs(
 mean(wtd.quantile(x[x$T==1,] [, covnames[i] ],probs=seq(0.01,1-0.01,0.01),weights=weights[x$T==1])-mean(wtd.quantile(x[x$T==0,] [, covnames[i] ],probs=seq(0.01,1-0.01,0.01),weights=weights[x$T==0]))))
	}

#mean and max differences of ecdf after weighting (nb: non si possono usare pesi con matchbalance per qq)
covnames          <- colnames(x[,-which(names(x)%in%c("Y","T"))])
vecqqmeanaw <- vecqqmaxaw<-vector() 
for(i in 1:length(covnames) ) {
	vals <- sort(unique(x[, covnames[i] ]))
	wt <- stepfun(wtd.Ecdf(x[x$T==1,][, covnames[i] ],weights=weights[x$T==1])$x[-1],wtd.Ecdf(x[x$T==1,][, covnames[i] ],weights=weights[x$T==1])$ecdf)
swt <- wt(vals)
wc <- stepfun(wtd.Ecdf(x[x$T==0,][, covnames[i] ],weights=weights[x$T==0])$x[-1],wtd.Ecdf(x[x$T==0,][, covnames[i] ],weights=weights[x$T==0])$ecdf)
swc <- wc(vals)
	vecqqmeanaw[i]        <- mean(abs(swt - swc))
	vecqqmaxaw[i] <- max(abs(swt - swc))
	}


# variance ratio
varratio_a =abs(log(wtd.var(ps[x$T==1],weights=weights[x$T==1])/wtd.var(ps[x$T==0],weights=weights[x$T==0])))  # log (var ratio) after weighting; reference in Inbens-Rubin Ch 10, eq. 14.5 and tobacco litigation




#~~ performance metrics: balance AFTER matching for w1 ...w10

bam = MatchBalance(T ~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10
+w1*w2+w1*w3+w1*w4+w1*w5+w1*w6+w1*w7+w1*w8+w1*w9+w1*w10
+w2*w3+w2*w4+w2*w5+w2*w6+w2*w7+w2*w8+w2*w9+w2*w10
+w3*w4+w3*w5+w3*w6+w3*w7+w3*w8+w3*w9+w3*w10
+w4*w5+w4*w6+w4*w7+w4*w8+w4*w9+w4*w10
+w5*w6+w5*w7+w5*w8+w5*w9+w5*w10
+w6*w7+w6*w8+w6*w9+w6*w10
+w7*w8+w7*w9+w7*w10
+w8*w9+w8*w10
+w9*w10, data=allmdata, ks = FALSE, nboots=0, print.level=0) 

# ASAM after matching
asb_am    <-vector()
for(i in 1:10){asb_am[[i]] <- bam$BeforeMatching[[i]]$sdiff}

# ASAM with interactions after matching
interasbam<-vector()
for(i in 1:(length(bam$BeforeMatching))){interasbam[[i]] <- bam$BeforeMatching[[i]]$sdiff}


# mean difference of empirical quantiles after matching. 

vecqqmeanrawam    <-vector()
for(i in 1:10){vecqqmeanrawam[[i]] <- bam$BeforeMatching[[i]]$qqsummary.raw$meandiff}

# mean diff in ecdf after matching

vecqqmeanam    <-vector()
for(i in 1:10){vecqqmeanam[[i]] <- bam$BeforeMatching[[i]]$qqsummary$meandiff}

interqqam<-vector()
for(i in 1:(length(bam$BeforeMatching))){interqqam[[i]] <- bam$BeforeMatching[[i]]$qqsummary$meandiff}


vecqqmaxam    <-vector()
for(i in 1:10){vecqqmaxam[[i]] <- bam$BeforeMatching[[i]]$qqsummary$maxdiff}

interqqmaxam<-vector()
for(i in 1:(length(bam$BeforeMatching))){interqqmaxam[[i]] <- bam$BeforeMatching[[i]]$qqsummary$maxdiff} 


 
 # variance ratio
 varratio_am =abs(log(var(ps[rr$index.treated])/var(ps[rr$index.control])))  # abs log (var ratio) after matching

 
#~~collect output


output<-list( 
"auc"=auc,
"hatgw"=hatg,"absrbiasw"=absrbias,"varhatgw"=varhatg,
"hatgm"=hatgm,"absrbiasm"=absrbiasm,"varhatgm"=varhatgm,
"hatgsem"=hatgsem,"hatgsew"=hatgsew,
"covw"=covw,"covm"=covm,
#  mean (over covariates) of balance summaries : 
# asam
"masb_b"=mean(abs(asb_b)),
"masb_aw"=mean(abs(asb_a)),
"masb_am"=mean(abs(asb_am)),
# mean asam>20
"over20masb_b"=sum(abs(asb_b[which(abs(asb_b)>20)])-20)/length(asb_b),
"over20masb_aw"=sum(abs(asb_a[which(abs(asb_a)>20)])-20)/length(asb_a),
"over20masb_am"=sum(abs(asb_am[which(abs(asb_am)>20)])-20)/length(asb_am),
# mean asam>10
"over10masb_b"=sum(abs(asb_b[which(abs(asb_b)>10)])-10)/length(asb_b),
"over10masb_aw"=sum(abs(asb_a[which(abs(asb_a)>10)])-10)/length(asb_a),
"over10masb_am"=sum(abs(asb_am[which(abs(asb_am)>10)])-10)/length(asb_am),
# asam with interactions
"masbinter_b"=mean(abs(interasbb)),
"masbinter_aw"=mean(abs(interasba)),
"masbinter_am"=mean(abs(interasbam)),
#"KLdivb"=KLdivb,"KLdivaw"=KLdivaw,"KLdivam"=KLdivam,
# 
"qqmeanraw_b"=mean(vecqqmeanrawb),
"qqmeanraw_aw"=mean(vecqqmeanrawaw),
"qqmeanraw_am"=mean(vecqqmeanam),
#
"qqmean_b"=mean(vecqqmeanb),
"qqmean_aw"=mean(vecqqmeanaw),
"qqmean_am"=mean(vecqqmeanam),
#"qqmeaninter_b"=mean(abs(interqqb)),"qqmeaninter_am"=mean(abs(interqqam)),
"qqmax_b"=mean(vecqqmaxb),
"qqmax_aw"=mean(vecqqmaxaw),
"qqmax_am"=mean(vecqqmaxam),
# variance ratio
"varratio_b"=varratio_b,
"varratio_aw"=varratio_a,
"varratio_am"=varratio_am,
#
"match.tnobs"=match.tnobs,"match.nobs"=match.nobs)
 return(output)
 }

 ####################################
# code funsim with additional argument:
# tp := tuning parameter
#####################################

 funsim_tp <- function(x,psmethod,par="ATT",tp=NULL){

#~~ estimate ps           
       
    
if (psmethod =="truelogit")
{   	
	  mod = glm(T~ w1 + w2 + w3 + w4, data=x, family=binomial)
	  ps = mod$fitted
	  
} else if (psmethod == "logit"){
     
    mod = glm(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, data=x ,family=binomial)
	 ps = mod$fitted
        
} else if (psmethod == "tree"){

    mod = rpart(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10,method="class",data=x, 
    cp=tp )
	 ps = predict(mod)[,2]     
        
}  else if (psmethod == "randomforest"){
	
      	mod = randomForest(factor(T)~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, data=x, ntree= 500,
      	mtry=tp)
 	   ps<-predict(mod , type="prob")[,2]
    
}  else if (psmethod == "gbm"){
              
     mod = gbm(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, data=x, distribution = "bernoulli", interaction.depth = 1, 
     n.trees = tp)  
	  ps = predict.gbm(mod,data=x,n.trees=100, type="response") 
	 
}  else if (psmethod == "gbmtwang"){
	
   mod = ps(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, data=x, n.trees = 10000,interaction.depth = 3,verbose=FALSE,shrinkage = 0.0005)
   ps = as.vector(mod$ps[,1]) 
     
} else if (psmethod == "bag"){
	
   mod = bagging(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, 
   cp=tp,data=x)
   ps =  predict(mod,newdata=x,type="prob")    
}
    else if (psmethod == "nnet")    {	
    	                         
   mod= nnet(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10,
   data=x, entropy=tp[1],
   size=tp[2],
   decay=tp[3],maxit=100,trace=F)
   ps<-as.numeric(mod$fitted.values)      
                               }
    else if (psmethod == "nb")     {
    
	 mod = naiveBayes(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, data=x,
	 laplace=tp)
      ps = predict(mod,newdata=x,type="raw")[,2] 
}

#### estimating ATT via propensity score weighting


weights     <- ifelse(x$T==1,1,ps/(1-ps))
if(par=="ATE") { weights     <- ifelse(x$T==1,1/ps,1/(1-ps)) }


modw     <- lm( Y ~ T , data=x, weights=weights)
hatgw <- summary(modw)$coefficients[c("T"),c("Estimate")]
sdhatgw <- summary(modw)$coefficients[c("T"),c("Std. Error")]


# estimating ATT via matching with caliper

rr = Match(Y=x$Y, Tr=x$T, X=ps, caliper=0.25, M=1, estimand="ATT",replace=TRUE,ties=FALSE)

if(par=="ATE"){
rr = Match(Y=x$Y, Tr=x$T, X=ps, caliper=0.25,M=1,estimand="ATE",replace=TRUE,ties=FALSE)


}


allmdata<-rbind(x[rr$index.treated,],x[rr$index.control,])
modm     <- lm( Y ~ T , data=allmdata )
hatgm <- summary(modm)$coefficients[c("T"),c("Estimate")]
sdhatgm <- summary(modm)$coefficients[c("T"),c("Std. Error")]

# ~  size of matched dataset

orig.nobs   <-rr$orig.nobs
orig.tnobs  <-rr$orig.treated.nobs
#match.nobs   <-rr$nobs # this is wrong! it is the nobs inthe unmatched dataset
match.nobs   <-length(rr$mdata$Tr)
match.tnobs  <-length(rr$mdata$Tr[rr$mdata$Tr==1])	      




# Balance BEFORE for w1 ... w10

bb = MatchBalance(T ~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10
+w1*w2+w1*w3+w1*w4+w1*w5+w1*w6+w1*w7+w1*w8+w1*w9+w1*w10
+w2*w3+w2*w4+w2*w5+w2*w6+w2*w7+w2*w8+w2*w9+w2*w10
+w3*w4+w3*w5+w3*w6+w3*w7+w3*w8+w3*w9+w3*w10
+w4*w5+w4*w6+w4*w7+w4*w8+w4*w9+w4*w10
+w5*w6+w5*w7+w5*w8+w5*w9+w5*w10
+w6*w7+w6*w8+w6*w9+w6*w10
+w7*w8+w7*w9+w7*w10
+w8*w9+w8*w10
+w9*w10, data=x, ks = FALSE, nboots=0, print.level=0) 


# ASAM
asb_b    <-vector();for(i in 1:10){asb_b[[i]] <- bb$BeforeMatching[[i]]$sdiff}

# ASAM with interactions
interasbb<-vector();for(i in 1:(length(bb$BeforeMatching))){interasbb[[i]] <- bb$BeforeMatching[[i]]$sdiff}

# # qq mean difference raw
# vecqqmeanrawb    <-vector();for(i in 1:10){vecqqmeanrawb[[i]] <- bb$BeforeMatching[[i]]$qqsummary.raw$meandiff}

# # qq mean diff
# vecqqmeanb<-vector();for(i in 1:10){vecqqmeanb[[i]] <- bb$BeforeMatching[[i]]$qqsummary$meandiff}
# interqqb<-vector();for(i in 1:length(bb$BeforeMatching)){interqqb[[i]] <- bb$BeforeMatching[[i]]$qqsummary$meandiff}

# # qq max difference
# vecqqmaxb<-vector();for(i in 1:10){vecqqmaxb[[i]] <- bb$BeforeMatching[[i]]$qqsummary$maxdiff}

# # qq max with difference
# interqqmaxb<-vector();for(i in 1:(length(bb$BeforeMatching))){interqqmaxb[[i]] <- bb$BeforeMatching[[i]]$qqsummary$maxdiff} 

# # variance ratio 
# # i.e. before log (var ratio) for ps; reference: Imbens-Rubin Ch 10, eq. 14.5 and tobacco litigation       
# varratio_b =abs(log(var(ps[x$T==1])/var(ps[x$T==0])))  

# # L1 distance
# #l1b=(imbalance(group = x$T, data = x[c("w1", "w2", "w3", "w4", "w5","w6", "w7", "w8", "w9", "w10")])$L1)[1]
 
    
#~~ performance metrics: balance AFTER weighting for w1 ...w10

ba = MatchBalance(T ~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10
+w1*w2+w1*w3+w1*w4+w1*w5+w1*w6+w1*w7+w1*w8+w1*w9+w1*w10
+w2*w3+w2*w4+w2*w5+w2*w6+w2*w7+w2*w8+w2*w9+w2*w10
+w3*w4+w3*w5+w3*w6+w3*w7+w3*w8+w3*w9+w3*w10
+w4*w5+w4*w6+w4*w7+w4*w8+w4*w9+w4*w10
+w5*w6+w5*w7+w5*w8+w5*w9+w5*w10
+w6*w7+w6*w8+w6*w9+w6*w10
+w7*w8+w7*w9+w7*w10
+w8*w9+w8*w10
+w9*w10, data=x,weights=weights, ks = FALSE, nboots=0, print.level=0)

# ASAM
asb_a    <-vector()
for(i in 1:10){asb_a[[i]] <- ba$BeforeMatching[[i]]$sdiff}


# ASAM with interactions 
interasba<-vector()
for(i in 1:(length(ba$BeforeMatching))){interasba[[i]] <- ba$BeforeMatching[[i]]$sdiff}



# # QQ mean raw
# # # mean difference of empirical quantiles for w1 with weights applied
# # # (cannot use weighths with matchBalance forr qq measure)
# covnames          <- colnames(x[,-which(names(x)%in%c("Y","T"))])
# vecqqmeanrawaw <- vector() 
# for(i in 1:length(covnames) ) {
	# vecqqmeanrawaw[i] <- abs(
 # mean(wtd.quantile(x[x$T==1,] [, covnames[i] ],probs=seq(0.01,1-0.01,0.01),weights=weights[x$T==1])-mean(wtd.quantile(x[x$T==0,] [, covnames[i] ],probs=seq(0.01,1-0.01,0.01),weights=weights[x$T==0]))))
	# }

# #mean and max differences of ecdf after weighting (nb: non si possono usare pesi con matchbalance per qq)
# covnames          <- colnames(x[,-which(names(x)%in%c("Y","T"))])
# vecqqmeanaw <- vecqqmaxaw<-vector() 
# for(i in 1:length(covnames) ) {
	# vals <- sort(unique(x[, covnames[i] ]))
	# wt <- stepfun(wtd.Ecdf(x[x$T==1,][, covnames[i] ],weights=weights[x$T==1])$x[-1],wtd.Ecdf(x[x$T==1,][, covnames[i] ],weights=weights[x$T==1])$ecdf)
# swt <- wt(vals)
# wc <- stepfun(wtd.Ecdf(x[x$T==0,][, covnames[i] ],weights=weights[x$T==0])$x[-1],wtd.Ecdf(x[x$T==0,][, covnames[i] ],weights=weights[x$T==0])$ecdf)
# swc <- wc(vals)
	# vecqqmeanaw[i]        <- mean(abs(swt - swc))
	# vecqqmaxaw[i] <- max(abs(swt - swc))
	# }


# # variance ratio
# varratio_a =abs(log(wtd.var(ps[x$T==1],weights=weights[x$T==1])/wtd.var(ps[x$T==0],weights=weights[x$T==0])))  # log (var ratio) after weighting; reference in Inbens-Rubin Ch 10, eq. 14.5 and tobacco litigation




#~~ performance metrics: balance AFTER matching for w1 ...w10

bam = MatchBalance(T ~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10
+w1*w2+w1*w3+w1*w4+w1*w5+w1*w6+w1*w7+w1*w8+w1*w9+w1*w10
+w2*w3+w2*w4+w2*w5+w2*w6+w2*w7+w2*w8+w2*w9+w2*w10
+w3*w4+w3*w5+w3*w6+w3*w7+w3*w8+w3*w9+w3*w10
+w4*w5+w4*w6+w4*w7+w4*w8+w4*w9+w4*w10
+w5*w6+w5*w7+w5*w8+w5*w9+w5*w10
+w6*w7+w6*w8+w6*w9+w6*w10
+w7*w8+w7*w9+w7*w10
+w8*w9+w8*w10
+w9*w10, data=allmdata, ks = FALSE, nboots=0, print.level=0) 

# ASAM after matching
asb_am    <-vector()
for(i in 1:10){asb_am[[i]] <- bam$BeforeMatching[[i]]$sdiff}

# ASAM with interactions after matching
interasbam<-vector()
for(i in 1:(length(bam$BeforeMatching))){interasbam[[i]] <- bam$BeforeMatching[[i]]$sdiff}


# # mean difference of empirical quantiles after matching. 

# vecqqmeanrawam    <-vector()
# for(i in 1:10){vecqqmeanrawam[[i]] <- bam$BeforeMatching[[i]]$qqsummary.raw$meandiff}

# # mean diff in ecdf after matching
# # eg for W1: abs(bam$BeforeMatching[[1]]$qqsummary$meandiff)
# vecqqmeanam    <-vector()
# for(i in 1:10){vecqqmeanam[[i]] <- bam$BeforeMatching[[i]]$qqsummary$meandiff}

# interqqam<-vector()
# for(i in 1:(length(bam$BeforeMatching))){interqqam[[i]] <- bam$BeforeMatching[[i]]$qqsummary$meandiff}

# # standardized max difference of ecdfs after matching

# vecqqmaxam    <-vector()
# for(i in 1:10){vecqqmaxam[[i]] <- bam$BeforeMatching[[i]]$qqsummary$maxdiff}

# interqqmaxam<-vector()
# for(i in 1:(length(bam$BeforeMatching))){interqqmaxam[[i]] <- bam$BeforeMatching[[i]]$qqsummary$maxdiff} 


 
 # # variance ratio
 # varratio_am =abs(log(var(ps[rr$index.treated])/var(ps[rr$index.control])))  # abs log (var ratio) after matching

 
#~~collect output


output<-list( 
#"auc"=auc,
 "hatgw"=hatgw, "sdhatgw"=sdhatgw,
 "hatgm"=hatgm, "sdhatgm"=sdhatgm,
#"covw"=covw,"covm"=covm,
#  mean (over covariates) of balance summaries : 
# asam
"masb_b"=mean(abs(asb_b)),
"masb_aw"=mean(abs(asb_a)),
"masb_am"=mean(abs(asb_am)),
# mean asam>20
"over20masb_b"=sum(abs(asb_b[which(abs(asb_b)>20)])-20)/length(asb_b),
"over20masb_aw"=sum(abs(asb_a[which(abs(asb_a)>20)])-20)/length(asb_a),
"over20masb_am"=sum(abs(asb_am[which(abs(asb_am)>20)])-20)/length(asb_am),
# mean asam>10
"over10masb_b"=sum(abs(asb_b[which(abs(asb_b)>10)])-10)/length(asb_b),
"over10masb_aw"=sum(abs(asb_a[which(abs(asb_a)>10)])-10)/length(asb_a),
"over10masb_am"=sum(abs(asb_am[which(abs(asb_am)>10)])-10)/length(asb_am),
# asam with interactions
"masbinter_b"=mean(abs(interasbb)),
"masbinter_aw"=mean(abs(interasba)),
"masbinter_am"=mean(abs(interasbam))
#"KLdivb"=KLdivb,"KLdivaw"=KLdivaw,"KLdivam"=KLdivam,
# # 
# "qqmeanraw_b"=mean(vecqqmeanrawb),
# "qqmeanraw_aw"=mean(vecqqmeanrawaw),
# "qqmeanraw_am"=mean(vecqqmeanam),
# #
# "qqmean_b"=mean(vecqqmeanb),
# "qqmean_aw"=mean(vecqqmeanaw),
# "qqmean_am"=mean(vecqqmeanam),
# #"qqmeaninter_b"=mean(abs(interqqb)),"qqmeaninter_am"=mean(abs(interqqam)),
# "qqmax_b"=mean(vecqqmaxb),
# "qqmax_aw"=mean(vecqqmaxaw),
# "qqmax_am"=mean(vecqqmaxam),
# #"qqmaxinter_b"=mean(abs(interqqmaxb)),
# #"qqmaxinter_am"=mean(abs(interqqmaxam)),
# # variance ratio
# "varratio_b"=varratio_b,
# "varratio_aw"=varratio_a,
# "varratio_am"=varratio_am,
# #
# "match.tnobs"=match.tnobs,"match.nobs"=match.nobs
)
 return(output)
 }
 #######################################
 #######     functions for CASE STUDY       ####
#######################################


###### install and load packages

# list.of.packages = c("Matching","MASS","arm","rpart","randomForest","gbm","twang","ipred","neuralnet",
# "nnet","e1071","klaR","xtable","flexmix","AUC","data.table","Hmisc") # replace xx and yy with package names
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages) > 0) {install.packages(new.packages)}
# lapply(list.of.packages, require, character.only=T)


 ##########################
 # function funcov_cs
 ###########################
 
 # generates dataset with real covariates and simulated T and Y
 # input: x: observed covariate from the case study
 # output T: binary treatment; 
 # output Y: a continuous or binary outcome variable.


 # load("cs_data/clean_data_10cov.Rdata")
 
 funcov_cs<-function(x,scenarioT,scenarioY,Ytype="disc"){
 x<-x[ , !(names(x) %in% c("Y","T"))]
 x<-data.frame(scale(x))
 # take covariates from real dataset
 w1<-x$w1;w2<-x$w2;w3<-x$w3;w4<-x$w4;w5<-x$w5
 w6<-x$w6;w7<-x$w7;w8<-x$w8;w9<-x$w9;w10<-x$w10 
 #T<- x$T #per prendere anche T
 
#~~ scenarios for data generation models 

# A: model with additivity and linearity
# B: mild non-linearity
# C: moderate non-linearity
# D: mild non-additivity
# E: mild non-additivity and non-linearity
# F: moderate non-additivity
# G: moderate non-additivity and non-linearity
# binary exposure modeling


#~~~~~~~~~~~~~~Global Variables~~~~~~~~~~~~~~~~~~~~#

#~~ coefficients for the treatment and potential outcomes equations

#~~ coefficients for treatment equation (roughly from model below)
# 	  modT = glm(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, data=data2, family=binomial)

  b0 <-  -4.5 
  b1 <-  0.01
  b2 <-  0.18
  b3 <-   0.014
  b4 <- 0 
  b5 <- 0.04 
  b6 <- 0.015
  b7 <-  0.1 
  b8 <- -0.04 
  b9 <- 0.13 
  b10 <-  0.05  
 
   #~~ scenarios for treatment assignment
 
 if (scenarioT == "T") {
#trueps <- (1 + exp( -(
linpred<-b0 + b1*w1 + b2*w2 + b3*w3 + b4*w4 + b5*w5 + b6*w6 + b7*w7 +  b8*w8 + b9*w9 + b10*w10 
#) ) ) ^ -1
} 

#~~ binary treatment Tnew 
unif1 <- runif(nrow(x),min(linpred),max(linpred))
 T <- ifelse(linpred>unif1,1,0) # there is a probability of unif1 that T=1
  # table(data2$T , T)

 # ~~ coefficients for outcome equation (roughly from model below)
 # ~~  modY = glm(Y~ Tnew + w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, data=data2, family=binomial)
 
  # ~~ coefficients for potential outcomes equation:
 # ~~ coefficients for Y(0)
 # ~~ coefficients for Y(1): same as Y(0) (zero effect)

 #~~ coefficients for Y0 equation from model (modY) above

 a00 <- -3
 a01 <-  0.05 
 a02 <- 1.89 
 a03 <- 0.025
 a04 <- 0.0004  
 a05 <-  0.07 
 a06 <- 0.033
 a07 <-  -0.03
 a08 <- -0.006
 a09 <- 0.42
 a10 <- -0.15

#~~ scenarios Y
if (scenarioY == "t") {# we assume no effect of T: Y1=Y0

Y0    <- 
#(1 + exp( -(  
a00 +       a01*w1 +       a02*w2 +       a03*w3 +       a04*w4       +
    a05*w5 +       a06*w6 +        a07*w7 + a08*w8 +       a09*w9 +        a10*w10
# ) ) )^-1             
Y1    <- #(1 + exp( -(
Y0-0.4
#) ) )^-1
}
 else
if (scenarioY == "a") {
Y0    <- #(1 + exp( -(  
    a00 +       a01*w1 +       a02*w2 +       a03*w3 +       a04*w4       +
    a05*w8 +       a06*w9 +        a07*w10   # ) ) )^-1             
Y1    <- #(1 + exp( -(
(a00-0.4) + (a01 + 0)*w1 + (a02 + 0)*w2 + (a03 + 0)*w3 + (a04 + 0)*w4 +      (a05+0)*w8 +  (a06+0)*w9 +     (a07+0)*w10 
#) ) )^-1
} 
if (Ytype=="cont"){
	
#~~ continuos outcome Y

   Y    <- T*Y1 + (1-T)*Y0
}
if (Ytype=="disc"){# discretize Y
	
#~~ discrete outcome Y

 pY1_0 <- (1+exp(-(Y0)))^-1  
# pY1_1 <- (1+exp(-(Y1)))^-1 


  unif1 <- runif(length(T),0,1)

 # Y1     <- ifelse( pY1_1> unif1, 1, 0)	# there is a probability of unif1 that Y=1
 
  Y0     <- ifelse( pY1_0 > unif1, 1, 0)
 
  Y1<-Y0
  Y1[Y1==1] <-  Y1[Y1==1]-rbinom(size=1, n=length(Y1[Y1==1]), prob=0.4)
 
  Y    <- T*Y1 + (1-T)*Y0
}
 # true individual effect of T on Y
if (Ytype=="disc"){
    indeff <- Y1-Y0
    }
 if (Ytype=="cont"){
    indeff <- Y1-Y0
    }
 
# create simulation dataset

 sim <- as.data.frame(cbind(w1, w2, w3 ,w4, w5, w6, w7, w8, w9, w10, T, Y,indeff ))
 return(sim)
}
 ###########################################
 ##############    funsim_cs (x, psmethod)
 ###########################################
 
 # funsim_cs uses psmethod to estimate ps; then use 1) weighting and 2) matching to estimate treatment effect. Finally, it calculates performance metrics.

#inputs: x = a simulated dataset created by funcov_cs, a psmethod
#output: a list containing performance metrics (bias and var of estimate, covariate balance diagnostics)

# note:
# the value of tuning parameters for random forest, tree and neural net are fixed using optimal tuning (caret package) and they are:
#mtry=3; asam aw =1.86 asam_am= 7
#cp=0.0015; asam am / aw less than 5
# for nn both size 7 decay 0.005 and size 19 decay 0.55 gave low asam am/ aw (less than 5)
 
 
 funsim_cs <- function(x,psmethod,Ytype="disc"){

#~~ estimate ps           
       
    
if (psmethod =="truelogit")
{   	
	  mod = glm(T~ w1 + w2 + w3 + w4, data=x, family=binomial)
	  ps = mod$fitted
	  
} else if (psmethod == "logit"){
     
    mod = glm(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, data=x ,family=binomial)
	 ps = mod$fitted
        
} else if (psmethod == "tree"){

    mod = rpart(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10,cp=0.0015,method="class",data=x)
	 ps = predict(mod)[,2]     
        
}  else if (psmethod == "randomforest"){
	
      	mod = randomForest(factor(T)~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, data=x, mtry=2,ntree= 500)
 	   ps<-predict(mod , type="prob")[,2]
    
}  else if (psmethod == "gbm"){
              
     mod = gbm(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, data=x, distribution = "bernoulli", n.trees = 100)  
	  ps = predict.gbm(mod,data=x,n.trees=100, type="response") 
	 
}  else if (psmethod == "gbmtwang"){
	
   mod = ps(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, data=x, n.trees = 10000,interaction.depth = 3,verbose=FALSE,shrinkage = 0.0005)
   ps = as.vector(mod$ps[,1]) 
     
} else if (psmethod == "bag"){
	
   mod = bagging(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, cp=0.001,data=x)
   ps =  predict(mod,newdata=x,type="prob")    
}
# else if (psmethod == "nn")    {
	                         
   # mod= neuralnet(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10,
   # data=x, hidden=10,lifesign="none", err.fct="sse",linear.output=FALSE)
   # ps = prediction(mod)$rep1[,"T"] #nb: anche predizioni fuori da [0,1]
   # ps=exp(ps)/(1+exp(ps))
   # #ps = ifelse(ps<0,0.01,ps);ps<-ifelse(ps>1,0.99,ps)               
                                # }
    else if (psmethod == "nnm")    {	
    	                         
   mod= nnet(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10,
   data=x, entropy=T,size=11,decay=1.3,maxit=2000,trace=F)
   ps = as.numeric(predict(mod, type='raw')) #nb: anche predizioni fuori da [0,1]
   #ps=exp(ps)/(1+exp(ps))          
                               }
   else if (psmethod == "nnw")    {	
    	                         
   mod= nnet(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10,
   data=x, entropy=T,size=13,decay=1.3,maxit=2000,trace=F)
   ps = as.numeric(predict(mod, type='raw')) #nb: anche predizioni fuori da [0,1]
   #ps=exp(ps)/(1+exp(ps))          
                               }
    else if (psmethod == "nb")     {
    
	 mod = naiveBayes(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, laplace=0,data=x)
      ps = predict(mod,newdata=x,type="raw")[,2] 
}

#### fit measure for ps and adjust extreme ps

auc<-auc(sensitivity(ps, factor(x$T), perc.rank = TRUE))

#### true (sample) att:

#g <- mean(x$indeff[x$T==1]) 

### true (sample) ate

g <- mean(x$indeff)

 
if(Ytype=="cont"){
#### estimating ATT via propensity score weighting

# weights     <- ifelse(x$T==1,1,ps/(1-ps))

#### estimating ATE via propensity score weighting

weights     <- ifelse(x$T==1,1/ps,1/(1-ps))

hatg          <- wtd.mean(x$Y[x$T==1],weights=weights[x$T==1])-wtd.mean(x$Y[x$T==0],weights=weights[x$T==0])
#equivalent to: sum(x$Y[x$T==1]*weights[x$T==1]/sum(weights[x$T==1]))-sum(x$Y[x$T==0]*weights[x$T==0]/sum(weights[x$T==0]))

modw     <- lm( Y ~   T, data=x, weights=weights)
hatgsew <- summary(modw)$coefficients[c("T"),c("Std. Error")]

absrbias    <- abs((hatg -g)/g)*100 

varhatg      <- (hatg-g)^2 

# estimating ATT via matching with caliper

rr = Match(Y=x$Y, Tr=x$T, X=ps, caliper=0.25, M=1, replace=TRUE,ties=FALSE)

#mY <-rr$mdata$Y ; mT <- rr$mdata$Tr 

allmdata<-rbind(x[rr$index.treated,],x[rr$index.control,])

hatgm          <-rr$est
 hatgsem          <-rr$se
# = mean(mdata$Y[mdata$Tr==1])-mean(mdata$Y[mdata$Tr==0])
# = mean(allmdata$Y[allmdata$T==1])-mean(allmdata$Y[allmdata$T==0])

absrbiasm    <- abs((hatgm -g)/g)*100 
varhatgm      <- (hatgm-g)^2 
}

if(Ytype=="disc"){
#### estimating ATT via propensity score weighting

#weights     <- ifelse(x$T==1,1,ps/(1-ps))

#### estimating ATE via propensity score weighting

weights     <- ifelse(x$T==1,1/ps,1/(1-ps))

#hatg          <- wtd.mean(x$Y[x$T==1],weights=weights[x$T==1])-wtd.mean(x$Y[x$T==0],weights=weights[x$T==0])
#equivalent to: sum(x$Y[x$T==1]*weights[x$T==1]/sum(weights[x$T==1]))-sum(x$Y[x$T==0]*weights[x$T==0]/sum(weights[x$T==0]))

modw     <- lm( Y ~  T , data=x, weight=weights)
hatg        <- summary(modw)$coefficients["T","Estimate"]#mow$coefficients["T"]
hatgsew <- summary(modw)$coefficients["T","Std. Error"]
covw <- ifelse(g > hatg-2*hatgsew  & g < hatg + 2*hatgsew , 1, 0)

absrbias    <- abs((hatg -g)/g)*100 

varhatg      <- (hatg-g)^2 

# estimating ATT via matching with caliper

rr = Match(Y=x$Y, Tr=x$T, X=ps, caliper=0.25, M=1, replace=TRUE,ties=FALSE)

#mY <-rr$mdata$Y ; mT <- rr$mdata$Tr 

allmdata<-rbind(x[rr$index.treated,],x[rr$index.control,])

modm     <- lm( Y ~   T , data=allmdata)
hatgm        <- summary(modm)$coefficients["T","Estimate"]
hatgsem <- summary(modm)$coefficients["T","Std. Error"]
 covm <- ifelse(g > hatgm-2*hatgsem  & g < hatgm + 2*hatgsem , 1, 0)

absrbiasm    <- abs((hatgm - g)/g)*100 
varhatgm      <- (hatgm-g)^2 
}

# ~  size of matched dataset

orig.nobs   <-rr$orig.nobs
orig.tnobs  <-rr$orig.treated.nobs
#match.nobs   <-rr$nobs # this is wrong! it is the nobs inthe unmatched dataset
match.nobs   <-length(rr$mdata$Tr)
match.tnobs  <-length(rr$mdata$Tr[rr$mdata$Tr==1])	      




# Balance BEFORE for w1 ... w10

bb = MatchBalance(T ~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10
+w1*w2+w1*w3+w1*w4+w1*w5+w1*w6+w1*w7+w1*w8+w1*w9+w1*w10
+w2*w3+w2*w4+w2*w5+w2*w6+w2*w7+w2*w8+w2*w9+w2*w10
+w3*w4+w3*w5+w3*w6+w3*w7+w3*w8+w3*w9+w3*w10
+w4*w5+w4*w6+w4*w7+w4*w8+w4*w9+w4*w10
+w5*w6+w5*w7+w5*w8+w5*w9+w5*w10
+w6*w7+w6*w8+w6*w9+w6*w10
+w7*w8+w7*w9+w7*w10
+w8*w9+w8*w10
+w9*w10, data=x, ks = FALSE, nboots=0, print.level=0) 
# e.g.  
#ASAM Before for W1: abs(bb$BeforeMatching[[1]]$sdiff) 
#ASAM Before for W1: abs(mean(data$w1[data$T==1])-mean(data$w1[data$T==0]))/sd(data$w1[data$T==1])*100

# ASAM
asb_b    <-vector();for(i in 1:10){asb_b[[i]] <- bb$BeforeMatching[[i]]$sdiff}

# ASAM with interactions
interasbb<-vector();for(i in 1:(length(bb$BeforeMatching))){interasbb[[i]] <- bb$BeforeMatching[[i]]$sdiff}

# qq mean difference raw
vecqqmeanrawb    <-vector();for(i in 1:10){vecqqmeanrawb[[i]] <- bb$BeforeMatching[[i]]$qqsummary.raw$meandiff}

# qq mean diff
vecqqmeanb<-vector();for(i in 1:10){vecqqmeanb[[i]] <- bb$BeforeMatching[[i]]$qqsummary$meandiff}
interqqb<-vector();for(i in 1:length(bb$BeforeMatching)){interqqb[[i]] <- bb$BeforeMatching[[i]]$qqsummary$meandiff}

# qq max difference
vecqqmaxb<-vector();for(i in 1:10){vecqqmaxb[[i]] <- bb$BeforeMatching[[i]]$qqsummary$maxdiff}

# qq max with difference
interqqmaxb<-vector();for(i in 1:(length(bb$BeforeMatching))){interqqmaxb[[i]] <- bb$BeforeMatching[[i]]$qqsummary$maxdiff} 

# variance ratio 
# i.e. before log (var ratio) for ps; reference: Imbens-Rubin Ch 10, eq. 14.5 and tobacco litigation       
varratio_b =abs(log(var(ps[x$T==1])/var(ps[x$T==0])))  

# L1 distance
#l1b=(imbalance(group = x$T, data = x[c("w1", "w2", "w3", "w4", "w5","w6", "w7", "w8", "w9", "w10")])$L1)[1]

#l1b_main<-mean((data.frame(imbalance(group = x$T, data = x[c("w1", "w2", "w3", "w4", "w5","w6", "w7", "w8", "w9", "w10")])[1]))[,3])      
    
#~~ performance metrics: balance AFTER weighting for w1 ...w10

ba = MatchBalance(T ~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10
+w1*w2+w1*w3+w1*w4+w1*w5+w1*w6+w1*w7+w1*w8+w1*w9+w1*w10
+w2*w3+w2*w4+w2*w5+w2*w6+w2*w7+w2*w8+w2*w9+w2*w10
+w3*w4+w3*w5+w3*w6+w3*w7+w3*w8+w3*w9+w3*w10
+w4*w5+w4*w6+w4*w7+w4*w8+w4*w9+w4*w10
+w5*w6+w5*w7+w5*w8+w5*w9+w5*w10
+w6*w7+w6*w8+w6*w9+w6*w10
+w7*w8+w7*w9+w7*w10
+w8*w9+w8*w10
+w9*w10, data=x,weights=weights, ks = FALSE, nboots=0, print.level=0)

# ASAM
asb_a    <-vector()
for(i in 1:10){asb_a[[i]] <- ba$BeforeMatching[[i]]$sdiff}

#asb_a for W1=abs(sum(x$w1[x$T==1]*weights[x$T==1])/sum(weights[x$T==1])-sum(x$w1[x$T==0]*weights[x$T==0])/sum(weights[x$T==0]))/sd(x$w1[x$T==1])*100
# equivalente a: abs(wtd.mean(x$w1[x$T==1],weights=weights[x$T==1])-wtd.mean(x$w1[x$T==0],weights=weights[x$T==0]))/sd(x$w1[x$T==1])*100
# equivalente a: ba$BeforeMatching[[1]]$sdiff

# ASAM with interactions 
interasba<-vector()
for(i in 1:(length(ba$BeforeMatching))){interasba[[i]] <- ba$BeforeMatching[[i]]$sdiff}

#vecKLdivaw<-vector();for(i in 1:10){vecKLdivaw[[i]] <-KLdiv(cbind(x[,i][x$T==1],rep(x[,i][x$T==0],weights[x$T==0])))[1,2]}

# QQ mean raw
# # mean difference of empirical quantiles for w1 with weights applied
# # (cannot use weighths with matchBalance forr qq measure)
covnames          <- colnames(x[,-which(names(x)%in%c("Y","T"))])
vecqqmeanrawaw <- vector() 
for(i in 1:length(covnames) ) {
	vecqqmeanrawaw[i] <- abs(
 mean(wtd.quantile(x[x$T==1,] [, covnames[i] ],probs=seq(0.01,1-0.01,0.01),weights=weights[x$T==1])-mean(wtd.quantile(x[x$T==0,] [, covnames[i] ],probs=seq(0.01,1-0.01,0.01),weights=weights[x$T==0]))))
	}

# # equivalent to:
# W1qqrawa=abs(mean(wtd.quantile(x[x$T==1,]$w1,probs=seq(0.01,1-0.01,0.01),weights=weights[x$T==1])-mean(wtd.quantile(x[x$T==0,]$w1,probs=seq(0.01,1-0.01,0.01),weights=weights[x$T==0]))))



#mean and max differences of ecdf after weighting (nb: non si possono usare pesi con matchbalance per qq)
covnames          <- colnames(x[,-which(names(x)%in%c("Y","T"))])
vecqqmeanaw <- vecqqmaxaw<-vector() 
for(i in 1:length(covnames) ) {
	vals <- sort(unique(x[, covnames[i] ]))
	wt <- stepfun(wtd.Ecdf(x[x$T==1,][, covnames[i] ],weights=weights[x$T==1])$x[-1],wtd.Ecdf(x[x$T==1,][, covnames[i] ],weights=weights[x$T==1])$ecdf)
swt <- wt(vals)
wc <- stepfun(wtd.Ecdf(x[x$T==0,][, covnames[i] ],weights=weights[x$T==0])$x[-1],wtd.Ecdf(x[x$T==0,][, covnames[i] ],weights=weights[x$T==0])$ecdf)
swc <- wc(vals)
	vecqqmeanaw[i]        <- mean(abs(swt - swc))
	vecqqmaxaw[i] <- max(abs(swt - swc))
	}

# # 
# vals <- sort(unique(x$w1))
# wt <- stepfun(wtd.Ecdf(x[x$T==1,]$w1,weights=weights[x$T==1])$x[-1],wtd.Ecdf(x[x$T==1,]$w1,weights=weights[x$T==1])$ecdf)
# swt <- wt(vals)
# wc <- stepfun(wtd.Ecdf(x[x$T==0,]$w1,weights=weights[x$T==0])$x[-1],wtd.Ecdf(x[x$T==0,]$w1,weights=weights[x$T==0])$ecdf)
# swc <- wc(vals)
# W1qqa<-mean(abs(swt - swc))#Standardized mean difference of ecdfs for w1 with weights applied: max(abs(swt-swc)) where swt = F_{W1 | T=1} (w1) and swc = F_{W1 | T=0} (w1); nb: normwt non efficace perché sono differenze tra frequenze cumulate nb: non si possono usare i pesi per qq con MatchBalance
# W1qqmaxa<-max(abs(swt - swc))



# L1 distance
#l1a=(imbalance(group = x$T, data = x[c("w1", "w2", "w3", "w4", "w5","w6", "w7", "w8", "w9", "w10")],weights=weights)$L1)[1]
#l1a_main<-mean((data.frame(imbalance(group = x$T, data = x[c("w1", "w2", "w3", "w4", "w5","w6", "w7", "w8", "w9", "w10")],weights=weights)[1]))[,3]) 

# variance ratio
varratio_a =abs(log(wtd.var(ps[x$T==1],weights=weights[x$T==1])/wtd.var(ps[x$T==0],weights=weights[x$T==0])))  # log (var ratio) after weighting; reference in Inbens-Rubin Ch 10, eq. 14.5 and tobacco litigation




#~~ performance metrics: balance AFTER matching for w1 ...w10

bam = MatchBalance(T ~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10
+w1*w2+w1*w3+w1*w4+w1*w5+w1*w6+w1*w7+w1*w8+w1*w9+w1*w10
+w2*w3+w2*w4+w2*w5+w2*w6+w2*w7+w2*w8+w2*w9+w2*w10
+w3*w4+w3*w5+w3*w6+w3*w7+w3*w8+w3*w9+w3*w10
+w4*w5+w4*w6+w4*w7+w4*w8+w4*w9+w4*w10
+w5*w6+w5*w7+w5*w8+w5*w9+w5*w10
+w6*w7+w6*w8+w6*w9+w6*w10
+w7*w8+w7*w9+w7*w10
+w8*w9+w8*w10
+w9*w10, data=allmdata, ks = FALSE, nboots=0, print.level=0) 

# ASAM after matching
asb_am    <-vector()
for(i in 1:10){asb_am[[i]] <- bam$BeforeMatching[[i]]$sdiff}

# ASAM with interactions after matching
interasbam<-vector()
for(i in 1:(length(bam$BeforeMatching))){interasbam[[i]] <- bam$BeforeMatching[[i]]$sdiff}

#vecKLdivam<-vector();for(i in 1:10){vecKLdivam[[i]] <-KLdiv(cbind(allmdata[,i][allmdata$T==1],allmdata[,i][allmdata$T==0]))[1,2]}

# mean difference of empirical quantiles after matching. 
# qqmeanrawam for W1=abs(bam$BeforeMatching[[1]]$qqsummary.raw$meandiff) 
#An approximation is: abs(mean(wtd.quantile(x[x$T==1,]$w1,probs=seq(0.01,1-0.01,0.01))-wtd.quantile(x[x$T==0,]$w1,probs=seq(0.01,1-0.01,0.01))))
vecqqmeanrawam    <-vector()
for(i in 1:10){vecqqmeanrawam[[i]] <- bam$BeforeMatching[[i]]$qqsummary.raw$meandiff}

# mean diff in ecdf after matching
# eg for W1: abs(bam$BeforeMatching[[1]]$qqsummary$meandiff)
vecqqmeanam    <-vector()
for(i in 1:10){vecqqmeanam[[i]] <- bam$BeforeMatching[[i]]$qqsummary$meandiff}

interqqam<-vector()
for(i in 1:(length(bam$BeforeMatching))){interqqam[[i]] <- bam$BeforeMatching[[i]]$qqsummary$meandiff}

# standardized max difference of ecdfs after matching
#qqmaxam for W1=abs(bam$BeforeMatching[[1]]$qqsummary$maxdiff) : max(abs(sxt-sxc)) where sxt = F_{W1 | T=1} (w1) and sxc = F_{W1 | T=0} (w1)
vecqqmaxam    <-vector()
for(i in 1:10){vecqqmaxam[[i]] <- bam$BeforeMatching[[i]]$qqsummary$maxdiff}

interqqmaxam<-vector()
for(i in 1:(length(bam$BeforeMatching))){interqqmaxam[[i]] <- bam$BeforeMatching[[i]]$qqsummary$maxdiff} 

#l1am=(imbalance(group = allmdata$T, data = allmdata[c("w1", "w2", "w3", "w4", "w5","w6", "w7", "w8", "w9", "w10")])$L1)[1]
 
 # variance ratio
 varratio_am =abs(log(var(ps[rr$index.treated])/var(ps[rr$index.control])))  # abs log (var ratio) after matching


#  mean (over simulated datasets) of balance summaries : 


masb_b<-mean(abs(asb_b))
masb_a<-mean(abs(asb_a))
masb_am<-mean(abs(asb_am))

masbinter_b<-mean(abs(interasbb))
masbinter_a<-mean(abs(interasba))
masbinter_am<-mean(abs(interasbam))

over20masb_b<-sum(abs(asb_b[which(abs(asb_b)>20)])-20)/length(asb_b)
over20masb_a<-sum(abs(asb_a[which(abs(asb_a)>20)])-20)/length(asb_a)
over20masb_am<-sum(abs(asb_am[which(abs(asb_am)>20)])-20)/length(asb_am)

over10masb_b<-sum(abs(asb_b[which(abs(asb_b)>10)])-10)/length(asb_b)
over10masb_a<-sum(abs(asb_a[which(abs(asb_a)>10)])-10)/length(asb_a)
over10masb_am<-sum(abs(asb_am[which(abs(asb_am)>10)])-10)/length(asb_am)

#KLdivb   <-mean(vecKLdivb)
#KLdivaw   <-mean(vecKLdivaw)
#KLdivam  <-mean(vecKLdivam)


qqmeanrawb<-mean(vecqqmeanrawb)
qqmeanrawaw<- mean(vecqqmeanrawaw)
qqmeanrawam<-mean(vecqqmeanrawam)

qqmeanb<-mean(vecqqmeanb)
qqmeanaw<-mean(vecqqmeanaw)
qqmeanam<-mean(vecqqmeanam)
qqmeaninterb<-mean(abs(interqqb))
qqmeaninteram<-mean(abs(interqqam))

qqmaxb<-mean(vecqqmaxb)
qqmaxaw<-mean(vecqqmaxaw)
qqmaxam<-mean(vecqqmaxam)
qqmaxinterb<-mean(abs(interqqmaxb))
qqmaxinteram<-mean(abs(interqqmaxam))

    
#~~collect output

output<-list( "auc"=auc,
"hatgw"=hatg,"absrbiasw"=absrbias,"varhatgw"=varhatg, #se di sim
"hatgm"=hatgm,"absrbiasm"=absrbiasm,"varhatgm"=varhatgm,
"hatgsem"=hatgsem,"hatgsew"=hatgsew,# se stimati
"covw"=covw,"covm"=covm,
"masb_b"=masb_b,"masb_aw"=masb_a,"masb_am"=masb_am,
"over20masb_b"=over20masb_b,"over20masb_aw"=over20masb_a,"over20masb_am"=over20masb_am,
"over10masb_b"=over10masb_b,"over10masb_aw"=over10masb_a,"over10masb_am"=over10masb_am,
"masbinter_b"=masbinter_b,"masbinter_aw"=masbinter_a,"masbinter_am"=masbinter_am,
#"KLdivb"=KLdivb,"KLdivaw"=KLdivaw,"KLdivam"=KLdivam,
"qqmeanraw_b"=qqmeanrawb,"qqmeanraw_aw"=qqmeanrawaw,"qqmeanraw_am"=qqmeanrawam,
"qqmean_b"=qqmeanb,"qqmean_aw"=qqmeanaw,"qqmean_am"=qqmeanam,
"qqmeaninter_b"=qqmeaninterb,"qqmeaninter_am"=qqmeaninteram,
"qqmax_b"=qqmaxb,"qqmax_aw"=qqmaxaw,"qqmax_am"=qqmaxam,
"qqmaxinter_b"=qqmaxinterb,"qqmaxinter_am"=qqmaxinteram,
"varratio_b"=varratio_b,"varratio_aw"=varratio_a,"varratio_am"=varratio_am,
#
"match.tnobs"=match.tnobs,"match.nobs"=match.nobs)
 return(output)
 }
 



 

  