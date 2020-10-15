
#########################################################################################

# In this file we analyzie weights distribution

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

############## function to generate covariates: funcov (size,scenario)

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


 ############## function to find weights: funsim_w (x, psmethod)
 
 ############## funsim is a simplified versione of funsim used in simulation having only weights in the output

# inputs: x = a simulated dataset created by funcov, a psmethod
# output: a list containing weights and treatment indicator
 
 
 
 funsim_w <- function(x,psmethod,par="ATE"){

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
    else if (psmethod == "nn")    {	
    	                         
   mod= nnet(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10,
   data=x, entropy=T,size=10,maxit=2000,trace=F)
   ps = as.numeric(predict(mod, type='raw')) #nb: anche predizioni fuori da [0,1]
   ps=exp(ps)/(1+exp(ps))          
                               }
    else if (psmethod == "nb")     {
    
	 mod = naiveBayes(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, data=x)
      ps = predict(mod,newdata=x,type="raw")[,2] 
}

#### separate ps

ps1<-ps[x$T==1]
ps0<-ps[x$T==0]

#### fit measure for ps

#auc<-auc(sensitivity(ps, factor(x$T), perc.rank = TRUE))

#### true (sample) att:

g <- mean(x$indeff[x$T==1]) 

### true (sample) ate

if(par=="ATE") {g <- mean(x$indeff) }

#### find weights for ATT

weights     <- ifelse(x$T==1,1,ps/(1-ps))
 weights0     <- weights[x$T==0]
 weights1     <- weights[x$T==1]

#### estimating weights for ATE

if(par=="ATE") {
 weights     <- ifelse(x$T==1,1/ps,1/(1-ps))
 weights0     <- weights[x$T==0]
 weights1     <- weights[x$T==1]
 }



 
#~~collect output


output <-list( 
"weights"=weights,"T"=x$T
)

 return(output)
 }
 
 
 
 ############ funsim_wo: funsim semplificata con soli weights in output 
 # e parametri di tuning ottimali per lo scenario Aa

  funsim_wo <- function(x,psmethod,par="ATT",tunopt=NULL){

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
    cp=0.002,
    method="class",data=x)
	 ps = predict(mod)[,2]     
        
}  else if (psmethod == "randomforest"){
	
      	mod = randomForest(factor(T)~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, 
      	mtry=2, 
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
   size=5,
   decay=2,maxit=2000,trace=F)
   ps = as.numeric(predict(mod, type='raw')) 
   ps=exp(ps)/(1+exp(ps))
   # mod for psw
      mod= nnet(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10,
   data=x, entropy=T,
   size=tunopt[5],
   decay=tunopt[6],maxit=2000,trace=F)
   psw = as.numeric(predict(mod, type='raw')) 
   psw=exp(psw)/(1+exp(psw))          
        
                               }
    else if (psmethod == "nb")     {
    
	 mod = naiveBayes(T~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, laplace=0,
	 data=x)
      ps = predict(mod,newdata=x,type="raw")[,2] 
}

#### fit measure for ps

#auc<-auc(sensitivity(ps, factor(x$T), perc.rank = TRUE))

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


output <-list( 
"weights"=weights,"T"=x$T
)

 return(output)
 }
 
 
 

  
  ############ simulation study: generates weights under different scenarios
  
  
  
  ########### scenarios
  scenarios<-c("Aa" ,"Ad", "Ae", "Ag", "Ca", "Cd", "Ce", "Cg","Fa", "Fd", "Fe", "Fg","Ga", "Gd", "Ge", "Gg")   
  
  ######## parameter
  par<-"ATE"
 
  ########### performance metrics

perfmetrics<-names(
funsim_w( funcov(500,"A","a"),"logit"))
 

# using default parameters

R   <- 1       #number of replications (keep 1 to separate plot of weights ofr T=1 and T=0)

size<- 500 

## ps model results
# logit

 lgresults <- matrix(0,nrow=size,ncol=length(perfmetrics))
 timestart<-Sys.time()
  for(i in 1:length(scenarios)){
	  	                   pb <- txtProgressBar(min = 0, max = length(scenarios), style = 3) 
  	                          partialresults<-matrix(0,nrow=size,ncol=length(perfmetrics))
  	                          for(j in 1:R)
                                             {
                partialresults <-partialresults+try(as.matrix(simplify2array(funsim_w(
                funcov(size,substr(scenarios[i],1,1),substr(scenarios[i],2,2)),"logit",par))))
   setTxtProgressBar(pb, i)  
                                              }
 lgresults<-cbind(lgresults,partialresults/R)
                               }
lgresults<-lgresults[,-c(1,2)]


###
# tree
 treeresults <- matrix(0,nrow=size,ncol=length(perfmetrics))
  for(i in 1:length(scenarios)){
	  	                   pb <- txtProgressBar(min = 0, max = length(scenarios), style = 3) 
  	                          partialresults<-matrix(0,nrow=size,ncol=length(perfmetrics))
  	                          for(j in 1:R)
                                             {
                partialresults <-partialresults+try(as.matrix(simplify2array(funsim_w(
                funcov(size,substr(scenarios[i],1,1),substr(scenarios[i],2,2)),"tree",par))))
   setTxtProgressBar(pb, i)  
                                              }
 treeresults<-cbind(treeresults,partialresults/R)
                               }
treeresults<-treeresults[,-c(1,2)]


###
# rf

rfresults <- matrix(0,nrow=size,ncol=length(perfmetrics))
 for(i in 1:length(scenarios)){
	  	                   pb <- txtProgressBar(min = 0, max = length(scenarios), style = 3) 
  	                          partialresults<-matrix(0,nrow=size,ncol=length(perfmetrics))
  	                          for(j in 1:R)
                                             {
                partialresults <-partialresults+try(as.matrix(simplify2array(funsim_w(
                funcov(size,substr(scenarios[i],1,1),substr(scenarios[i],2,2)),"tree",par))))
   setTxtProgressBar(pb, i)  
                                              }
 rfresults<-cbind(rfresults,partialresults/R)
                               }
rfresults<-rfresults[,-c(1,2)]


###
# bag
 bagresults <- matrix(0,nrow=size,ncol=length(perfmetrics))
  for(i in 1:length(scenarios)){
	  	                   pb <- txtProgressBar(min = 0, max = length(scenarios), style = 3) 
  	                          partialresults<-matrix(0,nrow=size,ncol=length(perfmetrics))
  	                          for(j in 1:R)
                                             {
                partialresults <-partialresults+try(as.matrix(simplify2array(funsim_w(
                funcov(size,substr(scenarios[i],1,1),substr(scenarios[i],2,2)),"bag",par))))
   setTxtProgressBar(pb, i)  
                                              }
 bagresults<-cbind(bagresults,partialresults/R)
                               }
bagresults<-bagresults[,-c(1,2)]



###
# nb
 nbresults <- matrix(0,nrow=size,ncol=length(perfmetrics))
  for(i in 1:length(scenarios)){
	  	                   pb <- txtProgressBar(min = 0, max = length(scenarios), style = 3) 
  	                          partialresults<-matrix(0,nrow=size,ncol=length(perfmetrics))
  	                          for(j in 1:R)
                                             {
                partialresults <-partialresults+try(as.matrix(simplify2array(funsim_w(
                funcov(size,substr(scenarios[i],1,1),substr(scenarios[i],2,2)),"nb",par))))
   setTxtProgressBar(pb, i)  
                                              }
 nbresults<-cbind(nbresults,partialresults/R)
                               }
nbresults<-nbresults[,-c(1,2)]


### 
# nn
 nnresults <- matrix(0,nrow=size,ncol=length(perfmetrics))
  for(i in 1:length(scenarios)){
	  	                   pb <- txtProgressBar(min = 0, max = length(scenarios), style = 3) 
  	                          partialresults<-matrix(0,nrow=size,ncol=length(perfmetrics))
  	                          for(j in 1:R)
                                             {
if(is.numeric(try(as.matrix(simplify2array(funsim_w(
    funcov(size,substr(scenarios[i],1,1),substr(scenarios[i],2,2)),"nn",par)))))==TRUE)
 {partialresults<- partialresults+try(as.matrix(simplify2array(funsim_w(
 	funcov(size,substr(scenarios[i],1,1),substr(scenarios[i],2,2)),"nn",par))))}else{partialresults<-partialresults}
    setTxtProgressBar(pb, i)  
                                              }
 nnresults<-cbind(nnresults,partialresults/R)
                               }
nnresults<-nnresults[,-c(1,2)]



###
# tw
 twresults <- matrix(0,nrow=size,ncol=length(perfmetrics))
  for(i in 1:length(scenarios)){
	  	                   pb <- txtProgressBar(min = 0, max = length(scenarios), style = 3) 
  	                          partialresults<-matrix(0,nrow=size,ncol=length(perfmetrics))
  	                          for(j in 1:1)
                                             {
                partialresults <-partialresults+try(as.matrix(simplify2array(
                funsim_w(
                funcov(size,substr(scenarios[i],1,1),substr(scenarios[i],2,2)),"gbmtwang",par))))
   setTxtProgressBar(pb, i)  
                                              }
 twresults<-cbind(twresults,partialresults/R)
                               }
twresults<-twresults[,-c(1,2)]

#### save workspace

save.image(file=file.path("results",paste("weights R",R,"size",size," par",par,".Rdata",sep="")))


# merge all results
longt<-rbind(lgresults,rfresults,twresults,bagresults,treeresults,nnresults,nbresults)
longt<-longt[,-seq(2,30,2)]
# column psmodel 
colnames(longt)<-c(scenarios,"T")
psmodel<-c("lg","rf","tw","bag","tree","nn","nb")
psmodel<-rep(psmodel,each=dim(lgresults)[1])
# add columns
simres<-data.frame(psmodel,longt)
# modify order of labels
simres$psmodel<-factor(simres$psmodel,levels=c("lg","rf","tw","bag","tree","nn","nb"))

# # set graphical parameters
# setEPS()
# postscript("results/Fig_weights.eps",width=13,height=10)
# pdf(paste("results/Fig_weights",size,".pdf",sep=""),width=13,height=10)
# set graphical window
par(mfrow=c(3,3),mai= c(0.3, 0.6, 0.5, 0.6)) #3 times 3 window (9 boxplots)
par(oma=c(0,4,4,0))#makes room for overall X and Y axis labels
#par(mar = c(2, 2, 1, 1))#makes the plot closer 
# loop on psmodel
for (i in 1:length(unique(psmodel))){
# select data for plot i
subdata <- subset(simres,psmodel==unique(psmodel)[i],select=c("Aa","T"))
head(subdata)
# plot
boxplot(subdata[ , "Aa" ] ~ subdata[,"T"],horizontal=TRUE,lwd=2,col=c("lightgreen","orange"),ylim=c(0,10),
ylim=c(0,max(subdata[ , "Aa" ])),#names.arg=levels(subdata$psmodel),
# treatment label on left side
ylab=if (isTRUE(i==1|i==4|i==7)){"treatment"}else{NULL},font.lab=2,
# scenario label on main
main=if (i<8) {unique(psmodel)[i]},legend=NULL)
abline(v=c(0,1,2,3,5,10,max(subdata[ , "Aa" ])),col="gray",lty="dotted")
abline(v=c( 
 mean(subdata[subdata$T==0,][,"Aa"]),
 mean(subdata[subdata$T==1,][,"Aa"])),col=c("lightgreen","orange"),lwd=2,lty="dotted")
abline(h=1.5,col="gray",lty="solid")
#mtext('....', side = 3, outer = TRUE, line = 2)
#mtext('....', side = 2, outer = TRUE, line = 2)
}# end loop graph
#dev.off()


# cycle as above but optimal tuning per tree, rf and nn
# c(mtry, cp, size, decay, size, decay)
  tunopt<-list( c(2,0.002,5,2,1,0.007) , c(5,0.02,1,0.0001,1,0.0001)  , c(1,0.0012,5,0.38,5,0.38),
 c(5,0.0015,11,3.17,1,0.001),
 c(2,0.037,3,0.52,3,0.52),c(2,0.001,7,2.77,13,0.43),c(2,0.034,1,3,15,0.83),
 c(2,0.034,3,0.52,3,0.52),c(3,0.002,1,0.02,9,0.35), c(3,0.002,1,0.02,9,0.35), c(3,0.0023,11,1.23,9,0.35),
 c(2,0.015,3,0.05,3,0.05),c(2,0.03,3,0.26,3,0.26), c(2,0.05,1,0.3,1,0.3), c(2,0.05,1,0.3,1,0.3),
 c(7,0.002,1,0.005,1,0.005))
## ps model results
# logit
 
 lgresults <- matrix(0,nrow=size,ncol=length(perfmetrics))
 timestart<-Sys.time()
  for(i in 1:length(scenarios)){
	  	                   pb <- txtProgressBar(min = 0, max = length(scenarios), style = 3) 
  	                          partialresults<-matrix(0,nrow=size,ncol=length(perfmetrics))
  	                          for(j in 1:R)
                                             {
                partialresults <-partialresults+try(as.matrix(simplify2array(funsim_wo(
                funcov(size,substr(scenarios[i],1,1),substr(scenarios[i],2,2)),"logit",par,tunopt[[i]]))))
   setTxtProgressBar(pb, i)  
                                              }
 lgresults<-cbind(lgresults,partialresults/R)
                               }
lgresults<-lgresults[,-c(1,2)]


###
# tree
 treeresults <- matrix(0,nrow=size,ncol=length(perfmetrics))
  for(i in 1:length(scenarios)){
	  	                   pb <- txtProgressBar(min = 0, max = length(scenarios), style = 3) 
  	                          partialresults<-matrix(0,nrow=size,ncol=length(perfmetrics))
  	                          for(j in 1:R)
                                             {
                partialresults <-partialresults+try(as.matrix(simplify2array(funsim_wo(
                funcov(size,substr(scenarios[i],1,1),substr(scenarios[i],2,2)),"tree",par,tunopt[[i]]))))
   setTxtProgressBar(pb, i)  
                                              }
 treeresults<-cbind(treeresults,partialresults/R)
                               }
treeresults<-treeresults[,-c(1,2)]


###
# rf

rfresults <- matrix(0,nrow=size,ncol=length(perfmetrics))
 for(i in 1:length(scenarios)){
	  	                   pb <- txtProgressBar(min = 0, max = length(scenarios), style = 3) 
  	                          partialresults<-matrix(0,nrow=size,ncol=length(perfmetrics))
  	                          for(j in 1:R)
                                             {
                partialresults <-partialresults+try(as.matrix(simplify2array(funsim_wo(
                funcov(size,substr(scenarios[i],1,1),substr(scenarios[i],2,2)),"tree",par,tunopt[[i]]))))
   setTxtProgressBar(pb, i)  
                                              }
 rfresults<-cbind(rfresults,partialresults/R)
                               }
rfresults<-rfresults[,-c(1,2)]


###
# bag
 bagresults <- matrix(0,nrow=size,ncol=length(perfmetrics))
  for(i in 1:length(scenarios)){
	  	                   pb <- txtProgressBar(min = 0, max = length(scenarios), style = 3) 
  	                          partialresults<-matrix(0,nrow=size,ncol=length(perfmetrics))
  	                          for(j in 1:R)
                                             {
                partialresults <-partialresults+try(as.matrix(simplify2array(funsim_wo(
                funcov(size,substr(scenarios[i],1,1),substr(scenarios[i],2,2)),"bag",par,tunopt[[i]]))))
   setTxtProgressBar(pb, i)  
                                              }
 bagresults<-cbind(bagresults,partialresults/R)
                               }
bagresults<-bagresults[,-c(1,2)]



###
# nb
 nbresults <- matrix(0,nrow=size,ncol=length(perfmetrics))
  for(i in 1:length(scenarios)){
	  	                   pb <- txtProgressBar(min = 0, max = length(scenarios), style = 3) 
  	                          partialresults<-matrix(0,nrow=size,ncol=length(perfmetrics))
  	                          for(j in 1:R)
                                             {
                partialresults <-partialresults+try(as.matrix(simplify2array(funsim_wo(
                funcov(size,substr(scenarios[i],1,1),substr(scenarios[i],2,2)),"nb",par,tunopt[[i]]))))
   setTxtProgressBar(pb, i)  
                                              }
 nbresults<-cbind(nbresults,partialresults/R)
                               }
nbresults<-nbresults[,-c(1,2)]


### 
# nn
 nnresults <- matrix(0,nrow=size,ncol=length(perfmetrics))
  for(i in 1:length(scenarios)){
	  	                   pb <- txtProgressBar(min = 0, max = length(scenarios), style = 3) 
  	                          partialresults<-matrix(0,nrow=size,ncol=length(perfmetrics))
  	                          for(j in 1:R)
                                             {
if(is.numeric(try(as.matrix(simplify2array(funsim_wo(
    funcov(size,substr(scenarios[i],1,1),substr(scenarios[i],2,2)),"nn",par,tunopt[[i]])))))==TRUE)
 {partialresults<- partialresults+try(as.matrix(simplify2array(funsim_wo(
 	funcov(size,substr(scenarios[i],1,1),substr(scenarios[i],2,2)),"nn",par,tunopt[[i]]))))}else{partialresults<-partialresults}
    setTxtProgressBar(pb, i)  
                                              }
 nnresults<-cbind(nnresults,partialresults/R)
                               }
nnresults<-nnresults[,-c(1,2)]



###
# tw
 twresults <- matrix(0,nrow=size,ncol=length(perfmetrics))
  for(i in 1:length(scenarios)){
	  	                   pb <- txtProgressBar(min = 0, max = length(scenarios), style = 3) 
  	                          partialresults<-matrix(0,nrow=size,ncol=length(perfmetrics))
  	                          for(j in 1:1)
                                             {
                partialresults <-partialresults+try(as.matrix(simplify2array(
                funsim_wo(
                funcov(size,substr(scenarios[i],1,1),substr(scenarios[i],2,2)),"gbmtwang",par,tunopt[[i]]))))
   setTxtProgressBar(pb, i)  
                                              }
 twresults<-cbind(twresults,partialresults/R)
                               }
twresults<-twresults[,-c(1,2)]

#### save workspace

save.image(file=file.path("results",paste("opt_weights R",R,"size",size," par",par,".Rdata",sep="")))



######################
# graph weights distribution
#####################

# merge all results
longt<-rbind(lgresults,rfresults,twresults,bagresults,treeresults,nnresults,nbresults)
longt<-longt[,-seq(2,30,2)]
# column psmodel 
colnames(longt)<-c(scenarios,"T")
psmodel<-c("lg","rf","tw","bag","tree","nn","nb")
psmodel<-rep(psmodel,each=dim(lgresults)[1])
# add columns
simres<-data.frame(psmodel,longt)
# modify order of labels
simres$psmodel<-factor(simres$psmodel,levels=c("lg","rf","tw","bag","tree","nn","nb"))

# # set graphical parameters
pdf(paste("results/","Fig_opt_weights",size,".pdf"),width=13,height=10)
# set graphical window
par(mfrow=c(3,3),mai= c(0.3, 0.6, 0.5, 0.3)) #3 times 3 window (9 boxplots)
par(oma=c(0,4,4,0))#makes room for overall X and Y axis labels
#par(mar = c(2, 2, 1, 1))#makes the plot closer 

# loop on psmodel
for (i in 1:length(unique(psmodel))){
# select data for plot i
subdata <- subset(simres,psmodel==unique(psmodel)[i],select=c("Aa","T"))
head(subdata)
# plot
boxplot(subdata[ , "Aa" ] ~ subdata[,"T"],horizontal=TRUE,lwd=2,col=c("lightgreen","orange"),ylim=c(0,10),
ylim=c(0,max(subdata[ , "Aa" ])),#names.arg=levels(subdata$psmodel),
# treatment label on left side
ylab=if (isTRUE(i==1|i==4|i==7)){"treatment"}else{NULL},font.lab=2,
# scenario label on main
main=if (i<8) {unique(psmodel)[i]},legend=NULL)
abline(v=c(0,1,2,3,5,10,max(subdata[ , "Aa" ])),col="gray",lty="dotted")
abline(v=c( 
 mean(subdata[subdata$T==0,][,"Aa"]),
 mean(subdata[subdata$T==1,][,"Aa"])),col=c("lightgreen","orange"),lwd=2,lty="dotted")
abline(h=1.5,col="gray",lty="solid")
#mtext('....', side = 3, outer = TRUE, line = 2)
#mtext('....', side = 2, outer = TRUE, line = 2)
}# end loop graph

dev.off()

 

  