## Writen for R version 3.2.2
library(lqa) # version 1.0-3
library(MASS) # version 3.3.1
  
# set information for simulating coviariates
mean_x = 0 
sig_x = 1
rho = 0
# sample size
n = 500
# total number of predictors
p = 200
pC = pI = pP = 2
pS = p - (pC+pI+pP)
var.list = c(paste("Xc",1:pC,sep=""),paste("Xp",1:pP,sep=""),paste("Xi",1:pI,sep=""),paste("Xs",1:pS,sep=""))
# Set strength of relationship between covariates and outcome
beta_v =  c( 0.6, 0.6, 0.6, 0.6, 0, 0, rep(0,p-6) )
# Set strength of relationship between covariates and treatment
alpha_v = c( 1.0, 1.0,   0,   0, 1, 1,  rep(0,p-6) )
names(beta_v) = names(alpha_v) = var.list
### set true average treatment effect
bA = 0

# set vector of possible lambda's to try
lambda_vec = c( -10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49)
names(lambda_vec) = as.character(lambda_vec)
# lambda_n (n)^(gamma/2 - 1) = n^(gamma_convergence_factor)
gamma_convergence_factor = 2
# get the gamma value for each value in the lambda vector that corresponds to convergence factor
gamma_vals = 2*( gamma_convergence_factor - lambda_vec + 1 )
names(gamma_vals) = names(lambda_vec)

### define some functions for generating data, ATE estimates, and the wAMD,
expit = function(x){ 
 pr = ( exp(x) / (1+exp(x)) ) 
 return(pr)
}
ATE_est = function(fY,fw,fA){
 t_ATE = fY*fw
 tt_ATE = ( ( sum(t_ATE[fA==1]) / sum(fw[fA==1]) ) - ( sum(t_ATE[fA==0]) /  sum(fw[fA==0]) ) )
 return(tt_ATE) 
}
create_weights = function(fp,fA,fw){
 fw = (fp)^(-1)
 fw[fA==0] = (1 - fp[fA==0])^(-1)
 return(fw)
}
wAMD_function = function(DataM,varlist,trt.var,wgt,beta){
 trt = untrt = diff_vec = rep(NA,length(beta)) 
 names(trt) = names(untrt) = names(diff_vec) = varlist
 for(jj in 1:length(varlist)){ 
  this.var = paste("w",varlist[jj],sep="") 
  DataM[,this.var] = DataM[,varlist[jj]] * DataM[,wgt] 
  trt[jj] = sum( DataM[DataM[,trt.var]==1, this.var ]) / sum(DataM[DataM[,trt.var]==1, wgt]) 
  untrt[jj] = sum(DataM[DataM[,trt.var]==0, this.var]) / sum(DataM[DataM[,trt.var]==0, wgt]) 
  diff_vec[jj] = abs( trt[jj] - untrt[jj] ) 
 } 
 wdiff_vec = diff_vec * abs(beta) 
 wAMD = c( sum(wdiff_vec))
 ret = list( diff_vec = diff_vec, wdiff_vec = wdiff_vec, wAMD = wAMD )
 return(ret) 
}

### simulate data
Sigma_x = matrix(rho*sig_x^2,nrow=length(var.list),ncol=length(var.list)) 
diag(Sigma_x) = sig_x^2
Mean_x = rep(mean_x,length(var.list))
Data = as.data.frame(mvrnorm(n = n,mu=Mean_x,Sigma = Sigma_x,empirical = FALSE))
names(Data) = var.list
gA_x = rowSums(Data[,var.list]*matrix(alpha_v,nrow=n,ncol=length(var.list),byrow=TRUE))
pA = expit( gA_x )
Data$A = as.numeric( runif(n=length(pA)) < pA) # simulate A 
gY_xA = rowSums(Data[,var.list]*matrix(beta_v,nrow=n,ncol=length(var.list),byrow=TRUE))   
Data$Y = gY_xA + rnorm(n=n,sd=sig_x)
Data$Y = Data$Y + Data$A*bA

# Normlize coviarates to have mean 0 and standard deviation 1
temp.mean = colMeans(Data[,var.list])
Temp.mean = matrix(temp.mean,ncol=length(var.list),nrow=nrow(Data),byrow=TRUE)
Data[,var.list] = Data[,var.list] - Temp.mean
temp.sd = apply(Data[var.list],FUN=sd,MARGIN=2)
Temp.sd = matrix(temp.sd,ncol=length(var.list),nrow=nrow(Data),byrow=TRUE)
Data[var.list] = Data[,var.list] / Temp.sd
rm(list=c("temp.mean","Temp.mean","temp.sd","Temp.sd"))

# estimate outcome model
y.form = formula(paste("Y~A+",paste(var.list,collapse="+")))
lm.Y = lm(y.form,data=Data)
betaXY = coef(lm.Y)[var.list] 

## Want to save ATE, wAMD and propensity score coefficients for each lambda value
ATE = wAMD_vec = rep(NA, length(lambda_vec))
names(ATE) = names(wAMD_vec) = names(lambda_vec)
coeff_XA = as.data.frame(matrix(NA,nrow=length(var.list),ncol=length(lambda_vec)))
names(coeff_XA) = names(lambda_vec)
rownames(coeff_XA) = var.list

######################################################################################
#####  Run outcome adaptive lasso for each lambda value 
######################################################################################
# weight model with all possible covariates included, this is passed into lasso function
w.full.form = formula(paste("A~",paste(var.list,collapse="+")))
for( lil in names(lambda_vec) ){
  il = lambda_vec[lil]
  ig = gamma_vals[lil]

  ### create the outcome adaptive lasso penalty with coefficient specific weights determined by outcome model
  oal_pen = adaptive.lasso(lambda=n^(il),al.weights = abs(betaXY)^(-ig) )
  ### run outcome-adaptive lasso model with appropriate penalty
  logit_oal = lqa( w.full.form, data=Data, penalty=oal_pen, family=binomial(logit) )
  # generate propensity score
  Data[,paste("f.pA",lil,sep="")] = predict(logit_oal)$mu.new
  # save propensity score coefficients
  coeff_XA[var.list,lil] = coef(logit_oal)[var.list]
  # create inverse probability of treatment weights
  Data[,paste("w",lil,sep="")] = create_weights(fp=Data[,paste("f.pA",lil,sep="")],fA=Data$A)
  # estimate weighted absolute mean different over all covaraites using this lambda to generate weights
  wAMD_vec[lil] = wAMD_function(DataM=Data,varlist=var.list,trt.var="A",
			wgt=paste("w",lil,sep=""),beta=betaXY)$wAMD
  # save ATE estimate for this lambda value
  ATE[lil] = ATE_est(fY=Data$Y,fw=Data[,paste("w",lil,sep="")],fA=Data$A)
} # close loop through lambda values

# print out wAMD for all the lambda values tried
wAMD_vec
# find the lambda value that creates the smallest wAMD
tt = which.min(wAMD_vec)
# print out ATE corresponding to smallest wAMD value
ATE[tt]
# print out the coefficients for the propensity score that corresponds with smalles wAMD value 
coeff_XA[,tt]