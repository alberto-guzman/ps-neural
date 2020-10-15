#load in libraries
library(MASS)
library(tidyverse)
library(lavaan)
library(MplusAutomation)
library(texreg)
library(furrr)
library(tictoc)

###########################################################
##########Structured Residual Generating Model#############
###########################################################

#tic() and toc() are simply tracking computing times for me
tic()
#set variable conditions,6 variables X 2 conditions
SxSyR<-c(0.1,0.5)
VSx<-c(3.44,13.76)
VSy<-c(4.22,16.88)
CLxy<-c(0.05,0.10)
ARx<-c(0.3,0.5)
ARy<-c(0.3,0.5)
#do 1000 replications
Niter<-1:1000

#fully crossing conditions for 12000 datasets
cond<-crossing(SxSyR,VSx,VSy,CLxy,ARx,ARy,Niter)	

#function for setting data generating models
SR_gen<-function(SxSyR,VSx,VSy,CLxy,ARx,ARy,Niter){

		#set mean structure
mu<-c(50.36,3.32,50.35,3.27,rep(0,12))

#set correlation matrix
R<-matrix(c(
  1,	0.3,	0.5,	0.1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
  0.3,	1,	0.1,	SxSyR,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
  0.5	, 0.1,	1,	0.3,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
  0.1	,SxSyR,	0.3,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
  0,	0	,0,	0,	1,	0,	0,	0,	0,	0,	0.2,	0,	0,	0,	0,	0,
  0,	0	,0,	0,	0,	1,	0,	0,	0,	0,	0,	0.2,	0,	0,	0,	0,
  0,	0,0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0.2,	0,	0,	0,
  0,	0	,0,	0	,0,	0,	0,	1,	0,	0,	0,	0,	0,	0.2,	0,	0,
  0,	0	,0,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0.2,	0,
  0,	0	,0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0.2,
  0,	0	,0,	0,	0.2,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0,
  0,	0	,0,	0,	0,	0.2,	0,	0	,0	,0	,0	,1	,0	,0	,0	,0,
  0,	0	,0,	0,	0,	0,	0.2,	0,	0,	0,	0,	0,	1,	0,	0,	0,
  0,	0	,0,	0,	0,	0,	0,	0.2,	0,	0,	0,	0,	0,	1,	0,	0,
  0,	0	,0,	0,	0,	0,	0,	0,	0.2,	0,	0,	0,	0,	0,	1,	0,
  0,	0	,0,	0,	0,	0,	0,	0,	0,	0.2,	0,	0,	0,	0,	0,	1
  
),16,16)

#create covariance from correlation
#Variances
V<-diag(c(
	72.08^.5,VSx^.5,78.90^.5,VSy^.5,
	28^.5,28^.5,28^.5,28^.5,28^.5,28^.5,
	25^.5,25^.5,25^.5,25^.5,25^.5,25^.5))

Sigma<-V%*%R%*%V

rownames(Sigma)<-c("Ix","Sx","Iy","Sy",
                   "Vx1","Vx2","Vx3","Vx4","Vx5","Vx6",
                   "Vy1","Vy2","Vy3","Vy4","Vy5","Vy6")
colnames(Sigma)<-c("Ix","Sx","Iy","Sy",
                   "Vx1","Vx2","Vx3","Vx4","Vx5","Vx6",
                   "Vy1","Vy2","Vy3","Vy4","Vy5","Vy6")

#raw dataframe from covariance matrix & mean structure
df<-mvrnorm(3109, mu, Sigma, tol = 1e-6, empirical = FALSE)

#create structured residuals dataframe
SR_df<-
df%>%as_tibble()%>%
	mutate(x1=Ix+Vx1) %>%
	mutate(y1=Iy+Vy1) %>%
	
	mutate(x2=(Ix+Sx)+(ARx*Vx1+CLxy*Vy1)+Vx2) %>%
	mutate(y2=(Iy+Sy)+(ARy*Vy1+0.05*Vx1)+Vy2) %>%
	
	mutate(x3=(Ix+2.2*Sx)+(ARx*Vx2+CLxy*Vy2)+Vx3) %>%
	mutate(y3=(Iy+2.3*Sy)+(ARy*Vy2+0.05*Vx2)+Vy3) %>%
	
	mutate(x4=(Ix+2.86*Sx)+(ARx*Vx3+CLxy*Vy3)+Vx4) %>%
	mutate(y4=(Iy+3.57*Sy)+(ARy*Vy3+0.05*Vx3)+Vy4) %>%
	
	mutate(x5=(Ix+3.49*Sx)+(ARx*Vx4+CLxy*Vy4)+Vx5) %>%
	mutate(y5=(Iy+4.35*Sy)+(ARy*Vy4+0.05*Vx4)+Vy5) %>%
	
	mutate(x6=(Ix+3.8*Sx)+(ARx*Vx5+CLxy*Vy5)+Vx6) %>%
	mutate(y6=(Iy+4.62*Sy)+(ARy*Vy5+0.05*Vx5)+Vy6) %>%

	subset(select=c(x1,x2,x3,x4,x5,x6,y1,y2,y3,y4,y5,y6))
return(SR_df)
}
toc()


#using pmap to apply data generating function across condtions
#to output simulated datasets into tibble collecting datasets
tic()
SR_simDats<-
  cond %>%
  group_by(SxSyR,VSx,VSy,CLxy,ARx,ARy,Niter) %>%
  mutate(
    datasets=pmap(list(SxSyR,VSx,VSy,CLxy,ARx,ARy,Niter),
                  possibly(SR_gen,NA))
  )
toc()

#better organize datasets with labels
tic()
tbl<-
SR_simDats%>%
  mutate(GenMod="LGCM_SR")%>%
	mutate(LSlpCov=factor(SxSyR,levels=c(0.1,0.5),labels=c("LoCov","HiCov")))%>%
	mutate(LVSx=factor(VSx,levels=c(3.44,13.76),labels=c("LoVSx","HiVSx")))%>%
	mutate(LVSy=factor(VSy,levels=c(4.22,16.88),labels=c("LoVSy","HiVSy")))%>%
	mutate(LCL=factor(CLxy,levels=c(0.05,0.10),labels=c("NoDom","Dom")))%>%		
	mutate(LARx=factor(ARx,levels=c(0.3,0.5),labels=c("AR_X_lo","AR_X_hi")))%>%
	mutate(LARy=factor(ARy,levels=c(0.3,0.5),labels=c("AR_Y_lo","AR_Y_hi")))%>%
	mutate(Cell=paste(LSlpCov,LVSx,LVSy,LCL,LARx,LARy))%>%
  group_by(Cell)

#allows individual datasets to be called into MPlus 
DataList<-group_split(tbl)
toc()

#######################################################################
##########Write out data and input files for MPlus to use##############
#######################################################################
tic()
for(k in 1:64){
  dir.create(paste0("C:/LGCM_SR/ConditionSet",k))}
for (j in 1:64){
      for (i in 1:1000){
        write.table(DataList[[j]]$datasets[i],
        paste0("C:/LGCM_SR/ConditionSet",j,"/","Rep",i,".dat"),
        sep="\t",row.names=FALSE,col.names=FALSE)
        write.table((paste0("Rep",i,".dat")),
                    paste0("C:/LGCM_SR/ConditionSet",j,"/Rep.dat"),
                    row.names=FALSE,col.names=FALSE,quote=FALSE)
      }
  write.table((paste0("Rep",1:1000,".dat")),
              paste0("C:/LGCM_SR/ConditionSet",j,"/Rep.dat"),
              row.names=FALSE,col.names=FALSE,quote=FALSE)

write.table((paste0(
    "TITLE: 
LGCM_SR fit to ",DataList[[j]]$GenMod[1],"  ",DataList[[j]]$Cell[1],
    
    "

DATA:
    FILE = C:/LGCM_SR/ConditionSet",j,"/Rep.dat;
  TYPE = MONTECARLO;
  
  VARIABLE:
    names=x1 x2 x3 x4 x5 x6 y1 y2 y3 y4 y5 y6;
  usevar=x1 x2 x3 x4 x5 x6 y1 y2 y3 y4 y5 y6;
  analysis: estimator=ml;
  model:
    !Set Intercept and Slope
  Xint Xslp | x1@0 x2@1 x3 x4 x5 x6; 
  Yint Yslp | y1@0 y2@1 y3 y4 y5 y6; 
  
  !set measurement error in observed vars to 0
  x1-x6@0; 
  y1-y6@0;
  
  !create residuals to structure
  Xres1 by x1@1; 
  Xres2 by x2@1; 
  Xres3 by x3@1; 
  Xres4 by x4@1;
  Xres5 by x5@1;
  Xres6 by x6@1;
  
  Yres1 by y1@1;
  Yres2 by y2@1;
  Yres3 by y3@1;
  Yres4 by y4@1;
  Yres5 by y5@1;
  Yres6 by y6@1;
  
  !AR structured resids
  Xres2-Xres6 pon Xres1-Xres5 (1);
  Yres2-Yres6 pon Yres1-Yres5 (2);
  
  !CL structured resids
  Xres2-Xres6 pon Yres1-Yres5 (3); 
  Yres2-Yres6 pon Xres1-Xres5 (4);
  
  !Innovation structured resids
  Xres1-Xres6 pwith Yres1-Yres6 (5);
  
  
  !Zero out Between & Within-level Exogenous Covariances
  Xint WITH Xres1@0 Yres1@0;
  Xslp WITH Xres1@0 Yres1@0;
  Yint WITH Xres1@0 Yres1@0;
  Yslp WITH Xres1@0 Yres1@0;
  
  output: 
    sampstat stdyx;
  "
)),paste0("C:/LGCM_SR/ConditionSet",j,"/LGCM_SR.inp"),
    row.names=FALSE,col.names=FALSE,quote=FALSE)

#<<<<<SNIPPING excess MPlus .inp templates>>>>>>>>>>

}
toc()

##################################################################
############Automate reading of batches from MPlus files##########
##################################################################
tic()
runModels("C:/LGCM_SR/", recursive=TRUE)
toc()

#################################################################
###########Bring in MPlus out file estimates into R##############
#################and prepare for analysis########################
#################################################################
tic()
for (m in 1:64){

###ALT###
read_alt<-readModels(paste0("C:/LGCM_SR/ConditionSet",m,"/alt.out"))
alt_params<-as.tibble(read_alt$parameters$unstandardized)
alt_ARCLest<-alt_params[c(11:12,21:22),c(1:2,4:6)]
alt_ARCLest<-alt_ARCLest%>%
            unite(Param,"paramHeader","param")%>%
            mutate(Param=factor(Param,
                                levels=c("X2.ON_X1","X2.ON_Y1","Y2.ON_Y1","Y2.ON_X1"),
                                labels=c("ARx","CLxy","ARy","CLyx")))%>%
            unite(pooledparams,"average","population_sd","average_se")%>%
          spread(Param,pooledparams)%>%
          separate(ARx,sep="_",into=c("ARx_est","ARx_sd","ARx_se"),convert=T)%>%
          separate(CLxy,sep="_",into=c("CLxy_est","CLxy_sd","CLxy_se"),convert=T)%>%
          separate(ARy,sep="_",into=c("ARy_est","ARy_sd","ARy_se"),convert=T)%>%
          separate(CLyx,sep="_",into=c("CLyx_est","CLyx_sd","CLyx_se"),convert=T)%>%
          mutate(FitMod="alt")%>%
          mutate(GenMod="LGCM_SR")
##do similar for fit stats (RMSEA,CFI,SRMR,BIC)
alt_FitInd<-as.tibble(read_alt$summaries[1,c(17,21,33,36,39)])

#cbind fit and "ARCL"
alt_FitInd_ARCL<-cbind(alt_FitInd,alt_ARCLest)

#cbind (fit&ARCL) with conditons & parameter values
alt_ModSum<-cbind(alt_FitInd_ARCL,
              DataList[[m]][1,c(4:6,9:14)])

#calculate Relative Bias for est & se to produce
#final 1 X 35 object per condition/loop
alt_ModSum<-
  alt_ModSum %>%
  mutate(bias_est_ARx=(ARx_est-ARx)/ARx,
        bias_est_ARy=(ARy_est-ARy)/ARy,
        bias_se_ARx=(ARx_se-ARx_sd)/ARx_sd,
        bias_se_ARy=(ARy_se-ARy_sd)/ARy_sd,
        bias_est_CLxy=(CLxy_est-CLxy)/CLxy,
        bias_est_CLyx=(CLyx_est-0.05)/0.05,
        bias_se_CLxy=(CLxy_se-CLxy_sd)/CLxy_sd,
        bias_se_CLyx=(CLyx_se-CLyx_sd)/CLyx_sd,
	  ConvRate=LL_NumComputations/10)

#<<<<SNIPPING same basic code repeated across 6 other models for bringing MPlus output into R>>>>

#bring output from different fitting models together
  #rbind all 7 model entries, resulting in  7 X 35 object
GeneralModSum<-rbind(
		lgcm_sr_ModSum,alt_ModSum,gclm_ModSum,
		ri_clpm_ModSum, alt_noslpvar_ModSum,gclm_meanstationary_ModSum,
		clpm_ModSum)		
write.table(GeneralModSum,
			paste0("C:/LGCM_SR/ConditionSet",m,"/Set4Analysis.dat"),
			sep="\t",row.names=FALSE)
}
toc()

#append tables from each loop for object with 7*64=448 rows X 37 columns
#use pmap to read external files into a tibble
file<-paste0("C:/LGCM_SR/ConditionSet",1:64,"/Set4Analysis.dat")
sep<-rep("\t",64)
header<-rep(TRUE,64)
args<-tibble(file,sep,header)
MastLGCM_SR<-
args%>%
	pmap(read.table)
#append all of the datasets together
MastDat<-do.call(rbind,MastLGCM_SR)

write.table(MastDat,"C:/LGCM_SR/SimResult_LGCM_SR.dat",
		sep="\t",row.names=FALSE)
#this object can then be rbinded to other generating model products

###############################################################
########################Analyze MCMC results###################
###############################################################

