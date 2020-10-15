library(tidyverse)
library(lavaan)
library(furrr)
library(tictoc)

###216 Conditions: 
#3(Vxs) X 3(Vys) X 3(CovS) X 2(ARx) X 2(ARy) X 2(CLdom)#
#slope conditions
x_slp_var<-c(1.8,5.4,9)
y_slp_var<-c(2.1,6.3,10.5)
xy_slp_corr<-c(0.4,0.6,0.8)
#AR conditions
ARx<-c(0.354,1)
ARy<-c(0.316,1)
#CL dom conditions
CLd<-c("0.05,0.05","0.05,0.1")
#number of iterations
N_iteration<-1:100

#setting conditions
cond<-crossing(
		x_slp_var,
		y_slp_var,
		xy_slp_corr,
		ARx,ARy,CLd,
		N_iteration
)
cond<-
cond %>%
	separate(CLd,c("CLyx","CLxy"),sep=",") %>%
	mutate(xy_slp_cov=xy_slp_corr*sqrt(x_slp_var*y_slp_var)) %>%
	mutate(larx=factor(ARx,levels=c(0.32,1),labels=c("Stationary","RandomWalk")))%>%
#<<<<SNIPPING out repetitive condition setting etc>>>>>>>>  
  mutate(YiXs_cov=0.1*sqrt(x_slp_var*73.1)) 
cond$CLxy<-as.numeric(cond$CLxy)
cond$CLyx<-as.numeric(cond$CLyx)
cond
#######Conditions Set###########

#############Function for Data Generating Models###########
##Structured Residuals
gen_sr<- 
function(ARx,ARy,CLxy,CLyx,xy_slp_cov,y_slp_var,x_slp_var,XiXs_cov,YiYs_cov,XiYs_cov,YiXs_cov,XiYi_cov,N_iteration){
model_sr<-  
'
<<<snipping details spelling out lavaan object to be used for generating>>>
	 '
#swap out param strings above with values from conditions
model_sr<-str_replace(model_sr,"ARx21",as.character(ARx))
#<<<<snipping out repetitive swap outs>>>>>>>>>>

#simulate datasets
data<-
	simulateData(model_sr,sample.nobs=3109, debug=FALSE)
 return(data)
}

#using pmap to apply data generating function to conditions
#output simulated datasets into a tibble
tic()
seeDatasets<-
  cond %>%
  group_by(ARx,ARy,CLxy,CLyx,y_slp_var,x_slp_var,xy_slp_cov,xy_slp_corr,N_iteration) %>%
  mutate(
    datasets=pmap(list(ARx,ARy,CLxy,CLyx,xy_slp_cov,y_slp_var,x_slp_var,XiXs_cov,YiYs_cov,XiYs_cov,YiXs_cov,XiYi_cov,N_iteration),
                  possibly(gen_sr,NA))
  )
toc()

#fitting function; validate by fitting LGCM-SR to LGCM-SR
validate_LGCM_SR<-function(df){
  #validation model
  lgcm_sr<- 
      '
<<<<snipping spell out of lavaan object for fitting
  '

    #fit the model lgcm-sr and compile results
fit_lgcm_sr<-sem(lgcm_sr,data=df,estimator="ML")
est_arx<-parameterEstimates(fit_lgcm_sr)[49,5]
se_arx<-parameterEstimates(fit_lgcm_sr)[49,6]
#<<<<<<<<snipping excess of identifying parameter estimates to compile for analysis>>>>

statNse<-
c(est_arx,se_arx,est_clyx,se_clyx,est_ary,se_ary,est_clxy,se_clxy,
#<<<<<<<<snip all named estimates to compile>>>>
est_My1,est_My2,est_My3,est_My4,est_My5,est_My6)
return(statNse)
}

#running simulation (automation piece)
plan(multiprocess)
tic()
results<-
	cond %>%
	group_by(ARx,ARy,CLxy,CLyx,y_slp_var,x_slp_var,xy_slp_cov,N_iteration) %>%
	mutate(
		datasets=pmap(list(ARx,ARy,CLxy,CLyx,xy_slp_cov,y_slp_var,x_slp_var,XiXs_cov,YiYs_cov,XiYs_cov,YiXs_cov,XiYi_cov,N_iteration),
		possibly(gen_sr,NA))
		)
results$values<-future_map(results$datasets,possibly(validate_LGCM_SR,NA))
toc()

#compile results in tidy way
results_tidy<-results %>%
	dplyr::select(-(datasets)) %>%
      unnest () %>%
	mutate(Estimates=c("est_arx","se_arx","est_clyx","se_clyx","est_ary","se_ary","est_clxy","se_clxy",
		#<<<<<<<<<<<<snip organizing estimates from results>>>>>>>>>>>>
		"est_My1","est_My2","est_My3","est_My4","est_My5","est_My6")) %>%
		separate(Estimates, c("Estimates","Parameter")) %>%
		spread(Estimates,values)
results_tidy	

summary<-results_tidy %>%
	group_by(ARx,ARy,CLxy,CLyx,xy_slp_cov,y_slp_var,x_slp_var,Parameter) %>%
	summarize(meanEST=mean(est),meanSE=mean(se),SEE=sd(est)) %>%
	mutate(biasSE=(meanSE-SEE)/SEE) 

MRBarx<-
summary %>%
  filter(Parameter=="arx") %>%
  mutate(bias=(meanEST-ARx)/ARx)

#<<<SNIPPING mean relative bias terms being created for analysis>>>>

MRBrescov<-
  summary %>%
  filter(Parameter=="ResCov") %>%
  mutate(bias=(meanEST-5.12)/5.12)

FullBias<-bind_rows(
		MRBarx,MRBary,MRBclyx,MRBclxy,
		MRBvix,MRBvsx,MRBviy,MRBvsy,
		MRBmx1,MRBmx2,MRBmx3,MRBmx4,MRBmx5,MRBmx6,
		MRBmy1,MRBmy2,MRBmy3,MRBmy4,MRBmy5,MRBmy6,
		MRBrescov,MRBslpcov 
			)
View(FullBias)

