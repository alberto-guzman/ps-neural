######################################################################
# This file generates dataset with
# a) real covariates  and
# b) simulated Treatment and Outcome using logit function
# Simulation is then carried on R replication of the datasets to examine performance of PSM and PSW

# load case study data
load("case_study/cs_data/clean_data_10cov.Rdata")
head(data)

# load functions 
source("functions.R")



######### SIMULATION using case-study covariates ############


# there is only one scenario:
 scenarios <-c("Tt") 

# performance metrics
   
   statistics<-names(
funsim_cs( funcov_cs(data,"T","t"),"logit"))

# number of replications
R   <- 500  

   
## ps models

### logit
 lgresults <- matrix(0,length(statistics),1)# mean statistics (over R)
 comp.lgresults<-list()# all statistics
 timestart<-Sys.time()
  for(i in 1:length(scenarios)){
	  	                   pb <- txtProgressBar(min = 0, max = R, style = 3) 
  	                          partialresults<-matrix(0,length(statistics),1)
  	                          comp.lgresults[[i]]<-matrix(0,length(statistics),R)
  	                          for(j in 1:R)
                                             {
                partialresults <-partialresults+try(as.matrix(simplify2array(funsim_cs( 
                funcov_cs(data,substr(scenarios[i],1,1),substr(scenarios[i],2,2))  ,"logit"))))
   setTxtProgressBar(pb, j)  
                comp.lgresults[[i]][,j]<-partialresults
                                              }
                lgresults<-cbind(lgresults,partialresults/R)
                           }
lgresults<-lgresults[,-1]#;colnames(lgresults)<-paste(scenarios)
#round(lgresults,2)
                               
timeend<-Sys.time();exectime<-timeend-timestart
exectime                              

#write.table(lgresults, file = file.path("cs_results",paste("logit","R",R,".txt",sep="")))


### tree

 treeresults <- matrix(0,length(statistics),1)
  for(i in 1:length(scenarios)){
	  	                   pb <- txtProgressBar(min = 0, max = R, style = 3) 
  	                          partialresults<-matrix(0,length(statistics),1)
  	                          for(j in 1:R)
                                             {
                         prova<-funcov_cs(data,substr(scenarios[i],1,1),substr(scenarios[i],2,2))
                partialresults <-partialresults+try(as.matrix(simplify2array(
                funsim_cs(prova, "tree"))))   
   setTxtProgressBar(pb, j)  
                                              }
 treeresults<-cbind(treeresults,partialresults/R)
                               }
treeresults<-treeresults[,-1]


#write.table(treeresults, file = file.path("cs_results",paste("tree","R",R,".txt",sep="")))

### random forest

rfresults <- matrix(0,length(statistics),1)
#timestart<-Sys.time()
 for(i in 1:length(scenarios)){
 	                      pb <- txtProgressBar(min = 0, max = R, style = 3)
 	                          #partialresults<-matrix(0,length(perfmetrics),1)
                                partialresults<-list()
	                          for(j in 1:R){
	                          	    setTxtProgressBar(pb, j)  
     prova<-funcov_cs(data,substr(scenarios[i],1,1),substr(scenarios[i],2,2))

 	partialresults[[j]]<- tryCatch(funsim_cs(prova,"randomforest"),
 	error=function(e) paste("rf error at ",j,sep=""))
                                             }
partialresults <- partialresults[lapply(partialresults,function(x) class(x))=="list"]
partialresults <- sapply(partialresults,simplify2array)
#partialresults <- sapply(partialresults,simplify2array)
     rfresults <- cbind(rfresults,apply(as.matrix(partialresults),1,mean))
                              }
rfresults<-rfresults[,-1]
#timeend<-Sys.time();exectime<i-timeend-timestart
#exectime

#write.table(rfresults, file = file.path("cs_results",paste("rf","R",R,".txt",sep="")))


### bag
 bagresults <- matrix(0,length(statistics),1)
  for(i in 1:length(scenarios)){
	  	                   pb <- txtProgressBar(min = 0, max = R, style = 3) 
  	                          partialresults<-matrix(0,length(statistics),1)
  	                          for(j in 1:R)
                                             {
                    prova<-funcov_cs(data,substr(scenarios[i],1,1),substr(scenarios[i],2,2))
                partialresults <-partialresults+try(as.matrix(simplify2array(
                funsim_cs(prova,"bag"))))
   setTxtProgressBar(pb, j)  
                                              }
 bagresults<-cbind(bagresults,partialresults/R)
                               }
bagresults<-bagresults[,-1]

#write.table(bagresults, file = file.path("cs_results",paste("bag","R",R,".txt",sep="")))

### nb
 nbresults <- matrix(0,length(statistics),1)
  for(i in 1:length(scenarios)){
	  	                   pb <- txtProgressBar(min = 0, max = R, style = 3) 
  	                          partialresults<-matrix(0,length(statistics),1)
  	                          for(j in 1:R)
                                             {
                                             	set.seed(7)
                         prova <- funcov_cs(data,substr(scenarios[i],1,1),substr(scenarios[i],2,2))
                partialresults <-partialresults+try(as.matrix(simplify2array(
                funsim_cs(prova,"nb"))))
   setTxtProgressBar(pb, j)  
                                              }
 nbresults<-cbind(nbresults,partialresults/R)
                               }
nbresults<-nbresults[,-1]

#write.table(nbresults, file = file.path("case_study/cs_results",paste("nb","R",R,".txt",sep="")))
### 
 nnresultsm <- matrix(0,length(statistics),1)
  for(i in 1:length(scenarios)){
	  	                   pb <- txtProgressBar(min = 0, max = R, style = 3) 
  	                          partialresults<-matrix(0,length(statistics),1)
  	                          for(j in 1:R)
                                             {
if(is.numeric(try(as.matrix(simplify2array(funsim_cs(
                funcov_cs(data,substr(scenarios[i],1,1),substr(scenarios[i],2,2)) ,"nnm")))))==TRUE)
 {partialresults<- partialresults+try(as.matrix(simplify2array(funsim_cs(
                funcov_cs(data,substr(scenarios[i],1,1),substr(scenarios[i],2,2))  ,"nnm"))))}
                else{partialresults<-partialresults}
    setTxtProgressBar(pb, j)  
                                              }
 nnresultsm<-cbind(nnresultsm,partialresults/R)
                               }
nnresultsm<-nnresultsm[,-1]
nnresults<-nnresultsm

###
### 
 nnresultsw <- matrix(0,length(statistics),1)
  for(i in 1:length(scenarios)){
	  	                   pb <- txtProgressBar(min = 0, max = R, style = 3) 
  	                          partialresults<-matrix(0,length(statistics),1)
  	                          for(j in 1:R)
                                             {
if(is.numeric(try(as.matrix(simplify2array(funsim_cs(
                funcov_cs(data,substr(scenarios[i],1,1),substr(scenarios[i],2,2)) ,"nnw")))))==TRUE)
 {partialresults<- partialresults+try(as.matrix(simplify2array(funsim_cs(
                funcov_cs(data,substr(scenarios[i],1,1),substr(scenarios[i],2,2))  ,"nnw"))))}
                else{partialresults<-partialresults}
    setTxtProgressBar(pb, j)  
                                              }
 nnresultsw<-cbind(nnresultsw,partialresults/R)
                               }
nnresults[c("hatgw","absrbiasw","varhatgw","hatgsew","covw","masb_aw","over20masb_aw","over10masb_aw","masbinter_aw","qqmeanraw_aw","qqmean_aw","qqmax_aw","varratio_aw")]<-nnresultsw[c("hatgw","absrbiasw","varhatgw","hatgsew","covw","masb_aw","over20masb_aw","over10masb_aw","masbinter_aw","qqmeanraw_aw","qqmean_aw","qqmax_aw","varratio_aw")]
#write.table(nnresults, file = file.path("case_study/cs_results",paste("nn","R",R,".txt",sep="")))

#### save workspace prima di twang

save.image(file=file.path("case_study/cs_results",paste("empsimR",R,".Rdata",sep="")))

###
 twresults <- matrix(0,length(statistics),1)
  for(i in 1:length(scenarios)){
	  	                   pb <- txtProgressBar(min = 0, max = R, style = 3) 
  	                          partialresults<-matrix(0,length(statistics),1)
  	                          for(j in 1:R)
                                             {
                partialresults <-partialresults+try(as.matrix(simplify2array(funsim_cs(
                funcov_cs(data,substr(scenarios[i],1,1),substr(scenarios[i],2,2))  ,"gbmtwang"))))
   setTxtProgressBar(pb, j)  
                                              }
 twresults<-cbind(twresults,partialresults/R)
                               }
twresults<-twresults[,-1]


#write.table(twresults, file = file.path("cs_results",paste("twang","R",R,".txt",sep="")))


#### save workspace after twang

save.image(file=file.path("case_study/cs_results",paste("empsim_R",R,".Rdata",sep="")))


