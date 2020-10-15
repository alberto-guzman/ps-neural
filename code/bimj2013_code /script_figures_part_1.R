#############
## This file merge all simulation results of all sizes in a single data.frame
############## 

# size 500
# load workspace with raw results for size 500
load("results/R500 size500 parATE.Rdata")

# merge all size 500 results
colnames(lgresults)<-colnames(rfresults)<-colnames(twresults)<-colnames(bagresults)<-colnames(nnresults)<-colnames(nbresults)<-scenarios
longt<-t(cbind(lgresults,rfresults,twresults,bagresults,treeresults,nnresults,nbresults))
# column psmodel 
psmodel<-c("logit","rf","tw","bag","tree","nn","nb")
psmodel<-rep(psmodel,rep(length(scenarios),length(psmodel)))
# column scenario
scenario<-rownames(longt)
# column size
size<-factor(rep(size,nrow(longt)))
# add columns
table500<-data.frame(size,scenario,psmodel,longt)

rm(list=setdiff(ls(), c("table500"))) 
# simres<-table500# decomment to analyze only size 500 results

# repeat for all sizes:

# size 1000
# load workspace with results for size 1000
load("results/R500 size1000 parATE.Rdata")

# merge all size 1000 results
colnames(lgresults)<-colnames(rfresults)<-colnames(twresults)<-colnames(bagresults)<-colnames(nnresults)<-colnames(nbresults)<-scenarios
longt<-t(cbind(lgresults,rfresults,twresults,bagresults,treeresults,nnresults,nbresults))
# column psmodel 
psmodel<-c("logit","rf","tw","bag","tree","nn","nb")
psmodel<-rep(psmodel,rep(length(scenarios),length(psmodel)))
# column scenario
scenario<-rownames(longt)
# column size
size<-factor(rep(size,nrow(longt)))
# add columns
table1000<-data.frame(size,scenario,psmodel,longt)
rm(list=setdiff(ls(), c("table500","table1000"))) 

# size 2000
# load workspace with raw results for size 2000
load("results/R500 size1000 parATE.Rdata")
# merge all size 2000 results as above:
colnames(lgresults)<-colnames(rfresults)<-colnames(twresults)<-colnames(bagresults)<-colnames(nnresults)<-colnames(nbresults)<-scenarios
longt<-t(cbind(lgresults,rfresults,twresults,bagresults,treeresults,nnresults,nbresults))
# column psmodel 
psmodel<-c("logit","rf","tw","bag","tree","nn","nb")
psmodel<-rep(psmodel,rep(length(scenarios),length(psmodel)))
# column scenario
scenario<-rownames(longt)
# column size
size<-factor(rep(size,nrow(longt)))
# add columns
table2000<-data.frame(size,scenario,psmodel,longt)
rm(list=setdiff(ls(), c("table500","table1000","table2000"))) 

#### Merge all simulation results in one single data.frame (all sizes together...)
simres<-rbind(table500,table1000,table2000)

# save table with all simulation results
rm(list=setdiff(ls(), c("simres")))
save.image(file="results/simresults.Rdata")


