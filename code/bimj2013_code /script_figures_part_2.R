
###### Code for figures  ######

# load simulation results
 load(file="results/simresults.Rdata")

### check
dim(simres)
table(simres$size);table(simres$scenario);table(simres$psmodel)

#### modify label order or labels for better plotting
simres$psmodel<-factor(simres$psmodel,levels=c("logit","rf","tw","bag","tree","nn","nb"))


############################
#### Fig 4: barplot bias of ATE estimators by scenario
#############################

# setEPS()
# postscript("results/Fig4_bias.eps",width=13,height=9)
pdf("results/Fig4_bias.pdf",width=13,height=9)
#tiff("results/Fig4_bias.tiff",width=12,height=8)
#tiff("Fig4.tiff",width=2000,height=2000,units="px",res=800)
# set graphical window
par(mfrow=c(4,4),mai= c(0.3, 0.6, 0.5, 0.3)) #4 times 4 window (16 barplots)
par(oma=c(0,4,4,0))#makes room for overall X and Y axis labels
# loop on scenarios
for (i in 1:nlevels(simres$scenario)){
# select data for plot i
subdata <- subset(simres, 
scenario==levels(scenario)[i],
select=c("psmodel","absrbiasw","absrbiasm"))
head(subdata)
# find statistics 
absrbiasw<-as.numeric(by(subdata$absrbiasw,subdata$psmodel,mean))# mean bias weight
absrbiasm<-as.numeric(by(subdata$absrbiasm,subdata$psmodel,mean))# meand bias match
# plot
barplot(as.matrix(rbind(absrbiasw,absrbiasm),nrow=2,ncol=nlevels(subdata$psmodel)),beside=TRUE,
ylim=c(0,50),names.arg=levels(subdata$psmodel),
# T scenario labels in first column
ylab=if (isTRUE(i==1|i==5|i==9|i==13)){substr(levels(simres$scenario)[i],1,1)} else{NULL},font.lab=2,
# forcedlabelT<-c("A","D","E","G")
# ylab=if (isTRUE(i==1|i==5|i==9|i==13)){forcedlabelT[i]} else{NULL},
# Y scenario labels in first row
main=if (i<5) {substr(levels(simres$scenario)[i],2,2)},legend=NULL)
#main=if (i<5) {newlabelY[i]} else{NULL},legend=NULL);
#grid(nx = 0, ny = 4, col = "gray", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)
abline(h=c(5,10,20,30,40,50),col="gray",lty="dotted")
mtext('outcome equation scenarios', side = 3, outer = TRUE, line = 2)
mtext('treatment equation scenarios', side = 2, outer = TRUE, line = 2)
}# end cycle
dev.off()

######################
# Fig 5 barplot RMSE of ATE estimators by scenario
######################

# setEPS()
# postscript("results/Fig5_rmse.eps",width=13,height=9)
pdf("results/Fig5_rmse.pdf",width=13,height=9)
# set graphical parameters
par(mfrow=c(4,4),mai= c(0.3, 0.6, 0.5, 0.3))
par(oma=c(0,4,4,0))
# select square distance of estimand for each ps model
for (i in 1:nlevels(simres$scenario)){
subdata <- subset(simres, 
scenario==levels(scenario)[i],
select=c("psmodel","size","varhatgw","varhatgm"))
#head(subdata)
# statistics (variance after weigthing and matching)
varhatgw<-as.numeric(by(subdata$varhatgw,subdata$psmodel,mean)) #by psmodel
varhatgm<-as.numeric(by(subdata$varhatgm,subdata$psmodel,mean))
# plot statistics
barplot(matrix(rbind(sqrt(varhatgw),sqrt(varhatgm)),nrow=2,ncol=nlevels(subdata$psmodel)),
beside=TRUE,ylim=c(0,0.3),names.arg=levels(subdata$psmodel),
ylab=if (isTRUE(i==1|i==5|i==9|i==13)){substr(levels(simres$scenario)[i],1,1)} else{NULL},font.lab=2,
main=if (i<5){substr(levels(simres$scenario)[i],2,2)} else{NULL},legend=NULL);
#main=levels(simres$scenario)[i],legend=NULL)
#;grid(nx = 0, ny = 4, col = "gray", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)
abline(h=c(0.05,0.1,0.15,0.2,0.25,0.3),col="gray",lty="dotted")
mtext('outcome equation scenarios', side = 3, outer = TRUE, line = 2)
mtext('treatment equation scenarios', side = 2, outer = TRUE, line = 2)
}# end cycle
#legenda
# #plot.new()#cambia plot
# par(xpd=TRUE)
# legend(x = "topleft",inset = 0,
        # legend = c("Weigthing","Matching"),col=c("black","grey"), lwd=7, cex=.7, horiz = TRUE)
# par(xpd=FALSE)
dev.off()

######################
# Fig 6 barplot coverage probability of 95% CI by scenario
###################### 

# setEPS()
# postscript("results/Fig6_coverage.eps",width=13,height=9)
pdf("results/Fig6_coverage.pdf",width=13,height=9)
# set window
# set graphical parameters
par(mfrow=c(4,4),mai= c(0.3, 0.6, 0.5, 0.3))
par(oma=c(0,4,4,0))

# select square distance of estimand for each ps model
for (i in 1:nlevels(simres$scenario)){
subdata <- subset(simres, 
scenario==levels(scenario)[i],
select=c("psmodel","size","covw","covm"))
#head(subdata)
# statistics (coverage of weigthing and matching CI)
covw<-as.numeric(by(subdata$covw,subdata$psmodel,mean)) #by psmodel
covm<-as.numeric(by(subdata$covm,subdata$psmodel,mean))
# plot statistics
barplot(matrix(rbind(covw,covm),nrow=2,ncol=nlevels(subdata$psmodel)),
beside=TRUE,ylim=c(0,1),names.arg=levels(subdata$psmodel),
ylab=if (isTRUE(i==1|i==5|i==9|i==13)){substr(levels(simres$scenario)[i],1,1)} else{NULL},font.lab=2,yaxt="n",
main=if (i<5){substr(levels(simres$scenario)[i],2,2)} else{NULL},legend=NULL);
ticks<-c(0,0.2,0.4,0.6,0.8,1)
axis(2,at=ticks,labels=c(0,0.2,0.4,0.6,0.8,1))
abline(h=c(0.2,0.4,0.6,0.8,1),col="gray",lty="dotted")
mtext('treatment equation scenarios', side = 3, outer = TRUE, line = 2)
mtext('outcome equation scenarios', side = 2, outer = TRUE, line = 2)
}# end cycle
#legenda
# par(xpd=TRUE)
# legend(x = "topleft",inset = 0,
        # legend = c("Weigthing","Matching"),col=c("black","grey"), lwd=7, cex=.7, horiz = TRUE)
# par(xpd=FALSE)
dev.off()

#### Graphs on correlations

# e.g.: corr(biasafterweigthing, asamafterweighting) for logit with size=500:

# subdata <- subset(data, size==500 & psmodel=="logit", select=c(absrbiasw, masb_aw))
# cor(subdata[,1],subdata[,2])
# Kendall(subdata[,1],subdata[,2])
# rcorr(as.matrix(subdata),type="spearman")

# eg: all correlations (bias and balancemetric[i]):

#select bias and balance metrics

awvars<-#aw: after weighting
c("size","scenario","psmodel","absrbiasw","auc","masb_aw","masbinter_aw","over20masb_aw","over10masb_aw","qqmean_aw","qqmax_aw","varratio_aw")
amvars<-#am: after matching
c("size","scenario","psmodel","absrbiasm","auc","masb_am","masbinter_am","over20masb_am","over10masb_am","qqmean_am","qqmax_am","varratio_am")


# select scenario: 

# only scenarios with "a" in Y (e.g: Stuart)
# dw <- subset(simres , scenario=="Aa" | scenario=="Ca" | scenario=="Fa" | scenario=="Ga", select=awvars)
# dm <- subset(simres , scenario=="Aa" | scenario=="Ca" | scenario=="Fa" | scenario=="Ga", select=amvars)

# all scenarios:
dw <- subset(simres , select=awvars)
dm <- subset(simres , select=amvars)
#dm<-na.omit(dm)
#dw<-na.omit(dw)

# colnames(dw)
# rename labels for better plotting 
namesvars<-
c("size","scenario","psmodel","abs rel bias",     "auc",          "asam"  , "asamint"  ,     "asam20" , "asam10" , "ecdfmean"   ,  "ecdfmax"    ,  "varratio"  )
colnames(dw)<-colnames(dm)<-namesvars
datawm<-rbind(dw,dm)
#datawm<-datawm[-which(datawm$varratio==Inf)];datawm<-na.omit(datawm)

# calculate correlation[ bias , metrics] for 1)weight 2) match 3) all data
# initialize correlation vectors
vector()->wcor->wtau->wspe->mcor->mtau->mspe->cor->tau->spe  
for(i in 5:length(datawm)){
wcor[i] <-cor(dw[,"abs rel bias"],dw[,i], use="complete.obs") 
wtau[i] <-cor(dw[,"abs rel bias"],dw[,i],method="kendall", use="complete.obs")
wspe[i] <-cor(dw[,"abs rel bias"],dw[,i],method="spearman", use="complete.obs")
mcor[i] <-cor(dm[,"abs rel bias"],dm[,i], use="complete.obs")
mtau[i] <-cor(dm[,"abs rel bias"],dm[,i],method="kendall", use="complete.obs")
mspe[i] <-cor(dm[,"abs rel bias"],dm[,i],method="spearman", use="complete.obs")
cor[i] <-cor(datawm[,"abs rel bias"],datawm[,i], use="complete.obs") 
tau[i] <-cor(datawm[,"abs rel bias"],datawm[,i],method="kendall", use="complete.obs")
spe[i] <-cor(datawm[,"abs rel bias"],datawm[,i],method="spearman", use="complete.obs")
}

#####   graphing correlations with scatterplots ########

#############################
# Fig 7 version. 1:   cor(bias , balance) for weighting and matching in the same plot
#############################

 # setEPS()
 # postscript("results/Fig7v1.eps",width=11,height=11)
 pdf("results/Fig7v1.pdf",width=11,height=11)
# opens the graphic window
par(mfrow=c(3,3),mai= c(0.8, 0.8, 0.3, 0.3)) 
# cycle for plotting cor (rel bias, i^th balancemetric)
for(i in 5:(length(datawm))){
plot(datawm[,i],datawm[,"abs rel bias"],
 main=paste("cor =",as.character(round(cor[i],digits=3)),
            "  tau =",as.character(round(tau[i],digits=3)),
            "  spe =",as.character(round(spe[i],digits=3)) ),
           #sub ="subtitle" ,
   xlab=namesvars[i], ylab=if (isTRUE(i==5|i==8|i==11)){"bias"}else{""}, pch=20,cex.lab=1.3,cex.main=0.8,font.main=1,col.sub="grey")
# add regression line
# abline(lm(datawm[,"abs rel bias"]~datawm[,i]),col="black")
# add lowess line
lines(lowess(datawm[,i],datawm[,"abs rel bias"]), col="grey")
} 
dev.off()
#############################
# Fig 7 version 2:   cor(bias , balance) for weighting and matching (same plot with separate colors)
#############################
 # setEPS()
 # postscript("results/Fig7v2.eps",width=11,height=11)
  pdf("results/Fig7v2.pdf",width=11,height=11)
# opens the graphic window 
par(mfrow=c(3,3),mai= c(0.8, 0.8, 0.4, 0.2))
# cycle for plotting:
for(i in 5:(length(datawm)-1)){
plot(dw[,i],dw[,"abs rel bias"],#plot (metric,bias aw)
 main=paste("Weighting: cor =",as.character(round(wcor[i],digits=3)),
            "   tau =",as.character(round(wtau[i],digits=3)),
            "   spe =",as.character(round(wspe[i],digits=3)),"\n","\n",
            "Matching: cor =",as.character(format(mcor[i],digits=3)),
            "   tau =",as.character(round(mtau[i],digits=3)),
            "   spe =",as.character(round(mspe[i],digits=3))),
           #sub ="subtitle" ,
   xlab=namesvars[i], ylab=if (isTRUE(i==5|i==8|i==11)){"bias"}else{""}, pch=20,cex.lab=1.3,cex.main=0.8,font.main=1,col.sub="grey")
#indna<-which(is.na(dw[,i]))
#lines(lowess(dw[,i][-indna],dw[,"abs rel bias"][-indna]), col="black")
lines(lowess(dw[,i],dw[,"abs rel bias"]), col="black") # lowess line for weighting
#indna<-which(is.na(dm[,i]))
#lines(lowess(dm[,i][-indna],dm[,"abs rel bias"][-indna]), col="gray")
points(dm[,i], dm[,"abs rel bias"], col='gray',pch=20) # plot (metric,bias am)
lines(lowess(dm[,i],dm[,"abs rel bias"]), col="grey")# lowess line for matching
} 
#legend
# plot.new()
# par(xpd=TRUE)
# legend(x = "top",inset = 0,
# legend = c("Weigthing","Matching"),col=c("black","grey"), lwd=11, cex=.7, horiz = TRUE)
# par(xpd=FALSE)
dev.off()

##########################
## Supplementary Figure 0: barplot ASAM by scenario
##########################

# # #setEPS()
# #postscript("Figbalance.eps",width=13,height=8)
#pdf("results/Figbalance.pdf",width=13,height=9)
par(mfrow=c(4,4),mai= c(0.3, 0.6, 0.5, 0.3))
#par(mfrow=c(4,4))
par(oma=c(0,4,5,0))#make room for overall X and Y axis labels
#i=7 # choose one scenario only
for (i in 1:nlevels(simres$scenario)){

# data for plot
subdata <- subset(simres, 
scenario==levels(scenario)[i] #& size==500
,select=c("psmodel","masb_aw","masb_am","masb_b"))
head(subdata)

# statistics (asam after weightin, asam after matching, asam before) by ps model
avg.masb_aw<-as.numeric(by(subdata$masb_aw,subdata$psmodel,mean))
avg.masb_am<-as.numeric(by(subdata$masb_am,subdata$psmodel,mean))
avg.masb_b<-as.numeric(by(subdata$masb_b,subdata$psmodel,mean))

#plot
barplot(as.matrix(rbind(avg.masb_aw,avg.masb_am),nrow=2,ncol=nlevels(subdata$psmodel)),beside=TRUE,
ylim=c(0,30),names.arg=levels(subdata$psmodel),
# T scenario labels in first column
ylab=if (isTRUE(i==1|i==5|i==9|i==13)){substr(levels(simres$scenario)[i],1,1)} else{NULL},font.lab=2,
# forcedlabelT<-c("A","D","E","G")
# ylab=if (isTRUE(i==1|i==5|i==9|i==13)){forcedlabelT[i]} else{NULL},
# Y scenario labels in first row
main=if (i<5) {substr(levels(simres$scenario)[i],2,2)},legend=NULL)
#grid(nx = 0, ny = 4, col = "gray", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)
abline(h=c(0,5,10,15,20,25,30),col="gray",lty="dotted")
abline(h=mean(avg.masb_b),col="red",lty="dashed")#asam before
if (i==1){text(3,mean(avg.masb_b+1),labels="asam before",cex=0.5)}
mtext(' treatment equation scenarios', side = 2, outer = TRUE, line = 2)
mtext('outcome equation scenarios', side = 3, outer = TRUE, line = 2)
}# end cycle
#dev.off()
#legenda
#plot.new()#cambia plot
par(xpd=TRUE)
legend(x = "topleft",inset = 0,
        legend = c("Weigthing","Matching"),col=c("black","grey"), lwd=7, cex=.7, horiz = TRUE)
par(xpd=FALSE)
dev.off()







