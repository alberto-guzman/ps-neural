
####################################
### Figure for results of empirical simulation (Fig 3)
####################################

# load data
load(file="case_study/cs_results/empsimR500.Rdata")
longt<-t(cbind(lgresults,rfresults,twresults,bagresults,treeresults,nnresults,nbresults))
#add column psmodel 
psmodel<-c("logit","rf","tw","bag","tree","nn","nb")
psmodel<-rep(psmodel,rep(length(scenarios),length(psmodel)))
#add column scenario
scenario<-rownames(longt)
empsimres<-data.frame(scenario,psmodel,longt)
empsimres$psmodel<-factor(empsimres$psmodel,levels=c("logit","rf","tw","bag","tree","nn","nb"))
# save plot in pdf 
pdf("results/Fig3_empsim_results_with_balance.pdf",width=11,height=9)
# graphical window
par(mfrow=c(2,2))
# plot BIAS
# select data
subdata <- subset(empsimres, select=c("psmodel","absrbiasw","absrbiasm"))
head(subdata)
# find statistics 
absrbiasw<-as.numeric(by(subdata$absrbiasw,subdata$psmodel,mean))# mean bias weight
absrbiasm<-as.numeric(by(subdata$absrbiasm,subdata$psmodel,mean))# meand bias match
# barplot
barplot(as.matrix(rbind(absrbiasw,absrbiasm),nrow=2,ncol=nlevels(subdata$psmodel)),beside=TRUE,main="abs relative bias",
ylim=c(0,50),names.arg=levels(subdata$psmodel),ylab="")
abline(h=c(10,20,30,40,50),col="gray",lty="dotted")

# gap barplot (to account for tall bars)
# require(plotrix)
# newdata<-(rbind(absrbiasw,absrbiasm))
# #newdata[newdata>60]<-newdata[newdata>60]-102
# barpos<-barplot(newdata,names.arg=levels(subdata$psmodel),
# ylim=c(0,90),beside=TRUE,main="abs relative bias",yaxt="n")
# axis(2,at=c(0,10,20,30,40,50,60,70,80,90,100),
 # #   labels=c(0,10,20,30,40,50,60,160,170,180,190))
# #box()
# #axis.break(2,60,style="gap")
# abline(h=c(10,20,30,40,50,60,70,80),col="gray",lty="dotted")

# plot RMSE
# select data
subdata <- subset(empsimres, select=c("psmodel","varhatgw","varhatgm"))
head(subdata)
# find statistics 
rmsew<-sqrt(as.numeric(by(subdata$varhatgw,subdata$psmodel,mean)))# mean rmse weight
rmsem<-sqrt(as.numeric(by(subdata$varhatgm,subdata$psmodel,mean)))# meand rmse match
# barplot
# gap barplot
newdata<-(rbind(rmsew,rmsem))
newdata[newdata>0.03]<-newdata[newdata>0.03]-0.055
barpos<-barplot(newdata,names.arg=levels(subdata$psmodel),
ylim=c(0,max(rmsem)+0.025),beside=TRUE,main="root mean square error",yaxt="n")
axis(2,at=c(0,0.01,0.02,0.03,0.04),
    labels=c(0,0.01,0.02,0.03,0.09))
box()
axis.break(2,0.03,style="gap")
abline(h=c(0,0.01,0.02,0.03,0.04),col="gray",lty="dotted")

#plot Coverage
subdata <- subset(empsimres, 
select=c("psmodel","covw","covm"))
head(subdata)
# find statistics 
covw<-as.numeric(by(subdata$covw,subdata$psmodel,mean))# mean bias weight
covm<-as.numeric(by(subdata$covm,subdata$psmodel,mean))# meand bias match
barplot(as.matrix(rbind(covw,covm),nrow=2,ncol=nlevels(subdata$psmodel)),beside=TRUE,main="Coverage of 95% c.i.",
ylim=c(0,1),names.arg=levels(subdata$psmodel),yaxt="n")
ticks<-c(0,0.2,0.4,0.6,0.8,1)
axis(2,at=ticks,labels=c(0,0.2,0.4,0.6,0.8,1))
abline(h=c(0.2,0.4,0.6,0.8,1),col="gray",lty="dotted")
dev.off()

# additional plot for ASAM
# plot ASAM version 2
# # subdata <- subset(empsimres, 
# select=c("psmodel","masb_aw","masb_am"))
# head(subdata)
# # find statistics 
# masb_aw<-as.numeric(by(subdata$masb_aw,subdata$psmodel,mean))# mean bias weight
# masb_am<-as.numeric(by(subdata$masb_am,subdata$psmodel,mean))# meand bias match
# newdata<-(rbind(masb_aw,masb_am))
# newdata[newdata>10]<-newdata[newdata>10]-240

# barplot(newdata,beside=TRUE,main="ASAM after psm/psw",
# ylim=c(0,40),names.arg=levels(subdata$psmodel),yaxt="n")
# ticks<-c(0,5,10,15,20,30)
# axis(2,at=ticks,labels=c(0,5,10,15,20,70))
# box()
# axis.break(2,20,style="gap")
# abline(h=c(0,5,10,15,20,30),col="gray",lty="dotted")
# #asam before
# abline(h=mean(empsimres$masb_b),col="red",lty="dashed",)
# text(3,mean(empsimres$masb_b)+1,labels="asam before",cex=0.5)

#plot legend
#plot.new()
# par(xpd=TRUE)
# legend(x = "topleft",inset = 0,
        # legend = c("Weigthing","Matching"),col=c("black","grey"), lwd=4, cex=0.6, horiz = TRUE)
# par(xpd=FALSE)




