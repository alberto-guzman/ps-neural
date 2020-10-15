####################################
# This scripts is for cross validating tuning parameters in ML algorithms for simulated data
# optimal tuning values are sought with covariate balance in sight

# load functions
source("functions.R")

# fix size of simulated data set


# generate dataset
set.seed(3)
scenario<-"Aa"
data<-funcov(1000,"A","a")



# add label for k-fold cross-validation
k=10
data$fold = cut(sample(1:nrow(data),nrow(data)), breaks=k, labels=F)
# oppure usa caret package
# data$fold <- createFolds(data$Y, k = 10,list=FALSE)

##################
## Selezione tuning parameter value con caret package
##################
require(caret)
metric <- "Accuracy"

# use caret with random search for random forest
# Random Search
control <- trainControl(method="repeatedcv", number=10, repeats=1, search="random")
set.seed(4)
#mtry <- sqrt(ncol(data))
rf_random <- train(factor(T)~w1+w2+w3+w4+w5+w6+w7+w8+w9+w10,
data=data, method="rf", metric=metric, tuneLength=15, trControl=control)
print(rf_random)
#plot(rf_random)
# find balance for random tuning values
mtry<-rf_random$results[,"mtry"]
rfacc<-rf_random$results[,"Accuracy"]
rfkap<-rf_random$results[,"Kappa"]
len<-length(mtry)
asam_am_rf<-asam_aw_rf<-vector(length=len)
for (j in 1:len){
masb_aw<-masb_am<-list()
for (i in 1:length(unique(data$fold))) {
 masb_aw[[i]] <- funsim_tp(data[data$fold!=i,],"randomforest", tp=mtry[j] )$masb_aw
 masb_am[[i]] <- funsim_tp(data[data$fold!=i,],"randomforest", tp=mtry[j] )$masb_am	
 }
 asam_am_rf[ j ]<-mean(unlist(masb_am))
 asam_aw_rf[ j ]<-mean(unlist(masb_aw))
 }
 # find optimal values of mtry for accuracy and balance
 mtry[which.max(rfacc)]#1
 mtry[which.min(asam_am_rf)]#3
 mtry[which.min(asam_aw_rf)]#3
 
# Using caret for random search for rpart 
# Random Search
control <- trainControl(method="repeatedcv", number=100, repeats=1, search="random")
set.seed(7)
rpart_random <- train(factor(T)~w1+w2+w3+w4+w5+w6+w7+w8+w9+w10,
data=data, method="rpart", metric=metric, tuneLength=7, trControl=control)
print(rpart_random)
# plot(rpart_random)
# find balance for random tuning values
cp<-c(rpart_random$results[,"cp"])#;cp<-cp[cp<0.03]
rpartacc<-rpart_random$results[,"Accuracy"]
rpartkap<-rpart_random$results[,"Kappa"]
 cp2<-c(cp,0.005,0.01)
 tunegrid <- treeGrid <-  expand.grid(cp=cp2)
  rpart_grid <- train(factor(T)~w1+w2+w3+w4+w5+w6+w7+w8+w9+w10,data=data, method="rpart", metric=metric, tuneGrid=tunegrid, trControl=control)
 rpartacc<-rpart_grid$results[,"Accuracy"]
len<-length(cp2)
asam_am_tree<-asam_aw_tree<-vector(length=len)
for (j in 1:len){
masb_aw<-masb_am<-list()
for (i in 1:length(unique(data$fold))) {
 masb_aw[[i]] <- funsim_tp(data[data$fold!=i,],"tree", tp=cp2[j] )$masb_aw
 masb_am[[i]] <- funsim_tp(data[data$fold!=i,],"tree", tp=cp2[j] )$masb_am	
 }
 asam_am_tree[ j ]<-mean(unlist(masb_am))
 asam_aw_tree[ j ]<-mean(unlist(masb_aw))
 }
 # find optimal value of cp for accuracy and balance
 cp2[which.max(rpartacc)]
 cp2[which.min(asam_am_tree)]
 cp2[which.min(asam_aw_tree)]
 
 # use caret for tuning net
 control <- trainControl(method="repeatedcv", number=100, repeats=1, search="random")
nnet_random <- train(factor(T)~w1+w2+w3+w4+w5+w6+w7+w8+w9+w10,
data=data, method="nnet", metric=metric, tuneLength=15, trControl=control)
#print(nnet_random)
#plot(nnet_random)
# find random tuning values
size<-nnet_random$results[,"size"]
decay<-nnet_random$results[,"decay"]
nnetacc<-nnet_random$results[,"Accuracy"]
#nnetkap<-nnet_random$results[,"Kappa"]

# # find balance for random tuning values
# len<-length(size)
# asam_am_net<-asam_aw_net<-vector(length=len)
# for (j in 1:len){ 
# masb_aw<-masb_am<-list()
# for (i in 1:10) {
 # masb_aw[[ i ]] <- funsim_tp(data[data$fold!=i,],"nnet", tp=c(TRUE,size[j],decay[j]))$masb_aw
 # masb_am[[ i ]] <- funsim_tp(data[data$fold!=i,],"nnet", tp=c(TRUE,size[j],decay[j]))$masb_am
                      # }
 # asam_aw_net[ j  ]<-mean(unlist(masb_aw))
 # asam_am_net[ j ]<-mean(unlist(masb_am))
 # }
# cbind(nnet_random$results,asam_am_net,asam_aw_net)
 # # find optimal value of size and decay for balance
 # size[which.max(nnetacc)]; decay[which.max(nnetacc)]
 # size[which.min(asam_am_net)]; decay[which.min(asam_am_net)]
 # size[which.min(asam_aw_net)]; decay[which.min(asam_aw_net)]

 # use caret with search on prespecified tuning values
 # grid values
 # #### tune nnet with tuning grid on number of hidden units and decay
 control <- trainControl(method="repeatedcv", number=10, repeats=1, search="grid")
 tunegrid <- nnetGrid <-  expand.grid(
                       #  entropy = c(TRUE), 
                        size = seq(5,19,2), 
                        decay =c(decay,0))
 nnet_grid <- train(factor(T)~w1+w2+w3+w4+w5+w6+w7+w8+w9+w10,data=data, method="nnet", metric=metric, tuneGrid=tunegrid, trControl=control)
#print(nnet_grid)
#plot(nnet_grid)
 # find balance for grid tuning values
 size<-nnet_grid$results[,"size"]
 decay<-nnet_grid$results[,"decay"]
 nnetacc<-nnet_grid$results[,"Accuracy"]
 size_opt<-size[which.max(nnetacc)]
 decay_opt<-decay[which.max(nnetacc)]
len<-length(size)
asam_am_net<-asam_aw_net<-vector(length=len)
for (j in 1:len){ 
masb_aw<-masb_am<-list()
for (i in 1:length(unique(data$fold))){ 
 masb_aw[[i]] <- tryCatch(funsim_tp(data[data$fold!=i,],"nnet", tp=c(TRUE,size[j],decay[j]))$masb_aw,error=function(e)paste("error at ",j,sep=""))
 masb_am[[i]] <- tryCatch(funsim_tp(data[data$fold!=i,],"nnet", tp=c(TRUE,size[j],decay[j]))$masb_am,error=function(e)paste("error at ",j,sep=""))	
}
 asam_am_net[ j ]<-mean(unlist(masb_am),na.rm=TRUE)
 asam_aw_net[ j ]<-mean(unlist(masb_aw),na.rm=TRUE)
 }

  # # find optimal value of size and decay for balance
  size[which.max(nnetacc)]; decay[which.max(nnetacc)]#13;1.3
  size[which.min(asam_am_net)]; decay[which.min(asam_am_net)]#11, 1.3
  size[which.min(asam_aw_net)]; decay[which.min(asam_aw_net)]#11, 1.3

 ##############################
 # plot together accuracy and balance
 ##############################

 #load("case_study/cs_results/tuningresults.Rdata")
 pdf(paste("results/Fig2_Tuning_",scenario,"v2.pdf",sep=""),width=13,height=9)
 # define window 
 par(mfrow=c(2,3))
 # plot rf
 par(mar = c(5,5,2,4))
 plot(mtry, rfacc, 
type = "l", lwd = 2, col ="black" , ylab = "Accuracy", ,lty="dashed",
 xlab = paste0("mtry: #randomly selected predictors"), 
main = paste0("random forest"), 
#ylim = c(min(rfacc)-0.01, max(rfacc)+0.01))
ylim=c(0.6,0.75))
points(mtry,rfacc,col="black",pch=17)
#add balance
par(new=T)
plot(mtry,asam_am_rf,col=gray(0.8),pch=19,axes=F,xlab=NA,ylab=NA,
#ylim=c(min(c(asam_am,asam_aw)),max(c(asam_am,asam_aw)))
ylim=c(0,15)
)
axis(side=4)
mtext(side=4,line=3,"ASAM")
#points(mtry,asam_am,col=gray(0.8),pch=17)
lines(mtry, asam_am_rf, lwd = 2, col = gray(0.8))#col = "darkgreen"
points(mtry,asam_aw_rf,col=gray(0.4),pch=1)
lines(mtry, asam_aw_rf, lwd = 2, col = gray(0.4))
#abline(v=0.01,col="red")#default value of cp
#text(0.01+0.01,asamb-2,labels="cp default=0.01",col="red",cex=0.8)
# abline(v=3,col="black",lty="dotted")#default value of mtry
segments(x0=3,y0=-01,x1=3,y1=15,lty="dotted")
#text(3.7,6,labels="mtry default = 3",col="black",cex=0.8)

# plot rpart
# plot together accuracy and balance
 par(mar = c(5,5,2,4))
  ord<-order(cp2)
  cp2<-cp2[ord]
  rpartacc<-rpart[ord]
  asam_am_tree<-asam_am_tree[ord]
  asam_aw_tree<-asam_aw_tree[ord]
 plot(cp2, rpartacc, 
type = "l", lwd = 1, col ="black" , ylab = "Accuracy",lty="dashed",
 xlab = paste0("cp: complexity parameter"), 
main = paste0("classification tree"), 
xlim=c(0, 0.011),
ylim=c(0.6,0.75)
#ylim = c(min(rpartacc)-0.01, max(rpartacc)+0.01)
)
points(cp2,rpartacc,col="black",pch=17)
#add balance
par(new=T)
plot(cp2,asam_am_tree,col=gray(0.8),pch=19,axes=F,xlab=NA,ylab=NA,
xlim=c(0, 0.011),
#yim<-c(6,17)
ylim=c(min(c(asam_am_tree,asam_aw_tree)),max(c(asam_am_tree,asam_aw_tree)))
)
axis(side=4)
mtext(side=4,line=3,"ASAM")
#abline(v=0.01,col="black",lty="dotted")#default value of cp
segments(x0=0.01,y0=012,x1=0.01,y1=19,lty="dotted")
#points(mtry,asam_am,col=gray(0.8),pch=17)
lines(cp2, asam_am_tree, lwd = 1, col = gray(0.8))#col = "darkgreen"
points(cp2,asam_aw_tree,col=gray(0.4),pch=1)
lines(cp2, asam_aw_tree, lwd = 1, col = gray(0.4))
#text(0.01+0.01,asamb-2,labels="cp default=0.01",col="red",cex=0.8)
cpdef <- 0.01
#text(3.7,6,labels="mtry default = 3",col="black",cex=0.8)
#add legend
plot.new()
legend(x = "bottomleft",inset=0.1, legend = c("Accuracy","ASAM after matching","ASAM after weighting"), 
 lty = c(2,1, 1), pch=c(17,19,1),lwd = rep(1, 3), col = c("black", gray(0.9),gray(0.4)), 
 text.width = 1.7, cex = 1,bty="n")

 # plot nnet
 # plot accuracy, asam_am and asam_aw for net
 nnet_tab<-cbind(nnet_grid$results,asam_am_net,asam_aw_net)
 lengt<-length(unique(nnet_tab$decay))
color=c(rep(gray(0.8),lengt-1),"black")
for (i in unique(nnet_tab$decay)){
	    dataplot<-nnet_tab[nnet_tab$decay==i,]
 plot(dataplot$size, dataplot$Accuracy, 
type = "l", lwd = 1, col=color[which(unique(nnet_tab$decay)==i)]#col =gray(1/(i+1))
,ylab = "Accuracy",lty="solid",
 xlab = paste0("size: # of hidden units"), 
#main = paste0("neural network"),
#xlim=c(0,0.05),
ylim = c(min(nnetacc)-0.00001, max(nnetacc)+0.0001)
)
points(dataplot$size, dataplot$Accuracy,col =color[which(unique(nnet_tab$decay)==i)] ,pch=17)
#points(dataplot$size, dataplot$Accuracy,col=gray(1/(i+1)),pch=10 )
par(new=T)
}
#asam_am
par(new=F)
# par(mar = c(4,4,4,3))
for (i in rev(unique(nnet_tab$decay))){
	    dataplot<-nnet_tab[nnet_tab$decay==i,]
 plot(dataplot$size, dataplot$asam_am_net, 
type = "l", lwd = 1, col =color[which(unique(nnet_tab$decay)==i)],
 ylab = "ASAM after matching",lty="solid",
 xlab = paste0("size: # of hidden units"), 
 #xlab = paste0("size: # of hidden units"), 
main = paste0("neural network"),
#xlim=c(0,0.05),
#ylim = c(min(asam_am_net)-0.00001, max(asam_am_net)+0.0001)
ylim=c(2,28)
)
points(dataplot$size, dataplot$asam_am_net,col =color[which(unique(nnet_tab$decay)==i)]
,pch=19)
par(new=T)
}

#asam_aw
par(new=F)
# par(mar = c(4,4,4,3))
for (i in rev(unique(nnet_tab$decay))){
	    dataplot<-nnet_tab[nnet_tab$decay==i,]
 plot(dataplot$size, dataplot$asam_aw_net, 
type = "l", lwd = 1, 
col=color[which(unique(nnet_tab$decay)==i)], 
ylab = "ASAM after weighting",lty="solid",
 xlab = paste0("size: # of hidden units"), 
#main = paste0("neural network"),
#xlim=c(0,0.05),
ylim=c(2,28)
#ylim = c(min(asam_am_net)-0.00001, max(asam_am_net)+0.0001)
)
points(dataplot$size, dataplot$asam_aw_net,
col =color[which(unique(nnet_tab$decay)==i)],pch=1)
par(new=T)
}
par(xpd=TRUE)
leglab<-rep(1,11)/unique(nnet_tab$decay)
legend(x = "bottomright",title="Decay", legend = round(unique(nnet_tab$decay),3), horiz=FALSE,inset=c(-0.32,0),
 lty = rep(1,length(leglab)), pch=unique(nnet_tab$decay)+1,
lwd = rep(1, length(leglab)), 
 #col =color[which(unique(nnet_tab$decay)==i)],
 col=c("black",rep("#CCCCCC",length(leglab)-1)),
 text.width = 1.6, cex = 0.75,bty="n")
dev.off()

