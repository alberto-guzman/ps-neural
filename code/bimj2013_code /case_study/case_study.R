
###############################
##########  Case study  ##########

##############################
# Using ps matching and weigthing to study the effect of
# labor induction versus expectant management
# in long-term, nulliparous, phisiological pregnancies
##############################

# load packages and functions
source("../functions.r")

# load data
 load("cs_data/clean_data_10cov.Rdata") 
# sample selection: alive, head, single newborns at 41 or later gestational weeks; no previous cs or pathological pregnancies
# variable selection:
# xvar<-c(
# "momage",
# #"Presentazione",
# "nullip",
# "Lunghezza",
# "Peso",
# "Circonferenza",
# #"Numero.parti.precedenti",
# #"Numero.tagli.cesarei",
# "Numero.ecografie.gravidanza",
# "Numero.visite.controllo.gravidanza",#, # solo in dati sardegna
# "Abitudine.fumo.madre",
# #"Metodica.antidolore",
# "Ricoveri",
# #"Amniocentesi"
# #"age3034",
# #"age35older",
# "Sesso"
# #"m_educ_2",
# #"m_educ_3"
# )

#  Balance before 
bbefore<-MatchBalance(T~.-Y,data=data)

# Balance after optimal tuning 
log    <- funsim_tp(data,"logit")
rf       <- funsim_tp(data,"randomforest", tp=2)
bag   <- funsim_tp(data,"tree", tp=0.0015)
tw      <- funsim_tp(data,"gbmtwang")
tree    <- funsim_tp(data,"tree", tp=0.0015)
nb      <- funsim_tp(data,"nb", tp=c(FALSE))
nn      <- funsim_tp(data,"nnet", tp=c(TRUE,15,2))

asam<-rbind(log$masb_am,log$masb_aw,rf$masb_am,rf$masb_aw,bag$masb_am,bag$masb_aw,tree$masb_am,tree$masb_aw,nb$masb_am,nb$masb_aw,nn$masb_am,nn$masb_aw)
asam_inter<-rbind(log$masbinter_am,log$masbinter_aw,rf$masbinter_am,rf$masbinter_aw,bag$masbinter_am,bag$masbinter_aw,tree$masbinter_am,tree$masbinter_aw,nb$masbinter_am,nb$masbinter_aw,nn$masbinter_am,nn$masbinter_aw)
asam20<-rbind(log$over20masb_am,log$over20masb_aw,rf$over20masb_am,rf$over20masb_aw,bag$over20masb_am,bag$over20masb_aw,tree$over20masb_am,tree$over20masb_aw,nb$over20masb_am,nb$over20masb_aw,nn$over20masb_am,nn$over20masb_aw)
asam10<-rbind(log$over10masb_am,log$over10masb_aw,rf$over10masb_am,rf$over10masb_aw,bag$over10masb_am,bag$over10masb_aw,tree$over10masb_am,tree$over10masb_aw,nb$over10masb_am,nb$over10masb_aw,nn$over10masb_am,nn$over10masb_aw)
hatg<-rbind(log$hatgw,log$hatgm,rf$hatgw,rf$hatgm,bag$hatgw,bag$hatgm,tree$hatgw,tree$hatgm,nb$hatgw,nb$hatgm,nn$hatgw,nn$hatgm)
sdhatg<-rbind(log$sdhatgw,log$sdhatgm,rf$sdhatgw,rf$sdhatgm,bag$sdhatgw,bag$sdhatgm,tree$sdhatgw,tree$sdhatgm,nb$sdhatgw,nb$sdhatgm,nn$sdhatgw,nn$sdhatgm)

baltab<-cbind(asam,asam_inter,asam20,asam10)

rownames(baltab)<-c("log_m","log_w","rf_m","rf_w","bag_m","bag_w","tree_m","tree_w","nb_m","nb_w","nn_m","nn_w")
colnames(baltab)<-c("asam","asam_inter","asam20","asam10");baltab

# outcome analysis on raw data
t.test(data$Y[data$T==1],data$Y[data$T==0],alternative="two.sided")
#sd

# outcome analysis for best balancing ps model (logit)

# final table
tab<-cbind(baltab,hatg,sdhatg)
colnames(tab)<-c(colnames(baltab),"hatg","sd.hatg");baltab

ci.left<-hatg-2*sdhatg
ci.right<-hatg+2*sdhatg





