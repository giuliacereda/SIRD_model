########codice per i grafici del paper SIRD semplice regioni#####

# LA lista VEC contenga i 20 vec_best (uno per ogni regione)
# La lista list.out contenga n.boot bootstrap (uno per ogni regione)

mainDir<-dirname(rstudioapi::getActiveDocumentContext()$path) #here put the path to the folder#load("100boot_0078p_1701.Rdata")
setwd(mainDir)
setwd(file.path(mainDir, "WS"))
load("100boot_0078p_3101_vncl.Rdata")

setwd(file.path(mainDir, "Figures"))

library(foreign)
library(parallel)
library(xtable())
library(optimr) 
library(splines)
library(future)

ordine_regioni<-c(c(1:11),c(13:17),12,c( 18:20))
REG[ordine_regioni]

q1<-function(x){return(quantile(x,0.05))}
q2<-function(x){return(quantile(x,0.95))}



schermo=1 #0 if you want the files to be saved in pdf
 

letalita=substr(as.character(p_seq), 3,7)
datafine=paste(substr(as.character(data.exit),9, 10), substr(as.character(data.exit),6, 7), sep="")
 

#########################
#   definisco index buoni (Serviva quando il codice produceva NA ora é solo un controllo)
##########################

index_buoni_list<-list() 
for(g in 1:20){
  
  index_buoni<-1:n.boot
ib<-  index_buoni[which(apply( is.na(par_boot[[g]]), 1, sum)==0)]

  index_buoni_list[[g]]<-ib
}

 
lapply(index_buoni_list, length) #controlla che siano tutti lunghi n.boot 
  
 #########################################################3################################################################3
#         FIT di morti positivi nuovipositivi e tabella a 4 date fisse
################################################################3################################################################3

schermo=0
tabella_I=list()
tabella_I2=list()
Icirco<-list()
Icircoq1<-list()
Icircoq2<-list()
 

 
par(mfrow=c(5,4))
 
REG[12]="Trentino"


#############################################################################
#plot decessi
##############################################################################


quale="Decessi"
schermo=0
if(schermo!=1 ){
  nomefile=paste("fit_", letalita , datafine, quale,".pdf", sep="")
  pdf(nomefile, width=9, height = 13)
}
par(mfrow=c(5,4))
for(g in ordine_regioni){
  print(REG[g])
  Regione<-REG[g]
  S.0<-S0.vec[g]
  print(S.0)
  data <- DATA[DATA$region==Regione,] 
  i0<-data$positivi[data$date==(as.Date(data.ini)-1)] # altrove  nel boot
  data$deceduti<-data$deceduti-data$deceduti[data$date==(as.Date(data.ini)-1)]
  data<- data[which(data$date<=data.exit & data$date>=data.ini),]
  nuovi_p<-data$nuovi_pos
  pos<-data$positivi
  obs.dead<-morti[[g]]
  i0=positivi.ini[[g]]
  print(i0)
  if(g==12){
    nuovi_p<-aggregate(data$nuovi_pos, by =list(data$date), sum)$x
    pos<-aggregate(data$positivi, by =list(data$date), sum)$x
    REG[g]="Trentino Alto Adige"
  }
  
  MAT.SPLINE<-list.out[[g]]
  min_spline<-apply(MAT.SPLINE, 1, min)
  index_buoni<- index_buoni_list[[g]]
  sim<- optim.parallel.ns(VEC[[g]], i0=i0, kp=kp, obs=obs.dead,  ####QUA###
                          last.day= data.exit, pesi=pesi,  per.opm=0)
  
  print(optim.parallel.ns(VEC[[g]], i0=i0, kp=kp, obs=obs.dead,  ####QUA###
                          last.day= data.exit, pesi=pesi,  per.opm=1))
  
  parameters.boot=par_boot[[g]]
  
  for(i in 1:n.boot){
    d<-optim.parallel.ns(parameters.boot[i,],i0=i0, kp=kp, obs=obs.dead, last.day= data.exit, pesi=pesi,  per.opm=0)
    d$D
  }
  
  f<-apply(parameters.boot,1,function(i) optim.parallel.ns(vec_parameters = i, i0=i0, kp=kp, obs=obs.dead,  ####QUA###
                                                           last.day= data.exit, pesi=pesi,  per.opm=0))
  
  
    mortiX<-lapply(f, function(x) x$D)
    morti_mat<-t(matrix(unlist(mortiX), nrow=t.o))
    q1_d<-apply(morti_mat[index_buoni,], 2, q1)
    q2_d<-apply(morti_mat[index_buoni,], 2, q2)
    plot(as.Date(1:length(obs.dead), origin="2020-08-01"),sim$D,type="l", ylab="Deaths", xlab="", cex=0.2, main=REG[g], ylim=c(0, max(q2_d)+10))
    lines(as.Date(1:length(obs.dead), origin="2020-08-01"),obs.dead,type="p",col="#267D5D", cex=0.5)
    #lines(as.Date(1:length(obs.dead), origin="2020-08-01"),q1_d , lty=2)
    #lines(as.Date(1:length(obs.dead), origin="2020-08-01") ,q2_d , lty=2)
    x=c(as.Date(1:length(obs.dead), origin="2020-08-01"), rev(as.Date(1:length(obs.dead), origin="2020-08-01")))
    y=c(q1_d, rev(q2_d))
    polygon(x,y,  col = "#d3d3d3", border=NA)
    lines(as.Date(1:length(obs.dead), origin="2020-08-01"),obs.dead,type="p",col="#267D5D", cex=0.5)
    lines(as.Date(1:length(obs.dead), origin="2020-08-01"),sim$D)
  
}

if(schermo!=1){  
  dev.off()
}



#############################################################################
#plot  nuovi positivi 
##############################################################################

REG[12]="Trentino"
quale="Nuovi Positivi"
schermo=0
if(schermo!=1 ){
  nomefile=paste("fit_", letalita , datafine, quale,".pdf", sep="")
  pdf(nomefile, width=9, height = 13)
}
par(mfrow=c(5,4))
g=12
for(g in ordine_regioni){
  print(REG[g])
  Regione<-REG[g]
  S.0<-S0.vec[g]
  print(S.0)
  data <- DATA[DATA$region==Regione,] 
  i0<-data$positivi[data$date==(as.Date(data.ini)-1)] # altrove  nel boot
  data$deceduti<-data$deceduti-data$deceduti[data$date==(as.Date(data.ini)-1)]
  data<- data[which(data$date<=data.exit & data$date>=data.ini),]
  nuovi_p<-data$nuovi_pos
  pos<-data$positivi
  obs.dead<-morti[[g]]
  i0=positivi.ini[[g]]
  print(i0)
  if(g==12){
    nuovi_p<-aggregate(data$nuovi_pos, by =list(data$date), sum)$x
    pos<-aggregate(data$positivi, by =list(data$date), sum)$x
    REG[g]="Trentino Alto Adige"
  }
  
  MAT.SPLINE<-list.out[[g]]
  min_spline<-apply(MAT.SPLINE, 1, min)
  index_buoni<- index_buoni_list[[g]]
  sim<- optim.parallel.ns(VEC[[g]], i0=i0, kp=kp, obs=obs.dead,  ####QUA###
                          last.day= data.exit, pesi=pesi,  per.opm=0)
  
  print(optim.parallel.ns(VEC[[g]], i0=i0, kp=kp, obs=obs.dead,  ####QUA###
                          last.day= data.exit, pesi=pesi,  per.opm=1))
  
  parameters.boot=par_boot[[g]]
 
    nuovi_posi<-lapply(f, function(x) x$Inew.cal)
    nuovi_posi_mat<-t(matrix(unlist(nuovi_posi), nrow=t.o))
    q1_nI<-apply(nuovi_posi_mat[index_buoni,], 2, q1)
    nuovi_posi_mat[which(nuovi_posi_mat<0)]=0
    q1_nI[which(q1_nI<0)]=0
    q2_nI<-apply(nuovi_posi_mat[index_buoni,], 2, q2)
    
    simulata<- sim$Inew.cal
    simulata[which(simulata<0)]=0
    plot(as.Date(1:length(obs.dead), origin="2020-08-01"),simulata,type="l",col=1, main=REG[g],  xlab="",ylab="new positive", ylim=c(0, max(q2_nI)+10))
    # lines(as.Date(1:length(obs.dead), origin="2020-08-01"),nuovi_p,type="p",col="#267D5D", cex=0.5)
    x=c(as.Date(1:length(obs.dead), origin="2020-08-01"), rev(as.Date(1:length(obs.dead), origin="2020-08-01")))
    y=c(q1_nI, rev(q2_nI))
    polygon(x,y,  col = "#d3d3d3", border=NA) 
    lines(as.Date(1:length(obs.dead), origin="2020-08-01"),nuovi_p,type="p",col="#267D5D", cex=0.5)
    lines(as.Date(1:length(obs.dead), origin="2020-08-01"),simulata,type="l", cex=0.5)
    
  
}

if(schermo!=1){  
  dev.off()
}



#############################################################################
#plot positivi circolanti
##############################################################################

quale="Positivi Circolanti"
schermo=0
if(schermo!=1 ){
  nomefile=paste("fit_", letalita , datafine, quale,".pdf", sep="")
  pdf(nomefile, width=9, height = 13)
}
par(mfrow=c(5,4))
REG[12]="Trentino"
for(g in ordine_regioni){
  print(REG[g])
  Regione<-REG[g]
  S.0<-S0.vec[g]
  print(S.0)
  data <- DATA[DATA$region==Regione,] 
  i0<-data$positivi[data$date==(as.Date(data.ini)-1)] # altrove  nel boot
  data$deceduti<-data$deceduti-data$deceduti[data$date==(as.Date(data.ini)-1)]
  data<- data[which(data$date<=data.exit & data$date>=data.ini),]
  nuovi_p<-data$nuovi_pos
  pos<-data$positivi
  obs.dead<-morti[[g]]
  i0=positivi.ini[[g]]
  print(i0)
  if(g==12){
    nuovi_p<-aggregate(data$nuovi_pos, by =list(data$date), sum)$x
    pos<-aggregate(data$positivi, by =list(data$date), sum)$x
    REG[g]="Trentino Alto Adige"
  }
  
  MAT.SPLINE<-list.out[[g]]
  min_spline<-apply(MAT.SPLINE, 1, min)
  index_buoni<- index_buoni_list[[g]]
  sim<- optim.parallel.ns(VEC[[g]], i0=i0, kp=kp, obs=obs.dead,  ####QUA###
                          last.day= data.exit, pesi=pesi,  per.opm=0)
  
  print(optim.parallel.ns(VEC[[g]], i0=i0, kp=kp, obs=obs.dead,  ####QUA###
                          last.day= data.exit, pesi=pesi,  per.opm=1))
  
  parameters.boot=par_boot[[g]]
  
 
    
    posi_circ<-lapply(f, function(x) x$I)
    posi_circ_mat<-t(matrix(unlist(posi_circ), nrow=t.o))
    q1_I<-apply(posi_circ_mat[index_buoni,], 2, q1)
    q2_I<-apply(posi_circ_mat[index_buoni,], 2, q2)
    
    Icirco[[g]]<-sim$I
    Icircoq1[[g]]<-q1_I
    Icircoq2[[g]]<-q2_I
    plot(as.Date(1:length(obs.dead), origin="2020-08-01"),sim$I,type="l",col=1, xlab="", ylab="all positive", main=REG[g], ylim=c(0, max( q2_I)+10))
    x=c(as.Date(1:length(obs.dead), origin="2020-08-01"), rev(as.Date(1:length(obs.dead), origin="2020-08-01")))
    y=c(q1_I, rev(q2_I))
    polygon(x,y,  col = "#d3d3d3", border=NA) 
    
    lines(as.Date(1:length(obs.dead), origin="2020-08-01"),pos,type="p",col="#267D5D", cex=0.5)
    lines(as.Date(1:length(obs.dead), origin="2020-08-01"),sim$I,type="l", cex=0.5)
    
    t1=32
    t2=30
    t3=31
    t4=30
    t5=31
    as.Date(data.ini)-1+cumsum(c(t1,t2,t3,t4,t5))
    I_inizio_mese<-sim$I[cumsum(c(t1,t2,t3,t4,t5))]
    q1_inizio_mese<-q1_I[cumsum(c(t1,t2,t3,t4,t5))]
    q2_inizio_mese<-q2_I[cumsum(c(t1,t2,t3,t4,t5))]
    
    matrix_tab<-cbind(I_inizio_mese, q1_inizio_mese,q2_inizio_mese)
    rownames(matrix_tab)<-c("1settembre", "1ottobre", "1novembre", "1dicembre", "1gennaio")
    colnames(matrix_tab)<-c("Positivi circolanti", "q1", "q2")
    
    tabella_I[[g]]<-matrix_tab
    
    t6<-t.o-30
    t7<-14
    t8<-7
    
    as.Date(data.ini)-1+cumsum(c(t6,t7,t8  ))
    I_inizio_meseb<-sim$I[cumsum(c(t6,t7,t8  ))]
    q1_inizio_meseb<-q1_I[cumsum(c(t6,t7,t8  ))]
    q2_inizio_meseb<-q2_I[cumsum(c(t6,t7,t8  ))]
    
    matrix_tab2<-cbind(I_inizio_meseb, q1_inizio_meseb,q2_inizio_meseb)
    rownames(matrix_tab2)<-c("one moth ago", "2 weeks ago", "1 week ago")
    colnames(matrix_tab2)<-c("Positivi circolanti", "q1", "q2")
    
    tabella_I2[[g]]<-matrix_tab2
    
    
    
    
  
  
}

if(schermo!=1){  
  dev.off()
}

################################################################3################################################################3
#        R0 e tabella a 5 date fisse
################################################################3################################################################3
schermo=0
if(schermo!=1){
  nomefile=paste("boot_Regioni_", letalita , datafine,".pdf", sep="")
  pdf(nomefile,width=9, height = 13)
}
 

tabella_R0=list()

tabella_R02=list()
par(mfrow=c(5,4))

#par(mar = rep(2, 4))

g=8
REG[12]="Trentino"
REG[12]

#par(mfrow=c(1,1))

#par(mar = rep(2, 4))
for(g in ordine_regioni){
 
  Regione<-REG[g]
  S.0<-S0.vec[g]
  data <- DATA[DATA$region==Regione,] 
  data$deceduti<-data$deceduti-data$deceduti[data$date==(as.Date(data.ini)-1)]
  data<- data[which(data$date<=data.exit & data$date>=data.ini),]
  obs.dead<-data$deceduti
  if(g==5){
    obs.dead<- obs.dead+c(rep(0, 14), rep(-154, t.o-14))
  }
  i0<-data$positivi[data$date==data.ini]
  if(g==12){
    obs.dead<-aggregate(data$deceduti, by =list(data$date), sum)$x
    i0<-sum(data$positivi[data$date==data.ini])
  }
  MAT.SPLINE<-list.out[[g]]
  min_spline<-apply(MAT.SPLINE, 1, min)
  index_buoni<- index_buoni_list[[g]]
    
 # index_buoni=1:100
  
  vec_best=VEC[[g]]  
  q1_par <- apply(MAT.SPLINE[index_buoni,], 2, q1)
  q2_par <- apply(MAT.SPLINE[index_buoni,], 2, q2)
   
  
  ######################################################
  # PLOT
  ######################################################
  plot_giorni_prima=14
  morti_a_Zero<-sum(obs.dead==0)
  fit_best<-nspline.r0(coeff=VEC[[g]][2:9], N=length(obs.dead))
  fit_best[1:max(1,(morti_a_Zero-plot_giorni_prima))]=NA
  q1_par[1:max(1,(morti_a_Zero-plot_giorni_prima))]=NA
  
  q2_par[1:max(1,(morti_a_Zero-plot_giorni_prima))]=NA
  
  fit_best[which(fit_best<0)]=0
  q1_par[which(q1_par<0)]=0
  fdate<-as.Date(1:length(obs.dead), origin="2020-08-01")
  plot(fdate, fit_best,main=REG[g], type="l", col=1, ylim=c(0,6), ylab="R0", xlab="")
  
  abline(v=fdate[cumsum(c(31, 30 , 31, 30, 31, 31))], col= "gray90" )
  x=c(fdate, rev(fdate))
  y=c(q1_par, rev(q2_par))
  polygon(x,y,  col = "#d3d3d3", border=NA)
  
  lines(fdate, fit_best)
  abline(h=1, col=2)
  
  t1=32
  t2=30
  t3=31
  t4=30
  t5=31
  as.Date(data.ini)-1+cumsum(c(t1,t2,t3,t4,t5))
  R0_inizio_mese<-fit_best[cumsum(c(t1, t2, t3,t4, t5))]
  q1_inizio_mese<-q1_par[cumsum(c(t1, t2, t3,t4,t5))]
  q2_inizio_mese<-q2_par[cumsum(c(t1, t2, t3,t4,t5))]
  
  
  t6<-t.o-30
  t7<-14
  t8<-7
  
  as.Date(data.ini)-1+cumsum(c(t6,t7,t8  ))
    
  R0_inizio_meseb<-fit_best[cumsum(c(t6,t7,t8  ))]
  q1_inizio_meseb<-q1_par[cumsum(c(t6,t7,t8  ))]
  q2_inizio_meseb<-q2_par[cumsum(c(t6,t7,t8  ))]
  
  
  matrix_tab<-cbind(R0_inizio_mese, q1_inizio_mese,q2_inizio_mese)
  rownames(matrix_tab)<-c("1settembre", "1ottobre", "1novembre", "1dicembre","1gennaio")
  colnames(matrix_tab)<-c("R0", "q1", "q2")
  
  tabella_R0[[g]]<-matrix_tab
  
  matrix_tab2<-cbind(R0_inizio_meseb, q1_inizio_meseb,q2_inizio_meseb)
  rownames(matrix_tab2)<-c("1gennaio", "15gennaio", "22gennaio")
  colnames(matrix_tab2)<-c("R0", "q1", "q2")
  
  tabella_R0[[g]]<-matrix_tab
  tabella_R02[[g]]<-matrix_tab2
  
}

if(schermo!=1){
dev.off()
}
 


#### ############################################################3################################################################3
#       stampa delle tabelle di interesse
################################################################3################################################################3

y=list()
for(i in 1:20){
  f1<-matrix(t(tabella_R0[[i]]), ncol=1)
  #f2<-matrix(t(sapply(tabella_I[[i]], FUN=round)), ncol=1)
  y0<-t(cbind(f1))
  colnames(y0)<-rep(c(NA, "low", "up"),5)
  y[[i]]<-t(y0)
}





y[[1]]

head(unlist(y))

y2<-t(matrix(t(unlist((y))), nrow=15))
rownames(y2)<-REG


save(y2, file=paste("R0", letalita, sep=""))



d<-xtable(y2,  align="l|ccc|ccc|ccc|ccc|ccc",caption=paste("R0\_", letalita, sep=""))
print(d)






for(i in 1:20){
  f1<-matrix(t(tabella_R0[[i]]), ncol=1)
  f2<-matrix(t(tabella_I[[i]]), ncol=1)
  y0<-t(cbind( f2))
  colnames(y0)<-rep(c(NA, "low", "up"),5)
  y[[i]]<-t(y0)
}


y[[1]]

head(unlist(y))

y2<-t(matrix(t(unlist((y))), nrow=15))
rownames(y2)<-REG
d<-xtable(y2,  align="l|ccc|ccc|ccc|ccc|ccc", digits=0, caption=paste("positivi_", letalita))
 
print(d)


#################################################################################
#               Plot Tassi
#################################################################################


cols= c(1:8, "orange", "brown")

if(schermo!=1){
  pdf(paste("tassi", letalita, datafine,".pdf", sep=""))
}
plot(as.Date(1:length(obs.dead), origin="2020-08-01"),1000*Icirco[[1]]/ S0.vec[1], type="l", ylim=c(0,140), col=cols[1], ylab="Infections/1000 inh", xlab="")


for(i in 1:10){
  lines(as.Date(1:length(obs.dead), origin="2020-08-01"),1000*Icirco[[i]]/ S0.vec[i], col=cols[i],lwd=2)
}

for(i in 11:20){
  lines(as.Date(1:length(obs.dead), origin="2020-08-01"),1000*Icirco[[i]]/ S0.vec[i], col=cols[i-10], lty=2, lwd=2)
}

legend("topleft", legend=REG, col=cols, lty=c(rep(1,10),rep(2,10)), cex = 0.5)

if(schermo!=1){
  dev.off()
}

#########################################3`
#     Tabella prevalenze date inizio mese
#######################################3`

yx<-list()
yx2<-list()
i=1

for(i in 1:20){ 
  f1<-matrix(t(1000*tabella_I[[i]]/S0.vec[i]))
  y0<-t(cbind(f1))
  colnames(y0)<-rep(c(NA, "low", "up"),5)
  yx[[i]]<-y0
  f2<-matrix(t(1000*tabella_I2[[i]]/S0.vec[i]))
  y02<-t(cbind(f2))
  colnames(y02)<-rep(c(NA, "low", "up"),3)
  yx2[[i]]<-y02
}
 


y2x<-t(matrix(unlist((yx)), nrow=15))
rownames(y2x)<-REG

save(y2x, file=paste("prevalenza", letalita, sep=""))

d<-xtable(y2x,  align="l|ccc|ccc|ccc|ccc|ccc", digits=2,caption=paste("Rapporto_", letalita, sep=""))



#########################################3#########################################3`
#     Tabella prevalenze ultimo mese e due settimane fa e una settimana fa
#######################################3#########################################3``


y2x2<-t(matrix(unlist((yx2)), nrow=9))
rownames(y2x2)<-REG

save(y2x2, file=paste("prevalenza2", letalita, sep=""))

d<-xtable(y2x2,  align="l|ccc|ccc|ccc|", digits=2,caption=paste("Rapporto_", letalita, sep=""))


# 
# d
# tab1<-get(load(("prevalenza0114")))
# tab2<-get(load(("prevalenza005")))
# tab3<-get(load(("prevalenza0078")))
# tab4<-get(load(("prevalenza0179")))
          
#prev_reg<-cbind(tab1[,10:12],tab2[,10:12],tab3[,10:12],tab4[,10:12])


#d<-xtable(prev_reg,  align="l|ccc|ccc|ccc|ccc", digits=1,caption="prev_reg")
#d



# tab1<-get(load(("R00114")))
# tab2<-get(load(("R0005")))
# tab3<-get(load(("R00078")))
# tab4<-get(load(("R00179")))
# 
# r0_reg<-cbind(tab1[,10:12],tab2[,10:12],tab3[,10:12],tab4[,10:12])
# 
# d<-xtable(r0_reg,  align="l|ccc|ccc|ccc|ccc", digits=2,caption="R0_reg")
# d

#plot barrette

name<-paste("prevalenze", letalita, datafine,".pdf", sep="")
pdf(name, width=9, height=16)

yx2[[1]][,1:3]

plot(yx2[[1]][1,1], 1, xlim=c(0,40), ylim=c(0,21), ylab="", las=1, yaxt="n")

axis(2, at=1:20, labels=REG, las=2)


segments(x1=yx2[[1]][1,3], x0=yx2[[1]][1,2], y1=1, y0=1, lty = 3)
for(i in 1:20){
  points(yx2[[i]][,1], i)
  segments(x1=yx2[[i]][1,3], x0=yx2[[i]][1,2], y1=i, y0=i, lty = 3)
}
 
 
dev.off()
 
 

DATA[DATA$region %in% c("Trentino"),]$lat<-mean(DATA[DATA$region %in% c("Trentino"),]$lat)
name<-paste("prevalenze_3", letalita, datafine,".pdf", sep="")

pdf(name, height=9, width=15 )
par(mfrow=c(1,3), oma=c(4,7,3,8))
 order.geo<-order(unique(DATA$lat))
plot(yx2[[order.geo[1]]][1,1], 1, xlim=c(0,40), ylim=c(0,21), ylab="", las=1, yaxt="n", xlab="prevalenza (x1000)",main="un mese fa")
axis(2, at=1:20, labels=REG_order, las=2)

#for(h in order(unique(DATA$lat))){
  segments(x1=yx2[[order.geo[1]]][1,3], x0=yx2[[order.geo[1]]][1,2], y1=1, y0=1, lty = 3)
  for(h in 1:20){
    points(yx2[[order.geo[h]]][,1], h)
    segments(x1=yx2[[order.geo[h]]][1,3], x0=yx2[[order.geo[h]]][1,2], y1=h, y0=h, lty = 3)
  }
   
 abline(v=c(10,20,30,40), col="Gray90")
 
 
  plot(yx2[[order.geo[1]]][,4], 1, xlim=c(0,40), ylim=c(0,21), ylab="", las=1, yaxt="n",xlab="prevalenza (x1000)",main="due settimane fa")
 #axis(2, at=1:20, labels=REG_order, las=2)
 
 #for(h in order(unique(DATA$lat))){
 segments(x1=yx2[[order.geo[1]]][1,6], x0=yx2[[order.geo[1]]][1,5], y1=1, y0=1, lty = 3)
 for(h in 1:20){
   points(yx2[[order.geo[h]]][,4], h)
   segments(x1=yx2[[order.geo[h]]][1,6], x0=yx2[[order.geo[h]]][1,5], y1=h, y0=h, lty = 3)
 }
 
 abline(v=c(10,20,30,40), col="Gray90")
 
 
 
 plot(yx2[[order.geo[1]]][,7], 1, xlim=c(0,40), ylim=c(0,21), ylab="", las=1, yaxt="n",xlab="prevalenza (x1000)",main="una settimana fa")
 axis(4, at=1:20, labels=REG_order, las=2)
 
 #for(h in order(unique(DATA$lat))){
 segments(x1=yx2[[order.geo[1]]][1,9], x0=yx2[[order.geo[1]]][1,8], y1=1, y0=1, lty = 3)
 for(h in 1:20){
   points(yx2[[order.geo[h]]][,7], h)
   segments(x1=yx2[[order.geo[h]]][1,9], x0=yx2[[order.geo[h]]][1,8], y1=h, y0=h, lty = 3)
 }
 
 abline(v=c(10,20,30,40), col="Gray90")
 
dev.off()

name<-paste("prevalenze_2", letalita, datafine,".pdf", sep="")

DATA[DATA$region %in% c("Trentino"),]$lat<-mean(DATA[DATA$region %in% c("Trentino"),]$lat)
pdf(name, height=9, width=15 )
par(mfrow=c(1,2), oma=c(4,7,3,8))
order.geo<-order(unique(DATA$lat))
plot(yx2[[order.geo[1]]][1,1], 1, xlim=c(0,40), ylim=c(0,21), ylab="", las=1, yaxt="n", xlab="prevalenza (x1000)",main="un mese fa")
axis(2, at=1:20, labels=REG_order, las=2)

#for(h in order(unique(DATA$lat))){
segments(x1=yx2[[order.geo[1]]][1,3], x0=yx2[[order.geo[1]]][1,2], y1=1, y0=1, lty = 3)
for(h in 1:20){
  points(yx2[[order.geo[h]]][,1], h)
  segments(x1=yx2[[order.geo[h]]][1,3], x0=yx2[[order.geo[h]]][1,2], y1=h, y0=h, lty = 3)
}

abline(v=c(10,20,30,40), col="Gray65")


#plot(yx2[[order.geo[1]]][,4], 1, xlim=c(0,40), ylim=c(0,21), ylab="", las=1, yaxt="n",xlab="prevalenza (x1000)",main="due settimane fa")
#axis(2, at=1:20, labels=REG_order, las=2)

#for(h in order(unique(DATA$lat))){
#segments(x1=yx2[[order.geo[1]]][1,6], x0=yx2[[order.geo[1]]][1,5], y1=1, y0=1, lty = 3)
#for(h in 1:20){
#  points(yx2[[order.geo[h]]][,4], h)
#  segments(x1=yx2[[order.geo[h]]][1,6], x0=yx2[[order.geo[h]]][1,5], y1=h, y0=h, lty = 3)
#}

#abline(v=c(10,20,30,40), col="Gray90")



plot(yx2[[order.geo[1]]][,7], 1, xlim=c(0,40), ylim=c(0,21), ylab="", las=1, yaxt="n",xlab="prevalenza (x1000)",main="una settimana fa")
axis(4, at=1:20, labels=REG_order, las=2)

#for(h in order(unique(DATA$lat))){
segments(x1=yx2[[order.geo[1]]][1,9], x0=yx2[[order.geo[1]]][1,8], y1=1, y0=1, lty = 3)
for(h in 1:20){
  points(yx2[[order.geo[h]]][,7], h)
  segments(x1=yx2[[order.geo[h]]][1,9], x0=yx2[[order.geo[h]]][1,8], y1=h, y0=h, lty = 3)
}

abline(v=c(10,20,30,40), col="Gray65")

dev.off()
 
