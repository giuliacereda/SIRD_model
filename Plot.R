########codice per i grafici del paper SIRD semplice regioni#####

# LA lista VEC contenga i 20 vec_best (uno per ogni regione)
# La lista list.out contenga i 1000 bootstrap (uno per ogni regione)

mainDir<-dirname(rstudioapi::getActiveDocumentContext()$path) #here put the path to the folder

setwd(mainDir)
setwd(file.path(mainDir, "WS"))

load("500boot_0114p_2712.Rdata")
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

schermo=1
schermo=0
letalita=substr(as.character(p_seq), 3,7)

datafine=paste(substr(as.character(data.exit),9, 10), substr(as.character(data.exit),6, 7), sep="")
 

#########################
#   definisco index buoni
##########################

index_buoni_list<-list()
for(g in 1:20){
  MAT.SPLINE<-list.out[[g]]
  min_spline<-apply(MAT.SPLINE, 1, min)
  index_buoni<-which(min_spline >=0)
   
  
  index_buoni_list[[g]]<-index_buoni
}
g=19
MAT.SPLINE<-list.out[[g]]
min_spline<-apply(MAT.SPLINE, 1, min)
index_buoni<-which(min_spline >=-2)
index_buoni_list[[g]]<-index_buoni

g=2
MAT.SPLINE<-list.out[[g]]
min_spline<-apply(MAT.SPLINE, 1, min)
index_buoni<-which(min_spline >=-1)
index_buoni_list[[g]]<-index_buoni

g=12
MAT.SPLINE<-list.out[[g]]
min_spline<-apply(MAT.SPLINE, 1, min)
index_buoni<-which(min_spline >=-2)
index_buoni_list[[g]]<-index_buoni


lapply(index_buoni_list, function(x) length(x))

REG[19]
 #########################################################3################################################################3
#         FIT di morti positivi nuovipositivi e tabella a 4 date fisse
################################################################3################################################################3

schermo=0
tabella_I=list()

Icirco<-list()
Icircoq1<-list()
Icircoq2<-list()

quale="Decessi"
quale="Nuovi_Positivi"
quale="Positivi_Circolanti"
schermo=1
if(schermo!=1 ){
  nomefile=paste("fit_", letalita , datafine, quale,".pdf", sep="")
  pdf(nomefile, width=9, height = 13)
}


#par(mfrow=c(1,1))
par(mfrow=c(5,4))
#par(mar = rep(2, 4))
 

g=17
REG[12]="Trentino"
for(g in ordine_regioni){
  #for(g in 1:20){
    
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
  
  
  
  
  

  if(quale=="Decessi"){
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

  
  
  
  if(quale=="Nuovi_Positivi"){
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
  
  
  
  
  if(quale=="Positivi_Circolanti"){
    
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
   t4=22
   as.Date(data.ini)-1+cumsum(c(t1,t2,t3,t4))
   I_inizio_mese<-sim$I[cumsum(c(t1, t2, t3,t4))]
   q1_inizio_mese<-q1_I[cumsum(c(t1, t2, t3,t4))]
   q2_inizio_mese<-q2_I[cumsum(c(t1, t2, t3,t4))]
   
   matrix_tab<-cbind(I_inizio_mese, q1_inizio_mese,q2_inizio_mese)
   rownames(matrix_tab)<-c("1settembre", "1ottobre", "1novembre", "16novembre")
   colnames(matrix_tab)<-c("Positivi circolanti", "q1", "q2")
   
   tabella_I[[g]]<-matrix_tab
   
      }
   }



if(schermo!=1){  
dev.off()
}


################################################################3################################################################3
#        R0 e tabella a 4 date fisse
################################################################3################################################################3
schermo=0
if(schermo!=1){
  nomefile=paste("boot_Regioni_", letalita , datafine,".pdf", sep="")
  pdf(nomefile,width=9, height = 13)
}
 

tabella_R0=list()
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
    
  
  
  vec_best=VEC[[g]]  
  q1_par <- apply(MAT.SPLINE[index_buoni,], 2, q1)
  q2_par <- apply(MAT.SPLINE[index_buoni,], 2, q2)
  
  
  #q1_par0 <- apply(MAT.SPLINE0[,], 2, q1)
  
  #q2_par0 <- apply(MAT.SPLINE0[,], 2, q2)
  
  ######################################################
  # PLOT
  ######################################################
  plot_giorni_prima=14
  morti_a_Zero<-sum(obs.dead==0)
  fit_best<-nspline.r0(coeff=VEC[[g]][2:8], N=length(obs.dead))
  fit_best[1:max(1,(morti_a_Zero-plot_giorni_prima))]=NA
  q1_par[1:max(1,(morti_a_Zero-plot_giorni_prima))]=NA
  
  
  
  q2_par[1:max(1,(morti_a_Zero-plot_giorni_prima))]=NA
  
  fit_best[which(fit_best<0)]=0
  q1_par[which(q1_par<0)]=0
  
  plot(as.Date(1:length(obs.dead), origin="2020-08-01"), fit_best,main=REG[g], type="l", col=1, ylim=c(0,6), ylab="R0", xlab="")
  
  x=c(as.Date(1:length(obs.dead), origin="2020-08-01"), rev(as.Date(1:length(obs.dead), origin="2020-08-01")))
  y=c(q1_par, rev(q2_par))
  polygon(x,y,  col = "#d3d3d3", border=NA)
  
  lines(as.Date(1:length(obs.dead), origin="2020-08-01"), fit_best)
  abline(h=1, col=2)
 
  
  t1=32
  t2=30
  t3=31
  t4=22
  as.Date(data.ini)-1+cumsum(c(t1,t2,t3,t4))
  R0_inizio_mese<-fit_best[cumsum(c(t1, t2, t3,t4))]
  q1_inizio_mese<-q1_par[cumsum(c(t1, t2, t3,t4))]
  q2_inizio_mese<-q2_par[cumsum(c(t1, t2, t3,t4))]
  
  matrix_tab<-cbind(R0_inizio_mese, q1_inizio_mese,q2_inizio_mese)
  rownames(matrix_tab)<-c("1settembre", "1ottobre", "1novembre", "16novembre")
  colnames(matrix_tab)<-c("R0", "q1", "q2")
  
  tabella_R0[[g]]<-matrix_tab
# plot(as.Date(1:length(obs.dead), origin="2020-08-01"), fit_best,main=REG[g], type="l", col=1, ylim=c(0,6), ylab="R0", xlab="")
  
   #for(i in 1:100){
  #   lines(as.Date(1:length(obs.dead), origin="2020-08-01"), MAT.SPLINE[i,])
  # }
  
  
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
  f2<-matrix(t(sapply(tabella_I[[i]], FUN=round)), ncol=1)
  y0<-t(cbind(f1))
  colnames(y0)<-rep(c(NA, "low", "up"),4)
  
  y[[i]]<-t(y0)
}





y[[1]]

head(unlist(y))

y2<-t(matrix(t(unlist((y))), nrow=12))
rownames(y2)<-REG


save(y2, file=paste("R0", letalita, sep=""))



d<-xtable(y2,  align="l|ccc|ccc|ccc|ccc",caption=paste("R0_", letalita))
d





i=1
for(i in 1:20){
  f1<-matrix(t(tabella_R0[[i]]), ncol=1)
  f2<-matrix(t(tabella_I[[i]]), ncol=1)
  y0<-t(cbind( f2))
  colnames(y0)<-rep(c(NA, "low", "up"),4)
 
  y[[i]]<-t(y0)
}


y[[1]]

head(unlist(y))

y2<-t(matrix(t(unlist((y))), nrow=12))
rownames(y2)<-REG
d<-xtable(y2,  align="l|ccc|ccc|ccc|ccc", digits=0, caption=paste("positivi_", letalita))
 
print(d)

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


plot(Icirco)


#tabella rapporti#
yx<-list()
i=1

for(i in 1:20){ 
  f2<-matrix(t(1000*tabella_I[[i]]/S0.vec[i]))
  y0<-t(cbind(f2))
  colnames(y0)<-rep(c(NA, "low", "up"),4)
  yx[[i]]<-y0
}




y2x<-t(matrix(unlist((yx)), nrow=12))
rownames(y2x)<-REG

save(y2x, file=paste("prevalenza", letalita, sep=""))

d<-xtable(y2x,  align="l|ccc|ccc|ccc|ccc", digits=2,caption=paste("Rapporto_", letalita))


tab1<-get(load(("prevalenza0114")))
tab2<-get(load(("prevalenza005")))
tab3<-get(load(("prevalenza0078")))
tab4<-get(load(("prevalenza0179")))
          
prev_reg<-cbind(tab1[,10:12],tab2[,10:12],tab3[,10:12],tab4[,10:12])


d<-xtable(prev_reg,  align="l|ccc|ccc|ccc|ccc", digits=1,caption="prev_reg")
d



tab1<-get(load(("R00114")))
tab2<-get(load(("R0005")))
tab3<-get(load(("R00078")))
tab4<-get(load(("R00179")))

r0_reg<-cbind(tab1[,10:12],tab2[,10:12],tab3[,10:12],tab4[,10:12])

d<-xtable(r0_reg,  align="l|ccc|ccc|ccc|ccc", digits=2,caption="R0_reg")
d




