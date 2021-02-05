
###########################################################
################# funzioni ################################
###########################################################

##################### funzione R0 bspline #####################

nspline.r0<-function(knots=NULL, coeff, N){
  
  ### b spline cubica per modellare R0
  # coeff vettore dei coefficienti (di lunghezza df)
  # N lunghezza della serie osservata
  # knots da specificare se si vogliono definire dei nodi 
  # specifici; i nodi devono essere length(knots)=df-2 ????? CONTROLLA
  
  df<-length(coeff)-1
  bmat<-ns(1:N, df=length(coeff)-1, knots=knots, 
           intercept = FALSE)
  r0<-coeff[1]+bmat%*%as.vector(coeff[-1])
  return(r0)
}


######### funzione che crea il SIRD con natural cubic spline #######################

optim.parallel.ns<- function(vec_parameters, i0, kp, obs=observed, last.day= data.exit, pesi=pesi, 
                             per.opm=1, knots=NULL){
  
  ####QUA###input i0, k
  
  
  # obs:       morti osservati (a partire dal primo morto, non dal tempo 0!)
  # obs_g:     dimessi guariti osservati (a partire dal giorno precedente al primo morto) 
  # obs_pos:   cumulata dei positivi circolanti a partire dal giorno precedente al primo morto (N.B. prima cumuli e poi tagli)
  
  l<-length(vec_parameters)
  num.r0<-(l-3)
  p<-vec_parameters[1]
  r0.vec<-vec_parameters[2:(num.r0+1)]
  dinf<-vec_parameters[num.r0+2]
  dmort<-vec_parameters[num.r0+3] 
  
  
  DEV.best<-NULL
  parameters<-NULL
  t.o<-length(obs)
  
  #p_v<-c(p,p+kp*(2:(t.o)))
  p_v<-rep(p,t.o)
  pstar<-p_v*dmort/(p_v*dmort+(1-p_v)*dinf)   #dmort=d tempo tampone-morte
  alpha.t<-(1-pstar)/dinf    #dtam=a+b tempo contagio-sintomi-tampone
  delta<-pstar/dmort
  
  r0<-nspline.r0(knots=knots, coeff=r0.vec, N=t.o)
  
  
  beta<-(r0*(alpha.t+delta))/S.0
  
  S<-c(S.0-i0,rep(NA, 200)) ####QUA### i0
  I<-c(i0,rep(NA, 200)) ####QUA### i0
  R.t<-c(0,rep(NA, 200))
  D<-c(0,rep(NA, 200))
  I.new<-c(0,rep(NA, 200))
  
  
  for (k in 1:199){ ####QUA### cambia il for , il while é un for controlls tuute le righe
    INF<- beta[k]*I[k]*S[k]    #nuovi infetti
    DEA<-I[k]*delta[k]
    REC<-I[k]*alpha.t[k]
    I.new[k+1] <- INF
    S[k+1]<-S[k]-INF
    I[k+1]<-I[k]+INF-REC-DEA
    R.t[k+1]<-R.t[k]+REC
    D[k+1]<-D[k]+DEA
  }
  
  # I.cal<-I[1:(t.o)]
  # D.cal<-D[1:(t.o)] #taglio la serie dei morti simulati dal giorno della prima morte
  # R.cal.t<-R.t[1:(t.o)]
  # I.new.cal<-I.new[1:(t.o)]
  # 
  # 
  I.cal<-I[2:(t.o+1)]
  D.cal<-D[2:(t.o+1)] #taglio la serie dei morti simulati dal giorno della prima morte
  R.cal.t<-R.t[2:(t.o+1)]
  I.new.cal<-I.new[2:(t.o+1)]

  
  ####### DEVIANZA  
  DEV.morti<-sqrt(mean((obs[1:length(obs)]-D.cal[1:length(obs)])^2,na.rm=TRUE))/mean(obs[1:length(obs)],na.rm=TRUE)
  
  ##################################
  #parameteri aggiuntivi che voglio vengano 
  #restituiti dalla funzione (in parameters)
  ##################################
  
  media_Devianze<-pesi[1]*DEV.morti
  # if(min(r0)<0){
  #   media_Devianze=1000
  #   }
  if(per.opm==1){
    return(media_Devianze)}
  if(per.opm==0){
    return(s=list(I=I.cal, D=D.cal, Rt=R.cal.t, Inew.cal=I.new.cal))   }
}

######################################################
##################################### versione di optim parallel da usare con lapply

optim.parallel.ns_dev<-function(par_all){
  return(optim.parallel.ns(par_all, i0=i0 ,kp=kp, obs=obs.dead, last.day= data.exit, ####QUA###
                           pesi=pesi, knots=NULL, per.opm=1))
}

########################################################
######################################### fop

fop.ns<-function(par_all){
  optim.parallel.ns(par_all, i0=i0, k=k, obs=obs.dead, last.day= data.exit, ####QUA###
                    pesi=pesi, knots=NULL, per.opm=1)
}


###########################################################
############################################# SIRD
######################################################################################

boot.SIRD<-function(N=100, D.est, i0, decessi.new.prima.morto=1, epsilon=epsilon.boot){####QUA###
  # N: numero di campioni bootstrap 
  # D.est: serie dei decessi stimata dal SIRD tramite la funzione calibrazione o optim. 
  
  # epsilon: margine per assicurare una minima variabilita anche nel caso in cui
  # la differenza D(t)-D(t-1)=0
  # la funzione restituisce in uscita una matrice che ha per colonne gli N campioni
  # bootstrap le cui prima n sono i morti
  
  n<-length(D.est)
  diff.D<-diff(D.est)
  diff.D[diff.D==0]<-epsilon
  
  mat.boot.D<-matrix(rep(NA,N*n),ncol=N)
  size.nb<-c(ceiling(decessi.new.prima.morto), ceiling(diff.D[1:(length(diff.D)-1)])) #### controlla!
  
  for(i in 1:N){
    sim.D.new<-c(D.est[1],rnbinom(n-1, size=size.nb,  mu=diff.D))
   # sim.D.new<-c(D.est[1],rpois(n-1, lambda = diff.D))
    sim.D<-cumsum(sim.D.new)
    mat.boot.D[,i]<-sim.D
  }   
  vec.boot.i0<-rpois(N,i0) ####QUA###
  
  return(list(mat.boot.D=mat.boot.D, vec.boot.i0=vec.boot.i0)) 
}

#########################################################################################
#########################################################################################
##CODICE per ARS
#########################################################################################

boot.calibr.par2<-function(output.boot.par_line){
  #prende un campione bootstrap e lo divide in tre tramite fx_boot e li fa diventare le osservazion
  #a quel punto ottimizza utilizzando optim 
  #restituisce i parametri stimati e il sird corrispondente
  #  k<-fx_boot(output.boot.par[,trial])
  morti.sim<-output.boot.par_line
  
  fop.boot<-function(vec_p){optim.parallel.ns(vec_p, obs=morti.sim,i0=i0, kp=kp,last.day= data.exit, pesi=pesi,  per.opm=1)}
  #x<-Rvmmin::Rvmmin(vec_parameters, fop.boot,  lower=lb, upper=ub, bdmsk=bdmsk.gc,control=list(maxfeval=3000, trace=TRUE))
   
 # x <- tryCatch(NlcOptim::solnl(vec_parameters,objfun =fop.boot,confun =  hin,lb=lower,ub=upper),
   #             error = function(e) return(list(par=rep(NA,length(vec_parameters)),fn='error')))
 #x<-NlcOptim::solnl(vec_parameters,objfun =fop.boot,confun =  hin2,lb=lower,ub=upper) 
  
 x<- nloptr::auglag(vec_parameters,fn =fop.boot,localtol = 1e-08,localsolver = 'LBFGS',control = list(maxeval = 10^4,xtol_rel = 1e-08), hin =  hin,lower=lower,upper=upper)
  #cambiato morti.sim invece che obs.dead
  c<-optim.parallel.ns(as.numeric(x$par[1:length(vec_parameters)]), obs=morti.sim,i0=i0, kp=kp,last.day= data.exit, pesi=pesi,  per.opm=0)    ###controlla qui
  return(list(parameters=as.numeric(x$par[1:length(vec_parameters)]),sird=c, devianza=x$value))
}


 
# 
###########################