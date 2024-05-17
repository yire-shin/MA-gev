library(ismev)
library(Rsolnp)
library(numDeriv)
library(lmomco)

library(scoringRules)



for(weight in c('gLd', 'med', 'like', 'cvt')) {
  for (pick in c(.90, .95)){
     for (start in c('mle', 'lme')){

     revise = MA.gev.BC_11(data= data, quant=c(.95,.99),
                       start=start, pick=pick, weight=weight, 
                       numk= 5, B=500, fig=T, numom=3,
                       BC=F, t3.scale= 1.0, cut = -0.3 )
     
     cat(start,pick,weight,"MA.zp= ", revise$zp_estim,"\n")
   }
  }
}

# t3.scale =0 means that no bias correction. 
# usual t3.scale is 1.0 -- t3.scale should be between 0 to 2.0

# choice of weight = 'gLd', 'med', 'like', 'cvt' (=conventional mle and aic weight)
# pick confidence inteval = .90 or .95
# start.est =  'mle' or 'lme'

data=wando

for(weight in c('gLd', 'med', 'like', 'cvt')) {
  for (pick in c(.90, .95)){
    for (start in c('mle','lme')){

numsim=200

sam.bootL=list()
nsample= length(data)
numq=2
est=matrix(NA, numsim, numq); qua.mle=matrix(NA, numsim, numq);
qua.lme=matrix(NA, numsim, numq);

for (ib in 1:numsim){
  
  sam.bootL[[ib]]= sample(data, nsample, replace=T)   # nonpara Boot

  revise= MA.gev.BC_11(data=sam.bootL[[ib]],  quant=c(.95,.99),
                        start=start, pick=pick, weight=weight, 
                        numk= 5, B=200, fig=F, numom=3,
                        BC=F, t3.scale= 1.0, cut = -0.3 )

  # cat("MA.zp= ", revise$zp_estim,"\n")
   est[ib,1:numq]= revise$zp_estim[1:numq]
  # qua.mle[ib,1:numq]= revise$qua.mle[1:numq]
  # qua.lme[ib,1:numq]= revise$qua.lme[1:numq]
  # print(ib)
  
}

  cat(start,pick,weight, sqrt(var(est[,1])), sqrt(var(est[,2])) , "\n")
#  print( sqrt(var(qua.mle[,1])) );  print( sqrt(var(qua.mle[,2])) );
#  print( sqrt(var(qua.lme[,1])) );  print( sqrt(var(qua.lme[,2])) )

# print(mean(est[,1]) ); print(mean(est[,2]) )
# print(mean(qua.mle[,1]) ); print(mean(qua.mle[,2]) )
# print(mean(qua.lme[,1]) ); print(mean(qua.lme[,2]) )

    }
  }
}

# 
#  hist(est[,2], xlim=c(300,700), ylim=c(0,35), nclass=13 ) ; abline(v=mean(est[,2]), lty=2); par(new=T)
#  hist(qua.lme[,2], xlim=c(300,700), ylim=c(0,35), nclass=13, col="pink"); abline(v=mean(qua.lme[,2]))



#----------------------------------------------------------------------------
MA.gev.BC_11 = function(data=NULL, quant=c(.95,.99),
                        start='mle', pick=.95, weight='gLd', 
                        numk=3, B=500, fig=T, numom= 3, 
                        BC = T, t3.scale= 1.0, cut = -0.3) 
{
  
  # start.est = 'mle' or 'lme'
  # pick: CI level between 0.5 and 0.99 (recommend .90 or .95)
  # choice of weight = 'gLd', 'med', 'like', 'cvt' (=conventional mle and aic weight)
 
  # t3.scale =0 means that no bias correction.
  # usual t3.scale is 1.0, applied for 'gLd' and 'med' weights
  # -- t3.scale should be between 0 to 2.0
  
  zx=list(); hosking=list()
  SMALL <- 1e-05
  numq=length(quant)
  
  names_quant=paste("q",as.character(quant[1:numq]), sep="")
  names_numk=rep(c("numk"),numk)
  method=paste(start,pick*100,weight,sep="")
  names_par=c("mu","sigma","shape")
  
  para3= matrix(NA,nrow=numk, ncol=3)
  wtgd= rep(NA, numk)
  zp= matrix(NA, numq, numk)   
  zpf=rep(NA, numq, dimnames=list(names_quant) )  
  
  # -------begin of action ----------------------------
  
  nsample=length(data)
  zx$type =list(method=method, num_pick=numk, B=B)
  
  # ------- mle and se by delta method-------------------
  
  delta= gev.rl.delta_new(data, ntry=10, quant)
  zx$qua.mle= delta$qua.mle
  zx$mle.hosking =  delta$mle
  
  # mle.crps = crps_gev( zx$qua.mle, loc=delta$mle[1], 
  #                      scale=delta$mle[2], shape= -delta$mle[3] )

  #-----  Lme and se by bootstrap --------------------
  
  hosking =lme.boots(data, B, quant)
  zx$qua.lme= hosking$qua.lme
  zx$lme = hosking$lme
  
  hosking$start = start; hosking$numk = numk; 
  hosking$BC = BC; hosking$cut = cut; hosking$weight = weight
  
  # lme.crps = crps_gev( zx$qua.lme, loc=zx$lme[1], 
  #                      scale=zx$lme[2], shape= -zx$lme[3] )

  #---  xi_k picking for cand submodels ---------------------------
  
  if(BC == F) t3.scale =0
  hosking$t3.scale = t3.scale
  
  kpar.list = cand.xi(pick, nint=128, delta, hosking, figure=fig)
  kpar = kpar.list$kpar

  #------- Zp and weights computing-----------------------------------
  
  mywt = weight.com(data, hosking, kpar, numom=numom)
  
  # # --- crps computing for MA  ------------------
  # 
  # para3 = mywt$prob.call$mle3
  # tcr= matrix(NA, numq, numk)
  # 
  #    if(BC== T) {
  #      for (iq in 1:numq){
  #        tcr[iq,1:numk] = crps_gev(mywt$zp[iq,1:numk], 
  #                                  loc =   mywt$mu.BC[1:numk, iq], 
  #                                  scale = para3[1:numk,2],  
  #                                  shape = -para3[1:numk,3])
  #      }
  #      
  #    }else if(BC == F){
  #      for (iq in 1:numq){
  #        tcr[iq,1:numk] = crps_gev(mywt$zp[iq,1:numk], 
  #                                  loc =   para3[1:numk,1], 
  #                                  scale = para3[1:numk,2],  
  #                                  shape = -para3[1:numk,3])
  #      }
  #    }   # end if
  # 
  #   MA.crps = tcr%*%mywt$wtgd

 #  crps_gev(mywt$zp, loc=100, scale=30, shape=0.45)%*%revise$weights
 #--------------------------------------------------------
  
  #------------final estimates -----------------------------------
  
  zp = mywt$zp
  wtgd= mywt$wtgd

  if( any(zp=='NaN') || any(wtgd=='NaN') ) {
    cat("Warning: zp or wtgd is NaN", "\n") }

  zpf= t(wtgd)%*%t(zp)
  zx$para_estim_MA = t(wtgd)%*% para3   
  
  if(fig==T & start=='mle'){
    lwd.w = wtgd*10
    lwd.w[which(lwd.w < 1.5)] = 1.5
    lwd.w[which(lwd.w > 5)] = 5

    abline( v= kpar, col="seagreen", lty=3, lwd=lwd.w)
  }
  
  # if(fig ==T & start=='mle'){
  #   abline(v = zx$para_estim_MA[3], col="pink", lwd=2, lty=6)
  #   text(x=  zx$para_estim_MA[3], y= kpar.list$ymin +0.4, "MA")
  # }
 
  # zpdiff=matrix(NA, numq, numk)  
  # for(ip in 1:numk){  
  #   zpdiff[1:numq,ip]=(zp[1:numq,ip]-zpf[1:numq])^2
  # }
  # mscm_zp = t(wtgd)%*%t(zpdiff)     # mean square of zp due to candidate models

  zx$zp= zp
  zx$zp_estim = zpf
  if(any(zpf=='NaN')) cat("Warning: MA estim is NaN", "\n")
  
  zx$weights = wtgd
  zx$pick_xi = kpar      
  if(weight=='gLd' | weight =='med'){
    zx$t3.scale = t3.scale
  }
  # zx$MA.crps = MA.crps
  # zx$lme.crps = lme.crps
  # zx$mle.crps = mle.crps
  
  return(zx)
}
# ---end of main program --------------------------------------- 
#------------------------------------------------------------
gev.max=function (xdat, ntry=20) 
{
  z <- list();  k =list()
  n=ntry
  
  nsample=length(xdat)
  z$nsample=nsample
  
  init= matrix(0, nrow=ntry, ncol=3)
  init <- ginit.max(xdat,ntry)
  
  #-------------------------------------------------  
  gev.lik.max <- function(a) {
    
    mu <- a[1]      #mulink(mumat %*% (a[1:npmu]))
    sc <- a[2]      #siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
    xi <- a[3]      #shlink(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))
    
    y <- (xdat - mu)/sc
    y <- 1 - xi * y        # park modify to negative, for xi in hosking
    
    for (i in 1:nsample){
      y[i] = max(0, y[i], na.rm=T) }
    
    if (any(y <= 0) || any(sc <= 0)) 
      return(10^6)
    
    if( abs(xi) >= 10^(-5) ) {ooxi= 1/xi
    }  else  {ooxi=sign(xi)*10^5}
    
    zz=nsample*(log(sc)) + sum( exp(ooxi *log(y)) ) + sum(log(y) * (1-(ooxi)) ) 
    
    return(zz)
  }
  #-------------------------------------------------------------
  tryCatch(
    for(i in 1:nrow(init)){
      
      value <- try(solnp(init[i,], fun=gev.lik.max, 
                         LB =c(-Inf,0,-5),UB =c(Inf,Inf,5), 
                         control=list(trace=0, outer.iter=10,
                                      delta=1.e-7, inner.iter=40, tol=1.e-5) ))
      
      if(is(value)[1]=="try-error"){
        k[[i]] <- list(value=10^6)
      }else{
        k[[i]] <- value
      }
      
    } #for
  ) #tryCatch
  
  optim_value  <-data.frame(num=1:n,value=sapply(k, function(x) x$value[which.min(x$value)]))

  optim_table1 <-optim_value[order(optim_value$value),]
  selc_num  <- optim_table1[1,"num"]
  
  x  <-k[[selc_num]]
  
#  mu <- x$par[1];  sc <- x$par[2];  xi <- x$par[3] 
  
  z$conv <- x$convergence
  z$nllh <- x$value[which.min(x$value)]
  z$mle <- x$par

  return(z)
}
#------------------------------------------------------
ginit.max <-function(data,ntry){
  
  n=ntry
  init <-matrix(rep(0,n*3),ncol=3)
  
  lmom_init = lmoms(data,nmom=5)
  lmom_est <- pargev(lmom_init)
  
  init[1,1]    <-lmom_est$para[1]
  init[1,2]    <-lmom_est$para[2]
  init[1,3]    <-lmom_est$para[3]
  
  maxm1=ntry; maxm2=maxm1-1
  init[2:maxm1,1] <- init[1,1]+rnorm(n=maxm2,mean=0,sd = 5)
  init[2:maxm1,2] <- init[1,2]+rnorm(n=maxm2,mean=3,sd = 3)
  init[2:maxm1,3] <- init[1,3]+runif(n=maxm2,min= -0.4,max=0.4)
  init[2:maxm1,2] = max(0.1, init[2:maxm1,2])
  
  return(init)
}
#------------------------------------------------
pargev.xifix= function (lmom, xifix= 0.1, checklmom = TRUE, ...) 
{

  para=rep(NA, 3)
  G=xifix
  para[3]= G
  
  if (length(lmom$L1) == 0) {
    lmom <- lmorph(lmom)
  }
 
   SMALL <- 1e-05
   if( abs(G) < SMALL) {
      para[1:2]= pargum(lmom)$para
      return(list(type = "gev", para = para, source = "pargev"))
   }

  GAM <- exp(lgamma(1 + G))
  para[2] <- lmom$L2 * G/(GAM * (1 - 2^(-G)))
  para[1] <- lmom$L1 - para[2] * (1 - GAM)/G

  return(list(type = "gev", para = para, source = "pargev"))
}
#---------------------------------------------------------
gev.profxi.mdfy=  function (data=NULL, mle=NULL, xlow, xup, 
                            pick = NULL, nint = 128, figure=T) 
{
  w <- list();  v <- numeric(nint)
  conf = pick
  
      xlow.t= min(xlow, xup);    xup.t = max(xlow, xup)
      xlow= max(-3.0, xlow.t, na.rm=T);  xup = min(3.0, xup.t, na.rm=T)

  x <- seq(xlow, xup, length.out = nint)
  sol <- c(mle[1], mle[2])      # mle = mle.coles
  
  #-------------------------------------------------------------   
  gev.plikxi <- function(a, data=data, xi=xi) {
    
    if (abs(xi) < 10^(-6)) {
      y <- (data - a[1])/a[2]
      if (a[2] <= 0) 
        l <- 10^6
      else l <- length(y) * log(a[2]) + sum(exp(-y)) + 
          sum(y)
    }
    else {
      y <- (data - a[1])/a[2]
      y <- 1 + xi * y
      if (a[2] <= 0 || any(y <= 0)) 
        l <- 10^6
      else l <- length(y) * log(a[2]) + sum(y^(-1/xi)) + 
        sum(log(y)) * (1/xi + 1)
    }
    l
  }
  #---------------------------------------------------------------   
  
  for (i in 1:nint) {
    xi <- x[i]         
    opt <- optim(sol, gev.plikxi, data=data, xi=xi)
    sol <- opt$par
    v[i] <- opt$value
  }
  
  for (i in 1:nint){
    if(v[i] >= 10^6) v[i]=NA }
  
  d <- data.frame(x=x,v=-v)  #%>% filter(v!=-10^6)
  d= na.omit(d)
  
  # --------------------------------------------------------- 
  irt=0
  idmax = which.max(d$v)
  if(idmax==1 ) {
    w$ci1 = xlow - 0.2
    w$ci2 = xup - 0.2
    irt=1
  }else if(idmax==length(d$x)){
    w$ci1 = xlow + 0.2
    w$ci2 = xup + 0.2
    irt=1               # irt=1 means the max is observed at the first or the last x
  }
  
  if(irt==1) {
    w$ci1_ci2 = c(w$ci1, w$ci2)
    w$ci_length <- w$ci2-w$ci1
    w$xmax= d$x[idmax]
    w$vmax= max(d$v)
    
    halfchi=  0.5 * qchisq(conf, 1)
    plot(-d$x, d$v, 
         xlim=c( -w$ci2-0.3, -w$ci1+0.3 ),
         ylim=c( w$vmax-(halfchi*2.7), w$vmax+halfchi),
         type = "l", xlab = "xi_hosking in GEV",  
         col="black",
         ylab = "Profile Log-likelihood", main=c("Profile CI for xi"))
    ma <-   w$vmax 
    abline(h = ma, col = 4)
    abline(h = ma - 0.5 * qchisq(conf, 1), col = 4)
    abline(v = c(-w$xmax, -w$ci1, -w$ci2), col=2, lty=2)
    
    w$ymin = w$vmax-(halfchi*2.7)+0.2
    
    return(w)
  }
  
  #-----------------------------------------------------------
  port=1.3
  x0= min(d$x); x9 = max(d$x) 
  bd0= max(0.2, abs(x0)) *sign(x0)
  bd9= max(0.2, abs(x9)) *sign(x9)
  
  xmin= min(-0.7, x0 - sign(x0)*port*bd0)
  xmax= max(0.7, x9 + sign(x9)*port*bd9)
  
#  if(xmax < x9) no need to extrapol
#  if(x0 < xmin) no need to extrapol

#  for (inter in 1:2){
  
  inter=1

#--- for linear interpolation ----------------------------------
#  if(inter==1){

  n.add= 32 
  yleft=rep(NA, n.add);  yright=rep(NA, n.add)
  
  if(xmax > x9) {
    idmost= which.max(d$x)
    b.over= (d$v[idmost]-d$v[idmost-1])/(d$x[idmost]-d$x[idmost-1])
    extra.over= seq( d$x[idmost], xmax, by=(xmax-d$x[idmost])/n.add )
    
    b.over = b.over * 1.2
    
    for (i in 1:length(extra.over)){
        yright[i] = d$v[idmost] -b.over*(d$x[idmost]- extra.over[i])
    }
  }
  
  if(x0 > xmin) {
    idlst= which.min(d$x)
    b.bef= (d$v[idlst+1]-d$v[idlst])/(d$x[idlst+1]-d$x[idlst])
    extra.bef= seq(xmin, d$x[idlst], by=( d$x[idlst]-xmin)/n.add )
    
    b.bef = b.bef * 1.2
    
    for (i in 1:(length(extra.bef) ) ){
      yleft[i] = d$v[idlst] -b.bef*(d$x[idlst]- extra.bef[i])
    }
  }
  
  msp1 = data.frame(x= c(extra.bef, d$x, extra.over), y= c(yleft, d$v, yright) )
  
  dv1 = msp1$y
  
  w1= comp.prof.ci(d=msp1, v= dv1, conf)
  
#  --- end of linear interpolation  --------------------------------
  # }else if(inter==2){
  # 
  #  nint.msp = round(nint*port+.5)
  # 
  #  msp.spl= spline(d$x, d$v, n=nint.msp, method=c("natural"),
  #                  xmin=xmin, xmax=xmax)
  # 
  #  msp2 = data.frame(x = msp.spl$x, y= msp.spl$y)
  #  
  #  dv2 = msp2$y
  #  
  #  w2= comp.prof.ci(d=msp2, v=dv2, conf)
   
#  }  # end if

# --- end of spline interpolation -------------------------------

#  }  # end for
   
#  if(w1$ci_length < w2$ci_length) {
    w = w1 
    inter = 1
    d= msp1
    d$v = msp1$y
  # }else{
  #   w = w2 
  #   inter = 2
  #   d= msp2
  #   d$v = msp2$y
  # }

  if(figure == T){
    
    halfchi=  0.5 * qchisq(conf, 1)
    # ci_max = w$ci2
    # ci_min = w$ci1
    
    xlim0 = -(w$ci2 + (w$ci_length)/4 ) 
    xlim1 = -(w$ci1 - (w$ci_length)/4 )
    
    if(inter==1) {co="green"} else if(inter==2) {co="blue"}
    
    plot(-d$x, d$v, 
         xlim=c(xlim0+0.02, xlim1-0.02),
         ylim=c( w$vmax-(halfchi*2.7), w$vmax+ 0.8*halfchi),
         type = "l", xlab = "xi_Hosking in GEV",  
         col= co, lwd= 2,
         ylab = "Profile Log-likelihood",main=c("Profile CI for xi"))
    ma <-   w$vmax  #-z$nllh
    abline(h = ma, col = 4)
    abline(h = ma - 0.5 * qchisq(conf, 1), col = 4)
    abline(v=c(-w$xmax, -w$ci1, -w$ci2), col=2, lty=2, lwd= 2)
#    abline(v = c(-mle[3]), lty=3, col="black", lwd=2)
    text( x=c(-mle[3]), y= w$vmax-(halfchi*2.7)+0.2, labels="mle")
    
    lm= lmoms(data)
    lme= pargev(lm)$para[3]
    
    abline( v = lme, lty=5, lwd=2, col="purple")
    text( x =lme, y= w$vmax-(halfchi*2.7)+0.2, labels="lme")
    
    w$ymin = w$vmax-(halfchi*2.7)+0.2
  }
  
  return(w)
}
#------------------------------------------------------------------------
comp.prof.ci = function(d, v, conf=NULL){
  
  w = list()
  d$v = v
  
  w$vmax= max(d$v)
  nmsp = length(d$v)
  
  w$nllh  <- w$vmax                                   # -z$nllh
  w$chisq <- w$vmax - 0.5 * qchisq(conf, 1)
  
  diff  <- d$v - w$chisq
  idmax = which.max(d$v)
  w$xmax = d$x[idmax] 
  
  w$ci1   <- d$x[which.min(abs(diff[1:idmax]))]
  w$ci2   <- d$x[(idmax+1):nmsp][which.min(abs(diff[(idmax+1):nmsp]))] 
  
  ci_min =min(w$ci1, w$ci2)
  ci_max =max(w$ci1, w$ci2)
  w$ci1 =ci_min
  w$ci2 =ci_max
  
  if (w$ci2 < w$xmax) {
    cat("Routine fails, try changing plotting interval", 
        fill = TRUE)
  }
  
  w$ci1_ci2 = c(w$ci1, w$ci2)
  w$ci_length <- w$ci2-w$ci1
  
  return(w)
}
# ----------------------------------------------------------------------  
#--------------------------------------------------------------  
gev.rl.delta_new <- function(data, ntry=10, quant){
  
  z=list()
  numq=length(quant)
  names_par= c("mu","sigma","shape")
  names_quant=paste("q",as.character(quant[1:numq]),sep="")
  
  kmle.h = gev.max(data, ntry=ntry)
  
  kmle.coles = gev.fit(data, show=F)
  
  if(kmle.h$nllh < kmle.coles$nllh) {
    better.mle= kmle.h$mle
  }else{
    better.mle= kmle.coles$mle
    better.mle[3]= -kmle.coles$mle[3]
  }
  
  savem= vec2par(better.mle,'gev')
  z$mle = better.mle
  z$qua.mle= quagev(quant[1:numq],savem)
  z$data = data
  
  names(z$qua.mle)=names_quant
  names(z$mle)=names_par

  z$quant=quant
  return(z)
}
#------------------------------------------------------------
cand.xi <- function(pick=.95, nint=128, delta, hosking, figure=T){
  
  SMALL= 1e-05
  BC = hosking$BC; t3.scale = hosking$t3.scale; numk= hosking$numk
  start = hosking$start; weight = hosking$weight; # cut = hosking$cut;
  
  if(start=='mle'){
    
    newse= max(0.2,abs(delta$mle[3]))*0.07

    mle.coles= delta$mle
    mle.coles[3]= -delta$mle[3]
    
    xlow= mle.coles[3]-8*newse
    xup=  mle.coles[3]+8*newse
    
    xlow= max(-4, min(0.5,xlow), na.rm=T)
    xup=  min(4, max(-0.5,xup), na.rm=T)
    xlow2 = min(xlow, xup)
    xup2 = max(xlow, xup)

    get.ci = gev.profxi.mdfy(data=delta$data, mle=mle.coles, 
                             xlow=xlow2, xup= xup2,
                             pick=pick, nint=nint, figure=figure)
    
    lo.lim =min(-get.ci$ci1, -get.ci$ci2)
    up.lim =max(-get.ci$ci1, -get.ci$ci2)
    
  }else if(start=='lme') {
    
    lo.lim= quantile( hosking$lme.boot[,3], prob=(1-pick)/2, na.rm=T)
    up.lim= quantile( hosking$lme.boot[,3], prob=1-(1-pick)/2, na.rm=T)

    lo.lim= max(-.9999, lo.lim)
  }
  
  lo.lim2 = min(lo.lim, up.lim);   up.lim2= max(lo.lim, up.lim)
  
  ra.lim = up.lim2 - lo.lim2
  fence= 0.1* pick
  dam =  1.2* pick
  mle.h = delta$mle[3]
  
  if(start=='mle') {

    if( lo.lim2 > mle.h | up.lim2 < mle.h ) {

      lo.lim = mle.h - sign(mle.h)*(fence) 
      up.lim = mle.h + sign(mle.h)*(fence)
      
      lo.lim2 = min(lo.lim, up.lim);  up.lim2= max(lo.lim, up.lim)
    }
    
    if( (lo.lim2 < mle.h) & (up.lim2 > mle.h) ){
      
      iwork=T
      if( up.lim2 - lo.lim2 < fence ) {
        iwork=F        # for the case when ra.lim is too small
        
        while(iwork==F){

          bd.lo= max(0.1, abs(lo.lim2)); 
          bd.up= max(0.1, abs(up.lim2)); 
          
          lo.lim2 = lo.lim2 - (fence/3)*bd.lo
          up.lim2 = up.lim2 + (fence/3)*bd.up

          if( up.lim2 - lo.lim2 >= fence ) iwork=T 
        } # end while
        
      }else if( up.lim2 - lo.lim2 > dam ){
        
        lo.lim = mle.h - sign(mle.h)*(dam/2) 
        up.lim = mle.h + sign(mle.h)*(dam/2)
        
        lo.lim2 = min(lo.lim, up.lim);  up.lim2= max(lo.lim, up.lim)

      } #end if-else if
      
    }  #end if
    
  } #end if start mle_prof
  
  up.lim2 = max(up.lim2, -0.45)
  lo.lim2 = min(lo.lim2,  0.45)
  
  iwork=T
  if( up.lim2 - lo.lim2 < fence ) {
    iwork=F   # for the case when ra.lim is too small
    
    while(iwork==F){
      bd.lo= max(0.1, abs(lo.lim2));
      bd.up= max(0.1, abs(up.lim2)); 
      
      lo.lim2 = lo.lim2 - (fence/3)*bd.lo
      up.lim2 = up.lim2 + (fence/3)*bd.up
      
      if( up.lim2 - lo.lim2 >= fence ) iwork=T 
    } # end while
  } #end if
  
  if(start=='mle' & BC==T ){
    lo.lim2 = lo.lim2 + abs(lo.lim2)*t3.scale/20
  }
  
  kpar = seq(from=lo.lim2, to=up.lim2, by=(up.lim2-lo.lim2)/(numk-1)) 
  
  if( any(abs(kpar) <= SMALL)) {
    kpar[ which(abs(kpar) <= SMALL) ] = SMALL*1.2 }
  
  # if(figure==T){
  #   abline( v=kpar, col="seagreen", lty=3, lwd=3) }
  
  wx = list()
  wx$kpar = kpar
  if(start=='mle') wx$ymin = get.ci$ymin
    
  return(wx)
}

#------- weight computing----------------------------------------
weight.com <- function(data, hosking, kpar, numom=3){
  
  z=list();  
  quant = hosking$quant; numk= hosking$numk; weight = hosking$weight
  
  numq= length(quant)
  later = rep(NA, numq)
  mu.BC = matrix(NA, numk, numq)
  
  B = hosking$B; BC = hosking$BC; t3.scale = hosking$t3.scale
  cut = hosking$cut;  start = hosking$start
  if(BC == F) t3.scale = 0
  
    if( t3.scale < 0  | t3.scale >= 3) {
      cat("Wrong specification for t3.scale ", t3.scale, "\n")
      stop()
    }

  zp= matrix(NA,nrow=numq, ncol=numk)
  wtgd=rep(NA,numk)
  
  if(weight == 'gLd' | weight == 'med' | weight =='cvt'){  

    prob.call = prob.dist.wt.sing.noboot(data, hosking, 
                                         kpar, numom, ntry=20)
    
    if(weight == 'cvt') {
       naic= prob.call$aic
       
       amin =naic -min(naic, na.rm=T)
       wtgd =exp(-amin)/sum( exp(-amin), na.rm=T)
    }
    
    if(weight == 'gLd' | weight == 'med'){
      if(weight=='gLd') {kcol=1} else if(weight=='med') {kcol=2}
    
      prob.mtx= prob.call$prob.mtx
      wtgd[1:numk]= prob.mtx[1:numk,kcol]/sum(prob.mtx[,kcol], na.rm=T) 

    }
    
    sigma = prob.call$mle3[,2]
    simple =  prob.call$mle3[,3]
      
    for (ip in 1:numk){
      if( any( is.na(prob.call$mle3[ip,1:3]) ) ) {
        zp[1:numq,ip] = 0.0
      }else{
        save= vec2par( prob.call$mle3[ip,1:3],'gev' )
        zp[1:numq,ip]= quagev(quant[1:numq],save)
      }
      
      if(BC == T){
      if(weight == 'gLd' | weight == 'med' ){

          zp[,ip] = zp[,ip]*(1 + prob.call$t3.diff[ip] ) 

          if(start=='lme'){
            
            if( t3.scale > 0 & t3.scale < 3 ){
            
              if( simple[ip] <  cut &  simple[ip] >  -0.6) {              
                 zp[,ip] = zp[,ip]*(1 - (simple[ip] -cut)/(4 -t3.scale) )  }
            
              if( simple[ip] <=  -0.6 ){
                 zp[,ip] = zp[,ip]*(1 - (-0.6 -cut)/(4 -t3.scale) )        }
            }
            
          }      # end if start==lme
          
        later[1:numq] = 1- (-log(quant[1:numq]))^(simple[ip])
        mu.BC[ip,1:numq] = zp[1:numq,ip] -later[1:numq]* sigma[ip]/simple[ip]
      }
      }          # end if BC=T
      
    }            # end for ip

  }else if(weight=='like'){
    
    prob.call= weight.lik.xifix(data, numk=numk, kpar)
    
    for (ip in 1:numk){
      if( any(prob.call$kfix[1:3,ip]=="NaN") ) {
        iNaN = which(prob.call$kfix[1:3,ip]=="NaN")
        prob.call$kfix[1:3,ip]= 0.1
        prob.call$wt[iNaN]=0.0
      } #end if
    } # end for
    
    wtgd = prob.call$wt
    
    prob.call$mle3 = t(prob.call$kfix)
    
    for (ip in 1:numk){
      save= vec2par(prob.call$kfix[1:3,ip],'gev')
      zp[1:numq,ip]= quagev(quant[1:numq],save)
    }   # end for ip

  }  # end of if weight= 'like'
  
  z$weight= weight
  z$wtgd= wtgd
  z$zp=zp
  z$prob.call=prob.call
  z$mu.BC = mu.BC
  # z$sigma = sigma
  # z$xi = simple

  return(z)
}
#--- getting weights based on bic -----
# --------------------------------------------------------------
weight.lik.xifix = function(data,numk,kpar){
  
  z=list(); work=list()
  kfix=matrix(NA,nrow=3,ncol=numk); 
  nbic=rep(NA,numk); bmin=rep(NA,numk)
  wt=rep(NA,numk)
  para=rep(NA,3); lm3=list()
  
  lm3=lmoms(data, nmom=5)
  for (ip in 1:numk){
    kfix[1:3,ip]= pargev.xifix(lm3, xifix=kpar[ip])$para[1:3]
  }
  
  for (ip in 1:numk){
    
    nbic[ip]= gev.xilik( kfix[1:2,ip], xifix= kpar[ip], xdat= data )
    
    if(nbic[ip] >= 10^6) nbic[ip]=NA
  }

  numpar=3
  nbic = 2*nbic + 2*numpar   # AIC
  bmin =nbic -min(nbic, na.rm=T)
  wt =exp(-bmin)/sum(exp(-bmin), na.rm=T)
  
  for (ip in 1:numk){ 
    if( is.na(wt[ip]) ) wt[ip]= 0 }

  z$wt= wt
  z$bic= nbic
  z$kfix= kfix
  
  return(z)
}
#--- getting the distance probability for weights calculation -----
# --------------------------------------------------------------
prob.dist.wt.sing.noboot = function(data, hosking,
                                    kpar, numom=NULL, ntry=15){
  
  numk= hosking$numk; weight = hosking$weight
  B = hosking$B;  BC = hosking$BC; 
  t3.scale = hosking$t3.scale;  cut = hosking$cut
  
  if(weight == 'gLd' | weight == 'med' | weight =='cvt') { 
  }else{ cat("wrong weight= ", weight, "\n")  }

  kfix=list();  zw=list()
 
  # dd= matrix(NA, nrow=numom, ncol=2)
  # gdd= matrix(NA, nrow=numk, ncol=2); 
  # prob.mtx= matrix(NA,nrow=numk, ncol=2); 
  mle3= matrix(NA, nrow=numk, ncol=3)
  lme.boot= matrix(NA, nrow=B, ncol=numom)
  lmB.med= rep(NA,B);   aic=rep(NA,numk)
  names_par= c("mu","sigma")

  nsample= length(data)
  lm= lmoms(data, nmom=5)
  medata= median(data)
  
  for (ip in 1:numk){
    
    kfix[[ip]]= gev.xifix.sing(data, xifix=kpar[ip], ntry=ntry)

    aic[ip]= 2*kfix[[ip]]$nllh + 6
    mle3[ip,1:3]= kfix[[ip]]$mle[1:3]
  }

  if(weight =='gLd'){
    
      cov.lm = cov( hosking$ratios[,1:numom])
      Vinv= solve(cov.lm)
      detV= det(cov.lm)

      zw$cov.lm=list(Vinv=Vinv, detV=detV)
      
  }   #end weight=gLd
  
  if(weight=='med'){
  
    lme.bb = hosking$ratios
    lme.bb[1:B,1]= hosking$med[1:B]

       cov.lm.med= cov( lme.bb[,1:numom])

    Vinv.med= solve(cov.lm.med)
    detV.med= det(cov.lm.med)
    
    Vinv= Vinv.med
    detV= detV.med
    
    zw$cov.lm=list(Vinv=Vinv.med, detV=detV.med)
    
  }   # end if med

  if( weight == 'gLd' | weight =='med'){

    normal= comp.prob.mtx(data, kfix=kfix, Vinv=Vinv, detV=detV, 
                          numom=numom, hosking)  
    
    # if(weight=='gLd') {kcol=1} else if(weight=='med') {kcol=2}
  
    # if( max( normal$prob.mtx[1:numk,kcol] ) < 1e-50  ) {
    # 
    #   cov.lm=matrix(0,numom,numom)
    #   diag(cov.lm)=1
    #   Vinv=cov.lm;  detV=1.0
    # 
    #   zw$cov.lm= list(Vinv=Vinv, detV=detV)
    # 
    #   normal= comp.prob.mtx(data, kfix=kfix, Vinv=Vinv, detV=detV, weight, 
    #                         numk, numom, t3.scale)
    # }

    for (ip in 1:numk){
      if( any( is.na(kfix[[ip]]$mle[1:3]) ) ) {
         mle3[ip,1:3]= NA
         normal$prob.mtx[ip,] = 0.0
      }else{
         mle3[ip,1:3]= kfix[[ip]]$mle[1:3]  }            # end if
    }    # end for
    
  }     #end if weight =gLd or med
  
  zw$aic= aic
  zw$mle3= mle3
  zw$kfix= kfix

  if( weight == 'gLd' | weight =='med'){
    zw$t3.diff = normal$t3.diff
    zw$prob.mtx= normal$prob.mtx
    zw$gdd= normal$gdd
    zw$t3.scale = normal$t3.scale
  }
  
  return(zw)
}
#----------------------------------------------------------
comp.prob.mtx = function(data, kfix=NULL, Vinv=NULL, detV=NULL, 
                         numom=NULL, hosking)
{
  numk= hosking$numk; weight = hosking$weight

  ww1=list(); kfix.lm=list()
  dd=matrix(NA, nrow=numom, ncol=2)
  gdd=matrix(NA, nrow=numk, ncol=2);  simple =rep(NA, numk)
  prob.mtx=matrix(NA, nrow=numk, ncol=2); t3.diff= rep(NA, numk) 
  
  BC = hosking$BC; t3.scale = hosking$t3.scale;  cut = hosking$cut
  if(BC == F) t3.scale = 0

  medata= median(data)
  lm = lmoms(data, nmom=5)

  for (ip in 1:numk){
    
    if(kfix[[ip]]$conv != 0) {
      prob.mtx[ip,]=0
      print(paste("No convergence in gev.xifix, ip=",ip))
      
    } else {

      if(kfix[[ip]]$mle[3] <=  -1.0) kfix[[ip]]$mle[3]= -0.999
      
      simple[ip] = kfix[[ip]]$mle[3]
      
      savek =vec2par(kfix[[ip]]$mle,'gev')
      kfix.lm[[ip]]= lmomgev(savek)       # mle[3] >= -1 

      if(BC == T){
        t3.diff[ip]= kfix.lm[[ip]]$ratios[3] - lm$ratios[3]
      }
    }
    
  }  # end for ip

  # -----------------------------Bias Correction ----------------------
  if(BC == F){ t3.diff[1:numk] = 0 
  
  }else if(BC == T){
    
      id = which(t3.diff > 0 )
      for (i in id){
         t3.diff[i] = min(0.4, t3.diff[i]) 
      }

       nid = which(t3.diff < 0)
       t3.diff[nid] = 0

       t3.diff = t3.scale* t3.diff
      
        bigm = which( simple <  cut )
        if(t3.scale > 0  &  t3.scale < 3 ){
           t3.diff[bigm] = t3.diff[bigm]*(1 - (simple[bigm]-cut)/(4-t3.scale) )
        }

      t3.diff[ which(t3.diff >= 1) ] = 1
  }  # end if BC
# -------------------------------------------------------------------      

    if(weight =='gLd'){
      
      for (ip in 1:numk){
        
        dd[1,1]= lm$lambdas[1]- kfix.lm[[ip]]$lambdas[1] 
    
        dd[2:numom,1]= lm$ratios[2:numom]- kfix.lm[[ip]]$ratios[2:numom]

        gdd[ip,1]= t(dd[1:numom,1])%*% Vinv %*% dd[1:numom,1]

        gdd[ip,1]=  gdd[ip,1]*(1- t3.diff[ip])
        
        prob.mtx[ip,1]= exp(-gdd[ip,1]/2 )/( (2*pi)^(numom/2) * sqrt(detV) )
        
        if(is.na(prob.mtx[ip,1]) )  prob.mtx[ip,1]= 0
      }
    }
      
      if(weight=='med') {
        
        for (ip in 1:numk){
          
        savek= vec2par(kfix[[ip]]$mle,'gev')
        dd[1,2]= medata -quagev(0.5, savek)
        
        dd[2:numom,2]= lm$ratios[2:numom]- kfix.lm[[ip]]$ratios[2:numom]
        
        gdd[ip,2]= t(dd[1:numom,2])%*% Vinv %*% dd[1:numom,2]
        
        gdd[ip,2]=  gdd[ip,2]*(1- t3.diff[ip])

        prob.mtx[ip,2]= exp(-gdd[ip,2]/2 ) /( (2*pi)^(numom/2) * sqrt(detV) )
        
        if(is.na( prob.mtx[ip,2]) )  prob.mtx[ip,2]=0
        
        }   #end for ip =1, numk
      }   #end if med

  ww1$prob.mtx = prob.mtx
  ww1$gdd= gdd
  ww1$t3.diff = t3.diff[1:numk]
  ww1$t3.scale = t3.scale
  
  return(ww1)
}
#------------------------------------------------------
ginit.xifix <-function(data,ntry){
  
  initx <-matrix(0,nrow=ntry,ncol=2)
  
  lmom_init = lmoms(data,nmom=5)
  lmom_est <- pargum(lmom_init, checklmom=T )
  
  initx[1,1]    <- lmom_est$para[1]
  initx[1,2]    <- lmom_est$para[2]
  
  maxm1=ntry; maxm2=maxm1-1
  initx[2:maxm1,1] <- initx[1,1]+rnorm(n=maxm2,mean=0,sd = 7)
  initx[2:maxm1,2] <- initx[1,2]+rnorm(n=maxm2,mean=2,sd = 2)
  initx[2:maxm1,2] = max(0.1, initx[2:maxm1,2])
  
  return(initx)
}
# gev estimation under xi fixed --single program
#------------------------------------------------------------
gev.xifix.sing =function (xdat, xifix= -0.1, ntry=15)   
{
  zx <- list();  kx =list()
  truepar=rep(0,2)

  init.xifix= matrix(0, nrow=ntry, ncol=2)
  
  init.xifix <- ginit.xifix(xdat, ntry)
  xi=xifix
  
  par.start=rep(1,2)
  
  tryCatch(
    for(itry in 1:nrow(init.xifix)){
      
      par.start[1:2]=init.xifix[itry,1:2]
      
      value <- try(solnp(par.start, fun=gev.xilik, xifix=xifix, xdat=xdat,
                         LB =c(-Inf,0),UB =c(Inf,Inf), 
                         control=list(trace=0, outer.iter=10, 
                                      inner.iter=20, tol=1.e-5,
                                      delta=1.e-6, rho=1.2) ))
      
      if(is(value)[1]=="try-error"){
        kx[[itry]] <- list(value=10^6)
      }else{
        kx[[itry]] <- value
      }
      
    } #for
  ) #tryCatch
  
  optim_value  <-data.frame(num=1:ntry,value=sapply(kx, function(x) x$value[which.min(x$value)]))

  optim_table1 <-optim_value[order(optim_value$value),]
  selc_num  <- optim_table1[1,"num"]
  
  x  <- kx[[selc_num]]

  zx$conv <- x$convergence
  zx$nllh <- x$value[which.min(x$value)]
  
  if(zx$conv != 0){
     zx$mle = NA
  }else{
     zx$mle <- x$par
     zx$mle[3]= xifix
     zx$grad <- grad(gev.xilik, x$par, xifix=xifix, xdat=xdat)
#     zx$cov <- solve(x$hessian)
  }
  
  class(zx) <- "gev.xifix"
  return(zx)
  invisible(zx)
}
#---------- end of gev.xifix program -----------------
#----------------------------------------------------------
gev.xilik <- function(a, xifix=xifix, xdat=xdat) {
  
  mu <- a[1]      #mulink(mumat %*% (a[1:npmu]))
  sc <- a[2]      #siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
                  # xi <- shlink(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))
  
  xi=xifix
  if(sc <= 0) return(10^6)
  
  nsam= length(xdat)
  y= rep(0,nsam)
  
  y <- (xdat - mu)/sc
  y <- 1 - xi * y      # park modify to negative, for xi in hosking
  
  for (i in 1:nsam){
    y[i] = max(0, y[i], na.rm=T) }
  
  if (any(y <= 0)) return(10^6)
  
  if( abs(xi) >= 10^(-5) ) {ooxi= 1/xi
  }  else  {ooxi=sign(xi)*10^5}
  
  nllh = ( nsam*(log(sc)) + sum( exp(ooxi *log(y)) )
           + (1-(ooxi))* sum(log(y)) )
  
  zz=nllh 
  return(zz)
}
#------------------------------------------------------------
#-----------------------------------------------------------  
lme.boots <-function(data, B=NULL, quant){

  z=list()
  numq=length(quant)
  names_par= c("mu","sigma","shape")
  names_quant=paste("q", as.character(quant[1:numq]), sep="")
  sam.bootL=list()
  lme.boot =matrix(NA,nrow=B,ncol=5);  ratios =matrix(NA,nrow=B,ncol=5)
  lme.seboot =matrix(NA,nrow=B,ncol=numq); med = rep(NA, B)

  lm=lmoms(data, nmom=5)
  klme=pargev(lm, checklmom=T); 
  
  z$lme = klme$para[1:3]
  savel= vec2par(klme$para[1:3],'gev')
  z$qua.lme= quagev(quant[1:numq],savel)
  
  names(z$qua.lme)=names_quant;   names(z$lme)=names_par
  nsample =length(data)
  
  for (ib in 1:B){
    
    sam.bootL[[ib]]= sample(data, nsample,replace=T)   # nonpara Boot

    lm.boot= lmoms(sam.bootL[[ib]])

    lme.boot[ib,1:3]= pargev(lm.boot)$para[1:3] 
    
    ratios[ib,1]= lm.boot$lambdas[1]
    ratios[ib,2:5] = lm.boot$ratios[2:5]
 
    med[ib]= median(sam.bootL[[ib]])
    
    save=vec2par(lme.boot[ib,1:3],'gev')
    lme.seboot[ib,1:numq]=quagev(quant[1:numq],save)
  }

  z$B= B
  z$lme.boot= lme.boot
  z$ratios = ratios
  z$med = med
  z$quant = quant

  qua.lme.se=rep(NA, numq)
  for (qi in 1:numq){
    qua.lme.se[qi]=sqrt(var(lme.seboot[1:B,qi],na.rm=T))
  }
  z$qua.lme.se= qua.lme.se  
  
  return(z)
}
#-------------------------------------------------------