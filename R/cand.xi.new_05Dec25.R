#------------------------------------------------------------
new.kpar = function(wtgd, numk, kpar, final){

wtgd[which(is.na(wtgd))] =0
zp=list()
aw=0
dist=0.03

nid = which(wtgd < 0.01); id = which(wtgd >= 0.01)
kpar = kpar[id]
wtgd = wtgd[id]
numk2 = length(kpar); if(numk2 < numk) aw=10 ; numk = numk2
kpar2= kpar

if(numk==0) {
  kpar2= -0.1
  numk=1
}

if(numk==1){
  ad=1.5
  kpar2= c(NA,NA, kpar2, NA,NA)
  kpar2[2]= max(kpar2[3] -dist/ad,   -0.99)
  kpar2[1]= max(kpar2[3] -2*dist/ad, -0.99)
  kpar2[4]= min(kpar2[3] +dist/ad,   0.49)
  kpar2[5]= min(kpar2[3] +2*dist/ad, 0.49)
  aw=11
  
}else{

tre1=0.1; tre8=0.8

work1=0; work2=0

if(wtgd[1] >= tre1 & wtgd[1] <= tre8) {
  kpar2= c(rep(NA,2), kpar)
  kpar2[2] = max(kpar[1]-dist, -0.99)
  kpar2[1] = max(kpar2[2]-dist, -0.99)
  work1=1
}else if(wtgd[1] > tre8) {
  kpar2= c(rep(NA,4), kpar)
  kpar2[4] = max(kpar[1]-dist/2, -0.99)
  kpar2[3] = max(kpar2[4]-dist/2, -0.99)
  kpar2[2] = max(kpar2[3]-dist/2, -0.99)
  kpar2[1] = max(kpar2[2]-dist/2, -0.99)
  work1=2
}

if(work1==0 & wtgd[numk] >= tre1 & wtgd[numk] <= tre8){
  kpar2= c(kpar, rep(NA,2))
  kpar2[numk+1]= min(kpar2[numk]+dist, 0.49)
  kpar2[numk+2]= min(kpar2[numk+1]+dist, 0.49)
  work2=1
}
if(work1==0 & wtgd[numk] > tre8){
  kpar2= c(kpar, rep(NA,4))
  kpar2[numk+1]= min(kpar2[numk]+dist/2, 0.49)
  kpar2[numk+2]= min(kpar2[numk+1]+dist/2, 0.49)
  kpar2[numk+3]= min(kpar2[numk+2]+dist/2, 0.49)
  kpar2[numk+4]= min(kpar2[numk+3]+dist/2, 0.49)
  work2=2
}

if(work1==1 & wtgd[numk] >= tre1){
  kpar2= c(kpar2,rep(NA,2))
  kpar2[numk+3]= min(kpar2[numk+2]+dist, 0.49)
  kpar2[numk+4]= min(kpar2[numk+3]+dist, 0.49)
  work2=2
}

aw= aw+ work1+work2

} #end if numk =1

numk=length(kpar2)

zp$kpar2=kpar2
zp$aw=aw
zp$numk=numk
return(zp)
}
#------------------------------------------------------------
cand.xi.new <- function(mle, pick=.95, nint=128, hosking=NULL, 
                        cov.lme, figure=T){
  
  wx = list()
  SMALL= 1e-05
  numk= hosking$numk
  start = hosking$start
  kpar=rep(NA, numk)
  data = hosking$data
  move=F
  
  while(start=='mle' | start=='mix'){
    
    newse= max(0.3,abs(mle[3]))*0.1
    
    mle.coles= mle
    mle.coles[3]= -mle[3]
    
    xlow= mle.coles[3]-8*newse; xup=  mle.coles[3]+8*newse
    xlow= max(-0.99, min(0.5,xlow), na.rm=T)
    xup=  min(0.99, max(-0.5,xup), na.rm=T)
    xlow2 = min(xlow, xup); xup2 = max(xlow, xup)
    
    get.ci = gev.profxi.mdfy.new(data, mle=mle.coles, 
                                 xlow=xlow2, xup= xup2,
                                 pick=pick, nint=nint, figure=figure)
    if(get.ci$fail==T) {
      start='lme'; move=T
      break
    }
    
    wx$get.ci = get.ci
    
    lo.lim =min(-get.ci$w1$ci1, -get.ci$w1$ci2)
    up.lim =max(-get.ci$w1$ci1, -get.ci$w1$ci2)
    
    kpar[1] = max(lo.lim, -0.99)
    kpar[numk] = min(up.lim, 0.99)

    if(is.odd(numk)) kpar[floor(numk/2)+1] = -get.ci$w1$xmax
    
    pilo= (1-pick)/2; piup =1-pilo
    eqpr = seq(pilo, piup, by=pick/(numk-1) )
    numpi= floor(numk/2)
    
    for (i in 2:numpi){
      
      pipick = 1- 2*eqpr[i] 
      
      get.ci = gev.profxi.mdfy.new(data, mle=mle.coles, 
                                   xlow=xlow2, xup= xup2,
                                   pick=pipick, nint=nint, figure=figure)
      if(get.ci$fail==T) {
        start='lme'; move=T
        break
      }
      
      lo.lim =min(-get.ci$w1$ci1, -get.ci$w1$ci2)
      up.lim =max(-get.ci$w1$ci1, -get.ci$w1$ci2)
      kpar[i]= max(lo.lim, -0.99)
      kpar[numk-i+1]= min(up.lim, 0.99)
    } #end for
    
    if(kpar[numk] < -0.3) {
      kpar2=rep(NA,numk+2)
      kpar2[numk+1]= max(-0.3, kpar[numk]+0.05)
      kpar2[numk+2]= max(-0.3, kpar2[numk+1]+0.05)
      kpar[1:numk] =c(kpar[3:numk], kpar2[numk+1], kpar2[numk+2])
    }
    
    wx$kpar.mle = kpar
    break
  } # end while
  
  #-------------------------------------------------------
  if(start=='lme' | start=='mix') {
    
    # lo.lim= quantile( hosking$lme.boot[,3], prob=(1-pick)/2, na.rm=T)
    # up.lim= quantile( hosking$lme.boot[,3], prob=1-(1-pick)/2, na.rm=T)
    # lo.lim= max(-.9999, lo.lim)
    
    pilo= (1-pick)/2; piup =1-pilo
    eqpr = seq(pilo, piup, by=pick/(numk-1) )
    
    if(!is.null(cov.lme)){
      hosking$lme.boot = cov.lme$lme.boot
      
    }else{
      
      if(hosking$boot.lme ==F & move==T){
        Bnew=hosking$B; qua= hosking$quant
        
        hosking = lme.boots.new(data=data, B=Bnew, quant=qua,
                                boot=T)
      } 
    } #end if

    kpar = quantile( hosking$lme.boot[,3], prob=eqpr, na.rm=T)
    
    if(kpar[numk] < -0.5) kpar[numk]= max(-0.5, kpar[numk]+0.05)
    wx$kpar.lme = kpar
    if(move==T) wx$kpar.mle = kpar
    hosking$start=start
  }
  
  if(start=='mix'){
    wx$kpar.mix = (wx$kpar.mle + wx$kpar.lme)/2
  }

  if(start=='mle') {kpar = wx$kpar.mle
  }else if(start=='lme') {kpar = wx$kpar.lme
  }else if(start=='mix') {kpar = wx$kpar.mix
  }
  if( any(abs(kpar) <= SMALL)) {
    kpar[ which(abs(kpar) <= SMALL) ] = SMALL*1.2 }
  
  if(kpar[numk] >= 0.5) kpar[numk]= 0.49
  wx$kpar = kpar
  wx$start = start
  if(start=='mle') wx$ymin = get.ci$ymin
  
  return(wx)
}
#---------------------------------------------------------
gev.profxi.mdfy.new=  function (data=NULL, mle=NULL, xlow, xup, 
                                pick = NULL, nint = 128, figure=F){
  
  w <- list();  v <- numeric(nint)
  
  w$fail=F; w$start='mle'
  
  conf = pick
  xlow.t= min(xlow, xup);    xup.t = max(xlow, xup)
  xlow= max(-0.9, xlow.t, na.rm=T);  
  xup = min(0.9, xup.t, na.rm=T)
  x <- seq(xlow, xup, length.out = nint)
  sol <- c(mle[1], mle[2])      # mle = mle.coles
  
  #-------------------------------------------------------------   
  gev.plikxi <- function(a, data=data, xi=xi) {
    
    if (a[2] <= 0) {
      l <- 10^6
      return(l)
    }
    if (abs(xi) < 10^(-6)) {
      y <- (data - a[1])/a[2]
      l <- length(y) * log(a[2]) + sum(exp(-y)) + 
        sum(y)
    } else {
      y <- (data - a[1])/a[2]
      y <- 1 + xi * y
      if (any(y <= 0)) { l <- 10^6
      }else{ l <- length(y) * log(a[2]) + sum(y^(-1/xi)) + 
        sum(log(y)) * (1/xi + 1)
      }
      l
    } #end if
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
 
    irt=1
  }else if(idmax==length(d$x)){
    
    irt=1               # irt=1 means the max is observed at the first or the last x
  }
  
  if(irt==1) {
    w$start='lme'
    w$fail=T
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
  
  inter=1
  
  #--- for linear interpolation ----------------------------------
  
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
  
  w$w1 = w1 
  inter = 1
  d= msp1
  d$v = msp1$y
  
  if(figure == T & pick >= .90){
    
    par(new=F)
#    par(mfrow=c(1,1))
    
    halfchi=  0.5 * qchisq(conf, 1)
    # ci_max = w$ci2
    # ci_min = w$ci1
    
    xlim0 = -(w$w1$ci2 + (w$w1$ci_length)/4 ) 
    xlim1 = -(w$w1$ci1 - (w$w1$ci_length)/4 )
    
    if(inter==1) {co="green"} else if(inter==2) {co="blue"}
    
   # --- temporary no action --------------------------------------
    plot(-d$x, d$v,
         xlim=c(xlim0-0.02, xlim1+0.02),
         ylim=c( w$w1$vmax-(halfchi*2.7), w$w1$vmax+ 0.8*halfchi),
         type = "l", xlab = "xi_Hosking in GEV",
         col= co, lwd= 2,
         ylab = "Profile Log-likelihood",main=c("Profile CI for xi"))
    ma <-   w$w1$vmax  #-z$nllh
    abline(h = ma, col = 4)
    abline(h = ma - 0.5 * qchisq(conf, 1), col = 4)
    abline(v=c(-w$w1$xmax, -w$w1$ci1, -w$w1$ci2), col=2, lty=2, lwd= 2)
    #    abline(v = c(-mle[3]), lty=3, col="black", lwd=2)
    text( x=c(-mle[3]), y= w$w1$vmax-(halfchi*2.7)+0.2, labels="mle")

    lm= lmoms(data)
    lme= pargev(lm)$para[3]

    abline( v = lme, lty=5, lwd=2, col="purple")
    text( x =lme, y= w$w1$vmax-(halfchi*2.7)+0.2, labels="lme")
    
    w$ymin = w$w1$vmax-(halfchi*2.7)+0.2
  }
  return(w)
}
#--------------------------------------------------------
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
  # --------------------------------------------------------------------