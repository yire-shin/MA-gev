# xdat=read.csv("C:/Users/user/Downloads/data_error.csv", header=T)
# as.vector(xdat)
# data= xdat$x
# xdat=data

  # # ----- min AD stat estimator -----------------------------
  # 
  #    xtry=3
  #    mup = abs(mean(data))*3; sigup = sqrt(var(data))*5
  #    init=matrix(NA,xtry,3)
  #    init[1,1:3]= hosking$lme
  #    init[2,1:3]= CD.mle
  #    init[3,1:3]= zx$remle
  #    for(i in 1:xtry){
  #      if(init[i,2] < 1.0) init[i,2] = 1.0
  #      if(abs(init[i,3]) < 0.01 ) init[i,3] = -0.01
  #    }
  # #   init[4,1:3]= delta$mle
  #    n=xtry; kstar=list(); k=list()
  # 
  #    for (itry in 1:xtry){
  #      work = try(mgedist(data, distr="gevad", gof="AD",
  #                    start= list(mu=init[itry,1], sigma=init[itry,2],
  #                                xi= init[itry,3]) )#,
  #                    # optim.method="L-BFGS-B",
  #                    # lower=c(-mup,0,-0.99), upper=c(mup, sigup,0.99))
  #      )
  # 
  #      if(is(work)[1]=="try-error"){
  #        k[[itry]] <- list(value=10^6)
  #      }else if(work$convergence != 0){
  #        k[[itry]] <- list(value=10^6)
  #      }else{
  #        k[[itry]] <- work
  #      }
  #    } #for
  # 
  #  optim_value  <-data.frame(num=1:n, value=sapply(k, function(x) x$value[which.min(x$value)]))
  #  optim_table1 <-optim_value[order(optim_value$value),]
  #  selc_num  <- optim_table1[1,"num"]
  #  kstar  <-k[[selc_num]]
  # 
  #  # z$conv <- x$convergence
  #   zx$minAD.value <-  kstar$value
  #   zx$minAD <- kstar$estimate    # hosking style xi
  #   if(abs(zx$minAD[3]) >= 1.0 | zx$minAD[2] <= 0) {
  #     zx$minAD = NA
  #     zx$qua.minAD= NA
  #   }else{
  #     zx$qua.minAD = quagev(quant,vec2par(zx$minAD,'gev'))
  #   }

  # 
  # #-------For minAD of gevd -----------------------------
  # dgevad= function(x, mu, sigma, xi){
  #   if(sigma <= 0) sigma= 0.01
  #   if(xi <= -1) xi= -0.98
  #   if(xi >= 1) xi= 0.98
  #   para1=vec2par(c(mu,sigma,xi),'gev')
  #   pdfgev(x,para1)
  # }
  # pgevad= function(q, mu, sigma, xi){
  #   if(sigma <= 0) {sigma= 0.01; return(5)}
  #   if(xi <= -1) {xi= -0.98; return(5)}
  #   if(xi >= 1) {xi= 0.98; return(5)}
  #   para2=vec2par(c(mu,sigma,xi),'gev')
  #   cdfgev(q,para2)
  # }
  # qgevad= function(p, mu, sigma, xi){
  #   if(sigma <= 0) sigma= 0.01
  #   if(xi <= -1) xi= -0.98
  #   if(xi >= 1) xi= 0.98
  #   para3=vec2par(c(mu,sigma,xi),'gev')
  #   quagev(p,para3)
  # }
  #----------------------------------------------------------
  #  quant=c(.99,.995)
   # gev.remle(xdat=xdat, ntry=5, rest='mean', quant, trim=0,
   #           CD.mle=NULL, mle=NULL, second=T, w.mpse=F) 
#------------------------------------------------------------
gev.remle =function (xdat=NULL, ntry=5, rest='mean', quant, trim=0,
                     CD.mle=NULL, mle=NULL, second=TRUE, w.mpse=FALSE) 
{
  z <- list();  k =list()
  data=xdat; # re=sort(data)
  if(rest=='med') rest='median'
  
  nsample=length(xdat)
  z$nsample=nsample
  
  if(is.null(CD.mle)){
    CD.mle= rep(NA,3)
    CD.mle= gev1.CD.hos(xdat=data, ntry=5)$mle  # hosking style xi
  }
  if(is.null(mle)){
    mle= rep(NA,3)
    mle= gev.max(xdat=data, ntry=5)$mle  # hosking style xi
  }
  
  if(ntry < 4) ntry =4
  init= matrix(NA, ntry,3)
  init <- ginit.max(data=xdat, ntry= ntry-2)
  init = rbind(init, matrix(NA,2,3))
  init[ntry-1,1:3]= mle[1:3]
  init[ntry,1:3]= CD.mle[1:3]
  
  init=na.omit(init);  ntry=nrow(init)
  
  upsig= sqrt(var(xdat))*3
  upmu=  abs(mean(xdat))*3
  
  ntrim=trim
  med.dat=median(data); mean.dat =mean(data)
  lmx=lmoms(data,nmom=3)

  nar=c('rest.Ex','rest.med')  #,'rest.tau3')
  if(rest=='mean'){sel=1
  }else if(rest=='median'){sel=2}
#  }else if(rest=='tau3'){sel=3}

  tryCatch(
    for(i in 1:ntry){
      
      work=list()
      work <- try(solnp(init[i,], fun=gev.lik.remax, 
                         eqfun= match.fun(nar[sel]), eqB=0,
                         LB =c(-upmu,0,-0.99),UB =c(upmu,upsig,0.99), 
                         control=list(trace=0, outer.iter=20,
                         delta=1.e-6, inner.iter=200, tol=1.e-5)
                        , xdat=xdat, lmx=lmx, med.dat=med.dat,
                        mean.dat=mean.dat) 
                  )
      
      if(is(work)[1]=="try-error"){
        k[[i]] <- list(values=10^6)
      }else if(work$convergence != 0){
        k[[i]] <- list(values=10^6)
      }else{
        k[[i]] <- work
      }
      
    } #for
  ) #tryCatch
  
  optim_value  <-data.frame(num=1:ntry,values=sapply(k, function(x) x$values[which.min(x$values)]))
  
  optim_table1 <-optim_value[order(optim_value$values),]
  selc_num  <- optim_table1[1,"num"]
  x= list()
  
  x  <- k[[selc_num]]

#  z$conv <- x$convergence
  z$remle1.value <- x$values[which.min(x$values)]
  z$remle1 <- x$pars

  if( is.null(x$pars) | any(is.na(x$pars)) | z$remle1.value==10^6) {
    x$pars= NA; z$remle1 <- NA; z$qua.remle1 = NA
  }else if( x$pars[2] <= 0 | abs(x$pars[3]) >= 1.0){
    x$pars= NA; z$remle1 <- NA; z$qua.remle1 = NA
  }else{ z$qua.remle1 = quagev(quant,vec2par(x$pars,'gev')) }
  
  z$rest.method= rest
  
  if(second==TRUE){
    
    k2=list()
    nar2=c('rest.Ex2','rest.med2') #,'rest.tau3')

    init2 = rbind(init, rep(NA,3))
    init2[nrow(init2),1:3]= z$remle1[1:3]
    init2=na.omit(init2); xtry= nrow(init2)

    tryCatch(
      for(i in 1:xtry){
        
        value2=list()
        value2 <- try(solnp(init2[i,], fun=gev.lik.remax, 
                           eqfun= match.fun(nar2[sel]), eqB=c(0,0),
                           LB =c(-upmu,0,-0.99),UB =c(upmu,upsig,0.99), 
                           control=list(trace=0, outer.iter=20,
                                   delta=1.e-6, inner.iter=200, tol=1.e-5)
                           , xdat=xdat, lmx=lmx, med.dat=med.dat,
                           mean.dat=mean.dat)
                   )
        
        if(is(value2)[1]=="try-error"){
          k2[[i]] <- list(values=10^6)
        }else if(value2$convergence != 0){
          k2[[i]] <- list(values=10^6)
        }else{
          k2[[i]] <- value2
        }
        
      } #for
    ) #tryCatch
    
    optim_value  <-data.frame(num=1:xtry,values=sapply(k2, function(x) x$values[which.min(x$values)]))
    
    optim_table1 <-optim_value[order(optim_value$values),]
    selc_num  <- optim_table1[1,"num"]
    xw=list()
    
    xw  <- k2[[selc_num]]
    
    #  z$conv <- x$convergence
    z$remle2.value <- xw$values[which.min(xw$values)]
    z$remle2 <- xw$pars
    
    if( is.null(xw$pars) | any(is.na(xw$pars)) | z$remle2.value==10^6) {
      xw$pars= NA; z$remle2 <- NA; z$qua.remle2 = NA
    }else if(xw$pars[2] <= 0 | abs(xw$pars[3]) >= 1.0){
      xw$pars= NA; z$remle2 <- NA; z$qua.remle2 = NA
    }else{ z$qua.remle2 = quagev(quant,vec2par(xw$pars,'gev')) }

  } # end if second

  #------ mpse ----------------------------------------------   
  if(w.mpse==TRUE){

    xtry=3
    mup = abs(mean(data))*3; sigup = sqrt(var(data))*3
    lm=lmoms(data,nmom=3)
    lme=pargev(lm)$para
  
  if(is.null(CD.mle)){
    CD.mle= rep(NA,3)
    CD.mle= gev1.CD.hos(xdat=data, ntry=5)$mle  # hosking style xi
  }
  
  init=matrix(NA,xtry,3)
  init[1,1:3]= lme
  init[2,1:3]= CD.mle
  init[3,1:3]= z$remle1
  init=na.omit(init);  xtry=nrow(init)
  
  for(i in 1:xtry){
    if(init[i,2] <= 0.0) init[i,2] = 0.1
    if(abs(init[i,3]) < 0.01 ) init[i,3] = -0.01
  }
  #   init[4,1:3]= delta$mle
  n=xtry; kstar=list(); k3=list()
  
  for (itry in 1:xtry){
    
    work2=list()
    work2 = try(mpsedist(data=xdat, distr="gevh",
                       start= list(mu=init[itry,1], sigma=init[itry,2], 
                                   xi= init[itry,3])) #,
                     #  optim.method="L-BFGS-B",
                     #   lower=c(-mup,0,-0.99), upper=c(mup, sigup,0.99)) 
    )
    
    if(is(work2)[1]=="try-error"){
      k3[[itry]] <- list(value=10^6)
    }else if(work2$convergence != 0){
      k3[[itry]] <- list(value=10^6)
    }else{
      k3[[itry]] <- work2
    }
  } #for
  
    optim_value  <-data.frame(num=1:n, value=sapply(k3, function(x) x$value[which.min(x$value)]))
    optim_table1 <-optim_value[order(optim_value$value),]
    selc_num  <- optim_table1[1,"num"]
  
    kstar  <- k3[[selc_num]]

    z$mpse.value= kstar$value
    z$mpse = kstar$estimate
    if(abs(z$mpse[3]) >= 1.0 | z$mpse[2] <= 0) {
      z$mpse= NA
      z$qua.mpse =NA
    }else{
      z$qua.mpse =quagev(quant,vec2par(z$mpse,'gev'))
    }
  
  }else if(w.mpse==F){
    z$mpse= NA
    z$qua.mpse =NA
  } # end if w.mpse
  
  return(z)
}
#----------------------------------------------------------------------
  #-------------------------------------------------  
  gev.lik.remax <- function(a, xdat=xdat, lmx=lmx, med.dat=med.dat,
                            mean.dat=mean.dat) {
    
    mu <- a[1]      #mulink(mumat %*% (a[1:npmu]))
    sc <- a[2]      #siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
    xi <- a[3]      #shlink(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))
    
    y <- (xdat - mu)/sc
    y <- 1 - xi * y        # park modify to negative, for xi in hosking
    
    nsample=length(xdat)
    for (i in 1:nsample){
      y[i] = max(0, y[i], na.rm=T) }
    
    if (any(y <= 0) || any(sc <= 0)) 
      return(10^6)
    
    if( abs(xi) >= 10^(-5) ) {ooxi= 1/xi
    }  else  {ooxi=sign(xi)*10^5}
    
    zz=nsample*(log(sc)) + sum( exp(ooxi *log(y)) ) + sum(log(y) * (1-(ooxi)) ) 
    
    return(zz)
  }
  #----------------------------------------------------------------------- 
  #-----------------------------------------------
  rest.Ex = function(a, xdat=xdat, lmx=lmx, med.dat=med.dat,
                     mean.dat=mean.dat){
    lam1= lmomgev(para=vec2par(a,'gev'))$lambdas[1]
    return(lam1-mean.dat )
  }
  rest.med =function(a,  xdat=xdat, lmx=lmx, med.dat=med.dat,
                     mean.dat=mean.dat){
    med= quagev(0.5,vec2par(a,'gev'))
    return(med-med.dat )
  }
  # rest.tau3 =function(a){
  #   tau3 = lmomgev(para=vec2par(a,'gev'))$ratios[3]
  #   t3=lmoms(data)$ratios[3]
  #   return(tau3-t3)
  # }
  #-------------------------------------------------------------
  
  #-----------------------------------------------
  rest.Ex2 = function(a,  xdat=xdat, lmx=lmx, med.dat=med.dat,
                      mean.dat=mean.dat){
    lamp= lmomgev(para=vec2par(a,'gev'))$lambdas[1:2]
    return( c(lamp[1]-lmx$lambdas[1],lamp[2]-lmx$lambdas[2]) )
  }
  rest.med2 =function(a,  xdat=xdat, lmx=lmx, med.dat=med.dat,
                      mean.dat=mean.dat){
    medx= quagev(0.5,vec2par(a,'gev'))
    lam2= lmomgev(para=vec2par(a,'gev'))$lambdas[2]
    return( c(medx-med.dat,lam2-lmx$lambdas[2]) )
  }
  #---------------------------------------------------
  
  #-------For mpse of gevd -----------------------------
  dgevh= function(x, mu, sigma, xi){
    if(sigma <= 0) {sigma= 0.01; return(500)}
    if(xi <= -1) {xi= -0.98; return(500)}
    if(xi >= 1) {xi= 0.98; return(500)}
    para1=vec2par(c(mu,sigma,xi),'gev')
    pdfgev(x,para1)
  }
  pgevh= function(q, mu, sigma, xi){
    if(sigma <= 0) {sigma= 0.01; return(5)}
    if(xi <= -1) {xi= -0.98; return(5)}
    if(xi >= 1) {xi= 0.98; return(5)}
    para2=vec2par(c(mu,sigma,xi),'gev')
    cdfgev(q,para2)
  }
  qgevh= function(p, mu, sigma, xi){
    if(sigma <= 0) sigma= 0.01
    if(xi <= -1) xi= -0.98
    if(xi >= 1) xi= 0.98
    para3=vec2par(c(mu,sigma,xi),'gev')
    quagev(p,para3)
  }
  #-----------------------------------------
