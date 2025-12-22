
library(ismev)
library(Rsolnp)
library(lmomco)
library(zoo)

#' @name ma.gev.open
#' @aliases ma.gev.open
#' @title model averaging with mixed criterion for estimating high quantile in the GEVD.
#' @description
#' The method in this code combines MLE and LME to construct candidate submodels and assign weights effectively, applied to generalized extreme value distribution.
#' @param data A numeric vector of data to be fitted.
#' @param quant The probabilities corresponding to high quantiles to be estimated. For example, quant=c(.90,.95,.98,.99,.995,.998)
#' @param weight The weighting method name; choice of weight = 'gLd', 'med', 'like', and 'cvt'(=conventional), and 'gLd1', 'gLd2', 'like0'(default), 'like1', 'med1', 'med2'
#' @param numk The number of submodels, K (default=12)
#' @param B The number of bootstrap samples (default=200).
#' @param varcom Whether or not to compute the variance of quantile estimation; It has value of T (default) or F.
#' @param trim The number of left trimming, usually 0 (default), 1, or 2.
#' @param start The method for picking initial K shape parameters (kpar) of submodels; choice of start= 'mle' (default) or 'lme'. 
#' @param fig TRUE or FALSE (default), for drawing a profile likelihood plot of the shape parameter
#' @param bma TRUE or FALSE (default), for the BMA procedure or not.
#' @param pen The penalty name; pen='norm' (default) or 'beta'.
#' @param CD TRUE or FALSE (default), for Coles and Dixon maximum penalized likelihood estimation.
#' @param remle TRUE or FALSE (default), for the restricted MLE (or mixed estimation)
#' @param cov.lme Only useful for simulation study. Default is NULL.
#' 
#' @return {
#'  A list containing the following components. Some of these components are printed after the fit, depending on the input options.
#' 
#'  \item{mle.hosking}{A vector containing the mle in hosking style (different sign for the shape parameter from Coles book).}
#'  \item{qua.mle}{The quantile estimates from mle.}
#'  \item{mle.cov3}{The covariance matrix of the mles of three parameters, obtained from the Fisher info matrix.}
#'  \item{qua.se.mle.delta}{A vector containing the standard errors of quantile estimate obtained from mle, using the delta method.}
#'  \item{lme}{A vector containing the lme.}
#'  \item{lme.cov3}{The covariance matrix of the lmes of three parameters, obtained by bootstrap.}
#'  \item{qua.lme}{The quantile estimate from the lme}
#'  \item{qua.se.lme.boots}{A vector containing the standard errors of quantile estimate obtained from lme, using bootstrap.}
#'  \item{leftrim}{Number of left trim.}
#'  \item{data}{The data that has been fitted. For non-stationary models, the data is standardized.}
#'  \item{fin.se.ma}{A vector containing the asymptotic SE of quantile estimates under fixed weights.}
#'  \item{adj.se.ma}{A vector containing the asymptotic SE of quantile estimates under random weights.}
#'  \item{pred.se.bma}{A vector containing the predictive SE of quantile estimates under BMA.}
#'  \item{surr$par}{Parameters for surrogate model.}
#'  \item{zp}{A matrix containing the quantile estimates from each submodel of K.}
#'  \item{zp.ma}{The quantile estimates from MA of K submodels.}
#'  \item{bma}{The quantile estimates from BMA of K submodels.}
#'  \item{w.ma}{A vector containing the weights used for MA.}
#'  \item{prior}{A vector containing hyperparameters of a data-adaptive prior.}
#'  \item{pick_xi}{A vector containing xi_k for K submodels, selected based on the 95% confidence interval of xi.}
#' }
#' @seealso \code{\link{lmoms}}
#' @export
#' 
#' @examples
#' \dontrun{
#' data(haenam)
#' ma.gev.open(data=haenam, quant=c(.98, .99, .995), weight='like', bma=TRUE)
#' }


#----------------------------------------------------------------------------
ma.gev.open= function(data=NULL, quant=c(.98,.99,.995), weight='like1',
                      numk=12, B=200, varcom=TRUE, trim=0, start= "mle",
                      fig=FALSE, bma=FALSE, pen="norm", CD=FALSE,
                      remle=FALSE, cov.lme=NULL) 
{
  
  # start.est = 'mle' (default) or 'lme'
  # pick: CI level between .90 or .95 (default)
  # choice of weight = 'gLd', 'med', 'like', and cvt' (=conventional mle and aic weight)
  #                    'gLd0', 'gLd1', 'gLd2', 'like0', 'like1', 'med1', 'med2'
  
  zx=list()
  SMALL <- 1e-05
  numq=length(quant)
  zx$weight = weight
  zx$start = start
  zx$quant=quant
  
  # defaults 
  pick=0.95; pertr=0.85; surr=TRUE; fcv.work=FALSE; 
  try.last=FALSE;  numom=3; order=3; 
  type='full'
  seek=FALSE; cov.type='lambda';  fix.kpar=FALSE; pick.xi=NULL; 
  
  if(!is.null(cov.lme)) {cov.lme.from.out = TRUE
  }else{cov.lme.from.out = FALSE }
  
  if(fix.kpar==FALSE) {pick.xi=NULL
  }else if(fix.kpar==TRUE) {numk= length(pick.xi)}
  
#  if(surr==FALSE) fig=FALSE
#  if(fig==TRUE) {surr=TRUE; varcom=TRUE}  #! boot.lme=TRUE;
  if(weight=='gld') weight='gLd'
  if(cov.type=='r') cov.type='ratio'
  if(cov.type=='l') cov.type='lambda'
  if(weight=='like0'){weight='like'; trim=0}
  if(weight=='like1'){weight='like'; trim=1}
  if(weight=='gld0' | weight=='gLd0'){weight='gLd'; trim=0}
  if(weight=='gld1' | weight=='gLd1'){weight='gLd'; trim=1}
  if(weight=='gld2' | weight=='gLd2'){weight='gLd'; trim=2}
  if(weight=='med1'){weight='med'; trim=1}
  if(weight=='med2'){weight='med'; trim=2}
  
  zx$data= data
  org.numk=numk
  
  if(fcv.work==TRUE) {    # using for the forward cross-validation for upper 100*pertr%
    ccv = cons.dte(data,pertr, type='full')
    data= ccv$dtr
    nsam= length(data)
  }

  para3= matrix(NA,nrow=numk, ncol=3)
  wtgd= rep(NA, numk)
  zp= matrix(NA, numq, numk)   
  zpf=rep(NA, numq)

  nsample=length(data)
  delta=list()
  
  # ------- mle and se by delta method-------------------
  
  delta= gev.rl.delta_new(data=data, ntry=10, quant=quant)
  
  zx$mle.hosking =  delta$mle
  zx$qua.mle= delta$qua.mle
  zx$mle.cov3 = delta$cov
  #  zx$mle.aic = delta$nllh*2 + 2*3
  
  if(varcom==TRUE & !is.null(delta$cov) ){
    qua.se.mle=rep(NA, numq)
    for(i in 1:numq){
      qua.se.mle[i]= delta.gev(gevf3= delta,quant=quant[i],
                               d3yes=TRUE)$v3
    }
    zx$qua.se.mle.delta = sqrt(qua.se.mle)
  }
  
  #-----  Lme and se by bootstrap --------------------
#  if(start == 'mle' | weight=="like") boot.lme= FALSE
  
  if(is.null(cov.lme)){
    if(start=='lme' | weight=='med' | start=='mix' | weight=='gLd'
       | weight=='lcv' | weight=='opt') boot.lme=TRUE
  }else if(!is.null(cov.lme)) boot.lme=FALSE
  
  if(weight=="like" & start =="mle") boot.lme= FALSE
  
  hosking=list()
  hosking = lme.boots.new(data=data, B=B, quant=quant, 
                          boot=boot.lme, trim=trim)
  
  if(cov.lme.from.out != TRUE) {zx$lme.cov3 = hosking$cov.par 
  }else if(cov.lme.from.out== TRUE) {zx$lme.cov3 = cov.lme$cov.par}
  
  zx$lme = hosking$lme
  zx$qua.lme= hosking$qua.lme; zx$leftrim=trim
  hosking$mle = delta$mle
  
  xi.hat = (delta$mle[3] + hosking$lme[3])/2
  #  if(xi.hat > 0.0) {cat("no need to do MA computing, xi.hat=",
  #                        xi.hat,"\n")}
  
# # ---- mple CD, re.mle, mpse  ------------------------------------------
  if(CD==TRUE){
    CD.mle= rep(NA,3)
    CD.mle= gev1.CD.hos(xdat=data, ntry=5)$mle   # hosking style xi

    if(any(is.na(CD.mle)) ) {
      zx$qua.CD=NA; zx$mle.CD =NA; CD.mle=NULL
    }else{
      zx$qua.CD = quagev(quant, vec2par(CD.mle,'gev'))
      zx$mle.CD = CD.mle
    }
  } #end if CD

  w.mpse=FALSE
  if(remle==TRUE){
    if(CD != TRUE) CD.mle= hosking$lme
    rem.try=list()
    rem.try = gev.remle(xdat=data, ntry=5, rest='mean', quant=quant,
                        trim=0, CD.mle=CD.mle, mle=delta$mle,
                        second=TRUE, w.mpse=w.mpse)

    zx$remle1 = rem.try$remle1
    zx$qua.remle1 = rem.try$qua.remle1
    zx$remle2 = rem.try$remle2
    zx$qua.remle2 = rem.try$qua.remle2
    zx$qua.rem.ave = (zx$qua.remle1 + zx$qua.remle2)/2
  }

    if(w.mpse==TRUE){
     zx$mpse = rem.try$mpse
     zx$qua.mpse = rem.try$qua.mpse
    }
# # --------------------------------------------------------------
  
  hosking$start = start; hosking$numk = numk; hosking$quant=quant
  hosking$weight = weight; hosking$B= B; hosking$trim=trim
  hosking$data = data; hosking$boot.lme = boot.lme
  
  #---  xi_k picking for cand submodels ---------------------------
  kpar=rep(NA, numk)
  
  if(fix.kpar==FALSE){
    kpar.list=list()
    
    kpar.list = cand.xi.new.paper(data, hosking=hosking, mle=zx$mle.hosking, 
                                  pick0=pick, nint=256, start=start, 
                                  numk=numk, figure=fig, bma=bma, pen=pen)

    kpar = kpar.list$kpar
    hosking$start= kpar.list$start
    
  }else if(fix.kpar==TRUE){
    kpar=pick.xi
    hosking$start=start
  }
  hosking$kpar = kpar; hosking$cov.type= cov.type
  
  #------- Zp and weights computing-----------------------------------
  run=1;

  xqa=c(0.5, .65, .8, .85, .9, .925, .95, .965,.98, 
        .985, .99, .9925, .995, .9965, .998, .999)
  mywt=list()
  hosking$pertr = pertr
  bmaw=rep(0, numk); wtgd=rep(NA, numk); #bmaw2=rep(NA, numk);
  
  mywt = weight.com.new(data, numk=numk, hosking=hosking, 
                        kpar=kpar, numom=numom,
                        xqa=xqa, varcom=varcom,  
                        cov.lme=cov.lme, surr=surr,
                        type=type, trim=trim, cov.type=cov.type,
                        bma=bma, pen=pen)
  
  para3 = mywt$prob.call$mle3; 
  npara3= na.omit(para3); numk=nrow(npara3); 
  if(mywt$numk != numk) print("trouble in numk after 1st mywt")
  notid= which(!is.na(para3[,1]))
  mywt$data = data; para3=npara3; mywt$prob.call$mle3 = npara3
  kpar= kpar[notid]; 
  
  zp = mywt$zp[,notid];  xzp = mywt$xzp[,notid]
  wtgd = mywt$wtgd[notid]; mywt$wtgd = wtgd
  bmaw = mywt$bmaw[notid]
  
  if(bma != TRUE){
    newk= new.kpar2(wtgd=wtgd, numk=numk, kpar=kpar)
  }else if(bma==TRUE){
    newk= new.kpar2(wtgd=bmaw, numk=numk, kpar=kpar)
  }
  
  numk= newk$numk
  aw= newk$aw
    
  while(aw >= 1){
    
    run=run+1
    kpar=newk$kpar2
    hosking$kpar = newk$kpar2; hosking$numk = numk
    
    mywt=list()
    bmaw=rep(0, numk); wtgd=rep(NA, numk)
    para3= matrix(NA, nrow=numk, ncol=3)
    zp= matrix(NA, numq, numk)
    
    mywt = weight.com.new(data, numk=numk, hosking=hosking, 
                          kpar=newk$kpar2, numom=numom,
                          xqa=xqa, varcom=varcom,  
                          cov.lme=cov.lme, surr=surr,
                          type=type, trim=trim, cov.type=cov.type,
                          bma=bma, pen=pen)
    
    para3 = mywt$prob.call$mle3
    npara3= na.omit(para3); numk=nrow(npara3);
    if(mywt$numk != numk) print("trouble in numk after 2nd mywt")
    notid= which(!is.na(para3[,1]))
    mywt$data = data; para3=npara3; mywt$prob.call$mle3 = npara3
    kpar= mywt$kpar
    kpar= kpar[notid] 
    
    zp = mywt$zp[,notid]
    xzp = mywt$xzp[,notid]
    wtgd = mywt$wtgd[notid]; mywt$wtgd = wtgd
    bmaw = mywt$bmaw[notid]; bmaw2 = mywt$bmaw2[notid]
    
    if(run >= 5) break
    mxid = which.max(wtgd)
    bmxid = which.max(bmaw)
    
    tre1=0.1; tre8=0.5
    
    if(bma != TRUE){
     if(wtgd[1] > tre1 | wtgd[numk] > tre1 | wtgd[mxid] > tre8){
      newk= new.kpar2(wtgd=wtgd, numk=numk, kpar=kpar)
      numk= newk$numk
      aw= newk$aw
     }else{ aw=0; break}
    }
    
    if(bma==TRUE){  
     if(bmaw[1] > tre1 | bmaw[numk] > tre1 | bmaw[bmxid] > tre8){
      newk= new.kpar2(wtgd=bmaw, numk=numk, kpar=kpar)
      numk= newk$numk
      aw= newk$aw
     }else{ aw=0; break}
    }
    
  } #end while

  idp3 = which( is.na(para3[,1]) )
  para3[idp3,]=0; wtgd[idp3]=0; bmaw[idp3]=0
  
  if(numq > 1 ) {
    zp[,which( is.na(t(zp)[,1]) )]=0 #zp[,id]=0
    wtgd[which( is.na(t(zp)[,1]) )]=0
    bmaw[which( is.na(t(zp)[,1]) )]=0
  }else if(numq==1) {
    zp[which(is.na(zp))]=0  #id = which( is.na(zp) )
    wtgd[which(is.na(zp))]=0
    bmaw[which(is.na(zp))]=0
  }
  
  if( all(wtgd==0) ) {
    cat("==== all RL are NA ====")
    zx$zp.ma = zx$qua.lme;  zx$zp.bma = zx$qua.lme; return(zx)
  }
  
  bmaw[which(is.na(bmaw))]=0.0; 
  wtgd[which(is.na(wtgd))]=0.0

  zpf=rep(NA, numq); zpf.bma=rep(NA, numq)
  if(numq > 1){
    for (iq in 1:numq){
      zpf[iq]= t(wtgd) %*% t(zp)[1:numk,iq]
      if(bma==TRUE) zpf.bma[iq] = t(bmaw) %*% t(zp)[1:numk,iq]
    }
  }else if(numq==1){
      zpf[1]= t(wtgd) %*% zp  #t(zp)[1:numk,iq]
      if(bma==TRUE) zpf.bma[1] = t(bmaw) %*% zp  #t(zp)[1:numk,iq]
  } # end if numq

  # ---- calculating asymptotic SE by delta method -----------------
  
  if(varcom==TRUE ){
    # ---- cov interpolation using spline, for NA cov cases ---
    covint =cov.interp.new(numk, para3, cov2= mywt$prob.call$cov2)
    
    #---asymptotic variance matrix MatC and asymp var of MA ---------
    
    avar = asymp.var(mywt, covint, qqq=quant, order)
    
  }      # end if varcom=T and if not 'like' weight
  
  #----- finding surrogate for MA using gev -----------
  if(surr==TRUE){
    zpf.surr = t(wtgd)%*%t(xzp)
    para.ma= t(wtgd)%*% para3
    surro= surrogate(zpf.surr, xqa, init.surr= para.ma)
    zx$surr$par = as.vector(surro$par)
  }

  # ----- Predictive posterior variance ----------------
  
  if(varcom==TRUE & bma==TRUE ){
    zpdiff= matrix(NA, numq, numk); loc.var= matrix(NA, numq, numk) 
    for(ip in 1:numk){  
      
     if(numq > 1) {
       zpdiff[1:numq,ip]= (zp[1:numq,ip]-zpf.bma[1:numq])^2
     }else if(numq==1){
       zpdiff[1:numq,ip]= (zp[ip]-zpf.bma[1:numq])^2
     }
      
      loc.var[1:numq,ip] = avar$MatC[ip,ip,1:numq]
    } #end for
    
    msm = t(bmaw) %*% t(zpdiff)  # mean square of zp due to candidate models
    mse = t(bmaw) %*% t(loc.var) # mean square error of zp due to internal variation
  } # end if
  
  if(varcom==TRUE ){
    zx$fin.se.ma =  avar$fin.se.MA.qua
    zx$adj.se.ma =  avar$adj.se.MA.qua
    if(bma==TRUE) { zx$pred.se.bma = as.vector(sqrt(msm + mse) )
                      #{pred.se.bma}{A vector containing the Bayesian 
                      #   predictive SE of quantile estimates.}
      zx$bma.se.between = sqrt(msm)
      zx$bma.se.within = sqrt(mse) }
  }
  if(cov.lme.from.out != TRUE){ zx$qua.se.lme.boots = hosking$qua.lme.se
  }else{zx$qua.se.lme.boots = cov.lme$qua.lme.se}
  
  # ------- finalization --------------------------- 
  if(bma != TRUE){
    zx$zp= zp
    zx$zp.ma = zpf 
    zx$w.ma = wtgd 
    zx$pick_xi_ma = kpar ; zx$numk_ma = numk
  }

  if(bma==TRUE) {
    zx$zp.ma = zx$bma= zpf.bma; zx$pen=pen 
    zx$w.bma= bmaw
    
    set= set.prior(pen=pen, numk=numk, xi_lme=hosking$lme[3],
                   kpar=kpar, weight=weight)
    
    if(pen=="norm") zx$prior_mu_std = set$prior_mu_std
    if(pen=="beta") zx$p_q_beta = set$p_q_beta
    
    zx$pick_xi_bma = kpar 
    zx$numk_bma = numk
  }
  zx$run_kpar = run
  zx$original.numk = org.numk

  return(zx)
}

# ---end of main program --------------------------------------- 
#------------------------------------------------------------
gev.max=function (xdat, ntry=5) 
{
  z <- list();  k =list()
  n=ntry
  
  nsample=length(xdat)
  z$nsample=nsample
  
  init= matrix(0, nrow=ntry, ncol=3)
  init <- ginit.max(xdat,ntry)
  upsig= sqrt(var(xdat))*5
  upmu=  abs(mean(xdat))*3
  
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
                         LB =c(-upmu,0,-1),UB =c(upmu,upsig,1), 
                         control=list(trace=0, outer.iter=10,
                                      delta=1.e-6, inner.iter=200, tol=1.e-5) ))
      
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
ginit.max <-function(data=NULL,ntry=NULL){
  
  if(ntry < 2) ntry=2
  n=ntry
  init <-matrix(NA,ntry,3)
  
  lmom_init = lmoms(data,nmom=3)
  lmom_est <- pargev(lmom_init)
  
  init[1,1]    <-lmom_est$para[1]
  init[1,2]    <-lmom_est$para[2]
  init[1,3]    <-lmom_est$para[3]
  
  maxm1=ntry; maxm2=maxm1-1
  init[2:maxm1,1] <- init[1,1]+rnorm(n=maxm2,mean=0,sd = 5)
  init[2:maxm1,2] <- init[1,2]+rnorm(n=maxm2,mean=3,sd = 3)
  init[2:maxm1,3] <- runif(n=maxm2,min= -0.45,max=0.4)
  init[2:maxm1,2] = max(0.1, init[2:maxm1,2])
  
  return(init)
}
#------------------------------------------------
pargev.xifix= function (lmom, xifix= 0.1, checklmom = TRUE, ...) 
{
  
  para=rep(NA, 3)
  G=xifix
  para[3]= G
  
  if (length(lmom$lambdas[1]) == 0) {
    lmom <- lmorph(lmom)
  }
  
  SMALL <- 1e-05
  if( abs(G) < SMALL) {
    para[1:2]= pargum(lmom)$para
    return(list(type = "gev", para = para, source = "pargev"))
  }
  
  GAM <- exp(lgamma(1 + G))
  para[2] <- lmom$lambdas[2] * G/(GAM * (1 - 2^(-G)))
  para[1] <- lmom$lambdas[1] - para[2] * (1 - GAM)/G
  
  return(list(type = "gev", para = para, source = "pargev"))
}
#--------------------------------------------------------------  
gev.rl.delta_new <- function(data, ntry=5, quant){
  
  zg=list()
  numq=length(quant)
  # names_par= c("mu","sigma","shape")
  # names_quant=paste("q",as.character(quant[1:numq]),sep="")
  
  kmle.h = gev.max(data, ntry=ntry)
  #   if(abs(kmle.h$mle[3]) >= 0.4 ) kmle.h$mle[3] = sign(kmle.h$mle[3])*0.4
  
  kmle.coles = gev.fit(data, show=F)
  
  if(kmle.h$nllh < kmle.coles$nllh) {
    
    if(abs(kmle.h$mle[3]) >= 0.49 ) kmle.h$mle[3] = sign(kmle.h$mle[3])*0.49
    better=  tryCatch( gev.fit(data, show=F, muinit=kmle.h$mle[1],
                               siginit=kmle.h$mle[2], shinit= -kmle.h$mle[3])
    )
    # we need the above computing to obtain covariance
    
    better.mle = better$mle
    better.mle[3]= -better$mle[3]
    
  }else{
    better = kmle.coles
    better.mle= kmle.coles$mle
    better.mle[3]= -kmle.coles$mle[3]
  }
  
  if(abs(better.mle[3]) >= 1) {
    #    cat("gev.rl.delta_new: mle less than -1.0 =", better.mle[3],"\n")
    better.mle[3]= sign(better.mle[3])*0.98
  }
  
  savem= vec2par(better.mle,'gev')
  zg$nllh = better$nllh
  zg$mle = better.mle
  zg$qua.mle= quagev(quant[1:numq],savem)
  zg$data = data
  zg$cov= better$cov
  
  # names(zg$qua.mle)=names_quant
  # names(zg$mle)=names_par
  zg$quant=quant
  return(zg)
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
gev.xifix.sing =function (xdat=NULL, xifix= -0.1, ntry=5, varcom=NULL)   
{
  zx <- list();  kx =list(); value=list()
  #  truepar=rep(0,2)
  
  init.xi= matrix(0, nrow=ntry, ncol=2)
  
  init.xi <- ginit.xifix(xdat, ntry)
  xi=xifix
  
  upsig= diff(quantile(xdat, probs = c(0.25, 0.75)))*5 #EnvStats::iqr(xdat)*5
  upmu=  abs(median(xdat))*3
  
  #  par.start=rep(1,2)
  
  tryCatch(
    for(itry in 1:nrow(init.xi)){
      
      #      par.start[1:2]=init.xifix[itry,1:2]
      
      value <- try(optim(init.xi[itry,], fn=gev.xilik2, 
                         lower =c(-upmu,0), upper =c(upmu,upsig), 
                         method="L-BFGS-B",
                         xifix=xifix, xdat=xdat) )
      # control=list(trace=1, outer.iter=40, 
      #              inner.iter=200, tol=1.e-5,
      #              delta=1.e-6) ))
      
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
    #     zx$grad <- numDeriv::grad(gev.xilik, x$par, xifix=xifix, xdat=xdat)
    if(varcom==TRUE){
      Hess = PrescottW(par2=x$par, xifix=xifix, nsam=length(xdat))
      #numDeriv::hessian(gev.xilik, x$par, xifix=xifix, xdat=xdat)
      if(is.na(Hess[1,1])) {
        zx$cov =matrix(NA,2,2)
      }else if(det(Hess) <= 0) {
        zx$cov =matrix(NA,2,2)
      }else{
        zx$cov <- solve(Hess)
      }
    } #end if varcom
  }
  class(zx) <- "gev.xifix"
  return(zx)
  invisible(zx)
}
#---------- end of gev.xifix program -----------------
# -- Fisher Info matrix for xi fixed by Prescott and Walden ----
PrescottW = function(par2=NULL, xifix=NULL, nsam=NULL){
  
  #   mu=par2[1]
  sig=par2[2]
  xi=xifix
  Hess=matrix(NA,2,2)
  
  if(xi < 0.5){
    pw= (1-xi)^2 * gamma(1-2*xi)
    Hess[1,1] = nsam*pw/(sig^2)
    Hess[2,2] = nsam*(1-2*gamma(2-xi)+pw)/(sig^2 * xi^2)
    Hess[1,2] = nsam*(pw-gamma(2-xi))/(sig^2 * xi)
    Hess[2,1] = Hess[1,2]
  }else{
    cat("FIM is not available because xi >= 0.5","\n")
  }
  return(Hess)
}
#--------------------------------------------------------------- 
cons.MatC = function(mywt, cov22, quant){
  
  numk= length(mywt$wtgd)
  numq= length(quant)  
  MatC = array(NA, c(numk, numk, numq))
  cov2i= matrix(NA,2,2); cov2j=matrix(NA,2,2)
  
  para3 = mywt$prob.call$mle3
  
  for (i in 1:numk){
    mle3i= para3[i,]
    cov2i= as.matrix(cov22[,,i])         #mywt$prob.call$cov2[[i]]
    
    for (j in i:numk){
      mle3j = para3[j,]
      cov2j = as.matrix(cov22[,,j])      #mywt$prob.call$cov2[[j]]
      
      if( any(is.na(cov2i)) | any(is.na(cov2j)) ) {
        MatC[i,j,1:numq]=NA
      }else{
        MatC[i,j,1:numq]=delta.gev(gevf3=NULL, mle3i=mle3i, cov2i=cov2i,
                                   mle3j=mle3j, cov2j=cov2j,
                                   quant=quant, d3yes=F)$covij[1:numq]
        MatC[j,i,1:numq]= MatC[i,j,1:numq]
      } #end if 
    } #end for j
  } #end for i
  
  wtgd = mywt$wtgd; bmaw = mywt$bmaw
  fin.se =rep(NA, numq); fin.se.bma =rep(NA, numq)
  
  for (iq in 1:numq){
    fin.se[iq] = sqrt(t(wtgd) %*% MatC[,,iq] %*% wtgd)
    fin.se.bma[iq] = sqrt(t(bmaw) %*% MatC[,,iq] %*% bmaw)
  }
  zz=list()
  zz$MatC = MatC 
  zz$fin.se = fin.se; zz$fin.se.bma = fin.se.bma
  zz$numk = numk; zz$numq=numq; zz$wtgd=wtgd; zz$bmaw=bmaw
  return(zz)
}

gev.xilik2 <- function(a, xifix=xifix, xdat=xdat) {
  gev.xilik(a, xifix=xifix, xdat=xdat)
}
#----------------------------------------------------------
gev.xilik <- function(a, xifix=NULL, xdat=NULL) {
  
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
  
  b2= sum(log(y)) * (1-ooxi)
  b3= sum(exp(ooxi * log(y)) )
  
  nllh =  nsam*(log(sc)) + b2 + b3
  
  zz=nllh 
  return(zz)
}

#------------------------------------------------------------
#-----------------------------------------------------------  
lme.boots.new <-function(data, B=NULL, quant, boot=TRUE, trim=NULL,
                         PB=FALSE, par=NULL){
  
  z=list();  numq=length(quant)
  # names_par= c("mu","sigma","shape")
  # names_quant=paste("q", as.character(quant[1:numq]), sep="")
  sam.bootL=list()
  lme.boot =matrix(NA,nrow=B,ncol=3);  ratios =matrix(NA,nrow=B,ncol=3)
  lme.seboot =matrix(NA,nrow=B,ncol=numq); med = rep(NA, B)
  lambdas =matrix(NA,nrow=B,ncol=3)
  
  nsam=length(data); re=sort(data)
  lmt= lmoms(data, nmom=3)
  klmet= pargev(lmt, checkmom=T)
  z$lme = klmet$para[1:3]
  savet= vec2par(klmet$para[1:3],'gev')
  z$qua.lme= quagev(quant[1:numq],savet)
  
  z$quant = quant
  
  # names(z$qua.lme)=names_quant; names(z$lme)=names_par
  
  if(boot==T | boot==TRUE){
    
    if(PB != TRUE){
      for (ib in 1:B){
        sam.bootL[[ib]]= sample(data, nsam,replace=T)   # nonpara Boot
      }
    }else if(PB== TRUE){
      for (ib in 1:B){
        sam.bootL[[ib]]= rlmomco(nsam, vec2par(par,'gev'))  # para Boot
      }
    } # end if PB
      
    for (ib in 1:B){
      lm.boot= lmoms(sam.bootL[[ib]], nmom=3)
      
      lme.boot[ib,1:3]= pargev(lm.boot, checklmom=F)$para[1:3] 
      save=vec2par(lme.boot[ib,1:3],'gev')
      lme.seboot[ib,1:numq]= quagev(quant[1:numq],save)
      
      # if(trim > 0){
      #   lm.boot= TLmoms(sam.bootL[[ib]], nmom=3, leftrim=trim)
      # }
      
      lambdas[ib,1:3] = lm.boot$lambdas[1:3]
      ratios[ib,1]= lm.boot$lambdas[1]
      ratios[ib,2:3] = lm.boot$ratios[2:3]
      med[ib]= median(sam.bootL[[ib]])
      
    } # end for ib
    
    lme.bb = ratios
    lme.bb[1:B,1]= med[1:B]
    z$cov.med.rat= cov(lme.bb[,1:3])
    
    lme.bb = lambdas
    lme.bb[1:B,1]= med[1:B]
    z$cov.med.lam= cov(lme.bb[,1:3])
    
    z$cov.lambda = cov(lambdas[,1:3])
    
    z$B= B; z$lme.boot= lme.boot
    #z$ratios = ratios;  z$med = med; z$lambdas = lambdas
    
    z$cov.rat= cov(ratios[,1:3])
    z$cov.par= cov(lme.boot[,1:3])
    
    qua.lme.se=rep(NA, numq)
    for (qi in 1:numq){
      qua.lme.se[qi]=sqrt(var(lme.seboot[1:B,qi],na.rm=T))
    }
    z$qua.lme.se= qua.lme.se  
    z$qua.lme.boot = lme.seboot
    
  }  #end if(boot==TRUE)
  
  return(z)
}
# #------------------------------------------------------------

#------- Reversed pf -----------------------------------
rcd = function(xi){

  minxi = -0.2
  nxi=length(xi)
  pfxi=rep(NA, nxi)

  for(i in 1:nxi){
    if( is.na(xi[i]) | is.null(xi[i]) ){
      pfxi[i]=0
    }else{
      if( -1.0 < xi[i] & xi[i] <= minxi) {
        b1= -1/(1+minxi)
        b0= 1-b1*minxi

        pfxi[i] = b0+b1*xi[i]
      }
    } #end if
    if( xi[i] <= -1.0 | xi[i] > minxi) pfxi[i]=1
  } #end for
  return(pfxi)
}
#-------------------------------------------------------
emp.prior = function(mle3, mle, lme){

  z=list()
  xi= mle3[,3]
  nxi=length(xi)
  pfxi=rep(NA, nxi);  pfmu=rep(NA, nxi);  pfsig=rep(NA, nxi)

  mnx= min(mle[3],lme[3]); mxx=max(mle[3],lme[3])
  u= -0.6; u0= -0.0;

  if(mxx <= u) {pfxi=1; return(pfxi)}
  if(mnx <= u) {mnx=u+0.05; b=0}
  else{ b= 1/(mnx-u-u0)}
  b2=1/(mnx-u-u0)

  for(i in 1:nxi){
    if( is.na(xi[i]) | is.null(xi[i]) ){
      pfxi[i]=0
    }else{
      if( u+u0 < xi[i] & xi[i] < mnx) {
        pfxi[i] = (xi[i] -u-u0)*b
      }else if(mnx <= xi[i] & xi[i] <= mxx) {
        pfxi[i] = 1
      }else if( mxx < xi[i] & xi[i] < mxx+mnx-u){
        pfxi[i] = 1 +(mxx-xi[i])*b2
      }else if( xi[i] <= u+u0 | xi[i] >= mxx+mnx-u) {pfxi[i]=0
      }
    } #end if
  } #end for

  return(pfxi)
}

#----------------------------------------------------   
movave= function(order,numk,numq, zp, wt=NULL){
  
  za=list()   
  q=order
  ez=matrix(NA,numq,numk);  ew=rep(NA,numk)
  
  if(numk >= (q+1)){
    for (ip in 1:numq){
      ez[ip,]=rollmean(zp[ip,1:numk], q, fill="extend",align="center")
    }
  }else if(numk <= q){
    for(ip in 1:numq){
      for(ik in 1:numk){
        ez[ip,ik] = mean(zp[ip,1:numk])
      }
    }
  } #end if
  
  if(numk >= (q+1)){
    ew = rollmean(wt,q,fill="extend",align="center")
    
  }else if(numk <= q){
    for(ik in 1:numk){
      ew[ik] = mean(wt[1:numk])
    }
  } #end if
  za$ez =ez; za$ew=ew
  return(za)
}
#-----------------------------------------------------------------
delta.gev = function(gevf3=NULL, mle3i=NULL, cov2i=NULL,
                     mle3j=NULL, cov2j=NULL,
                     quant=NULL, d3yes=FALSE){
  
  z=list()  
  numq=length(quant)
  # quant=.99 means T= 100 = 1/(1-quant) return period
  # gevf3  =list from gev.fit with mle.hosking
  # gevf2i,j =list from gev.fit(xi.fixed) with mle.hosking
  
  p= 1-quant
  yp= -log(1-p)
  
  if(d3yes==TRUE){                  # d3yes =T for computing delta3 cov mtx
    c3 = gevf3$cov
    mle = gevf3$mle  #mle.hosking
    xi =  mle[3] 
    sig = mle[2]
    d12= cbind(rep(1,numq), (1-yp^(xi))/xi )
    
    d3 = sig*(1-yp^(xi))/(xi^2) + sig * ( yp^(xi) ) * log(yp)/xi
    delta3 = cbind(d12, d3 )
    
    for(iq in 1:numq){
      z$v3[iq] = delta3[iq,] %*% c3 %*% delta3[iq,]
    }
    return(z)
    # z$v3 = asymptotic variance of zp from gev(theta hat):
    #                            without xi fixed
  }
  
  xi2i = mle3i[3] #gevf2i$mle[3]  # fixed xi_i 
  xi2j = mle3j[3] #gevf2j$mle[3]  # fixed xi_j 
  #    cov2i = cov2i #gevf2i$cov2[1:2,1:2]
  #    cov2j = cov2j #gevf2j$cov2[1:2,1:2] 
  
  # --- compute correlation between i and j ------------------
  
  qqc= c(seq(0.1,0.9,0.1), 0.95, 0.98, 0.99, 0.995, 0.998, 0.999)
  zp.cori= quagev(qqc,vec2par(mle3i,'gev'))
  zp.corj= quagev(qqc,vec2par(mle3j,'gev'))
  corr.rt = cor(zp.cori, zp.corj)
  
  # xi= (xi2i +1)/2; xj= (xi2j +1)/2
  # corr.rt2= exp(-(xi-xj)^2 /2)
  
  covij=rep(NA, length(quant)); vi=rep(NA,numq); vj=rep(NA,numq)
  deltai=matrix(NA, 2,numq); deltaj=matrix(NA, 2,numq)
  
  for (id in 1:length(quant)){
    deltai[1:2,id]= c(1, (1-yp[id]^(xi2i))/xi2i)
    deltaj[1:2,id]= c(1, (1-yp[id]^(xi2j))/xi2j)
    
    vi[id] = t(deltai[,id]) %*% cov2i %*% deltai[,id]
    vj[id] = t(deltaj[,id]) %*% cov2j %*% deltaj[,id]
    
    covij[id] = corr.rt* sqrt(vi[id] * vj[id])
  }
  
  z$covij = covij
  # z$vi=vi; z$vj =vj
  return(z)
}
#---------------------------------------------------------------------
#-----------------------------------------------------------
surrogate = function(zpf.surr=NULL,xqa, init.surr){
  
  press.fun = function(par){
    theta=par
    if(theta[2] <= 0)  return(10^10)
    if(theta[3] <= -1.0 | theta[3] >= 1)  return(10^10)
    mzp= quagev(xqa, vec2par(theta, 'gev'))
    fun = sum((zpf.surr - mzp)^2)
    #fun = sum(abs(zpf.surr - mzp))
    return(fun)
  }
  
  surrx=list()
  surrx = optim(par=init.surr, fn=press.fun )  #, method="L-BFGS-B", lower=c(-Inf, 0, -1), upper=c(Inf, Inf, 1))
  
  surrx$zp.surrmodel = quagev(xqa, vec2par(surrx$par, 'gev'))
  surrx$zp.MA = zpf.surr
  return(surrx)
}
#---------------------------------------------------
#------------------------------------------------------------
cov.interp.new = function(numk,para3,cov2=NULL){  
  
  cov22 = array(NA,c(2,2,numk))
  
  for(ik in 1:numk){ 
    cov22[,,ik] = cov2[[ik]]   #mywt$prob.call$cov2[[ik]] 
  }
  
  id= which(!is.na(cov22[1,1,]))
  id.na = which(is.na(cov22[1,1,]))
  if(length(id.na)==numk) stop("all cov22 are NA")
  
  if( length(id) < numk){
    
    cat("Need cov interpolation because some cov are NA",
        "id for cov interpol=",id.na,"\n")
    
    xsp = para3[id,3]; y1= cov22[1,1,id]
    y2= cov22[2,2,id]
    y12= cov22[1,2,id]
    
    covint = matrix(NA,3,numk)
    covint[1,id.na]= spline(x=xsp, y1, method="natural", 
                            xout=para3[id.na,3])$y
    covint[2,id.na]= spline(x=xsp, y2, method="natural", 
                            xout=para3[id.na,3])$y
    covint[3,id.na]= spline(x=xsp, y12, method="natural", 
                            xout=para3[id.na,3])$y
    
    cov22[1,1,id.na]= covint[1,id.na]
    cov22[2,2,id.na]= covint[2,id.na]
    cov22[1,2,id.na]= covint[3,id.na]
    cov22[2,1,id.na]= covint[3,id.na]
    
  } #end if 
  return(cov22)
} 
#-----------------------------------------------------
#---------------------------------------------------------
avar.ff.rl = function(mywt, covint, qqq=NULL,
                      hosking, numom=numom, consC.zp, order){
  
  zx=list()
  kpar= hosking$kpar; data =hosking$data; trim=hosking$trim
  wtgd= mywt$wtgd; bmaw = mywt$bmaw
  numq= length(qqq); numk = length(wtgd); cov.type=hosking$cov.type
  mywt2=list()
  
  mywt2 = weight.com.new(data,numk, hosking, kpar, numom=numom,
                         xqa=qqq, varcom=T, boot.lme, cov.lme,
                         surr=T, type="full", trim=trim, cov.type=cov.type)
  
  zp = mywt2$xzp
  consC.zp2 =list()
  
  consC.zp2 = cons.MatC(mywt2, covint, quant=qqq)
  
  # consC.zp2$fin.se
  # consC.zp2$fin.se.bma
  
  ezvar=rep(NA,numq); ewvar=rep(NA,numq); ewvar.bma=rep(NA,numq)
  trCD=rep(0, numq); trCD.bma=rep(0, numq)
  
  D = cov.dir(wtgd)
  Dbma = cov.dir(bmaw)
  
  move= movave(order, numk, numq, zp, wt=wtgd)
  ez=move$ez
  ew= move$ew
  
  moveb= movave(order, numk, numq, zp, wt=bmaw)
  ewbma= moveb$ew
  
  for(ip in 1:numq){
    ezvar[ip]= t(ez)[,ip] %*% D %*% ez[ip,] 
    ewvar[ip]= t(ew) %*% consC.zp2$MatC[,,ip] %*% ew
    ewvar.bma[ip]= t(ewbma) %*% consC.zp2$MatC[,,ip] %*% ewbma
    
    trCD[ip]= sum( diag(consC.zp2$MatC[,,ip] %*% D) ) + ezvar[ip]
    trCD.bma[ip]= sum( diag(consC.zp2$MatC[,,ip] %*% Dbma) ) + ezvar[ip]
  }
  zx$adj.se.MA.qua = sqrt(ewvar + trCD )
  zx$adj.se.bma.qua = sqrt(ewvar.bma + trCD.bma) 
  return(zx)
}
#------------------------------------------------
#------construction of asymptotic variance matrix MatC 
#---- and asymp var of MA and bma --------------------- 

asymp.var =  function(mywt, covint, qqq=NULL, order=2){
  
  zx=list(); consC=list()
  wtgd= mywt$wtgd; bmaw = mywt$bmaw
  numq=length(qqq); numk = length(wtgd)
  zp=mywt$zp
  
  consC = cons.MatC(mywt, covint, quant=qqq)
  
  zx$fin.se.MA.qua = consC$fin.se
  zx$fin.se.bma.qua = consC$fin.se.bma
  
  ezvar=rep(NA,numq); ewvar=rep(NA,numq); ewvar.bma=rep(NA,numq)
  trCD=rep(0, numq); trCD.bma=rep(0, numq)
  
  D = cov.dir(wtgd)
  Dbma = cov.dir(bmaw)
  
  move= movave(order, numk, numq, zp, wt=wtgd)
  ez=move$ez
  ew= move$ew
  
  moveb= movave(order, numk, numq, zp, wt=bmaw)
  ewbma= moveb$ew
  
  for(ip in 1:numq){
    ezvar[ip]= t(ez)[,ip] %*% D %*% ez[ip,]
    ewvar[ip]= t(ew) %*% consC$MatC[,,ip] %*% ew
    ewvar.bma[ip]= t(ewbma) %*% consC$MatC[,,ip] %*% ewbma
    
    trCD[ip]= sum( diag(consC$MatC[,,ip] %*% D) ) + ezvar[ip]
    trCD.bma[ip]= sum( diag(consC$MatC[,,ip] %*% Dbma) ) + ezvar[ip]
  }
  # zx$adj.se.MA.qua = sqrt(consC$fin.se^2 + trCD )
  # zx$adj.se.bma.qua = sqrt(consC$fin.se.bma^2 + trCD.bma) 
  
  zx$adj.se.MA.qua = sqrt(ewvar + trCD )
  zx$adj.se.bma.qua = sqrt(ewvar.bma + trCD.bma) 
  
  zx$MatC = consC$MatC
  
  # zx$ew.adj.se.MA.qua = sqrt(consC$ew.fin.se^2 + trCD )
  # zx$ew.adj.se.bma.qua = sqrt(consC$ew.fin.se.bma^2 + trCD.bma)
  return(zx)
}  
#--------------------------------------
#--------------------------------------------      
cov.dir = function(wtgd){
  
  numk=length(wtgd)
  D=matrix(NA, numk, numk)
  for(i in 1:numk) D[i,i] = wtgd[i]*(1-wtgd[i])/2
  for(i in 1:(numk-1)){
    for(j in (i+1):numk){
      D[i,j]= -wtgd[i]*wtgd[j]/2
      D[j,i]= D[i,j]
    }
  }
  return(D)
} 
# --------------------------------------------------------------------