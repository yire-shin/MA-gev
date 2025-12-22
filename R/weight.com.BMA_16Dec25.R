#------------------------------------------------  
prior.beta.stnary = function(x=NULL, xi_lme=NULL,
                             p=2, q=5, al, bl){
  
  # aa= max(-1, xi_lme-c0)
  # bb= min(0.5, xi_lme+c0)
  # al=min(aa,bb)
  # bl=max(aa,bb)
  # 
  # if(xi_lme <= 0) {
  #   qlim= min(abs(xi_lme)*c1, c2)
  # }else{ qlim =0.0 }
  # 
  # p=p; q= p+qlim
  
  Bef <- function(w) { ((-al+w)^(p-1)) * ((bl-w)^(q-1)) }
  Be  <- integrate(Bef, lower=al, upper=bl)[1]$value
  
  # x= xi 
  if( (-al+x) <= 0) {left.w =0
  }else{ left.w = (-al+x)}
  if( (bl-x) <= 0) {right.w=0
  }else{ right.w= (bl-x)}
  
  prior <- (left.w^(p-1))*(right.w^(q-1))/ Be 
  
  return(prior)
  
} # end function
#------------------------------------------------------------------
#------------------------------------------------------------------    
set.prior = function(pen=NULL, numk=NULL, xi_lme=NULL, kpar=NULL,
                     weight=NULL) {
  
  zp = list()
  prior=rep(NA, numk)
  lme=xi_lme  # just xi
  
  if(weight=='like'){
    if(pen=="beta"){
      p=2; c0=0.35; c1=40; c2=12
    }else if(pen=="norm"){
      multa=2.2; lme.cut= -0.5; std.cut= -0.45 #std.prior=0.12 
    }
    
  }else if(weight=='gLd' | weight=='med'){
    if(pen=="beta"){
      p=2; c0=0.35; c1=20; c2=7
    }else if(pen=="norm"){
      multa=1.5; lme.cut= -0.45; std.cut= -0.4 #std.prior=0.14
    }
  } # end if weight
  
  if(pen=="norm"){

    if(lme >= 0) {
      prior = rep(1/numk, numk)
      mu.prior=0; std.prior=1 }
    
    if(lme < 0){
      
      if(lme >= lme.cut){ mu.prior = multa* lme
      }else if(lme < lme.cut){ mu.prior = multa*lme.cut}
      
      if(weight=="like"){
        if(lme >= std.cut){ std.prior = -(std.cut-lme)/5 +0.11
        }else{ std.prior = 0.11 }  
      }
      if(weight=="gLd" | weight=="med"){
        if(lme >= std.cut){ std.prior = -(std.cut-lme)/4 +0.14
        }else{ std.prior = 0.14 }  
      }
      
      for (i in 1:numk){
        prior[i] = dnorm(x=kpar[i], mean= mu.prior, sd= std.prior) 
      }
    } # end if lme < 0
    
  }else if(pen=="beta"){
    
    aa= max(-1, xi_lme-c0)
    bb= min(0.5, xi_lme+c0)
    al=min(aa,bb)
    bl=max(aa,bb)
    
    if(xi_lme <= 0) {
      qlim= min(abs(xi_lme)*c1, c2)
    }else{ qlim =0.0 }
    
    p=p; q= p+qlim
    
    for (i in 1:numk){
      prior[i] = prior.beta.stnary(x=kpar[i], xi_lme=lme,
                       p=p, q=q, al=al, bl=bl)
    }
  } # end if pen
  
  zp$prior = prior
  if(pen=="norm") zp$prior_mu_std =c(mu.prior, std.prior)
  if(pen=="beta") zp$p_q_beta = round(c(p,q),2)
  return(zp)
}
#-------------------------------------------------
#------- weight computing----------------------------------------
weight.com.new <- function(data=NULL, numk=NULL, hosking=NULL, 
                           kpar=NULL, numom=3, xqa=NULL, varcom=TRUE, 
                           cov.lme=NULL, surr=FALSE,
                           type='full', trim=NULL, cov.type='ratio',
                           bma=TRUE, pen="norm"){
  
  z=list();  
  quant = hosking$quant; weight = hosking$weight
  B = hosking$B; 
  start = hosking$start

  if(bma==TRUE){

    prior = set.prior(pen=pen, numk=numk, xi_lme=hosking$lme[3],
                      kpar=kpar, weight=weight)$prior

    # id=which(!is.na(prior))
    # prior= na.omit(prior); numk= length(prior)
    # nkpar= kpar[id]; kpar= nkpar
    
    # pron= prior/sum(prior)
    # 
    # if(max(pron, na.rm=T) > 0.75) {
    #   cnpar = new.kpar2(wtgd=pron, numk=numk, kpar=kpar,
    #                     remove=1e-5, dist=0.03)
    #   kpar= cnpar$kpar2
    #   numk= cnpar$numk
    # 
    #   prior = set.prior(pen=pen, numk=numk, xi_lme=hosking$lme[3],
    #                   kpar=kpar, p.beta, c0.beta, c1.beta, c2.beta )
    # 
    # } # end if max(pron)
    # 
    # id=which(!is.na(prior))
    # prior= na.omit(prior); numk= length(prior)
    # nkpar= kpar[id]; kpar= nkpar

  }else if(bma != TRUE){
    prior = rep(1/numk, numk)
  } # end if bma
  
  prob.call=list()
  
  if(is.null(cov.lme)) {boot.lme=TRUE
  }else{boot.lme=FALSE}
  
  prob.call = dist.noboot(data=data, numk=numk, hosking=hosking, 
                          boot.lme= boot.lme,
                          kpar=kpar, numom=numom, ntry=5,varcom=varcom,
                          cov.lme=cov.lme, trim=trim, cov.type=cov.type)
  
  notid= which(!is.na(prob.call$mle3[,1]))
  prob.call$mle3= prob.call$mle3[notid,] ; mle3 = prob.call$mle3
  numk=length(notid)
  
  kpar=kpar[notid]
#  prior= rcd(kpar)

  if(surr==TRUE) { 
    nxq = length(xqa);
    xzp= matrix(NA, nrow= nxq, ncol=numk) } 
  numq= length(quant)

  zp= matrix(NA, nrow=numq, ncol=numk)
  wtgd=rep(NA,numk); bmaw=rep(NA,numk); #bmaw2=rep(NA,numk)
  
  lme=hosking$lme; mle= hosking$mle
  # eprior = emp.prior(mle3=mle3, mle=mle, lme=lme)
  # 
  # eprior[which(is.na(eprior))]=0
  # if(all(eprior==0)) eprior=1
  
  if(weight=='gLd' | weight=='med' | weight=='cvt' |weight=='fcv' ){
    
    if(varcom==TRUE | weight=='fcv'){
      cov2=list()
      for(i in 1:numk) cov2[[i]]= prob.call$cov2[[ notid[i] ]] 
    }
    
    if(weight == 'cvt') {
      
      naic= prob.call$aic[notid]  #-ibay*log(prior)
      id= which(is.na(naic));# nid = which(!is.na(naic))
      
      amin =naic -min(naic, na.rm=T)
      wtgd =exp(-amin)/sum( exp(-amin), na.rm=T)
      wtgd[id]=0
      
      bmaw = prior*exp(-amin)/sum( prior*exp(-amin), na.rm=T)
      bmaw[id]=0
      
      # bmaw2 = eprior*exp(-amin)/sum( eprior*exp(-amin), na.rm=T)
    }
    
    if(weight == 'gLd' | weight == 'med' | weight=='fcv'){
      
      if(weight=='gLd') {kcol=1} else if(weight=='med') {kcol=2}
      if(weight=='fcv') {kcol=3}
      
      prob.mtx= prob.call$prob.mtx[notid,] 
      wtgd[1:numk]= prob.mtx[1:numk,kcol]/sum(prob.mtx[,kcol], na.rm=T) 
      
      #      pr.ms = psigma(par=prob.call$mle3[,1:3])
      
      bmaw[1:numk] = prior*prob.mtx[1:numk,kcol]/sum(prior*prob.mtx[,kcol], na.rm=T)
      
      # bmaw2[1:numk] = eprior*prob.mtx[1:numk,kcol]/sum(eprior*prob.mtx[,kcol], na.rm=T)
      
    }  # end if
    
  } #end if weight

  if(weight=='lcv' | weight=='like' | weight=='opt'){
    
    pertr=hosking$pertr
    
    if(weight=='lcv' | weight=='like'){
      
      prob.call= wlik.xifix(data, numk=numk, kpar=kpar, weight=weight, 
                            pertr=pertr, varcom=varcom, type=type, 
                            prior=prior, trim= trim) # note that trim=0
      
      if(varcom==TRUE){ cov2=list()
        for(i in 1:numk) cov2[[i]]= prob.call$cov2[[ notid[i] ]]  
      } #end if varcom
    } 
    
    wtgd = prob.call$wt
    bmaw = prob.call$bmaw;  #bmaw2 = prob.call$bmaw2
    prob.call$mle3 = t(prob.call$kfix)

    for (ip in 1:numk){
      save= vec2par(prob.call$kfix[1:3,ip],'gev')
      zp[1:numq,ip]= quagev(quant[1:numq],save)
      
      if(surr==TRUE) xzp[1:nxq,ip]= quagev(xqa[1:nxq],save)
    }   # end for ip
    
  }else{      #if(weight=='gLd' | 'med' | 'fcv' | 'cvt')
    
    for (ip in 1:numk){
      if( any( is.na(prob.call$mle3[ip,1:3]) ) ) {
        zp[1:numq,ip] = NA
        if(surr==TRUE) xzp[1:nxq,ip] = NA
      }else{
        save= vec2par( prob.call$mle3[ip,1:3],'gev' )
        zp[1:numq,ip]= quagev(quant[1:numq],save)
        if(surr==TRUE) xzp[1:nxq,ip]= quagev(xqa[1:nxq],save)
      }
    }    # end for ip
    
  }  # end of if weight= 'lcv' or 'like'
  
  z$weight= weight ; z$numk=numk
  z$wtgd= wtgd; z$bmaw = bmaw ; #z$bmaw2= bmaw2
  z$zp=zp;  z$kpar=kpar
  z$prob.call= prob.call 
  
  if(varcom==TRUE) {z$prob.call$cov2 = cov2    
  }else{z$prob.call$cov2=NA}        # cov2[[ip]]: 2*2 cov.mat of (mu, sig)
  if(surr==TRUE) z$xzp = xzp           # for fitting a surrogate gev model
  return(z)
}
# --------------------------------------------------------------
wlik.xifix = function(data=NULL,numk=NULL,kpar=NULL,weight=NULL, 
                      pertr=0.85,varcom=FALSE, type=NULL, 
                      prior=NULL, trim=NULL){  
  
  z=list()
  eprior=prior
  kfix=matrix(NA,nrow=3,ncol=numk); 
  nbic=rep(NA,numk); bmin=rep(NA,numk)
  wt=rep(NA,numk); bmaw=rep(NA,numk); #bmaw2=rep(NA,numk)
  para=rep(NA,3); lm3=list()
  re=sort(data); nsam=length(data)
  
  ntrim= trim  #max(trim-1,0)
  lmtrim= lmoms(re[(ntrim+1):nsam], nmom=3)
  if(is.null(eprior)) eprior=1
  
  cov2.cv=list()
  
  for (ip in 1:numk){
    
    kfix[1:3,ip]= pargev.xifix(lmtrim, xifix=kpar[ip])$para[1:3] 
    
    #    cat(trim, quagev(quant,vec2par(kfix[1:3,ip],'gev')),"\n")
    # }
    #   There is NO sig. diff. between trim=0 and trim=1
    
    if(varcom==TRUE){
      Hess = PrescottW(par2=as.vector( kfix[1:2,ip]), xifix=kpar[ip],
                       nsam=length(re[(ntrim+1):nsam]))
      
      if(is.na(Hess[1,1])) {
        cov2.cv[[ip]] = matrix(NA,2,2)
      }else if(det(Hess) <= 0) {
        cov2.cv[[ip]] = matrix(NA,2,2)
      }else{
        cov2.cv[[ip]] <- solve(Hess)
      }
    }
    
  } #end for ip
  
#  prior = rcd(kpar)
  
  if(weight=='like'){
    for (ip in 1:numk){
      
      ntrim= trim
      dtrim=re[(ntrim+1):nsam]
      
      nbic[ip]= gev.xilik( a=kfix[1:2,ip], xifix= kpar[ip], xdat= dtrim )
      
      if(nbic[ip] >= 10^6)  nbic[ip]=NA
    } # end for ip
    
    numpar=2                      # Becoz xi fixed
    nbic = 2*nbic + 2*numpar #- ibay*prior     # AIC
    bmin = nbic -min(nbic, na.rm=T)
    wt =exp(-bmin)/sum(exp(-bmin), na.rm=T)
    #    for (ip in 1:numk){ if( is.na(wt[ip]) ) wt[ip]=0 }
    
    bmaw= prior*exp(-bmin)/sum(prior*exp(-bmin), na.rm=T)
#    bmaw2= eprior*exp(-bmin)/sum(eprior*exp(-bmin), na.rm=T)
    
  }else if(weight=='lcv'){
    
    prob.mtx=rep(NA, numk)
    
    cvp = fcv.kfix.new(data,pertr,mle3=NULL,numk,kpar, 
                       weight,type)
    
    prob.mtx= exp( cvp$prob )
    
    prob.mtx[which(is.na(prob.mtx))]=0
    if(all(prob.mtx==0)){
      prob.mtx= exp( cvp$prob/100 )/100
    }
    
    wt= prob.mtx[1:numk]/sum(prob.mtx, na.rm=T) 
    
    bmaw = prior*prob.mtx[1:numk]/sum(prior*prob.mtx, na.rm=T)
#    bmaw2 = eprior*prob.mtx[1:numk]/sum(eprior*prob.mtx, na.rm=T)
    
  }  #end if weight
  
  bmaw[which(is.na(bmaw))]=0; # bmaw2[which(is.na(bmaw2))]=0
  wt[which(is.na(wt))]=0
  
  z$bmaw=bmaw; z$wt= wt; #z$bmaw2=bmaw2;
  z$aic= nbic
  z$kfix= kfix
  if(varcom==TRUE){z$cov2 = cov2.cv}
  return(z)
}
#--- getting the distance probability for weights calculation -----
# --------------------------------------------------------------
dist.noboot = function(data=NULL, numk=NULL, hosking=NULL, boot.lme=TRUE,
                       kpar=NULL, numom=NULL, ntry=5, varcom=NULL, 
                       cov.lme=NULL, trim=NULL, cov.type='lambda'){
  
  cov2=list();  kfix=list();  zw=list()
  weight = hosking$weight
  B = hosking$B; #if(boot.lme==TRUE) boot.lme=TRUE
  
  mle3= matrix(NA, nrow=numk, ncol=3)
  lmB.med= rep(NA,B);   aic=rep(NA,numk)
  
  nsample= length(data)
  lm= lmoms(data, nmom=numom)
  medata= median(data)
  
  # re=sort(data); nsam=length(data)
  # dtrim=re[(trim+1):nsam]
  
  #  if(weight=='lcv') varcoma=T

  for (ip in 1:numk){
    
    if(is.na(kpar[ip])){
      aic[ip]=NA
      mle3[ip,1:3]=NA
      print(paste("NA in kpar, ip=",ip))
    }else{
      kfix[[ip]]= gev.xifix.sing(xdat=data, xifix=kpar[ip], ntry=ntry,
                               varcom=varcom)

    # cat(kpar[ip],quagev(c(.99,.995), vec2par(kfix[[ip]]$mle[1:3],'gev')),"\n")
    #   }
    
     if(kfix[[ip]]$conv != 0) {
      print(paste("No convergence in gev.xifix, ip=",ip))
      aic[ip]=NA
      mle3[ip,1:3]=NA
     }else{
      aic[ip]= 2*kfix[[ip]]$nllh + 2*2
      mle3[ip,1:3]= kfix[[ip]]$mle[1:3]
     }
     if(varcom==TRUE) cov2[[ip]] = kfix[[ip]]$cov
     
    } # end if(is.na..)
    
  } #end for
  
  
  if( all(is.na(mle3[1:numk,1]) ) ){
    lm3=lmoms(data, nmom=3)
    for(ip in 1:numk){
      mle3[ip,1:3] = pargev.xifix(lm3, xifix=kpar[ip])$para[1:3]
    }
  } #end if
  
  
  # if(weight=='lcv') {
  #   zw$cov2 = cov2
  #   zw$mle3= mle3
  #   return(zw)
  # }
  
  if(weight =='gLd'){
    
    if(cov.type=='ratio'){
      if(boot.lme==TRUE){
        cov.lm = hosking$cov.rat
      }else{
        cov.lm = cov.lme$cov.rat
      }
    }else if(cov.type=='lambda'){
      if(boot.lme==TRUE){
        cov.lm = hosking$cov.lambda
      }else{
        cov.lm = cov.lme$cov.lambda
      }
    }
    Vinv= solve(cov.lm)
    detV= det(cov.lm)
    
    zw$cov.lm=list(Vinv=Vinv, detV=detV)
    
  }   #end weight=gLd
  
  if(weight=='med'){
    
    if(cov.type=='ratio'){
      if(boot.lme==TRUE){
        cov.med= hosking$cov.med.rat
      }else{
        cov.med= cov.lme$cov.med.rat
      }
    }else if(cov.type=='lambda'){
      if(boot.lme==TRUE){
        cov.med= hosking$cov.med.lam
      }else{
        cov.med= cov.lme$cov.med.lam
      }
    }
    
    Vinv.med= solve(cov.med)
    detV.med= det(cov.med)
    
    Vinv=Vinv.med; detV=detV.med
    
    zw$cov.lm=list(Vinv=Vinv.med, detV=detV.med)
    
  }   # end if med
  
  if( weight == 'gLd' | weight =='med'){
    
    normal=list()
    normal$prob.mtx=matrix(NA, nrow=numk, ncol=2)
    
    normal= com.prdist(data=data, numk=numk, kfix=kfix, 
                       Vinv=Vinv, detV=detV, numom=numom, 
                       weight=weight, trim=trim, 
                       cov.type=cov.type)  
    
    for (ip in 1:numk){
      if( any( is.na(kfix[[ip]]$mle[1:3]) ) ) {
        mle3[ip,1:3]= NA
        normal$prob.mtx[ip,1:2] = 0.0
      }else{
        mle3[ip,1:3]= kfix[[ip]]$mle[1:3]  }            # end if
    }    # end for
    
  }     #end if weight =gLd or med
  
  if(weight=='fcv'){
    
    pertr=hosking$pertr  
    
    fcvp =fcv.kfix.new(data,pertr,mle3=mle3,numk,kpar, cov2, 
                       weight,type) 
  }
  
  zw$aic= aic
  zw$mle3= mle3
  zw$kfix= kfix
  if(varcom==TRUE) {zw$cov2 = cov2} else {zw$cov2 = NA}
  zw$prob.mtx= matrix(NA,numk, 3)
  
  if( weight == 'gLd' | weight =='med'){
    zw$prob.mtx= normal$prob.mtx
    zw$gdd= normal$gdd
  }else if(weight=='fcv') {
    zw$prob.mtx[1:numk,3]= exp( fcvp$prob[1:numk] )
    zw$cvp = fcvp$cvp
  }
  return(zw)
}
#----------------------------------------------------------
com.prdist = function(data=NULL, numk=NULL, kfix=NULL, 
                      Vinv=NULL, detV=NULL, numom=NULL, 
                      weight=NULL, trim=NULL, cov.type="lambda")
{
#  weight = hosking$weight
  
  ww1=list(); kfix.lm=list()
  dd=matrix(NA, nrow=numom, ncol=2)
  gdd=matrix(NA, nrow=numk, ncol=2);  simple =rep(NA, numk)
  prob.mtx=matrix(NA, nrow=numk, ncol=2);
  
  #medata= median(data); 
  re=sort(data); nsam=length(data)
  ntrim=trim
  
  lm = lmoms(data,nmom=numom)
  #  lmtrim= lmoms(re[(ntrim+1):nsam], nmom=3)
  #  lmtrim = TLmoms(data,leftrim=trim, nmom=3)
  
  for (ip in 1:numk){
    
    if(kfix[[ip]]$conv != 0) {
      prob.mtx[ip,1:2]=0
      print(paste("No convergence in gev.xifix, ip=",ip))
      
    } else {
      if(kfix[[ip]]$mle[3] <=  -1.0) kfix[[ip]]$mle[3]= -0.999
      
      simple[ip] = kfix[[ip]]$mle[3]
      
      savek =vec2par(kfix[[ip]]$mle,'gev')
      kfix.lm[[ip]]= lmomgev(savek)       # mle[3] > -1 
      
      #      kfix.lm[[ip]]= theoTLmoms(savek, leftrim=trim, nmom=3) 
    }
  }  # end for ip
  
  # lmomgev(vec2par(c(100,30,-0.3),'gev'), leftrim=1)
  
  if(weight =='gLd'){
    
    for (ip in 1:numk){ 
      if(kfix[[ip]]$conv == 0){
        
        # xipar= kfix[[ip]]$mle[3]
        # if(xipar <= -0.4) {ntrim=trim+2
        # }else if(xipar <= -0.3) {ntrim=trim+1
        # }else{ ntrim= trim}
        # 
        lmtrim.cvt= lmoms(re[(ntrim+1):nsam], nmom=3)
        
        dd[1,1]= lmtrim.cvt$lambdas[1]- kfix.lm[[ip]]$lambdas[1] 
        
        if(cov.type=='ratio'){
          dd[2:numom,1]= lmtrim.cvt$ratios[2:numom]- kfix.lm[[ip]]$ratios[2:numom]
          
        }else if(cov.type=='lambda'){
          dd[2,1]= lmtrim.cvt$lambdas[2]- kfix.lm[[ip]]$lambdas[2]
          dd[3,1]= lmtrim.cvt$lambdas[3]- kfix.lm[[ip]]$lambdas[3]
        }
        
        gdd[ip,1]= t(dd[1:numom,1])%*% Vinv %*% dd[1:numom,1]
        
        prob.mtx[ip,1]= exp(-gdd[ip,1]/2 )/( (2*pi)^(numom/2) * sqrt(detV) )
        
        if(is.na(prob.mtx[ip,1]) )  prob.mtx[ip,1]= 0
      }
    } #end for ip
  } #end if weight
  
  if(weight=='med') {
    
    for (ip in 1:numk){
      if( kfix[[ip]]$conv == 0){
        
        # xipar= kfix[[ip]]$mle[3]
        # if(xipar <= -0.4) {ntrim=trim+2
        # }else if(xipar <= -0.3) {ntrim=trim+1
        # }else{ ntrim= trim}
        # 
        #         lmtrim= lmoms(re[(ntrim+1):nsam], nmom=3)
        
        medata = median(re[(ntrim+1):nsam])
        
        lmtrim.cvt= lmoms(re[(ntrim+1):nsam], nmom=3)
        
        savek= vec2par(kfix[[ip]]$mle,'gev')
        dd[1,2]= medata -quagev(0.5, savek)
        
        if(cov.type=='ratio'){
          dd[2:numom,2]= lmtrim.cvt$ratios[2:numom]- kfix.lm[[ip]]$ratios[2:numom]
          
        }else if(cov.type=='lambda'){
          dd[2,2]= lmtrim.cvt$lambdas[2]- kfix.lm[[ip]]$lambdas[2]
          dd[3,2]= lmtrim.cvt$lambdas[3]- kfix.lm[[ip]]$lambdas[3]
        }
        
        gdd[ip,2]= t(dd[1:numom,2])%*% Vinv %*% dd[1:numom,2]
        
        prob.mtx[ip,2]= exp(-gdd[ip,2]/2 ) /( (2*pi)^(numom/2) * sqrt(detV) )
        
        if(is.na( prob.mtx[ip,2]) )  prob.mtx[ip,2]=0
      } #end if 
    }  #end for ip =1, numk
  }   #end if med
  
  ww1$prob.mtx = prob.mtx
  ww1$gdd= gdd
  return(ww1)
}
#-------------------------------------------------------------
#------------------------------------------------------------
new.kpar2 = function(wtgd=NULL, numk=NULL, kpar=NULL,
                    remove=0.004, dist=0.02, tre1=0.01, tre8=0.5){
  
  wtgd[which(is.na(wtgd))] =0
  zp=list()
  aw=0
#  dist=0.03
  
  nid = which(wtgd < remove); id = which(wtgd >= remove)
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
    
#    tre1=0.01; tre8=0.5
    
    work1= work2= work3= 0

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
      work1=20
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
      work2=20
    }
    
    if(work1==1 & wtgd[numk] >= tre1){
      kpar2= c(kpar2,rep(NA,2))
      kpar2[numk+3]= min(kpar2[numk+2]+dist, 0.49)
      kpar2[numk+4]= min(kpar2[numk+3]+dist, 0.49)
      work2=2
    }
    
    lek = length(kpar2)
    mxid= which.max(wtgd)
    if(wtgd[mxid] > tre8 & mxid != 1 & mxid != numk){
      kpar2 = c(kpar2,rep(NA,4))
      kpar2[lek+1]= max(kpar[mxid]-dist/3, -0.99)
      kpar2[lek+2]= max(kpar[mxid]-2*dist/3, -0.99)
      kpar2[lek+3]= min(kpar[mxid]+dist/3, 0.49)
      kpar2[lek+4]= min(kpar[mxid]+2*dist/3, 0.49)
      work3=4
    }
    aw= aw+ work1+work2 + work3
    
  } #end if numk =1
  
  numk=length(kpar2)
  
  zp$kpar2= sort(kpar2)
  zp$aw=aw
  zp$numk=numk
  return(zp)
}