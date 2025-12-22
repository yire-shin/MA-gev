
#-----------------------------------------------------------------
K_sensitivity_plot = function(data=NULL, q.cut=0.6, mink=4,
                              maxk=20, quant=c(0.99,0.995)){
  
  totk= maxk-mink+1
  test=list()
  numq=length(quant)
  diff=rt=se= matrix(NA,maxk,numq)
  
  qq= quant
  id12= seq(1,numq)
  
  bma= FALSE
  weight="like1"
  start="mle"
  
  for (k in mink:maxk){
    
    test[[k]] = ma.gev.open(data=data, quant=qq, weight=weight,
                            numk=k, B=200, varcom=TRUE, trim=0, start=start,
                            fig=FALSE, bma=bma, pen="norm", remle=FALSE)
    
    if(bma != TRUE){
      rt[k,id12]= test[[k]]$zp.ma[id12]
      se[k,id12]= test[[k]]$adj.se.ma[id12]
      
    }else if(bma == TRUE){
      rt[k,id12]= test[[k]]$bma[id12]
      se[k,id12]= test[[k]]$pred.se.bma[id12]
    }
    
    cat("k, zp.ma, se=", k, round(rt[k,id12],2), round(se[k,id12],2),"\n")
    
  } # end for k
  
  tr.diff=rep(NA,maxk)
  for (i in (mink+1):(maxk-1)){
    diff[i,]= abs( rt[i+1,]-rt[i,] ) + abs(rt[i,]- rt[i-1,]) 
  }
  
  nina = which(!is.na(rt[,1]))
  U =( se[,1]-min(se[,1],na.rm=T) )/( max(se[,1],na.rm=T)-min(se[,1],na.rm=T) )
  tr.se =U*(max(rt[,1],na.rm=T) -min(rt[,1],na.rm=T) ) +min(rt[,1],na.rm=T)
  
  dina= which(!is.na(diff[,1]))
  
  Ub =( diff[dina,1]-min(diff[dina,1],na.rm=T) )/
    ( max(diff[dina,1],na.rm=T)-min(diff[dina,1],na.rm=T) )
  tr.diff[(mink+1):(maxk-1)] =Ub*(max(rt[,1],na.rm=T) -min(rt[,1],na.rm=T) ) +min(rt[,1],na.rm=T)
  
  diff.cut = quantile(tr.diff[dina], probs=q.cut, na.rm=T) 
  id=  which(tr.diff <=  diff.cut )
  se.cut= quantile(tr.se[nina], probs=q.cut, na.rm=T)
  
  tid = which( tr.se <= se.cut ) 
  jd = sort(intersect(id,tid))
  selTk= min(jd)
  
  max.selTk = min(maxk, selTk+2); min.selTk = max(mink,selTk-2)
  iid= intersect( jd, seq(min.selTk, max.selTk) )
  star= iid[which.min(se[iid,1])]
  
  # png("D:Ksel_wando.png", width = 2000, height =1800, res = 350)
  
  my=min(rt[,1], na.rm=T)*.997; My=max(rt[,1], na.rm=T)*1.003
  plot(seq(mink,maxk),rt[nina,1], xlim=c(mink-2,maxk),
       ylim=c(my,My), 
       xlab="K",ylab="return level", pch=16, cex=1)
  
  axis(side=1, at=seq(mink-1,maxk,1))
  #abline(h=median(rt[nina,1],na.rm=T), lty=2, col="orange", lwd=1.5)
  
  par(new=T)
  plot(seq(mink,maxk),tr.se[nina], xlim=c(mink-2,maxk),ylim=c(my,My),
       type="l",xlab="K",ylab="return level",col="blue",
  )
  points(seq(mink,maxk),tr.se[nina], col="blue")
  axis(side=1, at=c(11,12,13,14,16,17,18,19))
  
  par(new=T)
  plot(seq(mink+1,maxk-1),tr.diff[(mink+1):(maxk-1)], 
       xlim=c(mink-2,maxk),ylim=c(my,My), 
       type="l",#lty=3,
       col="red",#lwd=2,
       xlab="K",ylab="return level")
  points(seq(mink+1,maxk-1),tr.diff[(mink+1):(maxk-1)],
         col="red", pch=2, cex=0.8)
  
  abline(h=diff.cut, lty=3,lwd=2, col="darkorange")
  abline(h=se.cut,lty=3,lwd=2, col="darkgreen")
  abline(v=star,col="purple")
  
  
  text(2.7, se.cut*1.0012 ,"q.6(SE_K)", cex=.8, col="darkgreen")
  text(2.5, diff.cut*1.0012, "q.6(d_K)", col="darkorange", cex=.8)
  axis(side=1, at=c(4,5,6,7,8,9,12,14))
  
  # dev.off()
}
#-------------------------------------------------