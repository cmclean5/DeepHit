##-----------------------------------------------------------------------------------------------------------------------------------
# get pseudo conditional survival probabilities 
# t is the survival time
# d is the censoring indicator
# qt is a vector of time points that are used to divide the time interval
# output has subject id (id) and time points (s) and pseudo conditional survival probabilities (pseudost) for subject=id and at time s
##------------------------------------------------------------------------------------------------------------------------------------
getPseudoConditional <- function(t, d, qt){
  ##browser()
  s <- c(0, qt)  
  n=length(t)
  ns=length(s)-1  ## the number of intervals

  ## Example i'th patient survival time is 8.3, and time intervals are:
  ## s=c(0,1,2,5,8,10)
  ## (0,1], (1,2], (2,5], (5,8], (8,10]
    
  ## For each uncensored patient, i.e. censor indicator == 1, define which time interval
  ## they belong too:
  ## D [number of patients, number of time points]
  ## D[i,] = [0,0,0,0,1]  
  D <- do.call(cbind, lapply(1:ns, function(j)  (s[j] < t) * (t <= s[j+1]) * (d == 1)))

  ## Define number of patients (censored or uncensored) at risk of experiencing
  ## an event at each time point:
  ## R [number of patients, number of time points]
  ## R[i,j] == 1 if patient 'i' still at risk at time point j,
  ## R[i,j] == 0 if patient 'i' already experienced the event before time point j.
  ## R[i,] = [1,1,1,1,1]
  R <- do.call(cbind, lapply(1:ns, function(j) ifelse(s[j] < t, 1, 0)))

  ## 'R*Delta' splits each patients survival-time into the 'ns' time intervals,
  ## i.e. running rowSum(R*Delta) will return each patients survival time.
  ## Delta [number of patients, number of time points]
  ## Delta[i,] = [1,1,3,3,0.3]  
  Delta   <- do.call(cbind, lapply(1:ns, function(j) pmin(t,s[j+1])-s[j]))
  Delta.t <- R*Delta
    
  ## format into long formate
  ##dd.tmp=cbind.data.frame(id=rep(1:n,ns),s=rep(c(0,qt[-length(qt)]), each=n),
  ##                        y=c(R*Delta),d=c(D))
  dd.tmp=cbind.data.frame(id=rep(1:n,ns),s=rep(c(0,qt[-length(qt)]), each=n),
                          y=c(Delta.t), d=c(D))
    
  dd=dd.tmp[dd.tmp$y>0,]
  pseudost=rep(NA, nrow(dd))
  for (j in 1:ns){
    index= (dd$s==s[j])
    dds=dd[index,]
    pseudost[index]=pseudosurv(time=dds$y, event=dds$d, tmax=s[j+1]-s[j])$pseudo
    print(j)
  }
  dd$pseudost=pseudost  
  
  return(dd[,c(1,2,5)])
}

##-----------------------------------------------------------------------------------------------------------------------------------
# get pseudo conditional survival probabilities for each competing event 
# t is the survival time
# d is the competeing event set
# qt is a vector of time points that are used to divide the time interval
# cencode is the censoring code using in d
# output has subject id (id) and time points (s) and pseudo conditional survival probabilities (pseudost) for subject=id and at time s
##------------------------------------------------------------------------------------------------------------------------------------
getPseudoConditionalCI <- function(t, d, qt,cencode=0){

 ##browser()
  s <- c(0, qt)  
  n=length(t)
  ns=length(s)-1  ## the number of intervals

  ## For each uncensored patient, i.e. censor indicator == 1, define which time interval
  ## they belong too:
  n.events   = get.comp.events(events=d, cencode=cencode)
  ne         = length(n.events)
  De <- array(0,dim=c(n,ns))
  D  <- array(0,dim=c(n,ns))
  for( e in 1:ne ){    
      De <- do.call(cbind, lapply(1:ns, function(j)
                                   (s[j] < t) * (t <= s[j+1]) * (d == n.events[e])))
      De[De==1] = n.events[e]
      D = D + De
  }
      
  ## Define number of patients (censored or uncensored) at risk of experiencing
  ## an event at each time point:
  R <- do.call(cbind, lapply(1:ns, function(j) ifelse(s[j] < t, 1, 0)))

  ## 'R*Delta' splits each patients survival-time into the 'ns' time intervals,
  ## i.e. running rowSum(R*Delta) will return each patients survival time.
  Delta   <- do.call(cbind, lapply(1:ns, function(j) pmin(t,s[j+1])-s[j]))
  Delta.t <- R*Delta
    
  ## format into long formate
  dd.tmp=cbind.data.frame(id=rep(1:n,ns),s=rep(c(0,qt[-length(qt)]), each=n),
                          y=c(Delta.t), d=c(D))
    
  dd=dd.tmp[dd.tmp$y>0,]
  pseudost = array(dim=c(nrow(dd),ne))
  colnames(pseudost) = n.events  
  for (j in 1:ns){
    index= (dd$s==s[j])
    dds=dd[index,]
    tmpci  = pseudoci(time=dds$y, event=dds$d, tmax=s[j+1]-s[j])
    causes = tmpci$cause
    for( e in 1:length(causes) ){
        pseudost[index,causes[e]]=tmpci$pseudo[[e]]
    }
    print(j)
  }
    
  dd = cbind(dd,pseudost)   
  return(dd[,c(1,2,match(n.events,colnames(dd)))])
}

check.tmin <- function(t.min,tau,print=FALSE){
    
    if( tau < t.min ){

        if(print){
            cat("warning: tau = ", tau, " < min obs.time = ", t.min,"\n")
            cat("       : tau set to ", t.min,"\n")
        }
        tau = t.min;
    }
    return(tau)
}


##-----------------------------------------------------------------------------------------------------------------------------------
## The three International Cancer Survival Standard (ICSS) weights used for
## age-standardization of net survival.
## Ref: https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-015-0057-3#Sec2
##-----------------------------------------------------------------------------------------------------------------------------------
compute_weight_stratified <- function(x,time,delta,
                                      x.strat=c(15,45,55,65,75,150),
                                      x.we=c("1"=0.07,"2"=0.12,"3"=0.23,"4"=0.29,"5"=0.29),
                                      method=c("custom")){

    method=match.arg(method);
    
    x.group=group.data(x=x,nbins=x.strat, method=method)
    x.group=prop.table(x.group$groups.tb)
    names(x.group) = names(x.we)

    n=length(time)
    pseudo <- data.frame(id=1:n,t=time,delta=delta)
 
    ## sort in time, if tied, put events before censoring
    ## order time lowest to highest, with events (i.e. delta==1) coming first if tied
    ## time, done by negating the values in vector delta  
    pseudo <- pseudo[order(pseudo$t,-pseudo$delta),]
  
    t=pseudo$t
    delta=pseudo$delta
  
    s=sort(unique(t[delta==1]))
 
    tmp  = do.call(rbind, lapply(1:n, function(j)  c(0,s) + x[j]))
    tmp2 = apply(tmp ,2,function(x) group.data(x=x,nbins=x.strat, method=method)[[2]])
    df   = apply(tmp2,2,function(x) x.we[x]/x.group[x])
    df   = df[,-1]

    return(df)
    
}

compute_unweight <- function(time,delta){                             

    n=length(time)
    pseudo <- data.frame(id=1:n,t=time,delta=delta)
 
    ## sort in time, if tied, put events before censoring
    ## order time lowest to highest, with events (i.e. delta==1) coming first if tied
    ## time, done by negating the values in vector delta  
    pseudo <- pseudo[order(pseudo$t,-pseudo$delta),]
  
    t=pseudo$t
    delta=pseudo$delta
  
    s=sort(unique(t[delta==1]))
    ns=length(s)
    
    df = matrix(1,nrow=n, ncol=ns)

    return(df)
    
}



##-----------------------------------------------------------------------------------------------------------------------------------
## Compute IPCW weights using random survival forest
## input are survival time (time), censoring indicator (delta) and covariate matrix
## output is a n*ns matrix containing IPCW weights ordered by unique event times, where n is the number of subjects 
##         and ns is the number of unique event times.
## For details on RandomForest see [1]: https://luminwin.github.io/randomForestSRC/reference/rfsrc.html
##-----------------------------------------------------------------------------------------------------------------------------------
compute_weight_rf_trunc <- function(time,delta,cov){
  
  n=length(time)
  ddw <- data.frame(t=time,d=1-delta,cov)
  res=rfsrc(Surv(t, d) ~ ., data=ddw)
  t.unique=res$time.interest      ## time.interest == Ordered unique death times.
  s=sort(unique(time[delta==1]))  ## unique event times
  ##survprob=t(sapply(1:n,function(j) approx(t.unique, foo$survival[j,], xout=s, method="constant", yleft=1)$y )) # n*ns dimension
  survprob=t(sapply(1:n,function(j) approx(t.unique, res$survival[j,], xout=s, method="constant", rule=2:2)$y )) # n*ns dimension
  trunc=max(1e-20,min(survprob[survprob>0])) # avoid close to zero prob
  weight=ifelse(survprob < trunc, 1/trunc, 1/survprob)                        
  return(weight)

}

##-----------------------------------------------------------------------------------------------------------------------------------
## Compute IPCW pseudo RMST at a single tau using Nelson-Aalen method; tau < the last event time
## input are survival time (time), censoring indicator (delta), 
## weight is a n*ns matrix, where n is the number of subjects and ns is the number of unique event times, and these weights
## will not change in leave-one-out part (see fourmula (9) in Binder:2014.)
## tau is the landmark time
## output is a vector of n pseudo values 
##-----------------------------------------------------------------------------------------------------------------------------------
pseudomean_ipcw <- function(time,delta,weight,tau,SMALL=1e-20){

  ## preparing the data
  n=length(time)
  pseudo <- data.frame(id=1:n,t=time,delta=delta,weight)

  tau   = check.tmin(t.min=min(time),tau=tau);
     
  ## sort in time, if tied, put events before censoring
  ## order time lowest to highest, with events (i.e. delat==1) coming first if tied
  ## time, done by negating the values in vector delta  
  pseudo <- pseudo[order(pseudo$t,-pseudo$delta),]
  
  t=pseudo$t
  delta=pseudo$delta
  weight=pseudo[,-c(1,2,3)]
  
  s=sort(unique(t[delta==1]))

  ns=length(s)  # the number of intervals
  D <- do.call(cbind, lapply(1:ns, function(j)  (s[j] == t)*(delta == 1)))## counting the number of any event, i.e. deaths, at each unique time point. 
  Y <- do.call(cbind, lapply(1:ns, function(j)  ifelse(s[j] <= t, 1, 0))) ## counting the number of patients at risk at each unique time point.

  inx=max(which(s<=tau))
  ttmp=c(0,s)
  tt=c(ttmp[ttmp<=tau],tau) ## add one extra column, may repeat, but diff=0
  dt=diff(tt) ## delta t (the interval of time) found by taking the difference of tt, i.e. diff[1] = tt[2] - tt[1],   

    ## xx=cuminc(ftime=survT$obs.time, fstatus=data$eventPrim, cencode=0)
    ## xx=timepoints(xx,times=s)

    ## ---------------------------
    ## Reference: https://mspace.lib.umanitoba.ca/bitstream/handle/1993/21696/galloway_katie.pdf?sequence=1&isAllowed=y
    ##  Estimated the survival function, S(t), using the Kaplan-Meier
    ##  estimator: S(t) = Pi_{tj <= t} (Nrisk_tj - Ndead_tj)/Nrisk_tj
    ##
    ##  The Nelson-Aalen estimator of the cumulative hazard function:
    ##  Lambda(t) = Sum_{tj <= t} Ndead_tj/Nrisk_tj
    ##  The Nelson-Aalen estimator (of survival) is related to the
    ##  Kaplan-Meier estimator by: Lambda(t) = -ln S(t)
    ##  Hence,  S(t) = exp( -Lambda(t) )
    ## Int S^W(t) dt = exp(- Lambda(t)^W ) * dt
    ##               = exp(- Sum_{tj <= t} Ndead_tj*W / Nrisk_tj*W) * dt
    ##               = exp(-cumsum(Ndead*W/Nrisk*W)) * dt
    ##               = sum( exp(-cumsum(Ndead*W/Nrisk*W)) * dt)
    ## ---------------------------
  denominator=colSums(Y*weight)
  ##denominator=apply(denominator, 2, function(x) ifelse(x==0, SMALL, x))  
  numerator=colSums(D*weight)
  IPCW_CH=cumsum(numerator/denominator)
  IPCW_surv=exp(-IPCW_CH)
  surv=c(IPCW_surv[1:inx],IPCW_surv[inx])
  IPCW_RM=sum(surv*dt)
  
  Yw=Y*weight
  Denominator=matrix(colSums(Yw),n,ns,byrow=TRUE)-Yw
  Denominator=apply(Denominator, 2, function(x) ifelse(x==0, SMALL, x))  
  Dw=D*weight
  Numerator=matrix(colSums(Dw),n,ns,byrow=TRUE)-Dw
  IPCW_CHi=t(apply(Numerator/Denominator,1,cumsum))
  IPCW_survi=exp(-IPCW_CHi)
  dt.mat=matrix(dt,nrow=n,ncol=length(dt),byrow=TRUE)
  survi=cbind(IPCW_survi[,1:inx],IPCW_survi[,inx])
  IPCW_RMi=rowSums(survi*dt.mat)

    ## weighted pseudo conditional probabilities survival for each patient
    ## at time t_j+1 = tau, conditioned on the mean number of patients at risk, i.e. n.
    ## Si^W(t_j+1 | R_j) = R_j * S^W(t_j+1 | R_j)
    ##                     - (R_j - 1) * S^W-i(t_j+1 | R_j)
    ##                   = n * IPCW_RM - (n-1) * IPCW_RMi
    
  pseudo_mean=n*IPCW_RM-(n-1)*IPCW_RMi
  ## back to original order
  pseudo <- as.vector(pseudo_mean[order(pseudo$id)])		##back to original order
  return(pseudo)
}


##-----------------------------------------------------------------------------------------------------------------------------------
## Compute pseudo RMST (the mean, or expected, survival time) at a single tau using Nelson-Aalen method; tau < the last event time
## input are survival time (time), censoring indicator (delta), 
## where n is the number of subjects and ns is the number of unique event times,
## and these weights will not change in leave-one-out part
## (see fourmula (9) in Binder:2014.)
## tau is the landmark time
## output is a vector of n pseudo values 
##-----------------------------------------------------------------------------------------------------------------------------------
pseudomean <- function(time,delta,tau,SMALL=1e-20){

  ## preparing the data
  n=length(time)
  pseudo <- data.frame(id=1:n,t=time,delta=delta)

  tau = check.tmin(t.min=min(time),tau=tau);
      
  ## sort in time, if tied, put events before censoring
  ## order time lowest to highest, with events (i.e. delat==1) coming first if tied
  ## time, done by negating the values in vector delta  
  pseudo <- pseudo[order(pseudo$t,-pseudo$delta),]
  
  t=pseudo$t
  delta=pseudo$delta
  
  s=sort(unique(t[delta==1]))

  ns=length(s)  ## the number of intervals
  D <- do.call(cbind, lapply(1:ns, function(j)  (s[j] == t)*(delta == 1)))
  Y <- do.call(cbind, lapply(1:ns, function(j)  ifelse(s[j] <= t, 1, 0)))
  inx=max(which(s<=tau))
  ttmp=c(0,s)
  tt=c(ttmp[ttmp<=tau],tau)
  dt=diff(tt)

  denominator=colSums(Y)
  numerator=colSums(D)
  CH=cumsum(numerator/denominator)
  surv=exp(-CH)
  surv=c(surv[1:inx],surv[inx])
  RM=sum(surv*dt)
  
  Denominator=matrix(colSums(Y),n,ns,byrow=TRUE)-Y
  Denominator=apply(Denominator, 2, function(x) ifelse(x==0, SMALL, x))  
  Numerator=matrix(colSums(D),n,ns,byrow=TRUE)-D
  CHi=t(apply(Numerator/Denominator,1,cumsum))
  survi=exp(-CHi)
  dt.mat=matrix(dt,nrow=n,ncol=length(dt),byrow=TRUE)
  survi=cbind(survi[,1:inx],survi[,inx])
  RMi=rowSums(survi*dt.mat)

  pseudo_mean=n*RM-(n-1)*RMi
  ## back to original order
  pseudo <- as.vector(pseudo_mean[order(pseudo$id)])		##back to original order
  return(pseudo)

}

##--------------------------------------------------------------------------------------
## Return order list of competeing event ids, given event set
##-------------------------------------------------------------------------------------
#get.comp.events <- function(events, cencode=0, print=FALSE){
#    ## number of competing events
#    n.events = unique(events)
#    n.events = n.events[order(n.events)]
#    n.events = n.events[-which(n.events==cencode)]
#    ne       = length(n.events)  
#    if( print ){ cat("no: competing events: ", ne,"\n"); }
#    return(n.events)
#}

##-----------------------------------------------------------------------------------------------------------------------------------
## Compute pseudo RMST for competing events at a single tau using Nelson-Aalen method;
## tau < the last event time
## input are survival time (time), censoring indicator (delta), competng events indicator (events) 
## where n is the number of subjects and ns is the number of unique event times,
## and these weights will not change in leave-one-out part
## (see fourmula (9) in Binder:2014.)
## tau is the landmark time
## output is a vector of n pseudo values 
## Reference: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3396320/
##-----------------------------------------------------------------------------------------------------------------------------------
pseudomean_comp <- function(time,delta,events,cencode=0,tau,SMALL=1e-20,
                            method=c("mean","risk")){

  method = match.arg(method)
    
  ## number of competing events
  n.events = get.comp.events(events=events, cencode=cencode)
  ne       = length(n.events)  
  #cat("no: competing events: ", ne,"\n")
    
  ## preparing the data
  n=length(time)
  pseudo <- data.frame(id=1:n,t=time,delta=delta,events=events)

  
  ## sort in time, if tied, put events before censoring
  ## order time lowest to highest, with events (i.e. delat==1) coming first if tied
  ## time, done by negating the values in vector delta  
  pseudo <- pseudo[order(pseudo$t,-pseudo$delta),]
  
  t=pseudo$t
  delta=pseudo$delta
  events=pseudo$events  
  
  s=sort(unique(t[delta==1]))

  tau = check.tmin(t.min=min(time),tau=tau);
  
    
  ns=length(s)  # the number of intervals
    
  Y    <- do.call(cbind, lapply(1:ns, function(j)  ifelse(s[j] <= t, 1, 0)))
  Ek   <- array(0,dim=c(n,ns,ne))
  Esum <- array(0,dim=c(n,ns))
  RMk  <- array(dim=c(ne)) 
  for( e in 1:ne ){
      Ek[,,e] = do.call(cbind, lapply(1:ns,
                                    function(j)  (s[j] == t)*(events == n.events[e])))
      Esum    = Esum + Ek[,,e]
  }

  inx=max(which(s<=tau))
  ttmp=c(0,s)
  tt=c(ttmp[ttmp<=tau],tau)
  dt=diff(tt)

  ## set the number of induviduals
  if( method == "mean" ){ Rj = n; }
  else {                  Rj = colSums(Y)[inx]; }  
  ##Rj = colSums(Y)[inx];
    
  CHkm = cumsum(colSums(Esum)/colSums(Y))
  surv = exp(-CHkm)
  KM   = shift(surv,1,1)
    
  for( e in 1:ne ){
      num    = colSums(Ek[,,e])
      dem    = colSums(Y)
      CHk    = cumsum( (num/dem)*KM )
      survk  = CHk
      survk  = c(survk[1:inx],survk[inx])
      RMk[e] = sum(survk*dt)
  }

  Numerator  = matrix(colSums(Esum),n,ns,byrow=TRUE)-Esum
  Denominator= matrix(colSums(Y),n,ns,byrow=TRUE)-Y
  Denominator= apply(Denominator, 2, function(x) ifelse(x==0, SMALL, x))  
 
  dt.mat = matrix(dt,nrow=n,ncol=length(dt),byrow=TRUE)  
  CHkmi  = t(apply(Numerator/Denominator,1,cumsum))
  survi  = exp(-CHkmi)
  KMi    = t(apply(survi,1,function(x) shift(x,1,1)))
           
  RMik        = array(dim=c(n,ne))
  pseudo_mean = array(dim=c(n,ne))
  colnames(pseudo_mean) = n.events
    
  for( e in 1:ne ){
      Numk     = matrix(colSums(Ek[,,e]),n,ns,byrow=TRUE)-Ek[,,e]
      CHik     = t(apply( (Numk/Denominator)*KMi,1,cumsum))
      survik   = CHik
      survik   = c(survik[,1:inx],survik[,inx])
      RMik[,e] = rowSums(survik*dt.mat)
      
      pseudo_mean[,e] = Rj*RMk[e]-(Rj-1)*RMik[,e]
      pseudo_mean[,e] = as.vector(pseudo_mean[order(pseudo$id),e])
  }
    
    return(pseudo_mean)

}

pseudomean_comp_weight<- function(time,delta,events,cencode=0,tau,weights,
                                  SMALL=1e-20,
                                  method=c("mean","risk")){                       

    method = match.arg(method)
    
  ## number of competing events
  n.events = get.comp.events(events=events, cencode=cencode)
  ne       = length(n.events)  
  #cat("no: competing events: ", ne,"\n")
    
  ## preparing the data
  n=length(time)
  pseudo <- data.frame(id=1:n,t=time,delta=delta,events=events)

  
  ## sort in time, if tied, put events before censoring
  ## order time lowest to highest, with events (i.e. delat==1) coming first if tied
  ## time, done by negating the values in vector delta  
  reorder <- order(pseudo$t,-pseudo$delta)
  pseudo  <- pseudo[reorder,]
  for( e in 1:ne ){
      weights[,,e] = weights[reorder,,e]
      weights[,,e] = apply(weights[,,e], 2, function(x) ifelse(x==0, SMALL, x))  
  }
  
  t=pseudo$t
  delta=pseudo$delta
  events=pseudo$events  
  
  s=sort(unique(t[delta==1]))

 tau = check.tmin(t.min=min(time),tau=tau)   
    
  ns=length(s)  # the number of intervals
    
  Y     <- do.call(cbind, lapply(1:ns, function(j)  ifelse(s[j] <= t, 1, 0)))
  Ek    <- array(dim=c(n,ns,ne))
  Esum  <- array(0,dim=c(n,ns))
  RMk   <- array(dim=c(ne)) 
  for( e in 1:ne ){
      Ek[,,e] = do.call(cbind, lapply(1:ns,
                                    function(j)  (s[j] == t)*(events == n.events[e])))
      Esum    = Esum + ( (Ek[,,e]*weights[,,e])/weights[,,e] )
  }

  inx=max(which(s<=tau))
  ttmp=c(0,s)
  tt=c(ttmp[ttmp<=tau],tau)
  dt=diff(tt)

  ## set the number of induviduals
  if( method == "mean" ){ Rj = n; }
  else {                  Rj = colSums(Y)[inx]; }  
  ##Rj = n;

  CHkm = cumsum(colSums(Esum)/colSums(Y))
  surv = exp(-CHkm)
  KM   = shift(surv,1,1)
    
  for( e in 1:ne ){
      num    = colSums(Ek[,,e]*weights[,,e])
      dem    = colSums(Y*weights[,,e])
      CHk    = cumsum( (num/dem)*KM )
      survk  = CHk
      survk  = c(survk[1:inx],survk[inx])
      RMk[e] = sum(survk*dt)
  }

  Numerator  = matrix(colSums(Esum),n,ns,byrow=TRUE)-Esum
  Denominator= matrix(colSums(Y),n,ns,byrow=TRUE)-Y
  Denominator= apply(Denominator, 2, function(x) ifelse(x==0, SMALL, x))  
    
  dt.mat = matrix(dt,nrow=n,ncol=length(dt),byrow=TRUE)  
  CHkmi  = t(apply(Numerator/Denominator,1,cumsum))
  survi  = exp(-CHkmi)
  KMi    = t(apply(survi,1,function(x) shift(x,1,1)))
           
  RMik        = array(dim=c(n,ne))
  pseudo_mean = array(dim=c(n,ne))
  colnames(pseudo_mean) = n.events

  for( e in 1:ne ){
      Ewei     = Ek[,,e]*weights[,,e]
      Ywei     = Y*weights[,,e]
      num      = matrix(colSums(Ewei),n,ns,byrow=TRUE)-Ewei 
      dem      = matrix(colSums(Ywei),n,ns,byrow=TRUE)-Ywei 
      CHik     = t(apply( (num/dem)*KMi,1,cumsum))
      survik   = CHik
      survik   = c(survik[,1:inx],survik[,inx])
      RMik[,e] = rowSums(survik*dt.mat)
      
      pseudo_mean[,e] = Rj*RMk[e]-(Rj-1)*RMik[,e]
      pseudo_mean[,e] = as.vector(pseudo_mean[order(pseudo$id),e])
  }
    
    return(pseudo_mean)

}


##---------------------------------------
## Calculate the expected probabilities for each competing risk type
##---------------------------------------

#---------------------------------------
## Kaplan-Meier estimator
##---------------------------------------
pseudoKM <- function(time,delta,events,cencode=0,tau,SMALL=1e-20,
                     method=c("mean","risk")){

  method = match.arg(method)
    
  ## number of competing events
  n.events = get.comp.events(events=events, cencode=cencode)
  ne       = length(n.events)  
    
  ## preparing the data
  n=length(time)
  pseudo <- data.frame(id=1:n,t=time,delta=delta,events=events)

  
  ## sort in time, if tied, put events before censoring
  ## order time lowest to highest, with events (i.e. delat==1) coming first if tied
  ## time, done by negating the values in vector delta  
  pseudo <- pseudo[order(pseudo$t,-pseudo$delta),]
  
  t=pseudo$t
  delta=pseudo$delta
  events=pseudo$events  
  
  s=sort(unique(t[delta==1]))

  tau = check.tmin(t.min=min(time),tau=tau);
      
  ns=length(s)  # the number of intervals

  D    <- do.call(cbind, lapply(1:ns, function(j)  (s[j] == t)*(delta == 1)))
  Y    <- do.call(cbind, lapply(1:ns, function(j)  ifelse(s[j] <= t, 1, 0)))

  inx  = max(which(ttmp<=tau))
  ttmp = c(0,s)  
  ds   = diff(ttmp)

  if( inx > length(ds) ){ inx = length(ds); }
    
  ## set the number of induviduals
  if( method == "mean" ){ Rj = n; }
  else {                  Rj = colSums(Y)[inx]; }  

  ## Reference: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3135581/
  ## See eqn 2.
  ## BMC Med Res Methodol. 2011; 11: 86.
  ## Published online 2011 Jun 3. doi: 10.1186/1471-2288-11-86
  ## Understanding competing risks: a simulation point of view.
  ## Arthur Allignol
  KM   = cumsum(colSums(D)/colSums(Y))
  surv = exp(-cumsum(KM*ds))
  RM   = surv[inx]
  ##RM   = sum(surv[1:inx]*ds[1:inx])

  Numerator   = matrix(colSums(D),n,ns,byrow=TRUE)-D
  Denominator = matrix(colSums(Y),n,ns,byrow=TRUE)-Y
  Denominator = apply(Denominator, 2, function(x) ifelse(x==0, SMALL, x))  
 
  ds.mat      = matrix(ds,nrow=n,ncol=length(ds),byrow=TRUE)  
  KMi         = t(apply(Numerator/Denominator,1,cumsum))
  survi       = exp(-t(apply(KMi*ds.mat,1,cumsum)))
  RMi         = survi[,inx]  
  ##RMi         = rowSums(survi[,1:inx]*ds.mat[,1:inx])  

  pseudo_mean = Rj*RM-(Rj-1)*RMi
  pseudo_mean = as.vector(pseudo_mean[order(pseudo$id)])    
    
  return(pseudo_mean)

}

##---------------------------------------
## Aalen-Johannsen estimator
##---------------------------------------
pseudoAJ_comp <- function(time,delta,events,cencode=0,tau,SMALL=1e-20,
                            method=c("mean","risk")){

  method = match.arg(method)
    
  ## number of competing events
  n.events = get.comp.events(events=events, cencode=cencode)
  ne       = length(n.events)  
  #cat("no: competing events: ", ne,"\n")
    
  ## preparing the data
  n=length(time)
  pseudo <- data.frame(id=1:n,t=time,delta=delta,events=events)

  
  ## sort in time, if tied, put events before censoring
  ## order time lowest to highest, with events (i.e. delat==1) coming first if tied
  ## time, done by negating the values in vector delta  
  pseudo <- pseudo[order(pseudo$t,-pseudo$delta),]
  
  t=pseudo$t
  delta=pseudo$delta
  events=pseudo$events  
  
  s=sort(unique(t[delta==1]))

  tau = check.tmin(t.min=min(time),tau=tau);
  
    
  ns=length(s)  # the number of intervals
    
  Y    <- do.call(cbind, lapply(1:ns, function(j)  ifelse(s[j] <= t, 1, 0)))
  Ek   <- array(0,dim=c(n,ns,ne))
  Esum <- array(0,dim=c(n,ns))
  AJk  <- array(dim=c(ns,ne))
  CIF  <- array(dim=c(ns,ne)) 
  RMk  <- array(dim=c(ne)) 
  for( e in 1:ne ){
      Ek[,,e] = do.call(cbind, lapply(1:ns,
                                    function(j)  (s[j] == t)*(events == n.events[e])))
      Esum    = Esum + Ek[,,e]
  }
    
  inx  = max(which(s<=tau))
  ttmp = c(0,s)
  ds   = diff(ttmp)

  if( inx > length(ds) ){ inx = length(ds); }
    
  ## set the number of induviduals
  if( method == "mean" ){ Rj = n; }
  else {                  Rj = colSums(Y)[inx]; }  
  ##Rj = colSums(Y)[inx];

  ## Reference: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3135581/
  ## See eqn 2.
  ## BMC Med Res Methodol. 2011; 11: 86.
  ## Published online 2011 Jun 3. doi: 10.1186/1471-2288-11-86
  ## Understanding competing risks: a simulation point of view.
  ## Arthur Allignol
  KM   = cumsum(colSums(Esum)/colSums(Y))
  surv = exp(-cumsum(KM*ds))
 
  for( e in 1:ne ){
      num     = colSums(Ek[,,e])
      dem     = colSums(Y)
      AJk[,e] = cumsum(num/dem)
      CIF[,e] = cumsum(surv*AJk[,e]*ds)
      RMk[e]  = CIF[,e][inx]
  }

  Numerator  = matrix(colSums(Esum),n,ns,byrow=TRUE)-Esum
  Denominator= matrix(colSums(Y),n,ns,byrow=TRUE)-Y
  Denominator= apply(Denominator, 2, function(x) ifelse(x==0, SMALL, x))  
 
  ds.mat     = matrix(ds,nrow=n,ncol=length(ds),byrow=TRUE)  
  KMi        = t(apply(Numerator/Denominator,1,cumsum))
  survi      = exp(-t(apply(KMi*ds.mat,1,cumsum))) 
            
  RMik        = array(dim=c(n,ne))
  pseudo_mean = array(dim=c(n,ne))
  colnames(pseudo_mean) = n.events

  for( e in 1:ne ){
      Numk     = matrix(colSums(Ek[,,e]),n,ns,byrow=TRUE)-Ek[,,e]
      AJki     = t(apply( (Numk/Denominator),1,cumsum))
      tmp      = do.call(rbind, lapply(1:n,
                                    function(j)
                                     cumsum(survi[j,1:inx]*AJki[j,1:inx]*ds.mat[j,1:inx]) ))
      ##tmp      = t(apply(survi[,1:inx]*AJki[,1:inx]*ds.mat[,1:inx],1,cumsum))
      RMik[,e] = tmp[,inx]
            
      pseudo_mean[,e] = Rj*RMk[e]-(Rj-1)*RMik[,e]
      pseudo_mean[,e] = as.vector(pseudo_mean[order(pseudo$id),e])
  }

    return(list(pseudo_mean=pseudo_mean, surv=surv, AJk=AJk, CIF=CIF, RMk=RMk, s=s))    
    #return(pseudo_mean)

}

##---------------------------------------
## Weighted Aalen-Johannsen estimator
##---------------------------------------
pseudoAJ_comp_weight<- function(time,delta,events,cencode=0,tau,weights,
                                SMALL=1e-20,
                                method=c("mean","risk")){                       

    method = match.arg(method)
    
  ## number of competing events
  n.events = get.comp.events(events=events, cencode=cencode)
  ne       = length(n.events)   
    
  ## preparing the data
  n=length(time)
  pseudo <- data.frame(id=1:n,t=time,delta=delta,events=events)

  
  ## sort in time, if tied, put events before censoring
  ## order time lowest to highest, with events (i.e. delat==1) coming first if tied
  ## time, done by negating the values in vector delta  
  reorder <- order(pseudo$t,-pseudo$delta)
  pseudo  <- pseudo[reorder,]
  weights <- weights[reorder,]
  weights <- apply(weights, 2, function(x) ifelse(x==0, SMALL, x))  
  
  t=pseudo$t
  delta=pseudo$delta
  events=pseudo$events  
  
  s=sort(unique(t[delta==1]))

  tau = check.tmin(t.min=min(time),tau=tau)   
    
  ns=length(s)  # the number of intervals

  Y     <- do.call(cbind, lapply(1:ns, function(j)  ifelse(s[j] <= t, weights[,j],0)))
  Ek    <- array(dim=c(n,ns,ne))
  Esum  <- array(0,dim=c(n,ns))
  AJk   <- array(dim=c(ns,ne))
  CIF  <- array(dim=c(ns,ne)) 
  RMk   <- array(dim=c(ne)) 
    for( e in 1:ne ){
        Ek[,,e] = do.call(cbind, lapply(1:ns,
                                        function(j)
                                            ifelse((s[j] == t)*(events == n.events[e]),
                                                   weights[,j], 0)
                                        ))
        Esum    = Esum + Ek[,,e]  
  }

  inx  = max(which(s<=tau))
  ttmp = c(0,s)
  ds   = diff(ttmp)

  if( inx > length(ds) ){ inx=length(ds); }
    
  ## set the number of induviduals
  if( method == "mean" ){ Rj = n; }
  else {                  Rj = colSums(Y)[inx]; }  
  

  KM   = cumsum(colSums(Esum)/colSums(Y))
  surv = exp(-cumsum(KM*ds))
 
  for( e in 1:ne ){
      num     = colSums(Ek[,,e])
      dem     = colSums(Y)
      AJk[,e] = cumsum(num/dem)
      CIF[,e] = cumsum(surv*AJk[,e]*ds)
      RMk[e]  = CIF[,e][inx]
  }

  Numerator  = matrix(colSums(Esum),n,ns,byrow=TRUE)-Esum
  Denominator= matrix(colSums(Y),n,ns,byrow=TRUE)-Y
  Denominator= apply(Denominator, 2, function(x) ifelse(x==0, SMALL, x))  
 
  ds.mat     = matrix(ds,nrow=n,ncol=length(ds),byrow=TRUE)  
  KMi        = t(apply(Numerator/Denominator,1,cumsum))
  survi      = exp(-t(apply(KMi*ds.mat,1,cumsum))) 
            
  RMik        = array(dim=c(n,ne))
  pseudo_mean = array(dim=c(n,ne))
  colnames(pseudo_mean) = n.events
    
  for( e in 1:ne ){
      Numk     = matrix(colSums(Ek[,,e]),n,ns,byrow=TRUE)-Ek[,,e]
      AJki     = t(apply( (Numk/Denominator),1,cumsum))
      tmp      = do.call(rbind, lapply(1:n,
                                    function(j)
                                     cumsum(survi[j,1:inx]*AJki[j,1:inx]*ds.mat[j,1:inx]) ))
      RMik[,e] = tmp[,inx]
      ##RMik[,e] = t(apply(survi*AJki*ds.mat,1,cumsum))[,inx]
            
      pseudo_mean[,e] = Rj*RMk[e]-(Rj-1)*RMik[,e]
      pseudo_mean[,e] = as.vector(pseudo_mean[order(pseudo$id),e])
  }
    
    return(list(pseudo_mean=pseudo_mean, surv=surv, AJk=AJk, CIF=CIF, RMk=RMk, s=s))   

}


##---------------------------------------
## Class weighted Aalen-Johannsen estimator
##---------------------------------------
pseudoAJ_comp_class_we <- function(time,delta,events,cencode=0,tau,
                                   SMALL=1e-20,
                                   method=c("mean","risk")){                       

  method = match.arg(method)
    
  ## number of competing events
  n.events = get.comp.events(events=events, cencode=cencode)
  ne       = length(n.events)   
    
  ## number of classes
  class          = c(cencode,n.events)  
  nclass         = length(class)
  weights        = rep(1,nclass)
  names(weights) = class  

    
  ## do we have class names matching event type names 
  #if( all(names(alpha) %in% c(cencode,n.events)) ){
  #    cat("> class names match event type names!\n")
  #} else {
  #    cat("> class names don't match event type names!\n")
  #} 

  ## patient weights
  alpha = 1-prop.table(table(events[time<=tau]))    
  weights[match(names(alpha),names(weights))] = alpha
  weights = weights[events]
    
  ## preparing the data
  n=length(time)
  pseudo <- data.frame(id=1:n,t=time,delta=delta,events=events)
  
  ## sort in time, if tied, put events before censoring
  ## order time lowest to highest, with events (i.e. delat==1) coming first if tied
  ## time, done by negating the values in vector delta  
  reorder <- order(pseudo$t,-pseudo$delta)
  pseudo  <- pseudo[reorder,]
  
  t=pseudo$t
  delta=pseudo$delta
  events=pseudo$events  
  
  s=sort(unique(t[delta==1]))

  tau = check.tmin(t.min=min(time),tau=tau)   
    
  ns=length(s)  # the number of intervals

  Y     <- do.call(cbind, lapply(1:ns, function(j)  ifelse(s[j] <= t, weights,0)))
  Ek    <- array(dim=c(n,ns,ne))
  Esum  <- array(0,dim=c(n,ns))
  AJk  <- array(dim=c(ns,ne))
  CIF  <- array(dim=c(ns,ne)) 
  RMk   <- array(dim=c(ne)) 
    for( e in 1:ne ){
        Ek[,,e] = do.call(cbind, lapply(1:ns,
                                        function(j)
                                            ifelse((s[j] == t)*(events == n.events[e]),
                                                   weights, 0)
                                        ))
        Esum    = Esum + Ek[,,e]  
  }
    
  inx  = max(which(s<=tau))
  ttmp = c(0,s)
  ds   = diff(ttmp)

  if( inx > length(ds) ){ inx=length(ds); }
    
  ## set the number of induviduals
  if( method == "mean" ){ Rj = n; }
  else {                  Rj = colSums(Y)[inx]; }  
  

  KM   = cumsum(colSums(Esum)/colSums(Y))
  surv = exp(-cumsum(KM*ds))
 
  for( e in 1:ne ){
      num     = colSums(Ek[,,e])
      dem     = colSums(Y)
      AJk[,e] = cumsum(num/dem)
      CIF[,e] = cumsum(surv*AJk[,e]*ds)
      RMk[e]  = CIF[,e][inx]
  }

  Numerator  = matrix(colSums(Esum),n,ns,byrow=TRUE)-Esum
  Denominator= matrix(colSums(Y),n,ns,byrow=TRUE)-Y
  Denominator= apply(Denominator, 2, function(x) ifelse(x==0, SMALL, x))  
 
  ds.mat     = matrix(ds,nrow=n,ncol=length(ds),byrow=TRUE)  
  KMi        = t(apply(Numerator/Denominator,1,cumsum))
  survi      = exp(-t(apply(KMi*ds.mat,1,cumsum))) 
            
  RMik        = array(dim=c(n,ne))
  pseudo_mean = array(dim=c(n,ne))
  colnames(pseudo_mean) = n.events
    
  for( e in 1:ne ){
      Numk     = matrix(colSums(Ek[,,e]),n,ns,byrow=TRUE)-Ek[,,e]
      AJki     = t(apply( (Numk/Denominator),1,cumsum))
      tmp      = do.call(rbind, lapply(1:n,
                                    function(j)
                                     cumsum(survi[j,1:inx]*AJki[j,1:inx]*ds.mat[j,1:inx]) ))
      RMik[,e] = tmp[,inx]
      ##RMik[,e] = t(apply(survi*AJki*ds.mat,1,cumsum))[,inx]
            
      pseudo_mean[,e] = Rj*RMk[e]-(Rj-1)*RMik[,e]
      pseudo_mean[,e] = as.vector(pseudo_mean[order(pseudo$id),e])
  }

    return(list(pseudo_mean=pseudo_mean, surv=surv, AJk=AJk, CIF=CIF, RMk=RMk, s=s))  

}


##---------------------------------------
## Fine & Gray estimator
##---------------------------------------
pseudoFG_comp <- function(time,delta,events,cencode=0,tau,SMALL=1e-20,
                           method=c("mean","risk")){

  ## use risk set as defined in fine & gray method
  method = match.arg(method)
    
  ## number of competing events
  n.events = get.comp.events(events=events, cencode=cencode)
  ne       = length(n.events)  
    
  ## preparing the data
  n=length(time)
  pseudo <- data.frame(id=1:n,t=time,delta=delta,events=events)

  
  ## sort in time, if tied, put events before censoring
  ## order time lowest to highest, with events (i.e. delat==1) coming first if tied
  ## time, done by negating the values in vector delta  
  pseudo <- pseudo[order(pseudo$t,-pseudo$delta),]
  
  t=pseudo$t
  delta=pseudo$delta
  events=pseudo$events  
  
  s=sort(unique(t[delta==1]))

  tau = check.tmin(t.min=min(time),tau=tau);  
    
  ns=length(s)  # the number of intervals

    
  Y    <- array(0,dim=c(n,ns,ne))
  Ysum <- do.call(cbind, lapply(1:ns, function(j)  ifelse(s[j] <= t, 1, 0)))
  Ek   <- array(0,dim=c(n,ns,ne))
  Esum <- array(0,dim=c(n,ns))
  FGk  <- array(dim=c(ns,ne))
  CIF  <- array(dim=c(ns,ne)) 
  RMk  <- array(dim=c(ne)) 
    for( e in 1:ne ){
        ## Define the fine & gray risk set
        ## count the number of patients, not including the risk type under study,
        ## at risk at each unique time point.
        Y[,,e] = do.call(cbind, lapply(1:ns,
                                       function(j) (s[j] <= t)*(events != n.events[e])))
        
        Ek[,,e] = do.call(cbind, lapply(1:ns,
                                        function(j)  (s[j] == t)*(events == n.events[e])))
        Esum    = Esum + Ek[,,e]
  }

  inx  = max(which(ttmp<=tau))
  ttmp = c(0,s)
  ds   = diff(ttmp)

  if( inx > length(ds) ){ inx = length(ds); }
    
  ## set the number of induviduals
  if( method == "mean" ){ Rj = n; }
  else {                  Rj = colSums(Y)[inx]; }  

  ## Reference: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3135581/
  ## See eqn 2.
  ## BMC Med Res Methodol. 2011; 11: 86.
  ## Published online 2011 Jun 3. doi: 10.1186/1471-2288-11-86
  ## Understanding competing risks: a simulation point of view.
  ## Arthur Allignol
  KM   = cumsum(colSums(Esum)/colSums(Ysum))
  surv = exp(-cumsum(KM*ds))
 
  for( e in 1:ne ){
      num     = colSums(Ek[,,e])
      dem     = colSums(Y[,,e])
      FGk[,e] = cumsum(num/dem)
      CIF[,e] = cumsum(surv*FGk[,e]*ds)
      RMk[e]  = CIF[,e][inx]
  }

  Numerator  = matrix(colSums(Esum),n,ns,byrow=TRUE)-Esum
  Denominator= matrix(colSums(Ysum),n,ns,byrow=TRUE)-Ysum
  Denominator= apply(Denominator, 2, function(x) ifelse(x==0, SMALL, x))  
 
  ds.mat     = matrix(ds,nrow=n,ncol=length(ds),byrow=TRUE)  
  KMi        = t(apply(Numerator/Denominator,1,cumsum))
  survi      = exp(-t(apply(KMi*ds.mat,1,cumsum))) 
            
  RMik        = array(dim=c(n,ne))
  pseudo_mean = array(dim=c(n,ne))
  colnames(pseudo_mean) = n.events

  for( e in 1:ne ){
      Numk     = matrix(colSums(Ek[,,e]),n,ns,byrow=TRUE)-Ek[,,e]
      Denk     = matrix(colSums( Y[,,e]),n,ns,byrow=TRUE)- Y[,,e]
      AJki     = t(apply( (Numk/Denk),1,cumsum))
      tmp      = do.call(rbind, lapply(1:n,
                                    function(j)
                                     cumsum(survi[j,1:inx]*AJki[j,1:inx]*ds.mat[j,1:inx]) ))
      ##tmp      = t(apply(survi[,1:inx]*AJki[,1:inx]*ds.mat[,1:inx],1,cumsum))
      RMik[,e] = tmp[,inx]
            
      pseudo_mean[,e] = Rj*RMk[e]-(Rj-1)*RMik[,e]
      pseudo_mean[,e] = as.vector(pseudo_mean[order(pseudo$id),e])
  }
    return(list(pseudo_mean=pseudo_mean, surv=surv, FGk=FGk, CIF=CIF, RMk=RMk, s=s))  
    ##return(pseudo_mean) 

}


##---------------------------------------
## Weighted Fine & Gray estimator
##---------------------------------------
pseudoFG_comp_weight <- function(time,delta,events,cencode=0,tau, weights,
                                 SMALL=1e-20,
                                 method=c("mean","risk")){

  ## use risk set as defined in fine & gray method
  method = match.arg(method)
    
  ## number of competing events
  n.events = get.comp.events(events=events, cencode=cencode)
  ne       = length(n.events)  

  ## preparing the data
  n=length(time)
  pseudo <- data.frame(id=1:n,t=time,delta=delta,events=events)
    
  ## sort in time, if tied, put events before censoring
  ## order time lowest to highest, with events (i.e. delat==1) coming first if tied
  ## time, done by negating the values in vector delta  
  reorder <- order(pseudo$t,-pseudo$delta)
  pseudo  <- pseudo[reorder,]
  weights <- weights[reorder,]
  weights <- apply(weights, 2, function(x) ifelse(x==0, SMALL, x))  
    
  ## sort in time, if tied, put events before censoring
  ## order time lowest to highest, with events (i.e. delat==1) coming first if tied
  ## time, done by negating the values in vector delta  
  pseudo <- pseudo[order(pseudo$t,-pseudo$delta),]
  
  t=pseudo$t
  delta=pseudo$delta
  events=pseudo$events  
  
  s=sort(unique(t[delta==1]))

  tau = check.tmin(t.min=min(time),tau=tau); 
    
  ns=length(s)  # the number of intervals
    
  Y    <- array(0,dim=c(n,ns,ne))
  Ysum <- do.call(cbind, lapply(1:ns, function(j)  ifelse(s[j] <= t, weights[j], 0)))
  Ek   <- array(0,dim=c(n,ns,ne))
  Esum <- array(0,dim=c(n,ns))
  FGk  <- array(dim=c(ns,ne))
  CIF  <- array(dim=c(ns,ne)) 
  RMk  <- array(dim=c(ne)) 
    for( e in 1:ne ){
        ## Define the fine & gray risk set
        ## count the number of patients, not including the risk type under study,
        ## at risk at each unique time point.
        Y[,,e] = do.call(cbind, lapply(1:ns,
                                       function(j)
                                           ifelse((s[j] <= t)*(events != n.events[e]),
                                                  weights[,j], 0)
                                       ))
        
        Ek[,,e] = do.call(cbind, lapply(1:ns,
                                        function(j)
                                            ifelse((s[j] == t)*(events == n.events[e]),
                                                   weights[,j], 0)
                                        ))
                                        
        Esum    = Esum + Ek[,,e]
  }

  inx  = max(which(s<=tau))
  ttmp = c(0,s)
  ds   = diff(ttmp)

  if( inx > length(ds) ){ inx = length(ds); }
    
  ## set the number of induviduals
  if( method == "mean" ){ Rj = n; }
  else {                  Rj = colSums(Y)[inx]; }  

  ## Reference: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3135581/
  ## See eqn 2.
  ## BMC Med Res Methodol. 2011; 11: 86.
  ## Published online 2011 Jun 3. doi: 10.1186/1471-2288-11-86
  ## Understanding competing risks: a simulation point of view.
  ## Arthur Allignol
  KM   = cumsum(colSums(Esum)/colSums(Ysum))
  surv = exp(-cumsum(KM*ds))
 
  for( e in 1:ne ){
      num    = colSums(Ek[,,e])
      dem    = colSums(Y[,,e])
      FGk[,e] = cumsum(num/dem)
      CIF[,e] = cumsum(surv*FGk[,e]*ds)
      RMk[e]  = CIF[,e][inx]
      ##AJk    = cumsum(num/dem)
      ##RMk[e] = cumsum(surv*AJk*ds)[inx]
  }

  Numerator  = matrix(colSums(Esum),n,ns,byrow=TRUE)-Esum
  Denominator= matrix(colSums(Ysum),n,ns,byrow=TRUE)-Ysum
  Denominator= apply(Denominator, 2, function(x) ifelse(x==0, SMALL, x))  
 
  ds.mat     = matrix(ds,nrow=n,ncol=length(ds),byrow=TRUE)  
  KMi        = t(apply(Numerator/Denominator,1,cumsum))
  survi      = exp(-t(apply(KMi*ds.mat,1,cumsum))) 
            
  RMik        = array(dim=c(n,ne))
  pseudo_mean = array(dim=c(n,ne))
  colnames(pseudo_mean) = n.events

  for( e in 1:ne ){
      Numk     = matrix(colSums(Ek[,,e]),n,ns,byrow=TRUE)-Ek[,,e]
      Denk     = matrix(colSums( Y[,,e]),n,ns,byrow=TRUE)- Y[,,e]
      AJki     = t(apply( (Numk/Denk),1,cumsum))
      tmp      = do.call(rbind, lapply(1:n,
                                    function(j)
                                     cumsum(survi[j,1:inx]*AJki[j,1:inx]*ds.mat[j,1:inx]) ))
      ##tmp      = t(apply(survi[,1:inx]*AJki[,1:inx]*ds.mat[,1:inx],1,cumsum))
      RMik[,e] = tmp[,inx]
            
      pseudo_mean[,e] = Rj*RMk[e]-(Rj-1)*RMik[,e]
      pseudo_mean[,e] = as.vector(pseudo_mean[order(pseudo$id),e])
  }
    return(list(pseudo_mean=pseudo_mean, surv=surv, FGk=FGk, CIF=CIF, RMk=RMk, s=s))  
    ##return(pseudo_mean) 

}

##---------------------------------------
## Class weighted Fine & Gray estimator
##---------------------------------------
pseudoFG_comp_class_we <- function(time,delta,events,cencode=0,tau, 
                                   SMALL=1e-20,
                                   method=c("mean","risk")){

  ## use risk set as defined in fine & gray method
  method = match.arg(method)
    
  ## number of competing events
  n.events = get.comp.events(events=events, cencode=cencode)
  ne       = length(n.events)  

  ## number of classes
  class          = c(cencode,n.events)  
  nclass         = length(class)
  weights        = rep(1,nclass)
  names(weights) = class  

  ## do we have class names matching event type names 
  #if( all(names(alpha) %in% c(cencode,n.events)) ){
  #    cat("> class names match event type names!\n")
  #} else {
  #    cat("> class names don't match event type names!\n")
  #} 

  ## patient weights
  alpha = 1-prop.table(table(events[time<=tau]))    
  weights[match(names(alpha),names(weights))] = alpha
  weights = weights[events]
    
  ## preparing the data
  n=length(time)
  pseudo <- data.frame(id=1:n,t=time,delta=delta,events=events)
    
  ## sort in time, if tied, put events before censoring
  ## order time lowest to highest, with events (i.e. delat==1) coming first if tied
  ## time, done by negating the values in vector delta  
  reorder <- order(pseudo$t,-pseudo$delta)
  pseudo  <- pseudo[reorder,]
    
  ## sort in time, if tied, put events before censoring
  ## order time lowest to highest, with events (i.e. delat==1) coming first if tied
  ## time, done by negating the values in vector delta  
  pseudo <- pseudo[order(pseudo$t,-pseudo$delta),]
  
  t=pseudo$t
  delta=pseudo$delta
  events=pseudo$events  
  
  s=sort(unique(t[delta==1]))

  tau = check.tmin(t.min=min(time),tau=tau); 
    
  ns=length(s)  # the number of intervals
    
  Y    <- array(0,dim=c(n,ns,ne))
  Ysum <- do.call(cbind, lapply(1:ns, function(j)  ifelse(s[j] <= t, weights, 0)))
  Ek   <- array(0,dim=c(n,ns,ne))
  Esum <- array(0,dim=c(n,ns))
  FGk  <- array(dim=c(ns,ne))
  CIF  <- array(dim=c(ns,ne))   
  RMk  <- array(dim=c(ne)) 
    for( e in 1:ne ){
        ## Define the fine & gray risk set
        ## count the number of patients, not including the risk type under study,
        ## at risk at each unique time point.
        Y[,,e] = do.call(cbind, lapply(1:ns,
                                       function(j)
                                           ifelse((s[j] <= t)*(events != n.events[e]),
                                                  weights, 0)
                                       ))
        
        Ek[,,e] = do.call(cbind, lapply(1:ns,
                                        function(j)
                                            ifelse((s[j] == t)*(events == n.events[e]),
                                                   weights, 0)
                                        ))
                                        
        Esum    = Esum + Ek[,,e]
  }

  inx  = max(which(s<=tau))
  ttmp = c(0,s)
  ds   = diff(ttmp)

  if( inx > length(ds) ){ inx = length(ds); }
    
  ## set the number of induviduals
  if( method == "mean" ){ Rj = n; }
  else {                  Rj = colSums(Y)[inx]; }  

  ## Reference: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3135581/
  ## See eqn 2.
  ## BMC Med Res Methodol. 2011; 11: 86.
  ## Published online 2011 Jun 3. doi: 10.1186/1471-2288-11-86
  ## Understanding competing risks: a simulation point of view.
  ## Arthur Allignol
  KM   = cumsum(colSums(Esum)/colSums(Ysum))
  surv = exp(-cumsum(KM*ds))
 
  for( e in 1:ne ){
      num     = colSums(Ek[,,e])
      dem     = colSums(Y[,,e])
      FGk[,e] = cumsum(num/dem)
      CIF[,e] = cumsum(surv*FGk[,e]*ds)
      RMk[e]  = CIF[,e][inx]
      #AJk    = cumsum(num/dem)
      #RMk[e] = cumsum(surv*AJk*ds)[inx]
  }

  Numerator  = matrix(colSums(Esum),n,ns,byrow=TRUE)-Esum
  Denominator= matrix(colSums(Ysum),n,ns,byrow=TRUE)-Ysum
  Denominator= apply(Denominator, 2, function(x) ifelse(x==0, SMALL, x))  
 
  ds.mat     = matrix(ds,nrow=n,ncol=length(ds),byrow=TRUE)  
  KMi        = t(apply(Numerator/Denominator,1,cumsum))
  survi      = exp(-t(apply(KMi*ds.mat,1,cumsum))) 
            
  RMik        = array(dim=c(n,ne))
  pseudo_mean = array(dim=c(n,ne))
  colnames(pseudo_mean) = n.events

  for( e in 1:ne ){
      Numk     = matrix(colSums(Ek[,,e]),n,ns,byrow=TRUE)-Ek[,,e]
      Denk     = matrix(colSums( Y[,,e]),n,ns,byrow=TRUE)- Y[,,e]
      AJki     = t(apply( (Numk/Denk),1,cumsum))
      tmp      = do.call(rbind, lapply(1:n,
                                    function(j)
                                     cumsum(survi[j,1:inx]*AJki[j,1:inx]*ds.mat[j,1:inx]) ))
      ##tmp      = t(apply(survi[,1:inx]*AJki[,1:inx]*ds.mat[,1:inx],1,cumsum))
      RMik[,e] = tmp[,inx]
            
      pseudo_mean[,e] = Rj*RMk[e]-(Rj-1)*RMik[,e]
      pseudo_mean[,e] = as.vector(pseudo_mean[order(pseudo$id),e])
  }

    return(list(pseudo_mean=pseudo_mean, surv=surv, FGk=FGk, CIF=CIF, RMk=RMk, s=s)) 
    ##return(pseudo_mean) 

}

##---------------------------------------
## use R package pseudo's pseudo value calculation
##---------------------------------------
Rpackage_pseudo <- function(time,delta,events,cencode=0,tau){

    ## number of competing events
    n.events = get.comp.events(events=events, cencode=cencode)
    ne       = length(n.events)  
    
    zero.indx  = which(tau<=0)
    nt.indx    = which(tau>0) 
    if( length(zero.indx) > 0 ){
        ## evalutation times
        eval.t = tau[-zero.indx]
    } else {
        eval.t = tau
    }

    res = pseudoci(event=events, time=time, tmax=eval.t)

    nr  = dim(res$pseudo[[1]])[1]
    nt  = length(tau)
    
    pseudo = array(0,dim=c(nr,nt,ne))

    for( e in 1:ne ){
        pseudo[,nt.indx,e] = res$pseudo[[e]]
    }

    return(pseudo)
}
    
eval.stats <- function( train=NULL, test=NULL, surv.prob.tr=NULL,
                       surv.prob.ts=NULL, qt=qt, time.points=NULL ){

    ## 
    
    if( !is.null(train) && !is.null(test) &&
        !is.null(surv.prob.tr) && !is.null(surv.prob.ts) &&
        !is.null(qt) && !is.null(time.points) ){

        N = length(qt)

        train   = as.data.frame(train)
        test    = as.data.frame(test)

        ## V1 ==> observed time
        ## V2 ==> censor indicator
        Surv.tr = Surv(train$V1, train$V2)
        Surv.ts = Surv(test$V1, test$V2)

        ## harrel cindex
        hc = matrix(NA, ncol=2, nrow=N)

        ## begg cindex
        bc = vector(length=N)

        ## uno cindex
        unoc = vector(length=N)

        ## auc
        auc     = matrix(NA, ncol=3, nrow=length(time.points))
        auc.tmp = matrix(NA, ncol=N, nrow=length(time.points))
        
        for( i in 1:N ){
            tmp.hc  = rcorr.cens(surv.prob.ts[,i], Surv.ts)
            hc[i,1] = tmp.hc["C Index"]
            hc[i,2] = tmp.hc["S.D."]
            bc[i]   = BeggC(Surv.tr, Surv.ts, -surv.prob.tr[,i],
                          -surv.prob.ts[,i])
            unoc[i] = UnoC(Surv.tr, Surv.ts, -surv.prob.ts[,i])

            AUC_CD  = AUC.uno(Surv.tr, Surv.ts, -surv.prob.ts[,i], time.points)

            auc.tmp[,i] = AUC_CD$auc

        }

        auc[,1] = time.points
        auc[,2] = apply(auc.tmp, 1, mean)
        auc[,3] = apply(auc.tmp, 1, sd)
        
    }

    return(list(hc=hc, bc=bc, unoc=unoc,auc=auc))
}
    
    #harrel cindex
    #harrelC1 <- rcorr.cens(pred_te[,1],with(test,Surv(time,status)))
    #hc<-harrelC1["C Index"]
    ##begg cindex
    #lp<- -pred_tr[,1]
    #lpnew <- -pred_te[,1]
    #Surv.rsp <- Surv(train$time, train$status)
    #Surv.rsp.new <- Surv(test$time, test$status)              
    #bc <- BeggC(Surv.rsp, Surv.rsp.new,lp, lpnew)
    #uno cindex
    #unoc<-UnoC(Surv.rsp, Surv.rsp.new, lpnew)
    
#}

## obtain the marginal survial prob. by multiplying series of conditional
## survival prob. (x) at each time.point (qt)
get.surv.prob <- function(x=NULL, qt=NULL){

    surv.prob = NULL

    if( !is.null(x) && !is.null(qt) ){
        ## obtain the marginal survial prob. by multiplying series of conditional survival prob.
        ## ypred.marg is a set of length(qt) vectors such that for patient 1:
        ## ypred.marg[[1]][1] = ypred.cond[1,1]
        ## ypred.marg[[2]][1] = ypred.cond[1,1] * ypred.cond[1,2]
        ## ypred.marg[[3]][1] = ypred.cond[1,1] * ypred.cond[1,2] * ypred.cond[1,3]
        ## ypred.marg[[n]][1] = ypred.cond[1,1] * ypred.cond[1,2] *...* ypred.cond[1,n] 
        ypred.marg = lapply(1:length(qt), function(i) apply(x[,1:i,
                                                            drop=FALSE], 1,
                                                            prod))
        surv.prob  = Reduce(cbind, ypred.marg)
    }
    
    return(surv.prob)

}


get.pred.multiState <- function(x, test.df=NULL, qt, nr){

    ## format predict multiState probability matrix 'x'.
    ##-------------------------------------------------

    if( is.null(test.df) ){ nr = nr; }
    else {
        nr = dim(test.df)[1] ## number of patients
    }
    nc = length(qt)      ## number of time points
    ns = dim(x)[2]       ## number of multiStates
    df = array(0,dim=c(nr,ns,nc))

    for( i in 1:nc ){
        if( i == 1 ){
            start=1
            end=nr
        } else {
            start=end+1
            end=i*nr            
        }        
        df[,,i] = x[start:end,]        
    }

    return(df)
    
}



get.comp.survival.time <- function(t,events,qt,cencode=0,
                                   method=c("aj","fg")){

    method=match.arg(method)
        
    ## get the nnumber of competing events
    n.events = get.comp.events(events=events, cencode=cencode)
    ne       = length(n.events)

    psedoComp =NULL
    pseudoNorm=NULL
    switch(method,
           "aj"={
               ### pseudoAJ_comp_weight calculates the expected probability of
               ### each patient experiencing each competing risk event at
               ### each landmark time.
                pseudoComp = do.call(cbind, lapply(qt, function(s)
                    pseudoAJ_comp(time=t$obs.time,
                                  delta=t$deltaC,
                                  events=events,
                                  method="mean",
                                  tau=s)))
           },
           "fg"={
               pseudoComp = do.call(cbind, lapply(qt, function(s)
                    pseudoFG_comp(time=t$obs.time,
                                  delta=t$deltaC,
                                  events=events,
                                  method="mean",
                                  tau=s)))
           },
           {
               print("> Error: set method to \"aj\" or \"fg\".\n")
           }
           )    

    if( !is.null(pseudoComp) ){
    
        ## format pseudo values per risk sets
        nr=dim(pseudoComp)[1]
        nc=length(qt)
        pseudoNorm = array(0,dim=c(nr,nc,ne))
        for( e in 1:ne ){
            indx = which(colnames(pseudoComp)==n.events[e])
            pseudoNorm[,,e] = pseudoComp[,indx]
        }    
    
        ## normalise pseudo values
        for( e in 1:ne ){
            tmp = exp(-pseudoNorm[,,e])
            tmp = t(sapply(1:nr, function(j) ifelse(tmp[j,]>1,1,tmp[j,])))
            tmp = 1-tmp
            pseudoNorm[,,e] = tmp
        }

    }        

    return(list(pseudoNorm=pseudoNorm, weights=we))

}
   

get.we.comp.survival.time <- function(t,events,qt,cencode=0,
                                      covar,x.indx,x.strat,x.we, alpha,
                                      type=c("age.strat","rf","unwe","class"),
                                      method=c("aj","fg","pseudo"),
                                      runNorm=TRUE){

    method = match.arg(method)
    type   = match.arg(type)
    
    ## get the nnumber of competing events
    n.events = get.comp.events(events=events, cencode=cencode)
    ne       = length(n.events)

    alpha.map   = alpha[match(events,names(alpha))]
    
    we=NULL
    switch(type,
           "age.strat"={
               we = compute_weight_stratified(x=covar[,x.indx],
                                              time=t$obs.time,
                                              delta=t$deltaC,
                                              x.strat=x.strat,
                                              x.we=x.we,
                                              method="custom")
           },
           "rf"={
               we = compute_weight_rf_trunc(time=t$obs.time,
                                            delta=t$deltaC,
                                            cov=covar)
           },
           "unwe"={
               we = compute_unweight(time=t$obs.time,
                                     delta=t$deltaC)
           },
           "class"={
               we = compute_unweight(time=t$obs.time,
                                     delta=t$deltaC)
               we = apply(we,2, function(x) x*alpha.map)
           },
           {
               print("> Error: no weight type given.\n") 
           }
           )
    

    pseudoComp = NULL
    pseudoNorm = NULL
   
    switch(method,
           "aj"={
               ### pseudoAJ_comp_weight calculates the expected probability of
               ### each patient experiencing each competing risk event at
               ### each landmark time.
                pseudoComp = do.call(cbind, lapply(qt, function(s)
                    pseudoAJ_comp_weight(time=t$obs.time,
                                         delta=t$deltaC,
                                         events=events,
                                         weights=we,
                                         method="mean",
                                         tau=s)[[1]]))
           },
           "fg"={
               pseudoComp = do.call(cbind, lapply(qt, function(s)
                    pseudoFG_comp_weight(time=t$obs.time,
                                         delta=t$deltaC,
                                         events=events,
                                         weights=we,
                                         method="mean",
                                         tau=s)[[1]]))
           },
           "pseudo"={
                pseudoNorm = Rpackage_pseudo(time=t$obs.time,
                                             delta=t$deltaC,
                                             events=events,
                                             tau=qt)
                runNorm=FALSE
           },
           
           {
               print("> Error: set method to \"aj\" or \"fg\".\n")
           }
           )    

    if( !is.null(pseudoComp) && runNorm==TRUE ){
    
        ## format pseudo values per risk sets
        nr=dim(pseudoComp)[1]
        nc=length(qt)
        pseudoNorm = array(0,dim=c(nr,nc,ne))
        for( e in 1:ne ){
            indx = which(colnames(pseudoComp)==n.events[e])
            pseudoNorm[,,e] = pseudoComp[,indx]
        }    
    
        ## normalise pseudo values
        for( e in 1:ne ){
            tmp = exp(-pseudoNorm[,,e])
            tmp = t(sapply(1:nr, function(j) ifelse(tmp[j,]>1,1,tmp[j,])))
            tmp = 1-tmp
            pseudoNorm[,,e] = tmp
        }

    }        

    return(list(pseudoNorm=pseudoNorm, weights=we))

}



##-------------------------------------------------------------------------------------
## compute each patients cause specific (i.e. competing event) probability of
## survival using survival boosted decision tree.
## Reference: https://luminwin.github.io/randomForestSRC/articles/competing.html
##-------------------------------------------------------------------------------------
compute_weight_rf.comp <- function(covs=NULL, time=NULL, event=NULL, cencode=0,
                                   method=c("prob","prob.wide",
                                            "weight","weight.wide",
                                            "ccp","ccp.wide") ){

    output=NULL

    method = match.arg(method)
    
    if( !is.null(covs) && !is.null(time) && !is.null(event) ){

        n=dim(covs)[1]
        data=data.frame(t=time,e=event,covs)
        res=rfsrc(Surv(t,e) ~., data)

        t.unique=res$time.interest
        s=sort(unique(time[event!=cencode]))
        ns=length(s)

        n.events = res$event.info$event.type
        ne       = length(n.events)
        cif      = array(,dim=c(n,ns,ne))
        sum.cif  = array(0,dim=c(n,ns))
        output   = array(,dim=c(n,ns,ne))
        
        for(e in 1:ne){
            indx = n.events[e]
            cif[,,indx]=t(sapply(1:n,
                                 function(j) approx(t.unique, res$cif[j,,indx],
                                                    xout=s, method="constant",
                                                    rule=2:2)$y) ) 
            sum.cif = sum.cif + cif[,,indx]
        }

        if( method=="weight" || method=="weight.wide" ){
            for( e in 1:ne ){
                indx = n.events[e]
                output[,,indx] = 1/(1-cif[,,indx])
            }
        }

        if( method=="ccp" || method=="ccp.wide" ){
            ## cumalative conditional probability (CCP):
            ## CCP the probability of
            ## experiencing an event c by time t, given that an
            ## individual has not experienced any of the other
            ## competing risks by time t.
            for( e in 1:ne ){
                indx = n.events[e]
                output[,,indx] = cif[,,indx]/(1-(sum.cif-cif[,,indx]))
            }
        }                   

        #if( grepl(".wide",method) ){
        #    
        #}
        
    }
    
    return(output)

}


norm.pseudo.values <- function( x, min, max ){    
    if( (min == 0   && max == 0) ||
        (is.na(min) && is.na(max)) ){ x; }
    else { (x-min)/(max-min); }
}

pseudo.values.long <- function(data, qt){
    colnames(data) = 1:dim(data)[2]
    long <- as.data.table(data) %>%
        rename(setNames(colnames(.),as.character(qt))) %>%
        rowid_to_column(var='id') %>%
        pivot_longer(cols=as.character(qt), names_to="s", values_to="pseudost") %>%
        mutate(id = as.numeric(id)) %>%
        mutate(s = as.numeric(s)) %>%
        mutate( pseudost = as.numeric(pseudost) )
    return(long);
}


 ## calculate the cumalative conditional probability (CCP):
 ## CCP is the probability of experiencing an event c by time t,
 ## given that an individual has not experienced any of the other
 ## competing risks by time t.
ccp.values <- function(x,n.events,norm.values=FALSE){

    nr = dim(x)[1]
    nc = length(which(colnames(x)==n.events[1]))
    ne = length(n.events)    
    
    if( norm.values ){ ## normalise pseudo values        
        for( e in 1:ne ){
            indx = which(colnames(x)==n.events[e])
            tmp  = x[,indx]
            min  = apply(tmp,2,min)
            max  = apply(tmp,2,max)
            tmp  = sapply(1:dim(tmp)[2], function(s) norm.pseudo.values(x=tmp[,s],
                                                                        min=min[s],
                                                                      max=max[s]))      
            x[,indx] = tmp
        }
    }

    tmp = array(0,dim=c(nr,nc,ne))
    ccp = array(0,dim=c(nr,nc,ne))
    tot = array(0,dim=c(nr,nc))
    for( e in 1:ne ){
        tmp[,,e] = x[,which(colnames(x)==n.events[e])]
        tot      = tot + tmp[,,e]
    }
    
    min  = apply(tot,2,min)
    max  = apply(tot,2,max)
    tot  = sapply(1:dim(tot)[2], function(s) (tot[,s]-min[s])/(max[s]-min[s]) )
    
    for( e in 1:ne ){
        ccp[,,e] = tmp[,,e]/(1-(tot-tmp[,,e]))
    }
    
    return(ccp)

}



risk.at.time.points <- function(Time_survival, Cause, Cause_int ){
    ## Ref 1: See eqn (3) in:
    ## Lee, C. et, al. DeepHit: A Deep Learning Approach to Survival Analysis with
    ##                 Competing Risks, AAAI Conference, 2018. 
    ##
    ## Reurn Akij object nedded for l2.norm function.

    n = length(Cause)
    ##A = matrix(0, nrow=n, ncol=n)
    
    ##for (i in  1:n){
    ##    A[i, intersect(which(Time_survival[i] < Time_survival),
    ##                   which(Cause==Cause_int) )]=1 ## eqn (3)
    ##}

    ##return(A)

    ## using rbind will give use exactly 'A' above,
    ## but if R is column-wise, use cbind and then t(A) at the end.
    A <- do.call(cbind, lapply(1:n, function(i)
                              (Time_survival[i] < Time_survival)*(Cause==Cause_int)))

    return(t(A))
}
