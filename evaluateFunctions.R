##---------------------------------------
## Aalen-Johannsen competing estimator
##---------------------------------------
## Inputs:
## time   = xtest$data.t$obs.time
## delta  = xtest$data.t$deltaC
## events = xtest$data.events
## tau    = qt[2]
AJ_comp <- function(time, delta, events, cencode=0, tau, SMALL=1e-20 ){

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

    names(RMk) = n.events
    
    return(RMk)
    
}

## Inputs:
## time   = xtrain$data.t$obs.time
## events = xtrain$data.t$deltaC
censor.probability <- function(time, events, interval=0.6){
    time.max = ceiling(max(time))
    km       = Surv(time, events)
    km_fit   = survfit(km ~ 1)
    km.prob  = summary(km_fit, times=seq(0,(time.max+1),interval))
    return(km.prob)
}
    
##### Cindex estimation for competeing events 
## Paper: KARTIK AHUJA, MIHAELA VAN DER SCHAAR:
##        Generalized Concordance for Competing Risks. (2021), 0, 0, pp. 1–30
##        doi:10.1093/joint˙concordance˙biostatistics˙oct˙23˙arxiv

## Ref:
## [1] http://medianetlab.ee.ucla.edu/papers/Kartikgenindex/New_Cindex_paper_upload.Rmd 
## [2] https://deepai.org/publication/generalized-concordance-for-competing-risks
## Input: Prediction vector for the event type of interest (risk predicted by the model ## for the event type of interest),
## Time_survival: time-to-event data for each subject,
## Censoring: censoring information for each subject (1==any event, 0==censored),
## Cause: vector of event types,
## Cause_int: event type of interest,
## Time: time horizon at which c-index is computed
##
## Output: Concordance index
## Example Inputs:
## Prediction=y_pred.dh[,1,2]
## Time=qt[2]
## Cause=xtest$data.events
## Cause_int=1 chose one of [1,2,3] for the competing risks
## Time_survival=xtest$data.t$obs.time
## Censoring=xtrain$data.t$deltaC
## method="hazard"
Cindex_estimator_efficient <- function(Prediction, Time_survival,
                                       Censoring, Cause, Cause_int, Time,
                                       method=c("survival","hazard")){

    ## A and B account for the number of risk order pairs.
    
    ## A == risk ordering of patients, small time means patient 'i' at higher risk
    ##      than patient 'j' experiencing event of interest 

    ## B == risk ordering of patients, large time for patient 'i' means lower risk 
    ##       than patient 'j' if not experienced the event of interest.
    
    ## Q == the risk ordering of the subjects, i.e., is subject i assigned a higher
    ##      risk by the model than the subject j, for event Ek until time t.

    ## N_t == number of subjects with survival time < time point and experience event
    ##        of interest

    ##          sum^N_i=1 sum^N_j=1 (A_ij + B_ij) * Q_ij * N^k_i(t)
    ## C_k(t) = -------------------------------------------------
    ##          sum^N_i=1 sum^N_j=1 (A_ij + B_ij) * N^k_i(t)

    method = match.arg(method)

    if( method=="survival" ){ Prediction=1-Prediction; }
    
    n = length(Prediction)
    A = matrix(0, nrow=n, ncol=n)
    B = matrix(0, nrow=n, ncol=n)
    Q = matrix(0, nrow=n, ncol=n) 
    N_t = matrix(0, nrow=n, ncol=n)
    Num_mat = matrix(0, nrow=n, ncol=n)
    Den_mat = matrix(0, nrow=n, ncol=n)
    
    Num=0
    Den=0
    for (i in  1:n){
        A[i,which(Time_survival[i] < Time_survival)] = 1 
        B[i, intersect(intersect(which((Time_survival[i] >= Time_survival)),
                                 which(Cause!=Cause_int) ), which(Censoring==1))] = 1
        Q[i,which(Prediction[i]>Prediction)]=1
    }

    for (i in 1:n){
        if(Time_survival[i]<=Time && Cause[i]==Cause_int && Censoring[i]==1){
            N_t[i,] = 1
        }
    }

    Num_mat = (A+B)*Q*N_t
    Den_mat = (A+B)*N_t

    Num = sum(Num_mat)
    Den = sum(Den_mat)

    return(Num/Den)

}

## Example Inputs:
## Prediction=y_pred.dh[,1,2]
## Time=qt[2]
## Cause=xtest$data.events
## Cause_int=1 chose one of [1,2,3] for the competing risks
## Time_survival=xtest$data.t$obs.time
Cindex_deephit <- function(Prediction, Time_survival,
                           Cause, Cause_int, Time){

   
    ## A == risk ordering of patients, small time means patient 'i' at higher risk
    ##      than patient 'j' experiencing event of interest 
    
    ## Q == the risk ordering of the subjects, i.e., is subject i assigned a higher
    ##      risk by the model than the subject j, for event Ek until time t.

    ## N_t == number of subjects with survival time < time point and experience event
    ##        of interest
     
    n = length(Prediction)
    A = matrix(0, nrow=n, ncol=n)
    Q = matrix(0, nrow=n, ncol=n) 
    N_t = matrix(0, nrow=n, ncol=n)
    Num_mat = matrix(0, nrow=n, ncol=n)
    Den_mat = matrix(0, nrow=n, ncol=n)
    
    Num=0
    Den=0
    for (i in  1:n){
        A[i,which(Time_survival[i] < Time_survival)] = 1         
        Q[i,which(Prediction[i]>Prediction)]=1
    }

    for (i in 1:n){
        if(Time_survival[i]<=Time && Cause[i]==Cause_int){
            N_t[i,] = 1
        }
    }

    Num_mat = A*Q*N_t
    Den_mat = A*N_t

    Num = sum(Num_mat)
    Den = sum(Den_mat)

    if( Num==0 && Den==0 ){
        result=-1;
    } else {
        result=Num/Den;
    }
    
    return(result)

}


####-------------------------------------------------
## Ref: https://onlinelibrary.wiley.com/doi/epdf/10.1002/sim.1802?saml_referrer
## Conservative estimate for ci Confindence Interval
## See eqn (14)
####-------------------------------------------------
## Example Inputs:
## Prediction=y_pred.dh[,1,2]
## Time_survival=xtest$data.t$obs.time
## cindex=Cindex_estimator_efficient(...)
## method="hazard"
Cindex_estimator_efficient_CI <- function(Prediction, Time_survival,
                                          cindex=1,alpha=0.05,
                                          method=c("survival","hazard")){

    upper=NA
    lower=NA

    if( !is.na(cindex) ){
    
    method = match.arg(method)

    if( method=="survival" ){ Prediction=1-Prediction; }    
    
    n = length(Prediction)
    A = matrix(0, nrow=n, ncol=n)
    B = matrix(0, nrow=n, ncol=n)

    for( i in 1:n ){
        A[i,intersect(which(Time_survival[i]<Time_survival),
                      which(Prediction[i]<Prediction))]=1
        B[i,intersect(which(Time_survival[i]>Time_survival),
                      which(Prediction[i]>Prediction))]=1
    }

    C=A+B
    Ch=rowSums(C)-diag(C)
    Dh=n-Ch
    P=1/(n*(n-1))*sum(Ch)
    D=1/(n*(n-1))*sum(Dh)

    Zalpha=qnorm(1-alpha/2)
    w=(2*Zalpha*Zalpha)/(n*(P+D))
    
    a=(w+2*cindex)/(2*(1+w))
    b=sqrt((w*w + 4*w*cindex*(1-cindex)))/(2*(1+w))

        upper=a+b
        lower=a-b
        
    }
        
    return(list(cindex=cindex,upper=upper,lower=lower))
    
}


## BRIER-SCORE
## Inputs:
## Prediction=y_pred.dh[,1,2]
## Time=qt[2]
## Cause=xtest$data.t$deltaC, or xtest$data.events
## Cause_int=1,               or [1,2,3] for the competing risks
## Time_survival=xtest$data.t$obs.time
brier_score <- function(Prediction, Time_survival, Cause, Cause_int, Time){

    ##n = length(Prediction)
    y_true = ((Time_survival <= Time) * (Cause==Cause_int))

    bscore = mean((Prediction-y_true)^2)
    
    return(bscore)
    
}

## Weighted BRIER-SCORE
## References:
## 1) https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-022-01679-6#Sec11
## 2): https://square.github.io/pysurvival/metrics/brier_score.html
## 
## Inputs:
## Prediction=y_pred.dh[,1,2]
## Time=qt[2]
## Cause=xtest$data.t$deltaC, or xtest$data.events
## Cause_int=1,               or [1,2,3] for the competing risks
## Time_survival=xtest$data.t$obs.time
## Train_DeltaC=xtrain$data.t$deltaC
## Train_Time=xtrain$data.t$obs.time
weighted_brier_score <- function(Train_Time, Train_DeltaC,
                                 Time, Time_survival,
                                 Cause, Cause_int,
                                 Prediction,
                                 calBS=TRUE){

    Gtmp = censor.probability(time=Train_Time, events=Train_DeltaC)
    G    = cbind(Gtmp$time, Gtmp$surv)

    BS   = NULL
    if( calBS ){
        BS   = brier_score(Prediction=Prediction, Time_survival=Time_survival, Cause=Cause, Cause_int=Cause_int, Time=Time)
    }
        
    n     = length(Prediction)
    W     = rep(0,n)
   
    Y_tilde = (Time_survival > Time)
    Y_true  = (Time_survival <= Time) * (Cause==Cause_int)    
    
    for( i in 1:n ){

        indx1 = which(G[,1] >= Time_survival[i])
        indx2 = which(G[,1] >= Time)

        if( length(indx1) > 0 ){
            G1 = G[indx1[1],2]
        } else {
            G1 = 1
        }

        if( length(indx2) > 0 ){
            G2 = G[indx2[1],2]
        } else {
            G2 = 1
        }

        W[i] = Y_tilde[i]/G2 + Y_true[i]/G1        
           
    }    

    S    = (1/n)*sum(Y_tilde)
    ##if( S>0 ){ 
    ##    norm = 1/(n*S)
    ##} else {
    ##norm = 1/n
    ##}

    norm        = 1/n
    BS_we_unorm = sum(W * (Y_true - Prediction)^2)
    BS_we       = norm * BS_we_unorm
    
    return(list(bs=BS, bs.we=BS_we, bs.we.unrom=BS_we_unorm, n=n, St=S, norm=norm))
    
}


## Scaled BRIER-SCORE
## Reference: https://arxiv.org/pdf/2212.05157.pdf
## K, .Monterrubio-Gómez et al. A REVIEW ON COMPETING RISKS METHODS FOR SURVIVAL ANALYSIS, 2022.
## doi: https://doi.org/10.48550/arXiv.2212.05157
## Inputs:
## Prediction=y_pred.dh[,1,2]
## qt=qt
## qt_indx=2
## Cause=xtest$data.t$deltaC, or xtest$data.events
## Cause_int=1,               or [1,2,3] for the competing risks
## Time_survival=xtest$data.t$obs.time
## Train_DeltaC=xtrain$data.t$deltaC
## Train_Time=xtrain$data.t$obs.time
scaled_brier_score <- function(Train_Time, Train_DeltaC,
                               qt, qt_indx, Time_survival,
                               Test_DeltaC,
                               Cause, #Cause_int,
                               Prediction,
                               SMALL=1e-5){


    ## set the evaluation time.point
    Tau = qt[qt_indx]

    ## Calculate the Brier score under the null model (no covariates) using the Aalen-Johansen estimator
    BS_null = AJ_comp(time=Time_survival, delta=Test_DeltaC, events=Cause, tau=Tau)


    BSwe_k =  vector("list", length(BS_null))
    names(BSwe_k) = names(BS_null)
    
    BSscaled_k        = rep(0,length(BS_null))
    names(BSscaled_k) = names(BS_null)
    
    for( i in 1:length(BSwe_k) ){
        e   = as.numeric(names(BSwe_k)[i])
        BSwe_k[[i]] = weighted_brier_score(Train_Time=Train_Time,
                                           Train_DeltaC=Train_DeltaC,
                                           Time=Tau,
                                           Time_survival=Time_survival,
                                           Cause=Cause,
                                           Cause_int=e,
                                           Prediction=Prediction[,e],##,qt_indx],
                                           calBS=TRUE)


        if( BS_null[i] < SMALL ){
            BSscaled_k[i] = 0
        } else {
            BSscaled_k[i] = 1-(BSwe_k[[i]][[2]]/BS_null[i])
        }
    }

    return(list(bs.null=BS_null, bs.we.k=BSwe_k, bs.scaled.k=BSscaled_k, tau=Tau))
}
