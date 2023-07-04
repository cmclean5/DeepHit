rm.zero.qt <- function(qt){
    if( qt[1] == 0 ){ return(qt[-1]) }

    return(qt)
}

time.points.by.subject <- function(n,qt){
    ## input: n == number of subjects,
    ##        qt == vector of time points
    ## return set of labels of time.point per subject
    ## For example:
    ## subject (n), time.point (qt)
    ## 1,t.1
    ## 2,t.1
    ## 3,t.1
    ## 1,t.2
    ## 2,t.2
    ## 3,t.2
    ## 1,t.3
    ## 2,t.3
    ## 3,t.3

    label.1 = rep(qt,each=n)
    label.2 = rep(seq(1,length(qt),1), each=n)

    return(list(label1=label.1,label2=label.2))
    
}

subject.by.time.points <- function(n,qt){
    ## input: n == number of subjects,
    ##        qt == vector of time points
    ## return set of labels of subjects per time.point
    ## For example:
    ## subject (n), time.point (qt)
    ## 1,t.1
    ## 1,t.2
    ## 1,t.3
    ## 2,t.1
    ## 2,t.2
    ## 2,t.3
    ## 3,t.1
    ## 3,t.2
    ## 3,t.3

    label.1 = rep(qt,times=n)
    label.2 = rep(seq(1,length(qt),1), times=n)

    return(list(label1=label.1,label2=label.2))
    
}


covariates.by.subject <- function(x,qt,rnd.col.indx=NULL){

    ## format subject covariate data 'x' by subject.
    ## x == [N subject, M covariates], qt = T time.points
    ## df == [N subject x T time.points, M covariates]
    ## For example, if we had 3 subjects, 3 time.points and 2 covariates,
    ## the output format by subject 'df' would look like this: 
    ##       c.1  c.2
    ## 1,t.1 .    .
    ## 2,t.1 .    .
    ## 3,t.1 .    .
    ## 1,t.2 .    .
    ## 2,t.2 .    .
    ## 3,t.2 .    .
    ## 1,t.3 .    .
    ## 2,t.3 .    .
    ## 3,t.3 .    .
    ##-------------------------------------------------
    ##if( !zero.qt ){ qt = rm.zero.qt(qt) }

    df=do.call(rbind,replicate(length(qt), x, simplify=FALSE))
    if( !is.null(rnd.col.indx) ){ df = randomiseCovariate(df,rnd.col.indx)}
    return(df)
}

covariates.by.time.point <- function(x,qt){##, zero.qt=FALSE){

    ## format subject covariate data 'x' by time.points.
    ## x == [N subject, M covariates], qt = T time.points
    ## df == [N subject x T time.points, M covariates]
    ## For example, if we had 3 subjects, 3 time.points and 2 covariates,
    ## the output format by time.points 'df' would look like this: 
    ##       c.1  c.2
    ## 1,t.1 .    .
    ## 1,t.2 .    .
    ## 1,t.3 .    .
    ## 2,t.1 .    .
    ## 2,t.2 .    .
    ## 2,t.3 .    .
    ## 3,t.1 .    .
    ## 3,t.2 .    .
    ## 3,t.3 .    .
    ##-------------------------------------------------
    ##if( !zero.qt ){ qt = rm.zero.qt(qt) }
    
    n  = dim(x)[1]  ## number of patients
    nc = dim(x)[2]  ## number of covariates
    ns = length(qt) ## number of time points
   
    ## Example, for first subject in 'x' 
    ## y    = replicate(ns,x[1,],simplify=FALSE)
    ## df.y = matrix(unlist(x),ncol=nc,byrow=TRUE)
    df = do.call(rbind, lapply(1:n, function(j) matrix(unlist(
                                                    replicate(ns,x[j,],simplify=FALSE)),
                                                    ncol=nc,
                                                    byrow=TRUE)))    
    if( !is.null(colnames(x)) ){
        colnames(df) = colnames(x);
    }
       
    return(df)
    
}


multiState.by.subject <- function(x){

    ## format multistate matrix 'x' by subject.
    ## x == [N subject, T time.points]
    ## df == [N subject x T time.points, K event types]
    ## For example, if we had 3 subjects, 3 time.points, and 3 event types  
    ## the output format by subject 'df' would look like this: 
    ##       K.1  K.2  K.3
    ## 1,t.1 .    .     .
    ## 2,t.1 .    .     .
    ## 3,t.1 .    .     .
    ## 1,t.2 .    1     .
    ## 2,t.2 .    .     .
    ## 3,t.2 .    .     .
    ## 1,t.3 .    1     .
    ## 2,t.3 .    .     1
    ## 3,t.3 .    .     .
    ##-------------------------------------------------
   
    nr    = dim(x)[1] ##number of patients
    nc    = dim(x)[2] ##number of time points
    max.s = max(apply(x,2,max)) ##max event type
    min.s = min(apply(x,2,min)) ##min event type
    ns    = length(min.s:max.s)
    df    = matrix(0,nrow=(nr*nc),ncol=ns)

    for( i in 1:nc ){
        if( i == 1 ){
            start=1
            end=nr
        } else {
            start=end+1
            end=i*nr            
        }
        tmp = do.call(cbind,lapply( 1:ns, function(j) ifelse(j==(x[,i]+1),1,0)))
        df[(start:end),] = tmp
    }

    colnames(df) = min.s:max.s
    
    return(df)
    
}

multiState.by.time.point <- function(x){

    ## format multistae matrix 'x' by time.point.
    ## x == [N subject, T time.points]
    ## df == [N subject x T time.points, K event types]
    ## For example, if we had 3 subjects, 3 time.points and 3 event types 
    ## the output format by subject 'df' would look like this:
    ##       K.1  K.2  K.3
    ## 1,t.1 1    .     .     
    ## 1,t.2 1    .     .
    ## 1,t.3 1    .     .
    ## 2,t.1 .    .     .
    ## 2,t.2 .    1     .
    ## 2,t.3 .    1     .
    ## 3,t.1 .    .     .
    ## 3,t.2 .    .     .
    ## 3,t.3 .    .     1 
    ##-------------------------------------------------
    nr    = dim(x)[1] ##number of patients
    nc    = dim(x)[2] ##number of time points
    max.s = max(apply(x,2,max)) ##max event type
    min.s = min(apply(x,2,min)) ##min event type
    ns    = length(min.s:max.s)
    
    ## Risk matrix
    df <- do.call(rbind, lapply(1:nr, function(j) multiStateMatrix(x[j,],ns)))

    colnames(df) = min.s:max.s
    
   return(df) 
}

##--------------------------------------------------------------------------------------
## Return order list of competeing event ids, given event set
##-------------------------------------------------------------------------------------
get.comp.events <- function(events, cencode=0, print=FALSE){
    ## number of competing events
    n.events = unique(events)
    n.events = n.events[order(n.events)]
    n.events = n.events[-which(n.events==cencode)]
    ne       = length(n.events)  
    if( print ){ cat("no: competing events: ", ne,"\n"); }
    return(n.events)
}

loss1.mask <- function(data.t, data.events, n.events, qt){


    ## get the number of competing events
    ne      = length(n.events)
   
    ## format covariate for dataset - only need to do this part once.
    ## find the risk sets per time point
    data.rs = riskSet(t=data.t$obs.time,qt=qt)
      
    ##------------------------------------------
    ## Mask to calculate LOSS_1 (log-likelihood loss)
    ##------------------------------------------
    ## mask is required to get the log-likelihood loss
    ##   mask size is [N, num_event, num_time.points]
    ##        if not censored : one element = 1 (0 elsewhere)
    ##        if censored     : fill elements with 1 after the censoring time (for all events)
    ##------------------------------------------
    n     = dim(data.rs)[1]
    nt    = dim(data.rs)[2]
    mask1 = array(0,dim=c(n,ne,nt))
    for( i in 1:n ){
        ##inx = max(which(data.rs[i,]==1))
        inx = which(data.rs[i,]==1)
        if(length(inx) > 0 ){
            inx = max(inx)
            if( data.events[i] != 0 ){
                mask1[i,data.events[i],inx]=1
            } else {
                if( inx < nt ){ inx=inx+1; }
                mask1[i,1:ne,inx:nt]=1
            }
        }    
    }
    mask1.by.subject = array(0,dim=c(n*nt,ne))
    for( e in 1:ne ){
        mask1.by.subject[,e] = stack(as.data.frame(mask1[,e,]))[[1]]
    }
    
    return(list(mask1=mask1, mask1.long=mask1.by.subject))
    
}

loss2.mask <- function(data.t, data.events, n.events, qt){

    ## get the number of competing events
    ne      = length(n.events)
   
    ## format covariate for dataset - only need to do this part once.
    ## find the risk sets per time point
    data.rs = riskSet(t=data.t$obs.time,qt=qt)
    
    ##------------------------------------------
    ## Mask to calculate LOSS_2 (ranking loss)
    ##------------------------------------------
    ## mask is required calculate the ranking loss (for pair-wise comparision)
    ##    mask size is [N, num_time.points].
    ##    - For longitudinal measurements:
    ##         1's from the last measurement to the event time (exclusive and inclusive, respectively)
    ##         denom is not needed since comparing is done over the same denom
    ##    - For single measurement:
    ##         1's from start to the event time(inclusive)
    ##------------------------------------------
    mask2 = data.rs
    mask2.by.subject = stack(as.data.frame(mask2))[[1]]

    return(list(mask2=mask2, mask2.long=mask2.by.subject))

}

multiStateMatrix <- function(x,ns){
    smat = do.call(cbind,lapply( 1:ns, function(j) ifelse(j==(x+1),1,0)))
    return(smat)
}


format.multiState <- function(x, events, cencode=0, censtate=FALSE ){

    ## format multistate matrix into a [No: subjects, time points] matrix,
    ## where each non-zero enter represents an event experienced by subject
    ## at that time point.
    ## x      ==> risk set
    ## events ==> event vector for each subject.

    N  = dim(x)[1]
    ns = dim(x)[2]

    n.events = get.comp.events(events=events, cencode=cencode)
    ne       = length(n.events)

    ## should we assign a unique state id for censoring events 
    if( censtate ){
        events[events==cencode] = ne+1
    }
    
    ## reverse the risk.set, i.e. convert every 1 to 0 and vice versa.
    rev.rs = do.call(cbind,lapply(1:ns, function(j) ifelse(x[,j]==1,0,1)))

    ## multistate matrix for each subject at each time point given the subject event
    ## vector.
    ms.mat = do.call(cbind,lapply(1:ns, function(j) ifelse(rev.rs[,j]==1,events,0)))
    
    return(ms.mat)
    
}

riskSet.by.subject <- function(x){

    ## format riskset matrix 'x' by subject.
    ## x == [N subject, T time.points]
    ## df == [N subject x T time.points, T time.points]
    ## For example, if we had 3 subjects, 3 time.points 
    ## the output format by subject 'df' would look like this: 
    ##       t.1  t.2  t.3
    ## 1,t.1 1    .     .
    ## 2,t.1 1    .     .
    ## 3,t.1 1    .     .
    ## 1,t.2 .    1     .
    ## 2,t.2 .    .     .
    ## 3,t.2 .    1     .
    ## 1,t.3 .    .     1
    ## 2,t.3 .    .     .
    ## 3,t.3 .    .     .
    ##-------------------------------------------------
   
    nr = dim(x)[1] ##number of patients
    nc = dim(x)[2] ##number of time points
    df = matrix(0,nrow=(nr*nc),ncol=nc)

    for( i in 1:nc ){
        if( i == 1 ){
            start=1
            end=nr
        } else {
            start=end+1
            end=i*nr            
        }        
        df[(start:end),i] = x[,i]        
    }

    indx = rowSums(df)
    indx = ifelse(indx==1,TRUE,FALSE)
    
    return(list(risk.mat=df,risk.indx=indx))
    
}


riskSet.by.time.point <- function(x,qt){

    ## format riskset matrix 'x' by time.point.
    ## x == [N subject, T time.points]
    ## df == [N subject x T time.points, T time.points]
    ## For example, if we had 3 subjects, 3 time.points 
    ## the output format by subject 'df' would look like this:
    ##       t.1  t.2  t.3
    ## 1,t.1 1    .     .     
    ## 1,t.2 .    1     .
    ## 1,t.3 .    .     1
    ## 2,t.1 1    .     .
    ## 2,t.2 .    .     .
    ## 2,t.3 .    .     .
    ## 3,t.1 1    .     .
    ## 3,t.2 .    1     .
    ## 3,t.3 .    .     . 
    ##-------------------------------------------------
    n = dim(x)[1]
    
    ## Risk matrix
    df <- do.call(rbind, lapply(1:n, function(j) riskMatrix(x[j,])))

    ## Risk matrix index
    indx <- x %>% as.data.table() %>% pivot_longer(cols=colnames(.))
    indx <- as.vector(indx$value)
    indx <- ifelse(indx!=0,TRUE,FALSE)

    return(list(risk.mat=df,risk.indx=indx))
    
}


riskMatrix <- function(x){
    n          <- length(x)
    rmat       <- matrix(0,n,n)
    diag(rmat) <- x
    return(rmat)
}

riskSet <- function(t,qt){

    ## return subject risk set, i.e. a indicator matrix [No: subjects, time points]
    ## highlighting those subjects still at risk at each time point, and risk matrix    
    ## [No: subjects * time points, time points]
    ##qt <- c(0,qt)
    qt <- unique(qt)
    s  <- qt
    s  <- c(qt,(qt[length(qt)]+1));
    n  <- length(t);
    ns <- length(s)-1;

    ## Define number of patients (censored or uncensored) at risk of experiencing
    ## an event at each time point:
    R <- do.call(cbind, lapply(1:ns, function(j) ifelse(s[j] < t, 1, 0)))

    ## 'R*Delta' splits each patients survival-time into the 'ns' time intervals,
    ## i.e. running rowSum(R*Delta) will return each patients survival time.
    Delta   <- do.call(cbind, lapply(1:ns, function(j) pmin(t,s[j+1])-s[j]))

    ## Risk set
    Delta.t <- R*Delta
    Delta.t[Delta.t!=0]=1
     
    return(Delta.t)
    
}

pseudoProb.by.subject <- function(x){

    ## format pseudo survival probability matrix 'x' by subject.
    ## x == [N subject, T time.points]
    ## df == [N subject x T time.points, T time.points]
    ## For example, if we had 3 subjects, 3 time.points 
    ## the output format by subject 'df' would look like this: 
    ##       t.1  t.2  t.3
    ## 1,t.1 1    .     .
    ## 2,t.1 1    .     .
    ## 3,t.1 1    .     .
    ## 1,t.2 .    1     .
    ## 2,t.2 .    .     .
    ## 3,t.2 .    1     .
    ## 1,t.3 .    .     1
    ## 2,t.3 .    .     .
    ## 3,t.3 .    .     .
    ##-------------------------------------------------
   
    nr = dim(x)[1] ##number of patients
    nc = dim(x)[2] ##number of time points
    df = matrix(0,nrow=(nr*nc),ncol=nc)

    for( i in 1:nc ){
        if( i == 1 ){
            start=1
            end=nr
        } else {
            start=end+1
            end=i*nr            
        }        
        df[(start:end),i] = x[,i]        
    }

    return(df)
    
}

pseudoProb.by.time.point <- function(x,qt){

    ## format pseudo probability matrix 'x' by time.point.
    ## x == [N subject, T time.points]
    ## df == [N subject x T time.points, T time.points]
    ## For example, if we had 3 subjects, 3 time.points 
    ## the output format by subject 'df' would look like this:
    ##       t.1  t.2  t.3
    ## 1,t.1 1    .     .     
    ## 1,t.2 .    1     .
    ## 1,t.3 .    .     1
    ## 2,t.1 1    .     .
    ## 2,t.2 .    .     .
    ## 2,t.3 .    .     .
    ## 3,t.1 1    .     .
    ## 3,t.2 .    1     .
    ## 3,t.3 .    .     . 
    ##-------------------------------------------------
    n = dim(x)[1]
    
    ## pseudo prob matrix
    df <- do.call(rbind, lapply(1:n, function(j) riskMatrix(x[j,])))

    return(df)
    
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

cum.prob.at.time.point <- function(x=NULL, qt=NULL, time.point=NULL){

    surv.prob = NULL

    nr = dim(x)[1]
    ne = dim(x)[2]
    
    surv.prob = array(NA,dim=c(nr,ne))

    if( !is.null(x) && !is.null(qt) && !is.null(time.point) ){
    
        inx = which(time.point>=qt)
        if( length(inx) > 0 ){
            inx = max(inx)
            for( e in 1:ne ){
                surv.prob[,e] = apply(y_pred.dh[,e,1:inx, drop=FALSE],1,sum)
            }
        }
    }
    
    return(surv.prob)
    
}

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

