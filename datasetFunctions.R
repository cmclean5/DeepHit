## Generalised function of:
## 1) 'generate.ms.training.data' and 'generate.ms.test.data',
generate.ms.data <- function(rnd_indx, covars, eventset, events,
                             survt, qt,
                             dataset.name=NULL){


    ## format patients obversed time
    t = as.data.frame(do.call(cbind, . %<-% survt))

    ## set data variables
    N           = dim(covars.)[1]
    data        = covars.[rnd_indx,]
    data.t      = t[rnd_indx,]
    data.events = eventset[rnd_indx]

    ## get the number of competing events
    n.events = events
    ne       = length(n.events)    

    ## generate DeepHit Loss1 and Loss2 masks
    m1 = loss1.mask(data.t=data.t, data.events=data.events, n.events=n.events, qt=qt)
    m2 = loss2.mask(data.t=data.t, data.events=data.events, n.events=n.events, qt=qt)    
    
    ## dataset normalisation
    Nsub    = dim(data)[1]
    Nfrac   = round(Nsub/N,3)
    mean    = apply(as.matrix(data), 2, mean)
    sd      = apply(as.matrix(data), 2, sd)
    data.df = scale(data, center=mean, scale=sd)

    ## format covariate for dataset - only need to do this part once.
    ## find the risk sets per time point
    data.rs = riskSet(t=data.t$obs.time,qt=qt)
    
    ## format riskset into multistate model
    data.ms = format.multiState(x=data.rs,events=data.events)

    ## format y.data for deephit, i.e. by competing risk (cr).
    ## [N patients, N events, N time.points]
    ## See Brier Score for details on why using this for y_true
    nn = dim(data)[1]
    nt = length(qt)
    y.data.by.cr = array(0, dim=c(nn,ne,nt))
    for( e in 1:ne ){
        for( t in 1:nt ){
            y.data.by.cr[,e,t] = (data.t$obs.time<=qt[t]) * (data.ms[,t]==e) ##as.numeric(data.ms==e);
        }
    }
                         
 
    return(list(data.t=data.t,data.events=data.events,                
                data.cov=data, data.cov.scale=data.df,
                y.data.by.cr=y.data.by.cr,
                mask1=m1, mask2=m2,
                n=N, nr=Nsub, nfrac=Nfrac,
                data.indices=rnd_indx,
                dataset=dataset.name))   
    
}


train.deephit.model <- function(train, val=NULL, events, qt, params, alpha=NULL, set.nc=TRUE){
  
    ## training data for model
    x_cov   = train$data.cov.scale
    x_t     = train$data.t$obs.time
    x_k     = train$data.events
    x_mask1 = train$mask1$mask1
    x_mask2 = train$mask2$mask2

    y_train = train$y.data.by.cr

    ## validation data for model
    if( !is.null(val) ){
        x_val_cov   = val$data.cov.scale
        x_val_t     = val$data.t$obs.time
        x_val_k     = val$data.events
        x_val_mask1 = val$mask1$mask1
        x_val_mask2 = val$mask2$mask2

        y_val       = val$y.data.by.cr
    } else {
        x_val_cov=NULL; x_val_t=NULL; x_val_k=NULL;
        x_val_mask1=NULL; x_val_mask2=NULL; y_val=NULL;
    }
    
    ## set additional parameters 
    params$n  = dim(x_mask1)[1]
    params$ne = dim(x_mask1)[2]
    params$nt = dim(x_mask1)[3]
    params$nc = dim(x_cov)[2]

    ## set no. of neurons to the number of covariates
    if( set.nc ){
        params$neurons = params$nc
        cat("> neurons = ", params$neurons, " = no: covaraite = ", params$nc,"\n")
    }
    
    model   = NULL   
    
    ##------------------------------------------------------------------
    ## Build & Train DeepHit model
    ##------------------------------------------------------------------
    model = deephit.model(x_cov       = x_cov,
                          x_k         = x_k,
                          x_t         = x_t,
                          x_mask1     = x_mask1,
                          x_mask2     = x_mask2,
                          y_train     = y_train,
                          x_val_cov   = x_val_cov,
                          x_val_k     = x_val_k,
                          x_val_t     = x_val_t,
                          x_val_mask1 = x_val_mask1,
                          x_val_mask2 = x_val_mask2,
                          y_val       = y_val,
                          qt          = qt,                          
                          params      = params,
                          alpha       = alpha
                          )
                            
        

    return(model) 
    
}


test.deephit.model <- function(model, test){
    
    ## test data for model
    x_cov   = test$data.cov.scale
    x_t     = test$data.t$obs.time
    x_k     = test$data.events
    x_mask1 = test$mask1$mask1
    x_mask2 = test$mask2$mask2

    test.events=test$data.events 
    
    ypred = deephit.predict(model   = model,
                            x_cov   = x_cov,
                            x_t     = x_t,
                            x_k     = x_k,
                            x_mask1 = x_mask1,
                            x_mask2 = x_mask2,
                            y_test  = test.events)   

    return(ypred)
        
}


eval.deephit.model <- function(model,
                               eval,
                               train=NULL,
                               qt,
                               eval.times,
                               n.events,
                               y_pred=NULL){

    ci=NULL
    ciCI=NULL
    br=NULL
    br.null=NULL
    br.we=NULL
    br.scaled=NULL
    br.tmp=list()
    
    if( is.null(y_pred) ){
        y_pred = test.deephit.model(model = model,
                                    test  = eval
                                    )
    } else {

        xdim = dim(y_pred)
        nr   = xdim[1]
        ne   = xdim[2]
        nt   = length(eval.times)##xdim[3]

        ci   = array(dim=c(ne,nt))
        ciCI = array(dim=c(nt,3,ne))
        dimnames(ciCI)[[2]] = c("ci","upper","lower")

        ci.dh   = array(dim=c(ne,nt))
        ci.dhCI = array(dim=c(nt,3,ne))
        dimnames(ci.dhCI)[[2]] = c("ci","upper","lower")
        
        br   = array(dim=c(ne,nt))
        dimnames(br)[[1]] = seq(1,ne,1)
        dimnames(br)[[2]] = eval.times

        br.we   = array(dim=c(ne,nt))
        dimnames(br.we)[[1]] = seq(1,ne,1)
        dimnames(br.we)[[2]] = eval.times

        br.null   = array(dim=c(ne,nt))
        dimnames(br.null)[[1]] = seq(1,ne,1)
        dimnames(br.null)[[2]] = eval.times

        br.scaled = array(dim=c(ne,nt))
        dimnames(br.scaled)[[1]] = seq(1,ne,1)
        dimnames(br.scaled)[[2]] = eval.times     

          
        for( s in 1:nt ){

            ## Calculate the Brier score under the null model (no covariates) using the
            ## Aalen-Johansen estimator
            br.null[,s] = AJ_comp(time=eval$data.t$obs.time,
                                  delta=eval$data.t$deltaC,
                                  events=eval$data.events,
                                  tau=eval.times[s])
            
            ## calculate cumulative probability at time point.
            prob.dh = cum.prob.at.time.point(x=y_pred, qt=qt, time.point=eval.times[s])

            for( e in 1:ne ){
                
                ci[e,s] = Cindex_estimator_efficient(
                    Prediction=prob.dh[,e],
                    Time_survival=eval$data.t$obs.time,
                    Censoring=eval$data.t$deltaC,
                    Cause=eval$data.events,
                    Cause_int=e,
                    Time=eval.times[s],
                    method="hazard")

                tmp = Cindex_estimator_efficient_CI(
                    Prediction=prob.dh[,e],
                    Time_survival=eval$data.t$obs.time,
                    cindex=ci[e,s],
                    method="hazard")
                
                for( i in 1:length(tmp) ){ ciCI[s,i,e] = tmp[[i]]; } 
                
                ci.dh[e,s] = Cindex_deephit(
                    Prediction=prob.dh[,e],
                    Time_survival=eval$data.t$obs.time,
                    Cause=eval$data.events,
                    Cause_int=e,
                    Time=eval.times[s]
                )
                
                tmp = Cindex_estimator_efficient_CI(
                    Prediction=prob.dh[,e],
                    Time_survival=eval$data.t$obs.time,
                    cindex=ci.dh[e,s],
                    method="hazard")
                
                for( i in 1:length(tmp) ){ ci.dhCI[s,i,e] = tmp[[i]]; } 
                
                br[e,s] = brier_score(
                    Prediction=prob.dh[,e],
                    Time_survival=eval$data.t$obs.time,
                    Cause=eval$data.events,
                    Cause_int=e,
                    Time=eval.times[s]
                )

                br.we[e,s] = weighted_brier_score(Train_Time=train$data.t$obs.time,
                                                  Train_DeltaC=train$data.t$deltaC,
                                                  Time=eval.times[s],
                                                  Time_survival=eval$data.t$obs.time,
                                                  Cause=eval$data.events,
                                                  Cause_int=e,
                                                  Prediction=prob.dh[,e],
                                                  calBS=FALSE)[[2]]
                 
            }
        }
                

    }          

    return(list(ci=ci, ciCI=ciCI, ci.dh=ci.dh, ci.dhCI=ci.dhCI, br=br, br.null=br.null,
                br.we=br.we, br.scaled=br.scaled))
    
}

