## Where the cumulative incidence function (CIF) method is either:
## 'aj' => Aalen-Johannsen estimator
## 'fg' => Fine & Gray estimator
generate.ds.data <- function(rnd_indx, covars, eventset, events,
                             survt, qt, x.indx, x.strat, x.we, alpha,
                             method=c("aj","fg","pseudo"),
                             type=c("unwe","age.strat","rf","class"),
                             dataset.name=NULL,
                             use.risk.set=FALSE,
                             use.pseudo.values=TRUE){

    ## method selection
    method = match.arg(method);
    type   = match.arg(type);

    desc      = c("aj"="Aalen-Johannsen",
                  "fg"="Fine&Gray",
                  "pseudo"="R's pseudo func",
                  "unkwn"="unknown")
    desc.indx = which(names(desc)==method)
    if( length(desc.indx) == 0 ){ desc.indx=4; }

    type.desc = c("unwe"="uncensored-weighted","age.strat"="age-stratified",
                  "rf"="random-forest","class"="class-weighted", "unkwn"="unkown")
    type.indx = which(names(type.desc)==type)
    if( length(type.indx) == 0 ){ type.indx=5; }
    
    t = as.data.frame(do.call(cbind, . %<-% survt))

    ## set data variables
    N           = dim(covars.)[1]
    data        = covars.[rnd_indx,]
    data.t      = t[rnd_indx,]
    data.events = eventset[rnd_indx]

    ## get the number of competing events
    n.events = events
    ne       = length(n.events)    
    
    ## compute un/weighted pseudo values for all landmark times    
    res = get.we.comp.survival.time(t=data.t,
                                    events=data.events,
                                    qt=qt,
                                    covar=data,
                                    x.indx=x.indx,
                                    x.strat=x.strat,
                                    x.we=x.we,
                                    alpha=alpha,
                                    type=type,
                                    method=method)

    pseudoComp_norm = res$pseudoNorm
   
    ## dataset normalisation
    Nsub    = dim(data)[1]
    Nfrac   = round(Nsub/N,3)
    mean    = apply(as.matrix(data), 2, mean)
    sd      = apply(as.matrix(data), 2, sd)
    data.df = scale(data, center=mean, scale=sd)

    ## format covariate for dataset - only need to do this part once.
    ## find the risk sets per time point
    data.rs         = riskSet(t=data.t$obs.time,qt=qt)
    data.rs.subject = riskSet.by.subject(data.rs)
    data.df.subject = covariates.by.subject(x=data.df,qt=qt)
    
    ## multistate
    data.ms            = format.multiState(x=data.rs,events=data.events)
    data.ms.by.subject = multiState.by.subject(x=data.ms)
    y.data.all         = data.ms.by.subject
    
    ## the multistate dataset
    nn = dim(data)[1]
    nt = length(qt)
    nr = dim(data.df.subject)[1]
    nc = dim(data.df.subject)[2]+nt#length(qt)
    np = dim(pseudoComp_norm)[2]
    x.data.all    = array(dim=c(nr,nc,ne))
    x.data.pseudo = array(dim=c(nr,np,ne))
     if( use.pseudo.values ){
         for( e in 1:ne ){
             data.pp.subject    = pseudoProb.by.subject(x=pseudoComp_norm[,,e])
             x.data.pseudo[,,e] = data.pp.subject
             x.data.all[,,e]    = cbind(data.df.subject, data.pp.subject)
         }
     } else {

         data.s  = rep(qt,each=nrow(data.df))
         smatrix = model.matrix(~as.factor(data.s)+0)
         
          for( e in 1:ne ){
              x.data.all[,,e] = cbind(data.df.subject, smatrix)
          }
     }
                      

    return(list(x.data=x.data.all, y.data=y.data.all,                
                x.data.cov=data.df.subject,
                x.data.pseudo=x.data.pseudo,
                data.t=data.t,data.events=data.events,                
                data.cov=data, data.cov.scale=data.df,
                pseudo.we=res$weights,               
                pseudo.probs=pseudoComp_norm,
                n=N, nr=Nsub, nfrac=Nfrac,
                method=desc[desc.indx], weight=type.desc[type.indx],
                dataset=dataset.name))
    
}


train.deepsurv.model <- function(train, val=NULL, events,
                                 qt,                           
                                 PARAMS,
                                 alpha=NULL,
                                 logdir=NULL){
    

    ## load training data
    x.train=train$x.data
    y.train=train$y.data
   
    ## validation data
    x.val = NULL
    y.val = NULL
    if( !is.null(val) ){
        x.val = val$x.data
        y.val = val$y.data
    }

    if( is.null(logdir) ){ logdir = tempdir(); }
    
    ## get the nnumber of competing events
    n.events = events
    ne       = length(n.events)

    model    = NULL
    
    ##------------------------------------------------------------------
    ## Train multistate model
    ##------------------------------------------------------------------
    model = multiState.func.model(x_train=x.train,
                                  y_train=y.train,
                                  x_val=x.val,
                                  y_val=y.val,
                                  params=PARAMS,
                                  alpha=alpha,
                                  logdir=logdir)    
        

    return(model)
}


test.deepsurv.model <- function(model, test, qt){
                      

    ## load test data
    x.test=test$x.data
    test.t=test$data.t
    test.events=test$data.events      
    nr = test$nr   

    ## get the number of competing events  
    nt = length(qt)   

    yprob = NULL
    ytemp = multiState.predict(model=model,
                               x_test=x.test)
    
    ypred = get.pred.multiState(x=ytemp, qt=qt, nr=nr)
    yprob = ypred

    return(yprob)

}


eval.deepsurv.model <- function(model=model$model,
                                eval.data,
                                eval.times,                                
                                y_pred,
                                acc.tests =
                                    c("Recall","Specificity","Precision",
                                      "F1","Balanced Accuracy","mcc")){

    
    if( length(acc.tests) == 1 ){
        test.metric = match.arg(acc.tests);
    } else {
        test.metric = acc.tests;
    }

    nt = length(eval.times)
    ns = dim(y_pred)[2]
    state.ids = seq(0,(ns-1),1)
    Ntests = length(test.metric)
    
    eval = array(dim=c(ns,Ntests,nt))
    dimnames(eval)[[1]] = state.ids
    dimnames(eval)[[2]] = test.metric
    dimnames(eval)[[3]] = eval.times

    events = eval.data$data.events
    
    ## using library caret 
    for( i in 1:nt ){            
        pred.events = apply(y_pred[,,i],1,which.max)-1
        cm=confusionMatrix(factor(pred.events),
                               factor(events), mode="everything")
        res=cm$byClass
        for( j in 1:Ntests ){
            indx=which(colnames(res)==test.metric[j])
            if( length(indx) != 0 ){
                eval[,j,i]=res[,indx]
            }
        }
    }
    
    return(eval)
    
}
