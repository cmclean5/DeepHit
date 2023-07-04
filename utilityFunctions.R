##-------------------------------------------
## load functions for preprocessing and analysis cardiac data
##-------------------------------------------


model.summary.stats <- function(model=NULL, data=NULL,
                                method=c("reReg","coxph"),
                                show.forest=TRUE, forest.title=NULL){

    method = match.arg(method)
    
    ## return the recurrent event 'r.event' object (if reReg model),
    ## and terminal event 't.event' object of coxph or reReg model
    r.event     = NULL;
    t.event     = NULL;
    forest.plot = NULL;
    
    if( !is.null(model) && method == "reReg" ){

        summ = summary(model)

        r.event = summ$coefficients.rec
        t.event = summ$coefficients.haz
        
    }

    if( !is.null(model) && method == "coxph" ){

        summ = summary(model)

        t.event = summ$coefficients

        if( show.forest ){
            if( is.null(forest.title) ){ forest.title="Hazard ratio"; }

            if( is.null(data) ){
                forest.plot = ggforest(model, main=forest.title)
            } else {
                forest.plot = ggforest(model, data=data, main=forest.title)
            }

        }
        
    }

    return(list(r.event=r.event,t.event=t.event,forest.plot=forest.plot))
    
}

format.age.covar <- function(data=NULL,age.label="age",nbins=5,factor.covar=TRUE,factor.int="1",plot.age=TRUE){

    binned.age = NULL
    p=NULL
    
    if( !is.null(data) ){

        age.indx   = which(colnames(data)==age.label)
        binned.age = group.data(x=data[,age.indx[1]],nbins=nbins, method="range")
        if( factor.covar == TRUE ){
            factor.age = factor2ind(binned.age$groups.int,factor.int)
            colnames(factor.age) = sprintf("%s_%s",age.label,binned.age$groups.map[2:nbins,1])
        
            ## remove unbinned age covariate and add binned covariate
            data = data[-age.indx[1]]
            data = cbind(data,factor.age)
            
        }

        if( plot.age == TRUE ){

            ##h=hist(binned.age$groups.int,breaks=nbins,plot=FALSE)

            df           = as.data.frame(binned.age$groups.int)
            labels       = binned.age$groups.map[1:nbins,1]
            colnames(df) = "years"
            gplot = ggplot(df, aes(years))+geom_bar()+scale_x_discrete(limit=labels) 
            
            #if( factor.covar == TRUE ){
            #    p=plot(h,xaxt="n",col="grey",main="",xlab="years",ylab="count")+
            #        axis(1,h$mids,labels=binned.age$groups.map[2:nbins,1],tick=TRUE,plot=FALSE)
            #} else {
            #    p=plot(h,xaxt="n",col="grey",main="",xlab="years",ylab="count")+
            #        axis(1,h$mids,labels=binned.age$groups.map[1:nbins,1],tick=TRUE,plot=FALSE)
            #}
        }
        
    }

    return(list(data=data,binned.age=binned.age,nbins=nbins,factor.int=factor.int,gplotp=gplot))
    
}

cardiac.time.count <- function(x=NULL,qt=NULL,startYear=NULL){

    groups.tb = NULL
    
    if( !is.null(x) && !is.null(qt) && !is.null(startYear) ){

        x = unlist(lapply(x,splitRecurrence))
        x = x[!is.na(x) & x!=""]
        x = as.numeric(x)
        x = x - as.numeric(startYear)

        breaks    = c(0,qt)
        breaks    = unique(breaks)
        groups    = cut((x-1),
                        breaks=breaks,
                        include.lowest=TRUE,
                        right=FALSE)
        groups.tb = table(groups)
        groups.tb = as.vector(groups.tb)
    }

    return(groups.tb)

}

group.data <- function(x=NULL, nbins=5, method=c("custom","range","centile") ){

    method     =  match.arg(method)
    
    groups     = NULL
    groups.int = NULL
    groups.tb  = NULL
    groups.map = NULL
    
    if( !is.null(x) ){

        x = as.numeric(x)
        
        Xmax=max(x,na.rm=T)
        Xmin=min(x,na.rm=T)

        if( method == "range" ){
            int    = (Xmax-Xmin)/nbins
            breaks = seq(Xmin,Xmax,int)
        }

        if( method == "centile" ){
            breaks = unique(c(0,nbins))
        }

        if( method == "custom" ){
            breaks = nbins
        }
        
        groups     = cut(x,breaks=breaks,include.lowest=TRUE,right=FALSE)
        groups.tb  = table(groups)
        labels     = names(groups.tb)
        groups     = as.vector(groups)
        groups.int = match(groups,labels)
        groups.map = cbind(labels,seq(1,length(labels),1))
    }

    return(list(groups=groups,groups.int=groups.int,
                groups.tb=groups.tb,groups.map=groups.map))
    
}

survivalTime <- function(start.time=NULL, end.time=NULL, c.time=NULL,
                         event.time=NULL, outcomes=NULL, outcome.e="2",
                         YEAR=365.25, ROUND=4, method="first"){

    obs.time = NULL
    deltaC   = NULL

    if( !is.null(start.time) && !is.null(end.time) &&
        !is.null(c.time)     && !is.null(event.time) &&
        !is.null(outcomes)  ){

        ## number of patients
        N = length(start.time)

        ## convert start date
        start.time = as.vector(sapply(start.time,getDate))

        ## convert end date
        end.time   = as.vector(sapply(end.time,getDate))

        ## convert multiple cardiac event dates
        event.time = as.vector(sapply(event.time, function(x)
            appendDates("",as.character(as.numeric(getDate(splitRecurrence(x)))))))        
        
        ## convert censor date
        c.time     = as.numeric(getDate(c.time))

        ## condition to make sure no negative time
        max.start  = max(start.time,na.rm=T)
        c.time     = ifelse( c.time < max.start, max.start, c.time)
        
        ## survival time in years
        ##Tstar = end.time - start.time
        Tstar = survival.time(start=start.time, end=end.time,
                              events=event.time, outcomes=outcomes,
                              outcome.e=outcome.e, YEAR=YEAR, method=method)
        
        ## censor time 
        C = round( (rep(c.time,N) - start.time)/YEAR, ROUND)

        ## observed survial time
        obs.time = apply(cbind(Tstar,C), 1, function(x) min(x,na.rm=T) )

        ## censor indicator
        deltaC = ifelse(Tstar <= C,1,0)
        deltaC = ifelse(is.na(deltaC),0,deltaC)
       
    }

    return(list(Tstar=Tstar, C=C, obs.time=obs.time,deltaC=deltaC))
    
}

## function will generate a one-hot encoding of vector 'x'
## with either the 'baseline' varaible removed,
## or the first variable removed. 
factor2ind <- function(x, baseline){

## 
## Given a factor variable x, create an indicator matrix of dimension 
## length(x) x (nlevels(x)-1) dropping the column corresponding to the 
## baseline level (by default the first level is used as baseline).
## Example:
## > x = gl(4, 2, labels = LETTERS[1:4])
## > factor2ind(x)
## > factor2ind(x, "C")
  xname <- deparse(substitute(x))
  n <- length(x)
  x <- as.factor(x)
  if(!missing(baseline)) x <- relevel(x, baseline)
  X <- matrix(0, n, length(levels(x)))
  X[(1:n) + n*(unclass(x)-1)] <- 1
  X[is.na(x),] <- NA
  dimnames(X) <- list(names(x), paste(xname, levels(x), sep = ":"))
  return(X[,-1,drop=FALSE])
}


check.missing.values <- function(data=NULL,id=NULL,pred.var=NULL){

    missing = NULL

    if( !is.null(data) && !is.null(id) && !is.null(pred.var) ){

        id.indx   = match(id,colnames(data))
        pred.indx = match(pred.var,colnames(data))
        pred.indx = pred.indx[!is.na(pred.indx)]

        if( !is.na(id.indx) && all(!is.na(pred.indx)) ){

            N      = length(pred.indx)
            cnames = c("variable","patient.id","row.indx")
            missing <- setNames(data.frame(matrix(ncol = length(cnames), nrow = 0)), cnames)
            
            for( i in 1:N ){

                val  = as.numeric(data[,pred.indx[i]])
                indx = which(is.na(val)==TRUE)
                Nt   = length(indx)
                ID   = as.character(data[indx,id.indx[1]])
                
                temp.table <- setNames(data.frame(matrix(ncol = length(cnames), nrow = Nt)), cnames)
                temp.table[,1] = rep(colnames(data)[pred.indx[i]],Nt)
                temp.table[,2] = ID
                temp.table[,3] = as.numeric(indx)

                missing = rbind(missing,temp.table)
            }
            
        }
            
    }

    return(missing)
    
}

cohort.summary.stats <- function(data=NULL, strat.var=NULL, strat.value=NULL,
                                 pred.var=NULL, pred.meth=NULL, ROUND=2, PREC=TRUE ){

    summ = NULL
    
    if( !is.null(data) && !is.null(strat.var) && !is.null(strat.value) &&
        !is.null(pred.var) && !is.null(pred.meth) ){

        if( length(strat.var) == length(strat.value) &&
            length(pred.var) == length(pred.meth) ){            
            
            strat.indx = match(strat.var,colnames(data))
            pred.indx  = match(pred.var,colnames(data))

            strat.indx = strat.indx[!is.na(strat.indx)]
            pred.indx  = pred.indx[!is.na(pred.indx)]

            if( all(!is.na(strat.indx)) && all(!is.na(pred.indx)) ){

                N      = dim(data)[1]
                
                Nc     = length(strat.indx)
                cnames = unlist(lapply(1:Nc,function(j) rep(colnames(data)[strat.indx][j],2)))
                cnames[seq(2,2*Nc,2)]=""

                Nr     = length(pred.var)
                rnames = "number"
                
                summ <- setNames(data.frame(matrix(ncol = length(cnames), nrow = 1)), cnames)

                ## patient number in each stratification
                k=1;
                for( i in 1:Nc ){
                    val = as.numeric(data[,strat.indx[i]])==strat.value[i]
                    val = sum(as.numeric(val),na.rm=T)
                    summ[1,k]     = val
                    summ[1,(k+1)] = ifelse( PREC, round(val/N,ROUND)*100, round(val/N,ROUND) )
                    k=k+2
                }


                for( i in 1:Nr ){

                    if( pred.meth[i] == "n" ){
                    
                        tmp    = table(as.numeric(data[,pred.indx[i]]))
                        ##total  = sum(tmp)
                        Nt     = length(tmp)
                        rnames = c(rnames,sprintf("%s_%s",colnames(data)[pred.indx[i]],names(tmp)))
                        
                        temp.table <- setNames(data.frame(matrix(ncol = length(cnames), nrow = Nt)), cnames)
                    
                        for( t in 1:Nt ){
                            k=1;
                            for( j in 1:Nc ){
                                val = data[,strat.indx[j]]==strat.value[j] & data[,pred.indx[i]]==names(tmp)[t]
                                val = sum(as.numeric(val),na.rm=T)
                                temp.table[t,k]     = val
                                temp.table[t,(k+1)] = val
                                k=k+2
                            }
                        }

                        ##normalise
                        k=2;
                        for( j in 1:Nc ){
                            val = temp.table[,k]/sum(temp.table[,k],na.rm=T)
                            if( PREC ){ val = val*100; }
                            temp.table[,k] = round(val,ROUND)
                            k=k+2;
                        }
                        
                    }

                    if( pred.meth[i] == "mean" ){

                        rnames = c(rnames,colnames(data)[pred.indx[i]])
                        
                        temp.table <- setNames(data.frame(matrix(ncol = length(cnames), nrow = 1)), cnames)

                        k=1;
                        for( j in 1:Nc ){
                            val  = as.numeric(data[,pred.indx[i]])
                            indx = data[,strat.indx[j]]==strat.value[j]                            
                            temp.table[1,k]     = round(mean(val[indx],na.rm=T),ROUND)
                            temp.table[1,(k+1)] = round(sd(val[indx],na.rm=T),ROUND)
                            k=k+2
                        }                                                
                    }
                        
                        summ = rbind(summ,temp.table)
                    
                }

                rownames(summ) = rnames

            }
            
        }

    }

    return(summ)
    
}

build.event.model <- function(data=NULL, id.cn="id", predictors=NULL, values=NULL,
                              methods=NULL, filter=c("found","id"),
                              end.cn="end", event.cn="event", death.cn="death"){

    set = NULL

    filter = match.arg(filter);
    
    if( !is.null(data) && !is.null(values) && !is.null(predictors) &&
       !is.null(methods) ){

        ## set of predictors to build event model. First predictor
        ## in set always name of the outcome variable
        N         = length(predictors)
        Nvals     = length(values)
        Nmeths    = length(methods)
        
        if( N == Nvals && N == Nmeths && N > 0 ){

            end.indx   = which(colnames(data)==end.cn)
            event.indx = which(colnames(data)==event.cn)
            death.indx = which(colnames(data)==death.cn)
            
            outcome = which(colnames(data)==predictors[1])
            id.indx = which(colnames(data)==id.cn)

            outcome.val = as.numeric(values[1])
            
            if( length(outcome) != 0 ){

                found     = data[,outcome[1]]==outcome.val
                
                pred.indx = match(predictors,colnames(data))
                pred.indx = pred.indx[-which(pred.indx==outcome)]
                values    = values[-1]
                methods   = methods[-1]
                
                if( !is.na(pred.indx) && length(pred.indx) != 0 ){
                    for( i in 1:length(pred.indx) ){            
                        
                        switch(methods[i],
                               "eq"={
                                   found = found &
                                       data[,pred.indx[i]]==as.numeric(values[i])   
                               },
                               "lt"={
                                   found = found &
                                       data[,pred.indx[i]]<as.numeric(values[i])   
                               },
                               "gt"={
                                   found = found &
                                       data[,pred.indx[i]]>as.numeric(values[i])   
                               },
                               "lteq"={
                                   found = found &
                                       data[,pred.indx[i]]<=as.numeric(values[i])   
                               },
                               "gteq"={
                                   found = found &
                                       data[,pred.indx[i]]>=as.numeric(values[i])  
                               }
                               )
                    }
                }

                switch(filter,
                       "found"={
                           set  = data[found,]
                           indx = !is.na(set[,end.indx[1]]) & set[,event.indx[1]]==1
                           ids  = unique(set[indx,id.indx[1]])

                           for( i in 1:length(ids) ){
                               tmp  = set[which(set[,id.indx[1]]==ids[i]),
                                          c(end.indx[1],event.indx[1],death.indx[1])]
                               Ntmp = dim(tmp)[1]
                               if( Ntmp > 0 ){
                                   if( !is.na(tmp[Ntmp,1]) && tmp[Ntmp,2] == 1 ){
                                       tmp.end = as.numeric(rownames(tmp)[Ntmp])
                                       indx    = which(rownames(set)==tmp.end)
                                       set[indx,end.indx[1]]   = NA
                                       set[indx,event.indx[1]] = 0
                                   }
                               }
                           }
                       },
                       "id"={
                           ids     = unique(data[found,id.indx[1]])
                           set.ids = match(data[,id.indx[1]],ids)
                           indx    = ifelse(is.na(set.ids),FALSE,TRUE)
                           set     = data[indx,]
                       },

                       {
                         ids     = unique(data[found,id.indx[1]])
                         set.ids = match(data[,id.indx[1]],ids)
                         indx    = ifelse(is.na(set.ids),FALSE,TRUE)
                         set     = data[indx,]   
                       }
                       )

                
            }
        }
    }

    return(set)
    
}

merge.event.dates <- function( data=NULL, COLLAPSE=";",
                              style=c("%Y-%m-%d","%d-%b-%Y","%d-%b-%y") ){

    result = NULL

    if( !is.null(data) && is.list(data) ){

        Nsets = length(data)
        N     = length(data[[1]])
        
        result = matrix(NA,ncol=3,nrow=N)

        for( i in 1:N ){
            events = ""; events.d = ""; events.y = "";
            for( j in 1:Nsets ){
                events = c(events,splitRecurrence(data[[j]][i]))
            }
            events   = events[!is.na(events)]
            events   = events[events!=""]
            if( length(events) > 0 ){
                events   = unique(events)
                events.d = as.numeric(getDate(events,style=style))
                events.y = getYear(events,style=style)
            
                if( length(events) > 1 ){
                    events   = paste(as.character(events),collapse=COLLAPSE)
                    events.d = paste(as.character(events.d),collapse=COLLAPSE)
                    events.y = paste(as.character(events.y),collapse=COLLAPSE)
                }
            
                result[i,1] = events
                result[i,2] = events.d            
                result[i,3] = events.y
            }
        }
    }

    return(result)
    
}

addPredictors <- function(data=NULL,rc.table=NULL,map=NULL,
                          predictors=NULL, data.id.name="id",
                          rc.table.id.name="id" ){

    ## add any predictor variables we want to recurrent events table

    res = rc.table
    
    if(!is.null(data) && !is.null(rc.table) && !is.null(map) &&
       !is.null(predictors) ){

        base.cn     = colnames(data)
        base.id     = which(base.cn==data.id.name)
        
        rc.cn   = colnames(rc.table)
        rc.id   = which(rc.cn==rc.table.id.name)
        
        to.add  = match(predictors,base.cn)
        to.add  = to.add[!is.na(to.add)]

        if( length(to.add) > 0  ){

            id.indx  = match(data[,base.id[1]],map[,1])
            map.indx = match(rc.table[,rc.id[1]],id.indx)

            append.cn = c()
            
            for( i in 1:length(to.add) ){
                if( base.cn[to.add[i]] != data.id.name ){
                    rc.table = cbind(rc.table,data[map.indx,to.add[i]])
                    append.cn = c(append.cn,base.cn[to.add[i]])
                }
            }

            colnames(rc.table) = c(rc.cn,append.cn)
            
        }
        
    }

    return(rc.table)
    
}

recurrentEventTable <- function(id=NULL,start=NULL,end=NULL,events=NULL,outcome=NULL,
                                style=c("%Y-%m-%d","%d-%b-%Y","%d-%b-%y"),
                                method=c("days","weeks","months","years") ){

    method =  match.arg(method);
    cnames = c("id","start","end","time","acc.time","enum","event","outcome","death")

    rc.table = data.frame( id=as.character(),start=as.numeric(),
                          end=as.numeric(),time=as.numeric(),
                          acc.time=as.numeric(), enum=as.numeric(),
                          event=as.numeric(),outcome=as.numeric(),
                          death=as.numeric() )

    id.map = NULL
    CENSOR = NULL
    YEAR   = 365.25
    DEM    = 1
    
    switch(method,
           "days"={ DEM=1;
           },
           "weeks"={ DEM=YEAR/52;
           },
           "months"={ DEM=YEAR/12;
           },
           "years"={ DEM=YEAR;
           },

           {
               DEM=1;
           }
    )
    
    if( !is.null(id) && !is.null(start) && !is.null(end) &&
        !is.null(events) && !is.null(outcome) ){

        N       = length(events)
        death   = ifelse(!is.na(end),1,0)

        id.map     = matrix(0,ncol=2,nrow=N)
        id.map[,1] = id
        id.map[,2] = seq(1,length(id),1)        
        id         = id.map[,2]
    
        start = as.numeric(sapply(start,occurrenceNumeric,style=style))                                    
        end   = as.numeric(sapply(end,occurrenceNumeric,style=style))                                    

        CENSOR = max(end,na.rm=T)
        
        for( i in 1:N ){
          
            e = splitRecurrence(events[i])

            if( is.na(e) || length(e) == 0 || e == "" ){

                tmp.table <- setNames(data.frame(matrix(ncol = length(cnames), nrow = 1)), cnames)
                
                tmp.table[,1] = id[i]
                tmp.table[,2] = start[i]
                tmp.table[,3] = end[i]
                tmp.table[,4] = abs(ifelse(!is.na(end[i]),((end[i] - start[i])/DEM), ((CENSOR-start[i])/DEM) ))
                tmp.table[,5] = tmp.table[,4]
                tmp.table[,6] = 1
                tmp.table[,7] = 0
                tmp.table[,8] = outcome[i]
                tmp.table[,9] = ifelse(!is.na(end[i]),1,0)
            } else {

                Ne = length(e)

                if( Ne == 1 ){

                    Ne=Ne+1

                    tmp.table <- setNames(data.frame(matrix(ncol = length(cnames), nrow = Ne)), cnames)
                
                    id.i    = rep(id[i],Ne)                
                    start.i = rep(0,Ne)
                    end.i   = rep(0,Ne)
                    time.i  = rep(0,Ne)                              
                    acc.i   = rep(0,Ne)
                    e.i     = as.vector(sapply(e,occurrenceNumeric,style=style))
                    e.i     = e.i[order(e.i)]
                    o.i     = rep(outcome[i],Ne)
                    death.i = rep(0,Ne)
                    if( !is.na(end[i]) ){ death.i[Ne] = 1; }
                    
                    start.i[1] = start[i]
                    end.i[1]   = e.i[1]
                    time.i[1]  = abs(((end.i[1]-start.i[1])/DEM))
                    acc.i[1]   = time.i[1]
                    
                    start.i[2] = e.i[1]
                    end.i[2]   = end[i]
                    time.i[2]  = abs(ifelse(!is.na(end.i[2]),((end.i[2] - start.i[2])/DEM), ((CENSOR-start.i[2])/DEM) ))
                    acc.i[2]   = acc.i[1] + time.i[2]
                    
                } else {
                    
                    tmp.table <- setNames(data.frame(matrix(ncol = length(cnames), nrow = Ne)), cnames)
                
                    id.i    = rep(id[i],Ne)                
                    start.i = rep(0,Ne)
                    end.i   = rep(0,Ne)
                    time.i  = rep(0,Ne)
                    acc.i   = rep(0,Ne)
                    e.i     = as.vector(sapply(e,occurrenceNumeric,style=style))
                    e.i     = e.i[order(e.i)]
                    o.i     = rep(outcome[i],Ne)
                    death.i = rep(0,Ne)
                    if( !is.na(end[i]) ){ death.i[Ne] = 1; }
                    
                
                    for( j in 1:(Ne-1) ){
                        if( j == 1 ){
                            start.i[j] = start[i]
                            end.i[j]   = e.i[j]
                        } else {
                            start.i[j] = e.i[j]
                            end.i[j]   = e.i[j+1]
                        }
                        time.i[j]  = abs(((end.i[j]-start.i[j])/DEM))
                    }
                    
                    start.i[Ne] = e.i[Ne]
                    end.i[Ne]   = end[i]
                    time.i[Ne]  = abs(ifelse(!is.na(end.i[Ne]),((end.i[Ne] - start.i[Ne])/DEM),((CENSOR-start.i[Ne])/DEM) ))

                    acc.i[1] = time.i[1]
                    for( j in 2:Ne ){ acc.i[j] = acc.i[j-1] + time.i[j]; }
                    
                }

                event.c = rep(1,Ne); event.c[Ne] = 0;
                
                tmp.table[,1] = id.i
                tmp.table[,2] = start.i
                tmp.table[,3] = end.i
                tmp.table[,4] = ifelse(is.na(time.i),0,time.i)
                tmp.table[,5] = acc.i
                tmp.table[,6] = seq(1,Ne,1);
                tmp.table[,7] = event.c;
                tmp.table[,8] = o.i
                tmp.table[,9] = death.i                
                
            }


            
            rc.table = rbind(rc.table,tmp.table)
            
        }

        
    }
    
    return(list(e.table=rc.table,map=id.map,censor=CENSOR))

}

splitRecurrence <- function(x,COLLAPSE=";"){

    x = as.character(x)

    if( !is.na(x) ){
        if( sum(grepl(COLLAPSE,x))>0 ){
            return(strsplit(x,COLLAPSE)[[1]])
        }
    }

    return(x)
    
}

countRecurrence <- function(x,COLLAPSE=";"){

    ## return number of events

    tally = 0

    x = splitRecurrence(x,COLLAPSE=COLLAPSE)

    x = x[x!=""]
    if( length(x) != 0 && is.na(x) == FALSE ){
        tally = length(x)
    }
    
    return(tally)

}


firstOccurrence <- function(x,COLLAPSE=";",method=c("first","last") ){

    ## find occurrence of first/last event/year

    event = ""

    method = match.arg(method)
    
    if( !is.na(x) && length(x) != 0 && !is.null(x) ){

        x = splitRecurrence(x,COLLAPSE=COLLAPSE)
        
        x = x[x!=""]
        if( length(x) != 0 ){

            switch(method,
                   "death.min" = { event = min(x,na.rm=T)
                    },
                   "death.max" = { event = max(x,na.rm=T)
                   },
                   "first"     = { event = min(x,na.rm=T)
                   },
                   "last"      = { event = max(x,na.rm=T)
                   },
                   {               ## default
                                   event = min(x,na.rm=T)
                   }
                   )
            
        }

    }
        
    return(event)

}


occurrenceNumeric <- function(x,COLLAPSE=";",get.year=FALSE, style="%Y-%m-%d"){

    ## return numerical value of each event/year occurrence

    events = ""

    x = splitRecurrence(x,COLLAPSE=COLLAPSE)    

    x = x[x!=""]
    if( length(x) != 0 && is.na(x) == FALSE ){
        events = as.numeric(getDate(x,style=style));
        if( get.year ){
            events = as.numeric(getYear(x,style=style));
        }
    }
    
    return(events)

}


appendDates <- function(oldStr,newStr,COLLAPSE=";",keep.all=TRUE){

    oldStr = as.character(oldStr)
    newStr = as.character(newStr)

    if( newStr == "" || is.na(newStr) ){ return(oldStr); }

    if( length(oldStr) == 0 ){ return(newStr); }
    
    if( sum(grepl(COLLAPSE,oldStr))>0 ){
        oldStr = strsplit(oldStr,COLLAPSE)[[1]]
        oldStr = oldStr[oldStr!=""]
        if( keep.all == FALSE ){ newStr = unique(c(newStr,oldStr)); }
        else {                    newStr = c(newStr,oldStr);         }
        newStr = paste(as.character(newStr),collapse=COLLAPSE)
    } else {
        if( keep.all == FALSE ){ newStr = unique(c(newStr,oldStr)); }
        else {                    newStr = c(newStr,oldStr);         }
        newStr = newStr[newStr!=""]
        newStr = paste(as.character(newStr),collapse=COLLAPSE)
    }
    
    return(newStr)

}


cardiacYearofEvent <- function(data=NULL,pIDs=NULL,date.field="ddisch",
                               style=c("%Y-%m-%d","%d-%b-%Y","%d-%b-%y"),
                               event.level=c("IP_Primary_Cardiac_All",
                                             "IP_Secondary_Cardiac_All"),
                               keep.all=TRUE){

    ## extract the years of a cardiac event, recording the year of first/last
    ## occurrence

    ## initialise the patient subgroup to be empty
    psubgroup = NULL;

    cnames    = c("pid","p.events","p.events.dec","p.events.year",
                        "s.events","s.events.dec","s.events.year");
    psubgroup = matrix(NA,ncol=length(cnames),nrow=length(pIDs))
    colnames(psubgroup) = cnames

    psubgroup[,1] = pIDs    
    
    if( !is.null(data) && !is.null(pIDs) ){

        date.indx  = which(colnames(data)==date.field)
        event.indx = match(event.level,colnames(data))
        event.indx = event.indx[!is.na(event.indx)]               
        
        
        pindx = match(data[,1],pIDs)
        summy = cbind(pindx,data[,c(date.indx,event.indx)])

        pindx = pindx[!is.na(pindx)]
        pindx = unique(pindx)

        t.indx = which(colnames(summy)==date.field)
        p.indx = which(colnames(summy)==event.level[1])
        s.indx = which(colnames(summy)==event.level[2])

        ## loop through all patients with cardiac event data
        for( i in 1:length(pindx) ){
            tmp = summy[which(summy[,1]==pindx[i]),]

            p.events = ""; p.events.d = ""; p.events.y = "";
            s.events = ""; s.events.d = ""; s.events.y = ""; 
            
            Ntmp = dim(tmp)[1]

            if( Ntmp > 0 ){
            
                for( t in 1:Ntmp ){

                    events     = as.character(tmp[t,t.indx[1]]);
                    t.events   = getDate(events,style=style);
                    t.events.d = as.numeric(getDate(events,style=style));
                    t.events.y = getYear(events,style=style);
                                        
                    ## found primary cardiac event(s) for patient
                    if( !is.na(tmp[t,p.indx[1]])  ){
                        p.events   = appendDates(oldStr=p.events,newStr=t.events,
                                                 keep.all=keep.all)
                        p.events.d = appendDates(oldStr=p.events.d,newStr=t.events.d,
                                                 keep.all=keep.all)
                        p.events.y = appendDates(oldStr=p.events.y,newStr=t.events.y,
                                                 keep.all=keep.all)
                    }

                    ## found secondary cardiac event(s) for patient
                    if( !is.na(tmp[t,s.indx[1]]) ){
                        s.events   = appendDates(oldStr=s.events,newStr=t.events,
                                                 keep.all=keep.all)
                        s.events.d = appendDates(oldStr=s.events.d,newStr=t.events.d,
                                                 keep.all=keep.all)
                        s.events.y = appendDates(oldStr=s.events.y,newStr=t.events.y,
                                                 keep.all=keep.all)
                    }
                
                }

            }          
            
            psubgroup[pindx[i],2] = p.events;
            psubgroup[pindx[i],3] = p.events.d;
            psubgroup[pindx[i],4] = p.events.y;
            
            psubgroup[pindx[i],5] = s.events;
            psubgroup[pindx[i],6] = s.events.d;
            psubgroup[pindx[i],7] = s.events.y;

        }
        
    }

    return(psubgroup)
    
}

cardiacSubGroups <- function(data=NULL,pIDs=NULL,binarise=FALSE,death.level=NULL){
    ## use IC10 codes to return patients in cohort 

    ## initialise the patient subgroup to be empty
    psubgroup = NULL;

    cardiac.groups = c("ischaemic","myopathy","other","notcardiac");
        
    group <- matrix(NA,ncol=length(cardiac.groups),nrow=27)
    colnames(group) = cardiac.groups
    group[1,1]  = "I219"; group[1,2] = "I500"; group[1,3]  = "I269"; group[1,4] = "I120";
    group[2,1]  = "I251"; group[2,2] = "I501"; group[2,3]  = "I48X"; group[2,4] = "I10X";   
    group[3,1]  = "I209"; group[3,2] = "I509"; group[3,3]  = "I471"; group[3,4] = "I119";
    group[4,1]  = "I259"; group[4,2] = "I428"; group[4,3]  = "I469"; group[4,4] = "I129"; 
    group[5,1]  = "I211"; group[5,2] = "I110"; group[5,3]  = "I350";
    group[6,1]  = "I210"; group[6,2] = "I420"; group[6,3]  = "I313";
    group[7,1]  = "I214"; group[7,2] = "I517"; group[7,3]  = "I442";
    group[8,1]  = "I200"; group[8,2] = "I255"; group[8,3]  = "I48"; 
    group[9,1]  = "I229";                      group[9,3]  = "I340"; 
    group[10,1] = "I249";                      group[10,3] = "I472"; 
    group[11,1] = "I201";                      group[11,3] = "I495"; 
    group[12,1] = "I208";                      group[12,3] = "I499"; 
    group[13,1] = "I221";                      group[13,3] = "I080";
    group[14,1] = "I228";                      group[14,3] = "I270";
    group[15,1] = "I250";                      group[15,3] = "I460"; 
    group[16,1] = "I258";                      group[16,3] = "I050";
    group[17,1] = "I212";                      group[17,3] = "I059";
    group[18,1] = "I220";                      group[18,3] = "I071";
    group[19,1] = "I248";                      group[19,3] = "I082";
                                               group[20,3] = "I139";
                                               group[21,3] = "I278";
                                               group[22,3] = "I319";
                                               group[23,3] = "I441";
                                               group[24,3] = "I480";
                                               group[25,3] = "I493";
                                               group[26,3] = "I516";
                                               group[27,3] = "I519";


    if( !is.null(data) && !is.null(pIDs) ){

        pindx     = match(pIDs,data[,1])
        psubgroup = matrix(0,ncol=(1+length(cardiac.groups)),nrow=length(pIDs))
        colnames(psubgroup) = c(colnames(data)[1],cardiac.groups)

        psubgroup[,1] = pIDs

        if( is.null(death.level) ){
            death.level = c("pcd",sprintf("scd%d",seq(0,9,1)))
        }
        
        level.indx  = match(death.level,colnames(data))
        level.indx  = level.indx[!is.na(level.indx)]

        if( length(level.indx) != 0 ){
        
            for( p in 1:length(pindx) ){
                if( !is.na(pindx[p]) ){
                    for( g in 1:length(cardiac.groups) ){
                        psubgroup[p,(1+g)] = sum(!is.na(match(data[pindx[p],level.indx],group[,g])))
                    }
                }
            }

            if( binarise ){
                for( p in 2:dim(psubgroup)[2] ){
                    psubgroup[,p] = ifelse(as.numeric(psubgroup[,p])==0,0,1);
                }
            }
        
        }
    }
        
    if( is.null(data) && !is.null(pIDs) ){
        psubgroup = matrix(0,ncol=(1+length(cardiac.groups)),nrow=length(pIDs))
        colnames(psubgroup) = c("",cardiac.groups) 
    }
    
    return(psubgroup)
    
}

getColumnIndx <- function(data=NULL,colNames=NULL,patient.rec=TRUE){

    indx = NULL; tally.tb = NULL;
    
    if( !is.null(data) && !is.null(colNames) ){

        N     = dim(data)[1]
        Nindx = length(colNames)

        if( length(Nindx) != 0 ){
            indx = vector(length=Nindx)
            for( i in 1:Nindx ){
                temp = which(colnames(data)==colNames[i]);
                indx[i] = ifelse(length(temp)!=0,temp,NA) 
            }                        

            if( patient.rec ){
                tally.tb = matrix(0,ncol=Nindx,nrow=N)
                colnames(tally.tb) = colNames
            }
            
        }
        
    }

    return(list(indx=indx,tally.tb=tally.tb,Nindx=Nindx))
    
}

tnm2stage <- function( tnm=NULL ){
    ## Find the correspondance between breast cancer staging and TNM scores can be
    ## found in Table 4 in: https://pubs.rsna.org/doi/full/10.1148/rg.2018180056

    stage = NULL

    if( !is.null(tnm) ){
    
        N     = dim(tnm)[1]
        Ncol  = dim(tnm)[2]

        stage = rep(NA,N)
    
        if( Ncol == 3 ){
            for( i in 1:N ){

                ## mapping works on truncated TNM scores
                TNM= truncateStage( x=as.character(tnm[i,]) )

                if( sum(is.na(TNM)) == 0 ){
        
                    if( TNM[3] == "1" ){
                        stage[i] = "4"; 
                    } else {
                        if( TNM[2] == "3" ){
                            stage[i] = "3C";
                        } else {
                            if(TNM[1] == "4" ){
                                stage[i] = "3B";
                            } else {
                                if( (TNM[1] == "3" && TNM[2] == "2") ||
                                    (TNM[1] == "3" && TNM[2] == "1") ||
                                    (TNM[1] == "2" && TNM[2] == "2") ||
                                    (TNM[1] == "1" && TNM[2] == "2") ||
                                    (TNM[1] == "0" && TNM[2] == "2") ){
                                    stage[i] = "3A";
                                } else {
                                    if( (TNM[1] == "3" && TNM[2] == "0") ||
                                        (TNM[1] == "2" && TNM[2] == "1") ){
                                        stage[i] = "2B";
                                    } else {
                                        if( (TNM[1] == "2" && TNM[2] == "0") ||
                                            (TNM[1] == "1" && TNM[2] == "1") ||
                                            (TNM[1] == "0" && TNM[2] == "1") ){
                                            stage[i] = "2A";
                                        } else {
                                            if( (TNM[1] == "1" && TNM[2] == "0")  ){
                                                stage[i] = "1A";
                                            }
                                        }
                                    }
                                }
                            }
                        }        
                    }
                }
            }
        }
    }  

    return(stage)
  
}

truncateStage <- function(x){

    ## Truncate patient's breast cancer stage measurement
    ## NOTES on breast cancer staging, the correspondance between staging and TNM can
    ## be found in Table 4 of ref: 4   
    ## 1) https://www.google.co.uk/amp/s/amp.cancer.org/cancer/breast-cancer/understanding-a-breast-cancer-diagnosis/stages-of-breast-cancer.html
    ## 2) https://www.cancerresearchuk.org/about-cancer/breast-cancer/stages-types-grades
    ## 3)  https://www.cancer.gov/about-cancer/diagnosis-staging/prognosis/tumor-grade-fact-sheet
    ## 4) https://pubs.rsna.org/doi/full/10.1148/rg.2018180056
    
    ## store truncated cancer stage values 
    N   = length(x)    
    res = rep(NA,N)

    ## truncate cancer stage values, i.e. if patients stage value is: "1",
    ## "1A", "1B" or "1C" truncate to "1"    
    stage <- c("0","1","2","3","4")

    for( i in 1:length(stage) ){
        res[grepl(stage[i],x)]=stage[i]
    }

    return(res)  
    
}


SIMDscore2percentile <- function( x, method=c("quintile","decile"),
                              rev.rank=FALSE ){

    ## Convert a SIMD score to percentile
    ## Ref: https://www.isdscotland.org/Products-and-Services/GPD-Support/Deprivation/SIMD/

    ## Dataset: ../SIMD/postcode_2006_2_simd2004.rds

    ## NOTE: SIMD
    ## Can find the 2004 SIMD scores, ranks, quintiles and decile for 
    ## all data zones from public health scotland:
    ## Ref: https://www.isdscotland.org/Products-and-Services/GPD-Support/Deprivation/SIMD/
    ## dataset: https://www.isdscotland.org/Products-and-Services/GPD-Support/Deprivation/SIMD/_docs/SIMD_2004/postcode_2006_2_simd2004.csv
    
    method =  match.arg(method);
    
    N   = length(x)
    res = rep(NA,N)
    
    SIMDrank = matrix(0,ncol=3,nrow=10)
    SIMDrank[1,1]  = 1; SIMDrank[1,2]  = 0.54;  SIMDrank[1,3]  = 5.33;
    SIMDrank[2,1]  = 2; SIMDrank[2,2]  = 5.33;  SIMDrank[2,3]  = 7.64;
    SIMDrank[3,1]  = 3; SIMDrank[3,2]  = 7.64;  SIMDrank[3,3]  = 10.42;
    SIMDrank[4,1]  = 4; SIMDrank[4,2]  = 10.42; SIMDrank[4,3]  = 13.49;
    SIMDrank[5,1]  = 5; SIMDrank[5,2]  = 13.49; SIMDrank[5,3]  = 17.02;
    SIMDrank[6,1]  = 6; SIMDrank[6,2]  = 17.02; SIMDrank[6,3]  = 21.18;
    SIMDrank[7,1]  = 7; SIMDrank[7,2]  = 21.18; SIMDrank[7,3]  = 26.36;
    SIMDrank[8,1]  = 8; SIMDrank[8,2]  = 26.36; SIMDrank[8,3]  = 33.9;
    SIMDrank[9,1]  = 9; SIMDrank[9,2]  = 33.9;  SIMDrank[9,3]  = 46.15;
    SIMDrank[10,1] = 10; SIMDrank[10,2] = 46.15; SIMDrank[10,3] = 87.57;

    if( method == "decile" ){
        
        if( rev.rank ){ SIMDrank[,1] = rev(SIMDrank[,1]); }
            
            for( i in 1:dim(SIMDrank)[1] ){
                res[x>=SIMDrank[i,2]&x<SIMDrank[i,3]] = SIMDrank[i,1]
            }
    }

    if( method == "quintile" ){

        if( rev.rank ){ SIMDrank[,1] = seq(5,1,-1); }
        
        k=1
        for( i in 1:dim(SIMDrank)[1] ){
            ff = k; ll = k+1;
            if( ll <= dim(SIMDrank)[1] ){
                res[x>=SIMDrank[ff,2]&x<SIMDrank[ll,3]] = SIMDrank[i,1]
                k=k+2
            }
        }
    }
    
   return(res)

}


getDate <- function(x,style=c("%Y-%m-%d","%d-%b-%Y","%d-%b-%y") ){
    if( x == " " || x == "" || is.na(x) ){ x = NA; }
    return(as.Date(x,tryFormats=style))
}

getYear <- function(x,style=c("%Y-%m-%d","%d-%b-%Y","%d-%b-%y") ){
    if( x == " " || x == "" || is.na(x) ){ x = NA; }
    return(format(as.Date(x,tryFormats=style),format="%Y"))
}

death.event <- function(x){ return( ifelse(is.na(x[1]),x[2],x[1]) ) }

first.event <- function(x){ return( ifelse(all(is.na(x)),NA,min(x,na.rm=T)) ) }

last.event  <- function(x){ return( ifelse(all(is.na(x)),NA,max(x,na.rm=T)) ) }

survival.time <- function(start=NULL, end=NULL, events=NULL,
                          outcomes=NULL, outcome.e=NULL, YEAR=c(365.25),
                          method=c("death.min","death.max","first","last") ){

    method = match.arg(method);
    
    s.time = NULL
    
    if( !is.null(start) && !is.null(end) ){

        if( is.null(events) || is.null(outcomes) || is.null(outcome.e) ){
            s.time = round( (as.numeric(end) - as.numeric(start))/YEAR,4);
        } else {

            ## default, find first occurrence of event, 
            ## if method="last" find last occurrence of event
            tmp = cbind(as.vector(end),
                        as.vector(sapply(events,firstOccurrence, method=method)));

            tmp[,1] = ifelse(tmp[,1]=="",NA,tmp[,1])
            tmp[,2] = ifelse(tmp[,2]=="",NA,tmp[,2])

            switch(method,
                   "death.min"={
                       event.tmp = as.vector(apply(tmp,1,death.event));
                   },
                   "death.max"={
                       event.tmp = as.vector(apply(tmp,1,death.event));
                   },
                   "first"={
                       event.tmp = as.vector(apply(tmp,1,first.event));
                   },
                   "last"={
                       event.tmp = as.vector(apply(tmp,1,last.event));
                   },
                   {
                       ## default
                       event.tmp = as.vector(apply(tmp,1,death.event));
                   }
                   )

            event.end = ifelse(outcomes==outcome.e, event.tmp,end)
            s.time    = round((as.numeric(event.end)-as.numeric(start))/YEAR,4)

        }
        
    }

    return(s.time);
    
}


every_nth <- function(x, nth, empty = TRUE, inverse = FALSE){
    ## https://stackoverflow.com/questions/34533472/insert-blanks-into-a-vector-for-e-g-minor-tick-labels-in-r

    if (!inverse) {
    if(empty) {
      x[1:nth == 1] <- ""
      x
      } else {
        x[1:nth != 1]
        }
    } else {
      if(empty) {
        x[1:nth != 1] <- ""
        x
        } else {
          x[1:nth == 1]
        }
    }

}
