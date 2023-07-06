##-------------------------------------------
## load plot functions for analysis cardiac data
##-------------------------------------------
##library(cowplot)
##library(gridExtra)

validation.plots <- function(hist=NULL,
                             metrics="loss",
                             point.size=2,
                             labels.short=c("training"="train","validation"="val"),
                             point.cols=c("darkorange","steelblue"),
                             rm.legend=FALSE){


    gp       = NULL
    cow.plot = NULL
    
    get_legend<-function(myggplot){
        tmp <- ggplot_gtable(ggplot_build(myggplot))
        leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
        legend <- tmp$grobs[[leg]]
        return(legend)
    }
    
    if( !is.null(hist) ){

        gp = list()
        
        hist = as.data.frame(hist)

        names(point.cols) = levels(hist$data)
        event.labs        = labels.short[match(names(labels.short),names(point.cols))]

        N = length(metrics)

        legend.names      = metrics
        title.names       = metrics
        if( rm.legend ){ legend.names = rep("",N); }
        
        for( i in 1:N ){
        
            df=hist[grepl(metrics[i],hist$metric),]
            
            gp[[i]] = ggplot(df,aes(x=epoch,y=value,color=data))+
                geom_point(size=rel(point.size),show.legend=TRUE)+
                geom_line(alpha=0.6,show.legend=FALSE)+
                {if(rm.legend)labs(title=title.names[i])}+
                scale_color_manual(values=point.cols,labels=event.labs,name=legend.names[i])+
                guides(colour = guide_legend(
                           override.aes = list(alpha=1,size=rel(5))))+
                theme_minimal()+
                theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
                      panel.grid.minor = element_line(colour="grey40",size=0.1),
                      panel.background = element_rect(fill = "white"))+
                theme(plot.title = element_text(face="bold",size=rel(1)),
                      axis.text.x = element_text(face="bold", size=rel(1.0), angle=40,),
                      axis.text.y = element_text(face="bold", size=rel(1.0)),
                      axis.title.x=element_text(face="bold",size=rel(1.0)),
                      axis.title.y=element_text(face="bold",size=rel(1.0)),
                      legend.title=element_text(face="bold",size=rel(1.0)),
                      legend.text=element_text(face="bold",size=rel(1.0)),
                      legend.position="top")

            names(gp)[i] = metrics[i]
        }

        if( rm.legend ){
        
            ## save legend
            legend <- get_legend(gp[[1]])

            ## remove legend from plots
            gp2 = list()
            for( i in 1:N ){ gp2[[i]] = gp[[i]] + theme(legend.position="none"); }

        } else {

            legend <- ggplot()+geom_blank()+cowplot::theme_nothing();

            gp2 = list()
            for( i in 1:N ){ gp2[[i]] = gp[[i]]; }
            
        }
            
        ## create blank plot
        blankPlot <- ggplot() + geom_blank(aes(1,1))+cowplot::theme_nothing()

        if( N == 1 ){ gp2[[2]] = ggplot()+geom_blank()+cowplot::theme_nothing(); }
        
        cow.plot = grid.arrange(legend, blankPlot,  gp2[[1]], gp2[[2]],
                                ncol=2, nrow = 2, 
                                widths = c(2.7, 2.7), heights = c(0.2, 2.5))
        
        ##cow.plot = plot_grid(plotlist=gp,labels=LETTERS[seq(1,N)],
        ##                     align="hv", axis="tblr")

        ##cow.plot = cowplot::plot_grid(
        ##                        plotlist=gp,
        ##                        labels=LETTERS[seq(1,N)],
        ##                        nrow = 2, 
        ##                        rel_heights = c(4, 1), 
        ##                        align = "v", 
        ##                        axis = "b"
        ##                    )
    }    

    return(list(gp=gp,cow.plot=cow.plot))
}

competing.risk.evalCI.plots <- function( ci=NULL, time.points=NULL,
                                        point.cols=NULL,
                                        event.labs=NULL,
                                        line.col="darkred",
                                        alpha.line=c(0.6),
                                        alpha.band=c(0.2),
                                        title=""){

    ## Plots to evaluate the survival probabilities made by the Deep learning algorithm
    ## ci.plot ==> concordance index plot
  
    ci.plot = NULL

    if( !is.null(ci) && !is.null(time.points) ){

        ##ci = as.data.frame(ci)

        qt = time.points
        ne = dim(ci)[3]
        nt = dim(ci)[1]

        if( is.null(point.cols) ){
            if( ne < 3 ){
                point.cols=c("darkorange","steelblue","black")[1:ne] 
            } else {                
                library(RColorBrewer)
                point.cols = brewer.pal(n=ne,name="Dark2")
            }
        }

        if( is.null(event.labs) ){
            event.labs = sprintf("event.type_%g",1:ne)
        }

        names(point.cols) = 1:ne
        names(event.labs) = 1:ne
        shapes = rep(15,ne)

        df = data.frame()
        
        for( e in 1:ne ){
            tmp = as.data.frame(ci[,,e])
            tmp = cbind(tmp,qt,rep(e,nt))
            df = rbind(df,tmp)
        }

        colnames(df) = c(dimnames(ci)[[2]],c("qt","event"))

        df       = df[!is.na(df$ci),]
        df$qt    = factor(df$qt,levels=unique(df$qt))
        df$event = as.factor(df$event)
        df       = droplevels(df)
        
        ci.plot = ggplot(df,aes(x=ci,y=qt,xmin=lower,xmax=upper,group=event))+
            geom_point(aes(col=event),size=rel(1.5),show.legend=TRUE)+
            geom_lineribbon(aes(col=event,fill=event),alpha=alpha.band,
                            show.legend=FALSE)+
            geom_path(aes(col=event),size=rel(1),alpha=alpha.line,
                      show.legend=FALSE)+
            scale_y_discrete(limits=rev(as.factor(qt)))+
            scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.1))+
            scale_color_manual(values=point.cols,labels=event.labs,name=NULL)+
            scale_fill_manual(values=point.cols,labels=event.labs,name=NULL)+
            guides(colour = guide_legend(
                       override.aes = list(shape = shapes, size=rel(7))))+               
            geom_vline(xintercept=0.5,colour=line.col,size=rel(1),show.legend=F)+
            theme_minimal()+
            labs(y="Years",x="Concordance Index",title=title)+
            theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
                  panel.grid.minor = element_line(colour="grey40",size=0.1),
                  panel.background = element_rect(fill = "white"))+
            theme(plot.title = element_text(face="bold",size=rel(1)),
                  axis.text.x = element_text(face="bold", size=rel(1.5)),
                  axis.text.y = element_text(face="bold", size=rel(1.5)),
                  axis.title.x=element_text(face="bold",size=rel(2.0)),
                  axis.title.y=element_text(face="bold",size=rel(2.0)),
                  legend.title=element_text(face="bold",size=rel(1.7)),
                  legend.text=element_text(face="bold",size=rel(1.7)),
                  legend.position="top")
               
    }

    return(list(pred.ci=df,ci.plot=ci.plot))
    
}


survival.plot <- function( fit=NULL, data=NULL, line.type=NULL, line.method="strata", colours=NULL,
                          legend.labs=NULL, legend.xy=NULL, fun=NULL, conf.int=FALSE, conf.style="ribbon",
                          censor=FALSE, pval=FALSE, rt=FALSE, rt.y.text.col=FALSE, rt.height=FALSE,
                          rt.y.text=FALSE, legend.title="", ylim=NULL, xlim=NULL, xlab="years" ){

    ## survival (ggplots) plots making use of the library survminer:
    sur.plot = NULL
    gplot    = NULL
    cow.plot = NULL
    
    ## legend.xy = c(0.2,0.2)  # bottom left
    ## legend.xy = c(0.8,0.87) # upper right
    ## legend.xy = c(0.2,0.87) # upper left
    
    if( !is.null(fit) && !is.null(data) &&
        !is.null(line.type) && !is.null(colours) &&
        !is.null(legend.labs) ){
        
        ##ggsurvplot_combine(fit,data)
        sur.plot <- ggsurvplot(fit=fit,data=data,combine=TRUE,
                               conf.int=conf.int,
                               conf.int.style=conf.style,
                               pval=pval,
                               censor=censor,
                               tables.theme=theme_cleantable(),
                               ggtheme = theme_bw(),
                               legend = legend.xy,
                               legend.title=legend.title,
                               legend.labs=legend.labs,
                               linetype=line.method,
                               xlab=xlab,
                               fun=fun,
                               risk.table=rt,
                               risk.table.y.text.col=rt.y.text.col,
                               risk.table.y.text=rt.y.text,
                               risk.table.height=rt.height,
                               xlim=xlim,
                               ylim=ylim#,
                               #palette=colours
                               )        
        
        gplot = sur.plot$plot +
            scale_linetype_manual(values=line.type)+
            scale_colour_manual(values=colours)+
            {if(conf.int)scale_fill_manual(values=colours)}+
            guides(colour = guide_legend(
                       override.aes = list(
                           colour=colours,
                           linetype=line.type,
                           fill="none",
                           alpha=c(1))               
                   ),
                   fill     = "none",
                   shape    = "none",
                   linetype = "none"
                   )

        if( rt ){
            ## join risk table to survival curve, with customised colours, using cowplot 
            
            risk.tab = sur.plot$table
            risk.tab$theme$axis.text.y$colour     = colours
            risk.tab$theme$axis.text.y$linetype   = line.type
            risk.tab$theme$axis.text.y$box.colour = "white"
            risk.tab$theme$axis.text.y$linewidth  = c(2)
            
            cow.plot = cowplot::plot_grid(
                                    gplot, 
                                    risk.tab + theme_cleantable(), 
                                    nrow = 2, 
                                    rel_heights = c(4, 1), 
                                    align = "v", 
                                    axis = "b"
                                )

        }
        
    }

    return(list(sur.plot=sur.plot,gplot=gplot,cow.plot=cow.plot))
    
}


plot.Fscore <- function(data, time.points, events, tests, colours, shapes,
                        shape=c(16,17,15),xlab="qt (years)",ylab="value",
                        text.y.size=2, text.x.size=2, legend.text.size=1.7,
                        x.tick.size=1, y.tick.size=1,
                        x.tick.face="plain", y.tick.face="plain", 
                        legend.box.size=5,
                        line.size=1, point.size=2,
                        width=480, height=480, x.scale=seq(0,1,0.1),
                        y.scale=seq(0,1,0.1),
                        line.x.pos = c(0.6),
                        zero.na=FALSE){
    gp1=NULL
    gp2=NULL
    
    if( !is.null(data) ){

        ## set NA values to zero 
        if( zero.na ){ data[is.na(data)]=0; }
        
        ## format data into data.frame for plotting
        df = data.frame()

        qt = time.points
        nt = length(qt)

        t.indx = match(names(tests),dimnames(data)[[2]])
        t.indx = t.indx[!is.na(t.indx)]
        
        for( s in 1:nt ){
            t.point=data[,t.indx,s]
            t.nr = dim(t.point)[1]
            t.nc = dim(t.point)[2]
            tmp = c()
            for( i in 1:t.nc ){
                if( i == 1 ){
                    tmp = cbind(rownames(t.point),rep(colnames(t.point)[i],t.nr),t.point[,i])
                } else {
                    tmp = rbind(tmp,cbind(rownames(t.point),rep(colnames(t.point)[i],t.nr),t.point[,i]))
                }
            }
            df = rbind(df,cbind(rep(qt[s],dim(tmp)[1]),tmp))               
        }

        rownames(df) = NULL
        colnames(df) = c("qt","event","score","value")        

        x.min = min(as.numeric(qt),na.rm=T)
        x.max = max(as.numeric(qt),na.rm=T)
        x.scale = as.numeric(qt)

        ## New facet label names for supp variable
        ##supp.labs        <- levels(df$class)
        ##new.labs         <- rep("",length(supp.labs))
        ##if( !is.null(names(class.labels)) ){
        ##    indx = match(names(class.labels),levels(df$class))
        ##    new.labs[indx] = as.vector(class.labels)
        ##} else {
        ##    new.labels = as.vector(class.labels)
        ##}   

        for( e in 1:length(events) ){
        
            df.tmp = df[df$event==events[e],]
            leg.tit = sprintf("%s Scores",names(events)[events==events[e]])
        
        gp1[[e]] <- ggplot(df.tmp,aes(x=as.numeric(as.vector(qt)),y=as.numeric(as.vector(value)),group=score,linetype=score,shape=score))+
            geom_line(aes(colour=score),size=line.size,alpha=0.6,show.legend=F)+
            geom_point(aes(colour=score),size=point.size)+
            #facet_wrap(~class, ncol=1)+
            labs(x=xlab,y=ylab)+
            scale_x_continuous(limits=c(x.min, x.max), breaks=x.scale)+
            scale_y_continuous(limits=c(0, 1), breaks=y.scale)+
            ##{if(plot.line)geom_vline(xintercept=line.x.pos,colour=line.colour,linetype="longdash", alpha=0.6, size=rel(1),show.legend=F);}+
            theme(            
                panel.grid.major = element_line(colour = "grey40",size=0.2),
                ## panel.grid.minor = element_line(colour="grey40",size=0.1),
                panel.background = element_rect(fill = "white"),
                ##legend.position="none",
                ##strip.text.x = element_text(size=rel(class.label.size), color="black",face="bold"),
                axis.text.x =
                    element_text(angle = 40, hjust = 1,
                                 face=x.tick.face, size=rel(x.tick.size)),
                axis.text.y = element_text(size=rel(y.tick.size),face=y.tick.face),
                axis.title.x=element_text(face="bold",size=rel(text.x.size)),
                axis.title.y=element_text(face="bold",size=rel(text.y.size)),
                legend.title=element_text(face="bold",size=rel(legend.text.size)),
                legend.text=element_text(face="bold",size=rel(legend.text.size)) )+
            guides(colour = guide_legend(override.aes = list(shape = shapes, size=rel(legend.box.size)),title=leg.tit),
                   size   = "none",##FALSE,
                   shape  = "none")+##FALSE)+
            scale_colour_manual(breaks=names(tests),values=tests,aesthetics = c("colour"))

            names(gp1)[e] = names(events)[e]
            
        }
            ##scale_colour_manual(breaks=levels(factor(df2$score)),values=tests,aesthetics = c("colour"))

        ##png(sprintf("%s/%s/%s.png",DIR,subDIR,shortTIT),width=width,height=height,units="px");
 ##       print(gplot)
   ##     dev.off()

        for( s in 1:length(tests) ){

            test=names(tests)[s]
            
            df.tmp = df[df$score==test,]
            fix.shape=shapes[names(shapes)==test]
            leg.tit=test
        gp2[[s]] <- ggplot(df.tmp,aes(x=as.numeric(as.vector(qt)),y=as.numeric(as.vector(value)),group=event,linetype=event),shape=fix.shape)+
            geom_line(aes(colour=event),size=line.size,alpha=0.6,show.legend=F)+
            geom_point(aes(colour=event),size=point.size)+
            #facet_wrap(~class, ncol=1)+
            labs(x=xlab,y=ylab)+
            scale_x_continuous(limits=c(x.min, x.max), breaks=x.scale)+
            scale_y_continuous(limits=c(0, 1), breaks=y.scale)+
            ##{if(plot.line)geom_vline(xintercept=line.x.pos,colour=line.colour,linetype="longdash", alpha=0.6, size=rel(1),show.legend=F);}+
            theme(            
                panel.grid.major = element_line(colour = "grey40",size=0.2),
                ## panel.grid.minor = element_line(colour="grey40",size=0.1),
                panel.background = element_rect(fill = "white"),
                ##legend.position="none",
                ##strip.text.x = element_text(size=rel(class.label.size), color="black",face="bold"),
                axis.text.x =
                    element_text(angle = 40, hjust = 1,
                                 face=x.tick.face, size=rel(x.tick.size)),
                axis.text.y = element_text(size=rel(y.tick.size),face=y.tick.face),
                axis.title.x=element_text(face="bold",size=rel(text.x.size)),
                axis.title.y=element_text(face="bold",size=rel(text.y.size)),
                legend.title=element_text(face="bold",size=rel(legend.text.size)),
                legend.text=element_text(face="bold",size=rel(legend.text.size)) )+
            guides(colour = guide_legend(override.aes = list(size=rel(legend.box.size)),title=leg.tit),
                   size   = "none",##FALSE,
                   shape  = "none")+##FALSE)+
            scale_colour_manual(breaks=names(colours),values=colours,
                                labels=names(events)[match(names(colours),events)],aesthetics = c("colour"))
           
        names(gp2)[s] = test

        }

    }

    
    return(list(gp1=gp1,gp2=gp2))

}
