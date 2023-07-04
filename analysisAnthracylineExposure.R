## References:
## Code : https://github.com/chl8856/DeepHit/blob/master/
## paper: http://medianetlab.ee.ucla.edu/papers/AAAI_2018_DeepHit

rm(list=ls())

## load all R libraries
source('Headers.R')

## check if we're working on mac M1 machine
#if( Sys.info()[["sysname"]] == "Darwin" && Sys.info()[["machine"]] == "arm64" ){
#    library(tensorflow)
#    tf_config()
#} else {
#    install_keras()                   
#}

#hp <- import("tensorboard.plugins.hparams.api")

## load utility functions
source('utilityFunctions.R')

## load plot functions
source('plotFunctions.R')

## functions to format datasets
source('formatFunctions.R')

## tensorflow deep learning model build using Keras API
source('model_architecture.R')
##source('Keras_model.R')

## function used in analysis
source('datasetFunctions.R')

## evaluation functions used in analysis
source('evaluateFunctions.R')

##  Examples define optimiser and learning rates
#1) https://github.com/tensorflow/tensorboard/issues/3688

#---Check or create summary stats dir
summdir = "stats"
if( !file_test("-d",summdir) ){
    dir.create(summdir)
}

#---Check or create plot dir
plotdir = "plots"
if( !file_test("-d",plotdir) ){
    dir.create(plotdir)
}

#---Check or create image dir
imgdir = "images"
if( !file_test("-d",imgdir) ){
    dir.create(imgdir)
}

#---Check or create keras model dir
modeldir = "Models"
if( !file_test("-d",modeldir) ){
    dir.create(modeldir)
}

## parameters
YEAR=365.25;

##WIDTH and HEIGHT for plots
WIDTH=480
HEIGHT=480

## Load hyperparameters for deephit model
source('Hyperparameters/Hyperparameter.R')

## load preprocessed data
dataFile <- vector(length=1)
dataFile[1] = "datasets/preprocess_anthracyline_dataset.csv"
data = read.delim(dataFile[1],sep="\t",header=T )

# Number of patients in cohort
N = dim(data)[1]

##------------------------- OBS. SURV. TIME --------------------------------
## Get observed survial time for each patients
##--------------------------------------------------------------------------

## set censore date
c.date = "2016-12-01";

## calculate survival time using first occurrence of cardiac event/death
survT = survivalTime(start.time=data$startDate,
                     end.time=data$endDate,
                     c.time=c.date,
                     event.time=data$cardiac.p.date,
                     outcomes=data$eventPrim,
                     outcome.e="2", YEAR=YEAR,
                     method="first");

## time points in years to train the model
qt = seq(0.5, 15, 0.5)

## The time points we want to evaluate predictions at
eval.times = c(0.6,1,5,10,15)

## Create a matrix consisting of the predictors.
covars. = cbind(data$chemo,data$anthra,data$notanthra,data$age,
               factor2ind(data$stage, "1"),
               factor2ind(data$simd, "1"),
               data$left_rt, data$er_pos,
               data$charlson.notcardiac,
               data$charlson.cardiac,
               factor2ind(data$grade,"1"),
               data$her2_pos,
               data$screen,
               factor2ind(data$startYear,"2000")
               )

colnames(covars.) <- c("chemo","anthra", "notanthra", "age",
                      "stage2", "stage3", "simdq2", "simdq3",
                      "simdq4", "simdq5", "left_rt",
                      "er_pos", "charlson.notMI","charlson",
                      "grade2","grade3","her2_pos","screen",
                      sprintf("200%d",seq(1,9,1)),"2010"
                      )

## get the nnumber of competing events
n.events = get.comp.events(events=data$eventPrim, cencode=0)
ne       = length(n.events)

## number of states (including censored)
n.states = c(0,n.events)
ns       = length(n.states)    
state.ids= seq(0,(ns-1),1)


##------------------------------------------------------------------
## sigmoid_focal_crossentropy
## Give high weights to the rare class and small weights
## to the dominating or common class. These weights are referred to as α.
##------------------------------------------------------------------
Y        = data$eventPrim
Yclasses = levels(factor(Y))
cr.frac  = prop.table(table(Y))
alpha    = 1-cr.frac ##prop.table(table(Y))

##------------------------- run deep learning model to calculate survival prob --------------------------------

## set random number
Seed=1234
set.seed(Seed)

# #number of random.runs
rnd.runs=5

## fraction of data set for training 
data_ids   = seq(1,N,1)
rnd_per    = floor(PARAMS.DH$ratio_train*N)
train_ids  = sample(data_ids,rnd_per,replace=F)
test_ids   = data_ids[!data_ids %in% train_ids]

if(PARAMS.DH$use.valData){
    Ntrain    = length(train_ids)
    rnd_per   = floor(PARAMS.DH$ratio_val*Ntrain)
    val_ids   = sample(train_ids,rnd_per,replace=F)
    train_ids = train_ids[!train_ids %in% val_ids] 
}

## Print some summary statistics
cat("> :== Summary ==:\n")
cat("> dataset: No patients = ", dim(covars.)[1], ", No cofounders: ", dim(covars.)[2],"\n")
cat("> evaluation time.points (years): ")
for( i in 1:length(eval.times) ){ cat(eval.times[i], ","); }
cat("\n")
cat("> competing risks: ")
for( i in 1:length(cr.frac) ){ cat( names(cr.frac)[i],"=",round(cr.frac[i],3)," "); }
cat("\n")

## Print dataset sizes
if(PARAMS.DH$use.valData){
    cat("> dataset split sizes:\n")
    cat("> data: ", length(data_ids),     "(", round(length(data_ids)/N,3)*100,
        "%), train: ", length(train_ids), "(", round(length(train_ids)/N,3)*100,
        "%), val: ", length(val_ids),     "(", round(length(val_ids)/N,3)*100,
        "%), test: ", length(test_ids),   "(", round(length(test_ids)/N,3)*100,
        "%).\n")
} else {
    cat("> data: ", length(data_ids),     "(", round(length(data_ids)/N,3)*100,
        "%), train: ", length(train_ids), "(", round(length(train_ids)/N,3)*100,
        "%), test: ", length(test_ids),   "(", round(length(test_ids)/N,3)*100,
        "%).\n")
}
cat("> :=============:\n")


## generate training data 
xtrain =  generate.ms.data(
    rnd_indx=train_ids,
    covars=covars.,
    eventset=data$eventPrim,
    events=n.events,
    survt=survT,
    qt=qt,
    dataset.name="training")    

xval = NULL
if(PARAMS.DH$use.valData){
    ## generate validation data 
    xval =  generate.ms.data(
        rnd_indx=val_ids,
        covars=covars.,
        eventset=data$eventPrim,
        events=n.events,
        survt=survT,
        qt=qt,
        dataset.name="validation")
}
    
## generate test data
xtest = generate.ms.data(
    rnd_indx=test_ids,
    covars=covars.,
    eventset=data$eventPrim,
    events=n.events,
    survt=survT,
    qt=qt,
    dataset.name="test")


##------------------------------------------------------------------
## Train deephit model
##------------------------------------------------------------------
model = NULL
model = train.deephit.model(train=xtrain,
                            val=xval,
                            events=n.events,
                            qt=qt,
                            params=PARAMS.DH
                            )

##------------------------------------------------------------------
## Test deephit model
##------------------------------------------------------------------
y_pred.dh = test.deephit.model(model = model$model,
                               test  = xtest
                               )
   
##------------------------------------------------------------------
## Evalutate deephit model
##------------------------------------------------------------------
eval.dh = eval.deephit.model(model=model$model,
                             eval=xtest,
                             train=xtrain,
                             qt=qt,
                             eval.times=eval.times,
                             n.events=n.events,
                             y_pred=y_pred.dh)


##------------------------------------------------------------------
## Plot of ci evaluation metric
##------------------------------------------------------------------
ci.plot = competing.risk.evalCI.plots(ci=eval.dh$ciCI,
                                      time.points=eval.times,##qt,
                                      point.cols=c("darkorange","steelblue","black"),
                                      event.labs=c("BrC","Cardiac","Other"),
                                      line.col="darkred",
                                      alpha.line=c(0.6),
                                      alpha.band=c(0.2),
                                      title="")

png(filename=sprintf("%s/deephit_ci_plot.png",plotdir),
    width=WIDTH, height=HEIGHT, units="px")
plot(ci.plot$ci.plot)
dev.off()

