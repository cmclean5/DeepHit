##--------------------------------------------------------------------------------------
## Suppose we want to combine each (e.g. 3) competing risk model together and predict
## each (e.g. 3) competing risk probability for each patient. 
## First, we need to use Keras functional API, instead of a Sequential model.
## Ref:
## [1] https://tomroth.com.au/keras/
## [2] https://theailearner.com/2019/01/25/multi-input-and-multi-output-models-in-keras/
## [3] https://datascience.stackexchange.com/questions/30537/keras-loss-function-for-multidimensional-regression-problem?rq=1
## [4] https://machinelearningmastery.com/keras-functional-api-deep-learning/
## [5] https://stackoverflow.com/questions/66845924/multi-input-multi-output-model-with-keras-functional-api
##-------------------------------------------------------------------------------------------

##-------------------------------------------------------------------------------------------
## You should tune the hyperparamters to find the best neural network for your own study.
##-------------------------------------------------------------------------------------------

##-------------------------------------------------------------------------------------------
## load custom loss functions
##-------------------------------------------------------------------------------------------
source("customLoss.R")


##-------------------------------------------------------------------------------------------
## utility functions for setting up model architecture
##-------------------------------------------------------------------------------------------

## set which initialiser to use
setInitializer <- function(params, seed=0){

    Initializer=NULL
    
    switch(params$initializer.name,
           "GlorotNorm"={
               Initializer = initializer_glorot_normal(seed=seed)
           },
            {
               print("Error: no initializer.name given.")
            }
           )

           return(Initializer)

}

## set which optimiser to use
setOptimizer <- function(params){

    Optimizer=NULL
    
    switch(params$optimizer.name,
           "AdamDyn"={
               Optimizer = optimizer_adam(lr_schedule(params))
           },
           "Adam"={
               Optimizer = optimizer_adam(learning_rate = params$learn.rate)
           },
           "sgd"={
               Optimizer = optimizer_sgd(learning_rate = params$learn.rate)
           },
           "rmsprop"={
               Optimizer = optimizer_rmsprop(learning_rate = params$learn.rate)
           },
           {
               print("Error: no optimizer.name given, 
                             use: \"Adam\",\"sgd\",\"rmsprop\".")
           }
           )      
    
    return(Optimizer)
    
}

## set which loss function to use.
## custom loss function can be found in script customLoss.R,
## which is loaded at the start.
setLoss <- function(params, alpha=NULL, num_classes=1){

    Loss=NULL
    switch(params$loss.name,
           "categorical_crossentropy"={
               Loss=loss_categorical_crossentropy()
           },           
           "mae"={
               Loss=loss_mean_absolute_error()
           },
           "mse"={
               Loss=loss_mean_squared_error()
           },
           "sigmoid_focal_crossentropy"={
               ##if( is.null(alpha) ){ alpha=0.25; }
               ##Loss=loss_sigmoid_focal_crossentropy(
               ##    alpha=alpha,
               ##    reduction=tf$keras$losses$Reduction$SUM_OVER_BATCH_SIZE)
           },
            "CCEwe"={
               if( is.null(alpha) ){ alpha = rep(1,num_classes); }
               Loss=CCEweLoss_wrapper(alpha=alpha)
           },          
            "focalLoss"={
               if( is.null(alpha) ){ alpha=rep(0.25,num_classes); }##alpha=0.25; } 
               Loss=focalLoss_wrapper(gamma=2.0, alpha=alpha)
            },
           "deephit_mse"={              
               Loss=deephitLoss_mse_wrapper(alpha = params$alpha, beta = params$beta, gamma = params$gamma,
                                            sigma1 = params$sigma1)               
           }, 
           "deephit"={              
               Loss=deephitLoss_wrapper(alpha = params$alpha, beta = params$beta, gamma = params$gamma,
                                        sigma1 = params$sigma1)               
           },          
           {
               print("Error: no loss.name given, 
                             use: \"...\".")
           }
           )    
    return(Loss)
}


## set which metrics to use
setMetric <- function(metric.name, num_classes=1){

    multi_label = as.integer(0)
    if( num_classes > 1 ){ multi_label = as.integer(1); }
    num_classes = as.integer(num_classes);
    
    Metric=NULL
    switch(metric.name,
           "accuracy"={
               Metric=metric_accuracy()
           },           
           "auc"={
               Metric=metric_auc(multi_label=multi_label,
                                 num_labels=num_classes)
           },
           "categorical_accuracy"={
               Metric=metric_categorical_accuracy()
           },
           "categorical_crossentropy"={
               Metric=metric_categorical_crossentropy()
           },
           "f1score"={
               Metric=metric_fbetascore(num_classes=num_classes)
           },
           "mae"={
               Metric=metric_mean_absolute_error()
           },
           "mse"={
               Metric=metric_mean_squared_error()
           },
           {
               print("Warning: no metric.name given.")
           }
           )    
    return(Metric)
}


## set which regulariser to use
setRegularizer <- function(params,name,value){

    Regularizer=NULL

    if( params$use_regularizer ){
        ##if dropout layer is not true, use regularizer  
        
    switch(name,
           "l1"={
               Regularizer = regularizer_l1(value)
               },
           "l2"={
               Regularizer = regularizer_l2(value)
           },
           {
               print("Error: no optimizer.name given, 
                             use: \"l2\".")   
           })

    }
        
    return(Regularizer)
}

setMax.Weight.Const <- function(cond,value){

    Mx.We.Con=NULL

    if(cond){
        ##if dropout layer is true, use use max.weight.constraint        
        Mx.We.Con = value
    }

    return(Mx.We.Con)
}

## Many models train better if you gradually reduce the learning rate during training:
## Ref: https://tensorflow.rstudio.com/tutorials/keras/overfit_and_underfit
## For the anthracycline dataset, i don't see an imporvement.
lr_schedule <- function(params){
    learning_rate_schedule_inverse_time_decay(
        0.01,
        decay_steps = params$steps.per.epoch*4,
        decay_rate = 1,
        staircase = FALSE
    )
}

##----------------------------------------------------------------------------
## Use callback_tensorboard() to generate TensorBoard logs for the training.
## Ref: https://tensorflow.rstudio.com/tutorials/keras/overfit_and_underfit
##----------------------------------------------------------------------------
get_callbacks <- function(model.name, dir.name) {
    list(
        ##callback_early_stopping(monitor = 'categorical_crossentropy', patience = 200),
        callback_early_stopping(monitor = 'val_loss', patience = 50),
        callback_tensorboard(file.path(dir.name, model.name))
  )
}


##-------------------------------------------------------------------------------------------
## END LOAD UTILITY FUNCTIONS
##-------------------------------------------------------------------------------------------


##-------------------------------------------------------------------------------------------
## MODEL ARCHITECTURES
##-------------------------------------------------------------------------------------------

##-------------------------
## DeepSurv: covariate Feed-forwards network model
##-------------------------
##1) https://livebook.manning.com/book/deep-learning-with-r/chapter-7/44
##2) https://keras.rstudio.com/articles/guide_keras.html
##3) https://tensorflow.rstudio.com/guide/keras/functional_api/
multiState.func.model <- function(x_train, y_train, x_val=NULL, y_val=NULL, params,
                                  alpha=NULL, logdir=NULL, img.dir="images",
                                  model.dir="Models"){


    model.name  = params$model.name
    ## number of event types (including censored event tpes)
    num_classes = dim(y_train)[[2]]
    Loss        = setLoss(params=params, alpha=alpha, num_classes=num_classes)
    model       = NULL
    history     = NULL
    
    if( params$load.model ){
        model = load_model_tf(sprintf("%s/%s.tf",model.dir,model.name),
                              custom_objects = {'loss': Loss}, compile = TRUE) 
        } else {
       
    type      = list()
    pred_type = list()
    data      = list()
    ## number of non-censored event types
    ne        = dim(x_train)[[3]]

    ##use.callback
    callback = NULL
    if( params$use.callback ){ callback = get_callbacks(model.name=model.name,
                                                        dir.name=logdir); }
    
    Regularizer.1 = setRegularizer(params,
                                   params$regularizer1.name,
                                   params$reg.1)

    Regularizer.2 = setRegularizer(params,
                                   params$regularizer2.name,
                                   params$reg.2)

    Mx.We.Con.1   = setMax.Weight.Const(params$dropout.1, params$mx.we.con)
    Mx.We.Con.2   = setMax.Weight.Const(params$dropout.2, params$mx.we.con)
    
    Optimizer = setOptimizer(params)
    ##Loss      = setLoss(params=params,   alpha=alpha, num_classes=num_classes)

    Metric      = list()
    Metric[[1]] = setMetric(metric.name=params$metric.name,  num_classes=num_classes)
    Metric[[2]] = setMetric(metric.name=params$metric.name2, num_classes=num_classes)
    
    ## is a validation dataset provided
    val.data = NULL
    
    for( e in 1:ne ){

        input.name=sprintf("type%g",e)
        
        type[[e]] = layer_input(shape=c(dim(x_train)[[2]]),name=input.name)
    
        pred_type[[e]] = type[[e]] %>% layer_dense(units=params$neurons.1,
                                         kernel_regularizer=Regularizer.1,
                                         kernel_constraint=constraint_maxnorm(max_value=Mx.We.Con.1),
                                         activation=params$act.fun.1) %>%
            {if(params$dropout.1) layer_dropout(.,rate=params$drop.rate.1) else .} %>%
                                       layer_dense(units=params$neurons.2,
                                         kernel_regularizer=Regularizer.2,
                                         kernel_constraint=constraint_maxnorm(max_value=Mx.We.Con.2),
                                         activation=params$act.fun.2) %>%
            {if(params$dropout.2) layer_dropout(.,rate=params$drop.rate.2) else .} %>%
                                       layer_dense(units=dim(y_train)[[2]],
                                         activation=params$act.fun.out)

        data[[e]]      = x_train[,,e]
        names(data)[e] = input.name

        ## construct the validation dataset
        if( params$use.valData ){ val.data[[e]] = x_val[,,e]; names(val.data)[e] = input.name; }
        
    }

    ## set validation dataset to use in fit
    val.dataset = NULL;     
    if( params$use.valData ){ val.dataset = list(val.data,y_val); }
    
    ## concatenate the predicted event types
    ## softmax layer to convert each subjects state value into a probability
    ## which activation function to use?
    ## Ref: https://stats.stackexchange.com/questions/218542/which-activation-function-for-output-layer
    concatenated <- layer_concatenate(pred_type)
    ms.probs <- concatenated %>% 
        layer_dense(units=dim(y_train)[[2]], activation = "softmax")
    
    ## Instantiate the model given inputs and outputs.
    model <- keras_model(inputs=type,
                         outputs=ms.probs)
    
    ## The compile step, i.e. specify the training configuration.
    model %>% compile(
                  optimizer=Optimizer,
                  loss=Loss,
                  metrics=Metric)

    
    ## Train model & evaluate
    history <- model %>% fit(
                             x=data,
                             y=y_train,
                             batch_size=params$batch,
                             epochs=params$epochs,
                             ##shuffle=params$shuffle,
                             callbacks = callback,
                             verbose=params$verbose,
                             validation_data=val.dataset,
                             validation_split=params$fit.valSplit)

    ## plot model
    if( params$plot.model ){
        keras$utils$plot_model(
                        model, sprintf("%s/%s.png", img.dir, model.name),
                        show_shapes = TRUE
                    )
    }

    ## save model
    if( params$save.model ){
        save_model_tf(model,sprintf("%s/%s.tf",model.dir, model.name))
    }

  }
    
    return(list(model=model,history=history))

}

#----------------------------------------------------------------------------
# multistate prediction based on keras functional multistate model
#----------------------------------------------------------------------------
multiState.predict <- function(model, x_test){
    data = list()

    ## number of event types
    ne = dim(x_test)[[3]]
    for( e in 1:ne ){
        data[[e]] = x_test[,,e]
        names(data)[e] = model$input_names[e]
    }

    ypred <- model %>% predict(data)

    ypred
}

