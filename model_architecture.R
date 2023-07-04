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
## load custom  Keras models, which may depend in custom loss functions above 
##-------------------------------------------------------------------------------------------
source("Keras_model.R")


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
## DeepHit: covariate Feed-forwards network model
##-------------------------
## x_dim     : number columns in data (e.g. 12)
## input     : input (tensor) (e.g. data)
## num_layers: number of layers in FCNet, i.e. hidden layers  
## h_dim     : number of hidden units, i.e. nodes.
## h_fn      : activation function for hidden layers (default: tf.nn.relu)
## o_dim     : number of output units, i.e. nodes.
## o_fn      : activation function for output layers (defalut: None)
## w_init    : initialization for weight matrix (defalut: Xavier)
##           : Notes:
##           : Xavier and Glorot are 2 names for the same initializer algorithm,
##           : in TF version2 use,
##           : tf$initializers$glorot_normal()
## keep_prob : keep probabilty [0, 1]  (if None, dropout is not employed)
##-------------------------
## Hyperparameters (DeepHit/main_RandomSearch.py)
##-------------------------
##    SET_BATCH_SIZE    = [32, 64, 128] #mb_size
## 
##    SET_LAYERS        = [1,2,3,5] #number of layers
##    SET_NODES         = [50, 100, 200, 300] #number of nodes
##
##    SET_ACTIVATION_FN = ['relu', 'elu', 'tanh'] #non-linear activation functions
##
##    SET_ALPHA         = [0.1, 0.5, 1.0, 3.0, 5.0] #alpha values -> log-likelihood loss 
##    SET_BETA          = [0.1, 0.5, 1.0, 3.0, 5.0] #beta values -> ranking loss
##    SET_GAMMA         = [0.1, 0.5, 1.0, 3.0, 5.0] #gamma values -> calibration loss
##    'iteration': 50000,
##
##    'keep_prob': 0.6,
##    'Adam optimizer',
##    'lr_train': 1e-4,
##-------------------------
## Page 6 in in paper describe model as follows:
## 1 fully-connected layer for covariates,
##    - (No of nodes 3 x No: of covariates)
## 1 fully-connected layer for each competing risk event,
##   - (No of nodes 5 x No: of covariates)
## 1 fully-connected layer for output layer.
##   - (No of nodes 3 x No: of covariates)
## So 3 layers for one competing risk event,
##    4 layers for two competing risk events.
##
##---------------------------------------
## Ref: https://tensorflow.rstudio.com/guides/keras/writing_a_training_loop_from_scratch
##---------------------------------------
deephit.model <- function(x_cov, x_t, x_k, x_mask1, x_mask2, y_train,
                          x_val_cov=NULL, x_val_k=NULL, x_val_t=NULL,
                          x_val_mask1=NULL, x_val_mask2=NULL,
                          y_val=NULL,                          
                          qt, params, alpha=NULL, logdir=NULL, img.dir="images",
                          model.dir="Models"){
    ## Set params values
    model.name = params$model.name
    Loss       = setLoss(params)
    model      = NULL
    history    = NULL

    ## If we've save a model we can load it...
    if( params$load.model ){
        model = load_model_tf(sprintf("%s/%s.tf",model.dir,model.name),
                              custom_objects = {'loss': Loss}, compile = TRUE) 
    } else {
        ## ...else can build our model

        ## for training dataset
        data.cov = list()
        input.cr = list()
        data.cr  = list()
        cs_out   = list()
        
        ## for evaluation dataset
        val.cov  = list()
        val.cr   = list()

        ## number of patients
        n         = params$n

        ## number of non-censored event types
        ne        = params$ne
        
        ## number of time.points
        nt        = params$nt
    
        ## number of covariates
        nc        = params$nc

        ## cast values to type int32 
        N         = tf$cast(n,  dtype=tf$int32)
        NE        = tf$cast(ne, dtype=tf$int32)
        NT        = tf$cast(nt, dtype=tf$int32)

        ## number of neurons
        cov.neurons  = params$neurons.cov  * params$neurons ##nc
        type.neurons = params$neurons.type * params$neurons ##nc
        
        Initializer       = setInitializer(params, seed=params$seed)
        Regularizer.base  = setRegularizer(params,
                                           name=params$regularizer.base,
                                           value=params$reg.base
                                           )
        Regularizer.out   = setRegularizer(params,
                                           name=params$regularizer.out,
                                           value=params$reg.out
                                           )
        Optimizer         = setOptimizer(params)
        Metric            = list()
        Metric[[1]]       = setMetric(params$metric.name)
            
        ## placeholder for dimensions of covariates
        input.layer.name  = "input_covars"
        input.data.name   = "input_covars"

        #data.cov      = x_cov
        data.cov[[1]]      = x_cov
        names(data.cov)[1] = input.layer.name

        if( params$use.valData ){
            val.cov[[1]]       = x_val_cov
            names(val.cov)[1]  = input.layer.name
        }
        
        input    = layer_input(shape=c(dim(x_cov)[[2]]), name=input.layer.name)
    
        ## Shared subnetwork for covariates
        shared_out = input %>% layer_dense(units=cov.neurons,
                                           activation=params$act.fun.base,
                                           kernel_initializer=Initializer,
                                           kernel_regularizer=Regularizer.base,
                                           name="shared_network_out")

        ## Cause-specific subnetworks
        for( e in 1:ne ){

            input.layer.name  = sprintf("input_cr_type%g",e)
            input.data.name   = sprintf("input_cr_type%g",e)
            concat.layer.name = sprintf("concat_type%g",e)
            output.layer.name = sprintf("output_cr_type%g",e)
        
            input.cr[[e]]     = layer_input(shape=c(dim(x_cov)[[2]]),
                                            name=input.layer.name)

            data.cr[[e]]      = x_cov
            names(data.cr)[e] = input.data.name

            if( params$use.valData ){
                val.cr[[e]]       = x_val_cov
                names(val.cr)[e]  = input.data.name
            }
            
            h <- list(shared_out, input.cr[[e]]) %>%
                 layer_concatenate(axis=1L, name=concat.layer.name)
              
            cs_out[[e]] = h %>% layer_dense(units=type.neurons,
                                            activation=params$act.fun.cr,
                                            kernel_initializer=Initializer,         
                                            kernel_regularizer=Regularizer.base,
                                            name=output.layer.name)

        }                       
    
        out.neurons = tf$cast(ne*type.neurons, dtype=tf$int32)
        sm.neurons  = tf$cast(ne*nt, dtype=tf$int32)
    
        out = tf$stack(cs_out, axis=1L) ## stack referenced on subject
        out = tf$reshape(out, c(-1L, out.neurons) )    

        out <- out %>%
            layer_dropout(rate=params$drop.rate,
                          name="out_layer_w_dropout") %>%
            layer_dense(units=sm.neurons,
                        activation=params$act.fun.out,
                        kernel_initializer=Initializer,         
                        kernel_regularizer=Regularizer.out,
                        name="softmax_layer")

        out = tf$reshape(out, c(-1L, NE, NT))    
   
        ## Instantiate a DeepHit model given inputs and outputs.
        model <- DeepHitModel(inputs=list(input,input.cr),
                              outputs=list(out))

        ## compile model 
        model %>% compile(optimizer=Optimizer, 
                          loss=Loss,
                          metrics=Metric,
                          alpha=params$alpha,
                          beta=params$beta,
                          gamma=params$gamma,
                          sigma1=params$sigma1
                          )

    ## plot model
    if( params$plot.model ){
        keras$utils$plot_model(
                        model, sprintf("%s/%s.png",img.dir, model.name),
                        show_shapes = TRUE
                    )
    }
    

    ## prepare training data to pass to our model 
    train_dataset <- list(data.cov, data.cr,
                          x_t, x_k,
                          x_mask1, x_mask2,
                          y_train) %>%
        tensor_slices_dataset() %>% 
        dataset_shuffle(params$buffer) %>% 
        dataset_batch(params$batch)

    ## Prepare evaluation data
    eval_dataset = NULL
    if( params$use.valData ){    
        eval_dataset <- list(val.cov, val.cr,
                             x_val_t, x_val_k,
                             x_val_mask1, x_val_mask2,
                             y_val) %>%
            tensor_slices_dataset() %>%
            dataset_batch(params$batch)
    }

    ## fit & evaluate model
    history <- model %>% fit(train_dataset,
                             epochs = params$epochs,
                             validation_data = eval_dataset
                             )

        ## save model
        if( params$save.model ){
            save_model_tf(model, sprintf("%s/%s.tf",model.dir, model.name))
        }

    }
    
    return(list(model=model,history=history))

}


#----------------------------------------------------------------------------
#DeepHit prediction 
#----------------------------------------------------------------------------
deephit.predict <- function(model, x_cov, x_t, x_k, x_mask1, x_mask2, y_test){

    ## set parameters 
    n  = dim(x_mask1)[1]
    ne = dim(x_mask1)[2]
    nt = dim(x_mask1)[3]
    nc = dim(x_cov)[2]   
    
    data.cov = list()
    data.cr  = list()
    
    input.data.name    = "input_covars"    
    data.cov[[1]]      = x_cov
    names(data.cov)[1] = input.data.name

    ## Cause-specific subnetworks
    for( e in 1:ne ){    
        input.data.name   = sprintf("input_cr_type%g",e)
        data.cr[[e]]      = x_cov
        names(data.cr)[e] = input.data.name
    }

    ## Prepare test dataset.
    data <- list(data.cov, data.cr, x_t, x_k, x_mask1, x_mask2, y_test)
    
    ypred <- model %>% predict(data)

    return(ypred)
    
}
