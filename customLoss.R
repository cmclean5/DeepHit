##------------------------------------------
## CUSTOM LOSS FUNCTIONS
##------------------------------------------   

## Weighted Categorical Crossentropy
CCEweLoss_wrapper <- function(alpha){
    list(alpha)
    fn <- function(y_true, y_pred){
        CCEwe(y_true, y_pred,
               alpha)
        }
    custom_metric("CCEwe", fn)
}



CCEwe <- function(y_true, y_pred, alpha){

    y_true   = tf$cast(y_true, tf$float32)
    y_pred   = tf$cast(y_pred, tf$float32)
    class_we = tf$cast(alpha,  tf$float32)
    
    ## row.dim
    row.dim = tf$cast(1,dtype=tf$int32)

    ## smallest number
    epsilon = tf$keras$backend$epsilon()    
     
    ## categorical crossentropy-entropy
    Z    =  tf$math$multiply(y_true,tf$math$log(y_pred))   

    ## weighted CCE
    Zwe   = tf$math$multiply(Z,class_we)
    CCEWe = -tf$reduce_sum(Zwe,axis=row.dim)
    loss  = tf$reduce_mean(CCEWe)
    
    return(loss)
    
}

##----------------------------------------------------------------------------------------
## Example testLoss function
## Ref: https://github.com/rstudio/keras/issues/1271
##-----------------------------------------------------------------------------------------
testLoss_wrapper <- function(gamma, alpha){
    list(gamma, alpha)
    fn <- function(y_true, y_pred){
        testLoss(y_true, y_pred,
                 gamma=gamma, alpha=alpha)
        }
    custom_metric("testLoss", fn)
}



testLoss <- function(y_true, y_pred, gamma=2.0, alpha=0.25){

    y_true = tf$cast(y_true, tf$float32)
    y_true = tf$keras$backend$flatten(y_true)
    y_pred = tf$keras$backend$flatten(y_pred)    
    
    BCE        = tf$keras$losses$binary_crossentropy(y_true,y_pred)
    BCE_EXP    = tf$math$exp(-BCE)
    focal_loss = tf$reduce_mean((alpha * tf$math$pow((1-BCE_EXP), gamma) * BCE))
    
    return(focal_loss)
    
}
##-------------------


##----------------------------------------------------------------------------------------
## Example focalLoss function
## Ref: 1 https://www.programmersought.com/article/60001511310/
## Ref: 2 https://github.com/umbertogriffo/focal-loss-keras/blob/master/src/loss_function/losses.py
## Ref 3: https://gist.github.com/PsycheShaman/ea39081d9f549ac410a3a8ea942a072b
## Ref 4: https://amaarora.github.io/2020/06/29/FocalLoss.html
## Ref 5: https://towardsdatascience.com/quirky-keras-custom-and-asymmetric-loss-functions-for-keras-in-r-a8b5271171fe
##-----------------------------------------------------------------------------------------
focalLoss_wrapper <- function(gamma, alpha){
    list(gamma, alpha)
    fn <- function(y_true, y_pred){
        focalLoss(y_true, y_pred,
                  gamma=gamma, alpha=alpha)
    }
    custom_metric("focalLoss", fn)
}

focalLoss <- function(y_true, y_pred, gamma=2.0, alpha){

    alpha   <- tf$constant(alpha, dtype=tf$float32)
    gamma   <- tf$constant(gamma, dtype=tf$float32)
    epsilon <- tf$keras$backend$epsilon()

    y_true <- tf$cast(y_true,dtype=tf$float32)    

    ##clip to prevent NaNs and Infs
    y_pred <- tf$keras$backend$clip(y_pred,epsilon,1.-epsilon)   
    y_pred <- tf$cast(y_pred,dtype=tf$float32)

    alpha_t    <- y_true * alpha + (tf$ones_like(y_true)-y_true) * (1-alpha)
    y_t        <-  tf$multiply(y_true,y_pred) + tf$multiply(1-y_true,1-y_pred)

    ce         <- -tf$math$log(y_t)
    weight     <- tf$pow(tf$subtract(1.0,y_t),gamma)
    focal_loss <- tf$multiply(tf$multiply(weight,ce),alpha_t)
    loss       <- tf$reduce_mean(focal_loss)
    
    return(loss)   
        
}
