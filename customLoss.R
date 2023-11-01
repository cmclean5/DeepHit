##------------------------------------------
## CUSTOM LOSS FUNCTIONS
##------------------------------------------   

## Define wrapper for our custom loss function
deephitLoss_mse_wrapper <- function(alpha, beta, gamma, sigma1) {
    fn <- function(y_true, y_pred, events, time, mask1, mask2) {
      deephitLoss_mse(y_true, y_pred,
                      events = events, time = time, mask1 = mask1, mask2 = mask2,
                      alpha = alpha, beta = beta, gamma = gamma, sigma1 = sigma1)
  }
    custom_metric("deephitLoss_mse", fn)
}

## Define custom MSE loss function
deephitLoss_mse <- function(y_true, y_pred,
                            events=NULL, time=NULL, mask1=NULL, mask2=NULL,
                            alpha = 1, beta = 1, gamma = 1, sigma1 = 0.1) {

    y_true <- tf$cast(y_true, dtype = tf$float32)
    y_pred <- tf$cast(y_pred, dtype = tf$float32)
    
    loss = tf$reduce_mean(tf$math$squared_difference(y_true,y_pred))
    return(loss)
}


## Define wrapper for DeepHit custom loss function
deephitLoss_wrapper <- function(alpha, beta, gamma, sigma1) {
  fn <- function(y_true, y_pred, events, time, mask1, mask2) {
      deephitLoss(y_true, y_pred,
                  events = events, time = time, mask1 = mask1, mask2 = mask2,
                  alpha = alpha, beta = beta, gamma = gamma, sigma1 = sigma1)
  }
  #return(fn)
  custom_metric("deephitLoss", fn)
}


## Define custom DeepHit loss function...
## this is loss function described in paper.
deephitLoss <- function(y_true, y_pred,
                        events=NULL, time=NULL, mask1=NULL, mask2=NULL,
                        alpha = 1, beta = 1, gamma = 1, sigma1 = 0.1) {

    epsilon <- tf$keras$backend$epsilon()
    x_dims  <- dim(mask1)

    n  <- x_dims[1]
    ne <- x_dims[2]
    nt <- x_dims[3]
    N  <- tf$cast(n, dtype = tf$int32)
    NE <- tf$cast(ne, dtype = tf$int32)
    NT <- tf$cast(nt, dtype = tf$int32)

    y_true <- tf$cast(y_true, dtype = tf$float32)
    y_pred <- tf$cast(y_pred, dtype = tf$float32)
    mask1  <- tf$cast(mask1, dtype = tf$float32)
    mask2  <- tf$cast(mask2, dtype = tf$float32)
    events <- tf$cast(events, dtype = tf$int32)
    time   <- as_tensor(time, shape = c(n, 1), dtype = tf$float32)

    ##------------------------------------------
    ## LOSS_1: Log-likelihood loss
    ## As described in deephit's github code 
    ##------------------------------------------    
    I1    <- tf$sign(events)
    I1    <- tf$cast(I1, dtype = tf$float32)

    ## Hard coded for the moment, but choose
    ## which implementation of loss1 to run?
    run.authors = 1 ## run authors implementation of loss1
    #run.authors = 0 ## run our     implementation of loss1  

    if( run.authors ){

        tmp1  <- tf$reduce_sum(tf$reduce_sum(mask1 * y_pred, axis = 2L), axis = 1L, keepdims = TRUE)
        tmp1  <- tf$keras$backend$clip(tmp1, epsilon, 1 - epsilon)
        ## 
        ## for uncensored: log P(T=t,K=k|x)   'I1*log(tmp1)'
        ## for censored:   log \sum P(T>t|x)  '(1.0-I1)*log(tmp1)'
        loss1 <- -tf$reduce_sum(I1 * log(tmp1) + (1 - I1) * log(tmp1))
        ##------------------------------------------

    } else {
        
        ##------------------------------------------
        ## TEST... code below should match eqn 2 in the deepHit paper,
        ## as unsure about the definition of loss1 using deephit's github code
        ## above?
        ##------------------------------------------
        ## Gather uncensored patients
        crLoss=list()
        for( e in 1:ne ){
            E      = tf$cast(e, dtype=tf$int32)
            crTemp = tf$equal(events,E)
            s_mask = tf$boolean_mask(mask1[,(E-1L),],crTemp)
            y_mask = tf$boolean_mask(y_pred[,(E-1L),],crTemp)
            crSum  = tf$reduce_sum(s_mask*y_mask, axis=1L)
            crSum  = tf$keras$backend$clip(crSum, epsilon, 1-epsilon)
            crLoss[[e]] = tf$reduce_sum(log(crSum))
        }
        tmp1 = -tf$reduce_sum(tf$stack(crLoss))    
        ##
        ## Gather censored patients
        cenLoss=list()
        for( e in 1:ne ){
            E        = tf$cast(e, tf$int32)
            cenTemp  = (1-I1)
            cen_mask = tf$boolean_mask(mask1[,(E-1L),],cenTemp)
            y_mask   = tf$boolean_mask(y_pred[,(E-1L),],cenTemp)
            cenLoss[[e]] = tf$reduce_sum(cen_mask*y_mask, axis=1L)
        }
        ## tmp2 ==> dims(cr, No: patients)
        tmp2 = tf$reduce_sum(tf$stack(cenLoss, axis=1L), axis=1L)
        tmp2 = tf$keras$backend$clip(tmp2, epsilon, 1-epsilon)
        tmp2 = -tf$reduce_sum(log(1-tmp2))
        loss1 = (tmp1+tmp2)
        ##------------------------------------------
    }
    
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
    ##one_vector = tf$ones_like(time, dtype=tf$float32)    
    eta       = list()
    for( e in 1:ne ){
        E     = tf$cast(e,dtype=tf$int32)
        one_vector = tf$ones_like(time, dtype=tf$float32)
        I2    = tf$cast(tf$equal(events,E),dtype=tf$float32)
        I2    = tf$linalg$diag(tf$squeeze(I2))
        tmp_e = tf$reshape(tf$slice(y_pred, c(0L,(E-1L),0L), c(-1L,1L,-1L)), c(-1L, NT))
        tmp_e = tf$cast(tmp_e, dtype=tf$float32)
        R     = tf$matmul(tmp_e, tf$transpose(mask2))
        DIAGr = tf$reshape(tf$linalg$diag_part(R), c(-1L,1L))
        R     = tf$matmul(one_vector, tf$transpose(DIAGr)) - R
        R     = tf$transpose(R)
        T     = tf$nn$relu(tf$sign(
                                  tf$matmul(one_vector, tf$transpose(time)) -
                                  tf$matmul(time, tf$transpose(one_vector))
                              ))
        T        = tf$matmul(I2, T)
        tmp_eta  = tf$reduce_mean(T * tf$exp(-R/sigma1), axis=1L, keepdims=TRUE)
        eta[[e]] = tmp_eta
    }
    eta   = tf$stack(eta, axis=1L)
    eta   = tf$reduce_mean( tf$reshape(eta, c(-1L,NE)), axis=1L, keepdims=TRUE)
    loss2 = tf$reduce_sum(eta)
    ##------------------------------------------

    
    ##------------------------------------------
    ## LOSS_3: Calibration Loss
    ##------------------------------------------    
    eta = list()
    for( e in 1:ne ){
        E     = tf$cast(e,dtype=tf$int32)      
        I2    = tf$cast(tf$equal(events,E),dtype=tf$float32)
        tmp_e = tf$reshape(tf$slice(y_pred, c(0L,E-1L,0L), c(-1L,1L,-1L)), c(-1L, NT))
        tmp_e = tf$cast(tmp_e, dtype=tf$float32)
        ## Writing in author's code, but does not seem to work?
        ##R     = tf$reduce_sum(tmp_e * mask2, axis=0L)
        ##tmp_eta  = tf$reduce_mean( (R-I2)**2, axis=1L, keepdims=TRUE)
        ## Below code does work.
        R     = tf$reduce_sum(tmp_e * mask2, axis=1L)
        tmp_eta  = tf$reduce_mean( (R-I2)**2, axis=0L, keepdims=TRUE)
        eta[[e]] = tmp_eta
    }

    eta = tf$stack(eta, axis=1L)
    eta = tf$reduce_sum( tf$reshape( eta, c(-1L,NE)), axis=1L, keepdims=TRUE)

    loss3 = tf$reduce_sum(eta)   
    ##------------------------------------------    
   
    loss = loss1*alpha + loss2*beta + loss3*gamma
    
    return(loss)

}
