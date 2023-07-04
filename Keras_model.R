##---------------------------------------------------------------------------------------
## Custom Keras Models
##-------------------------------------------------------------------------------------------

## Define the DeepHit model class
DeepHitModel <- keras::new_model_class(
  classname = "DeepHitModel",    
  
  ## define our compile function to our custom loss function and parameter values
  compile = function(optimizer = "adam",
                     loss = NULL, metrics="mae",                     
                     alpha = 1, beta = 1, gamma = 1, sigma1 = 0.1) {      

      ## Use custom loss if provided
      if (is.null(loss)) {
          loss <- deephitLoss_wrapper(alpha = alpha, beta = beta, gamma = gamma,
                                      sigma1 = sigma1)
      }

      if (is.character(metrics)) {
          metrics <- list(metrics)
      }

      ## We need this function to update our custom loss metric
      ## We can define the loss tracker inside our compile function
      self$loss_tracker <- metric_mean(name = "loss")       

      ## define your custom metric result retrieval function
      self$get_metric_results <- function() {
          metric_results <- list()
          for (metric in self$metrics) {
              metric_results[[metric$name]] <- metric$result()
          }
          return(metric_results)
      }
      
      ## set parameters in our model
      super$compile(
                optimizer = optimizer,
                loss = loss,
                metrics = metrics
            )
  },

    ## define training step for our model
    train_step = function(data) {
      c(x1, x2, x_t, x_k, x_mask1, x_mask2, y_true) %<-% data

      with(tf$GradientTape() %as% tape, {
        y_pred <- self(list(x1, x2), training = TRUE)

        ## Compute our own loss, i.e. the metric loss
        loss <- self$loss(
            y_true = y_true,
            y_pred = y_pred,
            events = x_k,
            time = x_t,
            mask1 = x_mask1,
            mask2 = x_mask2
        )       

        ## Compute the regularization loss
        regularization_loss <- tf$reduce_sum(self$losses)

        ## total loss, i.e. our metric loss plus regularization loss
        loss <- loss + regularization_loss
        
      })

      ## Compute gradients
      trainable_vars <- self$trainable_variables
      gradients <- tape$gradient(loss, trainable_vars)

      ## Update weights
      self$optimizer$apply_gradients(zip_lists(gradients, trainable_vars))

      ## Update the custom loss metric
      ## This ensures that the loss_tracker object is updated
      ## with the current loss value during each training step.
      self$loss_tracker$update_state(loss)
  
      ## Update the other metrics
      self$compiled_metrics$update_state(y_true, y_pred)
  
      ## Retrieve the metric results
      ##metric_results <- list()
      ##for (metric in self$metrics) {
      ##    metric_results[[metric$name]] <- metric$result()
      ##}

      ## Retrieve the metric results using the helper function
      metric_results <- self$get_metric_results()
    
      
      ##Finally, return the loss and metric values together
      return(c(loss = self$loss_tracker$result(), metric_results))      
    
    },
 
  ## Define here the metrics that will be tracked during training.
  ## It is separate from the metrics used in the train_step method to update
  ## and retrieve the metric results.
  metrics = mark_active(function() {
    list(self$loss_tracker)
  }),

  ## add function test_step to evaluate model on test data.
  test_step = function(data){
      ## Unpack the data. Its structure depends on your model and
      ## on what you pass to `fit()`.
      c(x1, x2, x_t, x_k, x_mask1, x_mask2, y_true) %<-% data

      ## Compute predictions
      y_pred <- self(list(x1,x2), training=FALSE)

      ## Updates the metrics tracking the loss
      loss <- self$loss(
            y_true = y_true,
            y_pred = y_pred,
            events = x_k,
            time = x_t,
            mask1 = x_mask1,
            mask2 = x_mask2
        )       

      ## Compute the regularization loss
      regularization_loss <- tf$reduce_sum(self$losses)

      ## total loss, i.e. our metric loss plus regularization loss
      loss <- loss + regularization_loss
       
      ## update loss
      self$loss_tracker$update_state(loss)
      
      ## Update the metrics.
      self$compiled_metrics$update_state(y_true, y_pred)

      ## Retrieve the metric results using the helper function
      metric_results <- self$get_metric_results()
      
      ## Return a named list mapping metric names to current value.
      ## Note that it will include the loss (tracked in self$metrics).
      return(c(loss = self$loss_tracker$result(), metric_results))  
      
  },

  ## add function evaluation_step to evaluate model on test data.
  evaluation_step = function(data){
      ## Unpack the data. Its structure depends on your model and
      ## on what you pass to `fit()`.
      c(x1, x2, x_t, x_k, x_mask1, x_mask2, y_true) %<-% data

      ## Compute predictions
      y_pred <- self(list(x1,x2), training=FALSE)

      ## Updates the metrics tracking the loss
      loss <- self$loss(
            y_true = y_true,
            y_pred = y_pred,
            events = x_k,
            time = x_t,
            mask1 = x_mask1,
            mask2 = x_mask2
        )       

      ## Compute the regularization loss
      regularization_loss <- tf$reduce_sum(self$losses)

      ## total loss, i.e. our metric loss plus regularization loss
      loss <- loss + regularization_loss
       
      ## update loss
      self$loss_tracker$update_state(loss)
      
      ## Update the metrics.
      self$compiled_metrics$update_state(y_true, y_pred)

      ## Retrieve the metric results using the helper function
      metric_results <- self$get_metric_results()
      
      ## Return a named list mapping metric names to current value.
      ## Note that it will include the loss (tracked in self$metrics).
      return(c(loss = self$loss_tracker$result(), metric_results))  
      
  }


)
