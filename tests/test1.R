library(reticulate)
library(tensorflow)
library(tfdatasets) ## format training dataset in tensorflow format
library(keras)      ## Deep learning packages, Keras API and tensorflow

CustomModel <- new_model_class(
  classname = "CustomModel",
  train_step = function(data) {
    # Unpack the data. Its structure depends on your model and on what you pass
    # to `fit()`.  A third element in `data` is optional, but if present it's
    # assigned to sample_weight. If a thrid element is missing, sample_weight
    # defaults to NULL
    c(x, y, sample_weight = NULL) %<-% data
    
    with(tf$GradientTape() %as% tape, {
      y_pred <- self(x, training = TRUE)  # Forward pass
      # Compute the loss value.
      # The loss function is configured in `compile()`.
      loss <- self$compiled_loss(y,
                                 y_pred,
                                 sample_weight = sample_weight,
                                 regularization_losses = self$losses)
    })
    
    # Compute gradients
    trainable_vars <- self$trainable_variables
    gradients <- tape$gradient(loss, trainable_vars)
    
    # Update weights
    self$optimizer$apply_gradients(zip_lists(gradients, trainable_vars))
    
    # Update the metrics.
    # Metrics are configured in `compile()`.
    self$compiled_metrics$update_state(y, y_pred, sample_weight = sample_weight)
    
    # Return a named list mapping metric names to current value.
    # Note that it will include the loss (tracked in self$metrics).
    results <- list()
    for (m in self$metrics)
      results[[m$name]] <- m$result()
    results
  }
)


# Construct and compile an instance of CustomModel

inputs <- layer_input(shape(32))
outputs <- inputs %>% layer_dense(1)
model <- CustomModel(inputs, outputs)
model %>% compile(optimizer = "adam",
                  loss = "mse",
                  metrics = "mae")

# You can now use sample_weight argument

x <- k_random_uniform(c(1000, 32))
y <- k_random_uniform(c(1000, 1))
sw <- k_random_uniform(c(1000, 1))
model %>% fit(x, y, sample_weight = sw, epochs = 3)
