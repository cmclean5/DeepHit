## Hyperparameter flags for deephit model
PARAMS.DH <- flags(
    flag_integer('neurons', 0), ## we'll set this to 'nc', i.e. the number of covariates
    flag_integer('neurons.cov', 3),
    flag_integer('neurons.type', 5),
    flag_integer('neurons.out', 1),
    flag_string('initializer.name',"GlorotNorm"),
    flag_numeric("seed",0.0),
    flag_boolean("use_regularizer",1),
    flag_string('regularizer.base',"l2"),    
    flag_numeric("reg.base", 1e-4),
    flag_string('regularizer.out',"l1"),    
    flag_numeric("reg.out", 1e-4),
    flag_boolean("dropout",1),
    flag_numeric("drop.rate",0.6),
    flag_string("optimizer.name", "Adam"),
    flag_numeric("learn.rate", 1e-4),    
    flag_string("act.fun.base", "relu"),
    flag_string("act.fun.cr", "relu"),
    flag_string("act.fun.out", "softmax"),
    ##flag_string("loss.name", "deephit_mse"),
    flag_string("loss.name", "deephit"),
    flag_numeric("alpha", 0.001), ## 0.1, default value = 1.0
    flag_numeric("beta",  1.0), ## 2.0
    flag_numeric("gamma", 0.0), ## 1.0, default value = 1.0
    flag_numeric("sigma1", 0.1),
    flag_numeric("n", 0),  ## no: of patients
    flag_numeric("ne", 0), ## no: of non-censored events
    flag_numeric("nt", 0), ## no: of time.points 
    flag_numeric("nc", 0), ## no: of covariates, set in function 'train.deephit.model'
    flag_string("metric.name",  "mse"),
    #flag_string("metric.name2", ""),
    flag_integer("batch", 64),
    flag_integer("buffer", 1024),
    flag_integer("epochs", 25),#50),##200),
    flag_boolean("verbose", 0),
    #flag_boolean("use.valSplit",0),
    flag_boolean("use.valData",1),
    flag_numeric("ratio_train",0.8),
    flag_numeric("ratio_val",0.2),
    flag_string("model.name", "deephit.model"),
    flag_boolean("use.callback",0),
    flag_boolean("save.model",0),
    flag_boolean("load.model",0),
    flag_boolean("plot.model",1)
)
