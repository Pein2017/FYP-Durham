
getTrainingInfo = function(x) {
  attr(x$learner.model, "mlr.train.info") # attr(x, "mlr.train.info") 
}

attachTrainingInfo = function(x, info) {
  attr(x, "mlr.train.info") = info
  x
}

getFixDataInfo = function(data, restore.levels = FALSE, factors.to.dummies = FALSE, ordered.to.int = FALSE) {
  
  cl = vcapply(data, getClass1)
  factors = lapply(data[cl == "factor"], levels)
  ordered = lapply(data[cl == "ordered"], levels)
  
  makeS3Obj("FixDataInfo",
            factors = factors,
            ordered = ordered,
            restore.levels = restore.levels,
            factors.to.dummies = factors.to.dummies && length(factors) > 0L,
            ordered.to.int = ordered.to.int && length(ordered) > 0L
  )
}

fixDataForLearner = function(data, info) {
  
  cn = c(names(info$factors), names(info$ordered))
  not.found = which.first(cn %nin% names(data))
  if (length(not.found) > 0L) {
    stopf("Column '%s' found in info, but not in new data", cn[not.found])
  }
  
  if (info$restore.levels) {
    if (!info$factors.to.dummies && length(info$factors) > 0L) {
      cols = names(info$factors)
      data[cols] = Map(factor, x = data[cols], levels = info$factors)
    }
    if (!info$ordered.to.int && length(info$ordered) > 0L) {
      cols = names(info$ordered)
      data[cols] = Map(factor, x = data[cols], levels = info$ordered, ordered = TRUE)
    }
  }
  
  if (info$factors.to.dummies) {
    cols = names(info$factors)
    new.cols = Map(function(x, lvls) {
      as.data.frame(setNames(lapply(lvls, "==", x), lvls))
    }, x = data[cols], lvls = info$factors)
    data = cbind(dropNamed(data, cols), do.call(cbind, new.cols))
  }
  
  if (info$ordered.to.int) {
    cols = names(info$ordered)
    data[cols] = lapply(data[cols], as.integer)
  }
  
  data
}
## new learner

makeRLearner.regr.trunc.lasso = function() {
  makeRLearnerRegr(
    cl = "regr.trunc.lasso",
    package = "glmnet",
    par.set = makeParamSet(
      makeDiscreteLearnerParam(id = "family", values = c("gaussian", "poisson"), default ="gaussian"),
      makeNumericLearnerParam(id = "alpha", default = 1, lower = 0, upper = 1),
      makeNumericLearnerParam(id = "s", lower = 0, when = "predict") #,
      # makeLogicalLearnerParam(id = "exact", default = FALSE, when = "predict"),
      # makeIntegerLearnerParam(id = "nlambda", default = 100L, lower = 1L),
      # makeNumericLearnerParam(id = "lambda.min.ratio", lower = 0, upper = 1),
      # makeNumericVectorLearnerParam(id = "lambda", lower = 0),
      # makeLogicalLearnerParam(id = "standardize", default = TRUE),
      # makeLogicalLearnerParam(id = "intercept", default = TRUE),
      # makeNumericLearnerParam(id = "thresh", default = 1e-07, lower = 0),
      # makeIntegerLearnerParam(id = "dfmax", lower = 0L),
      # makeIntegerLearnerParam(id = "pmax", lower = 0L),
      # makeIntegerVectorLearnerParam(id = "exclude", lower = 1L),
      # makeNumericVectorLearnerParam(id = "penalty.factor", lower = 0, upper = 1),
      # makeNumericVectorLearnerParam(id = "lower.limits", upper = 0),
      # makeNumericVectorLearnerParam(id = "upper.limits", lower = 0),
      # makeIntegerLearnerParam(id = "maxit", default = 100000L, lower = 1L),
      # makeLogicalLearnerParam(id = "standardize.response", default = FALSE),
      # makeNumericLearnerParam(id = "fdev", default = 1.0e-5, lower = 0, upper = 1),
      # makeNumericLearnerParam(id = "devmax", default = 0.999, lower = 0, upper = 1),
      # makeNumericLearnerParam(id = "eps", default = 1.0e-6, lower = 0, upper = 1),
      # makeNumericLearnerParam(id = "big", default = 9.9e35),
      # makeIntegerLearnerParam(id = "mnlam", default = 5, lower = 1),
      # makeNumericLearnerParam(id = "pmin", default = 1.0e-9, lower = 0, upper = 1),
      # makeNumericLearnerParam(id = "exmx", default = 250),
      # makeNumericLearnerParam(id = "prec", default = 1e-10),
      # makeIntegerLearnerParam(id = "mxit", default = 100L, lower = 1L),
      # makeUntypedLearnerParam(id = "offset", default = NULL),
      # makeDiscreteLearnerParam(id = "type.gaussian", values = c("covariance", "naive"), requires = quote(family == "gaussian")),
      # makeLogicalLearnerParam(id = "relax", default = FALSE)
    ),
    properties = c("numerics"), #  , "factors", "ordered", "weights"
    par.vals = list(s = 0.01),
    name = "GLM with Lasso or Elasticnet Regularization",
    short.name = "glmnet",
    callees = c("glmnet", "glmnet.control", "predict.glmnet")
  )
}



trainLearner.regr.trunc.lasso = function (.learner, .task, .subset,  ...) 
{
  d = getTaskData(.task, .subset, target.extra = TRUE)
  info = getFixDataInfo(d$data, factors.to.dummies = TRUE, ordered.to.int = TRUE)
  args = c(list(x = as.matrix(fixDataForLearner(d$data, info)), y = d$target), list(...))
  glmnet::glmnet.control(factory = TRUE)
  saved.ctrl = glmnet::glmnet.control()
  is.ctrl.arg = names(args) %in% names(saved.ctrl)
  if (any(is.ctrl.arg)) {
    on.exit(glmnet::glmnet.control(factory = TRUE))
    do.call(glmnet::glmnet.control, args[is.ctrl.arg])
    args = args[!is.ctrl.arg]
  }
  attachTrainingInfo(do.call(glmnet::cv.glmnet, args), info)
}

predictLearner.regr.trunc.lasso = function (.learner, .model, .newdata , ...) 
{
  info = getTrainingInfo(.model)
  .newdata = as.matrix(fixDataForLearner(.newdata, info))
  mins = .model$learner.model$lambda.min
  print(mins)
  p = predict(.model$learner.model, newx = .newdata, s = mins)
  index = which( p < 0)
  p[index] = 0
  return ( as.vector(p) )
}

