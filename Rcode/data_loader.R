### load the packages
#library(Rgraphviz)
#library(rJava)
library(AppliedPredictiveModeling)
transparentTheme(trans = .4)
library(caret)
library(gbm)
library(mboost)
library( dplyr )
library(KRLS) # Polynomial Kernel Regularized Least Squares
library(monomvn) # Bayesian Ridge Regression (Model Averaged)
library(doParallel)
library( e1071 )
library( xgboost )
library( kernlab )
library( KRLS )
library(randomForest)
#library(RWeka)
library( partykit )
library( mlbench)
library(stats)

### parallel function 
removeParallel = function() {
  env = foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

### denoise functioon
mean_xt = function( data  , center , windows_size )
{
  data = c( data[  (center - windows_size) : (center-1)  ] , data[(center+2) : (center + windows_size+1)] )
  radial = length( data ) / 2
  weight = c(  1:radial , radial:1 )
  weight = weight / sum(weight)
  res = weighted.mean( data , w = weight )  
  #res = mean(data)
  return ( res )
}

### Pre-process function 
remove_columns = function( data_train , data_test = NULL)
{
  remove_col = c()
  for (j in 1:(ncol(data_train)-2) ) 
  {
    if (  mean(data_train[,j]) < 0.01 | length( table( data_train[,j] ) ) == 1  )
      remove_col = c(remove_col , j)
    
  }
  return ( remove_col )
}

correct = function( data_train , thres = 0.01 )
{
  data_train =  data_train[complete.cases(data_train) ,]
  for ( j in 1:ncol( data_train ) )
  {
    index = which( data_train[,j] < thres)
    data_train[index,j] = 0
  }
  return (data_train)
}

### Counting classes function
count_dataframe = function( data , thres =5 )
{
  index = c()
  collections = c()
  for ( j in 1:ncol(data) )
  {
    cont = table( data[,j] )
    collection = length( cont )
    if ( collection < thres ) 
    {
      index = c(index , j )
      collections = c( collections , collection )
      
    }
  }
  return ( list( columns = index  , count = collections ) )
}

### Ploting function
plot_data = function(  data , index , counts=-1 )
{
  par( mfrow = c(2,1 ) )
  name = names(data)[index]
  main = paste( name , "with counts: ", counts)
  plot( data[,index] , data$interpolated_acetic , main = main , xlab = name , ylab = 'acid' ) 
  plot( data$EFT_H , data.matrix( data[,index] ) , xlab = 'EFT ' , ylab = name)
}

### Merge by interpolation
cb_merge = function( ac19 , data  )
{
  x = ac19$EFT
  y = ac19$acetic
  x_out = data$EFT_H 
  interpolated__acetic = spline(x = x ,y = y , xout = x_out )$y
  interpolated__acetic[ which(interpolated__acetic < 0) ] = 0
  ####### Setting thresholds
  ## disolved 16~20
  ## methane 21~25
  ## oxygen flow 26~30
  cols = c(16:30)
  #beforeMean = apply(data[,cols] , 2 , mean)
  for (j in cols)
  {
    index = which( data[,j] < 0.5 )
    data[index,j] = 0
  }
  #afterMean = apply(data[,cols] , 2 , mean)
  
  ############ Ading ratio
  ratioCols = c(26:30)
  k=1
  newNames = c()
  new = c()
  for (j in ratioCols)
  {
    name = paste( "O2_Ch4_ratio_staion",k,sep="")
    newNames = c(newNames , name )
    k = k+1
    column = data[,j] / data[,(j-5) ]
    index = which( column == Inf | is.na(column) )
    column[index] = 0
    new = cbind(new ,  column)
  }
  colnames(new) = newNames
  data = cbind(data , new , interpolated_acetic = interpolated__acetic ) 
  return (data)
}

### Adding the ratio

### TPP19
data_loader = function( load_19 = True, Hour = 361)
{
  data = read.csv("C:/Users/85212/Desktop/Pro-data/TPP19_online.csv", header=TRUE)
  ac19 = read.csv( "C:/Users/85212/Desktop/Pro-data/TPP19_ac.csv" ,header=TRUE)
  #data = read.csv("/home/gmqx85/FYP/Pro-Data/TPP19_online.csv", header=TRUE)
  #ac19 = read.csv( "/home/gmqx85/FYP/Pro-Data/TPP19_ac.csv" ,header=TRUE)
  data= data[,-c(2,3)]
  ac19 = dplyr:: filter( ac19 , EFT_H < Hour)
  data = dplyr:: filter( data , EFT_H < Hour)
  index = c( which( names(data) == 'headspace_pressure') , which(names(data) == 'pressure_at_position_4') )
  data = data[,-index]  ## remove constant headspac19e.pressure
  data = correct(data)
  remove_index = remove_columns( data )
  data = data[, - remove_index ]
  
  return  (   list( data_19 = cb_merge( ac19, data) , ac_19 = ac19)   )
}


### inputs data, cols need to be removed, selected subsets of row that
### is specified, threshold of VIF value that need to be removed
vifSelection = function( data, removeCols = c(1) , vifThres =30 , subsets = NULL )
{
  library(car)
  # if ( is.null(removeCols) ) removeCols = c()
  nm = names(data)
  p = ncol( data )
  run = 1
  while ( run < p )
  {
    if ( is.null( subsets) )
    {
      reg1 = lm( interpolated_acetic ~., data = data[,-removeCols])
    }
    else
    {
      reg1 = lm( interpolated_acetic ~., data = data[subsets, -removeCols])
    }
    vif1 =  vif(reg1) 
    vif1 =  sort( vif1 , decreasing = TRUE)
    removeTarget = vif1[1]
    if ( removeTarget > vifThres )
    {
      addCol = which( nm == names(removeTarget) )
      removeCols = c( removeCols,  addCol )
      run = run + 1
    }
    else break
  }
  return ( list(name = nm[removeCols] , columns = removeCols )    )   
}
## Inputs data and preProcess
## return preprocessed data and processing collector for further prediction 
process = function( data , preProcess = c("center" , "scale") , normalization = NULL  )
{
  p = ncol(data)
  if (is.null( normalization ) )
  { normalization = preProcess(data[,-p] , method = preProcess) }
  
  x =  as.data.frame( predict(normalization, data[,-p]) )
  preData = cbind(x, interpolated_acetic = data[,p])
  return ( list( data = preData, normalization = normalization) )
}


### inputs methods with 'model','RMSE' and 'MErrorPercentage'
### outputs the matrix/vector of importance following the orders of 
#   'names_of_predictors'
getVarImp = function( method , names_of_predictors ) 
{
  current_importance = varImp( method$model)$importance
  if ( length( current_importance$Overall ) == length( names_of_predictors) )
  {
    variable_importance = current_importance$Overall
    return ( variable_importance )
  }
  else
  {
    cat(method$model$method,"has few predictors, needs modification")
    current_variables_names = rownames( current_importance )
    index = c()
    for ( current_name in current_variables_names)
    {
      index = c(index , which( names_of_predictors == current_name) )
    }
    emptyVector = rep( 0 ,length( names_of_predictors ) )
    emptyVector[index] = current_importance$Overall
    return ( emptyVector )
  }
}

### Wrapper training function
# parallel is set to be off
# Test data will be off 
modelTrain = function( data, method, trControl = NULL, tuneGrid = NULL, weights = NULL, tuneLength=2 , preProcess = NULL , cpu = 0 , test = NULL)
{
  begin = Sys.time()
  if ( cpu != 0  )
  {
    cl = makeCluster( cpu )
    registerDoParallel(cl)
  }
  model = train( interpolated_acetic~. ,
                 data = data,
                 method = method,
                 trControl = trControl,
                 tuneGrid = tuneGrid,
                 tuneLength = tuneLength,
                 preProcess = preProcess,
                 weights = weights
  )
  if ( cpu != 0  )
  {
    stopCluster( cl )
    removeParallel()
  }
  if ( is.null( test ) ) 
  {
    end = Sys.time()
    cat("model:",method,"takes:",(end-begin),"to finish with",cpu,"cpus \n")
    return (model) 
  }
  else
  {
    Predict = predict( model, newdata = test )
    y = test$interpolated_acetic
    RMSE =  mean( (y - Predict )**2 )**0.5
    MErrorPercentage = sum( abs(y - Predict ) ) / sum( y )
    end = Sys.time()
    cat("model:",method,"takes:",(end-begin),"to finish with",cpu,"cpus")
    return ( list( model = model , RMSE = RMSE , MErrorPercentage = MErrorPercentage) )
  }
}


models_train = function( data , models_wrapper , train_control = NULL, tune_length = 2, test = NULL, weights = NULL,cpu = 0)
  ### inputs a 'models_wrapper' with a list of models with class ('method','tune_grid')
  # lm = list( method = "lm" , tune_grid = NULL)
  # mars = list( method = 'earth' , tune_grid = marsGrid)
  # lasso = list( method = 'glmnet' , tune_grid = lassoGrid)
  # models_wrapper = list( lm = lm , lasso = lasso )
  
{
  models_names = names(wrappers)
  order_of_current_model = 0
  if ( cpu != 0 )
  {
    cl = makeCluster( cpu )
    registerDoParallel(cl)
  }
  total_start = Sys.time()
  collectors = list()
  for ( model in models_wrapper)
  {
    order_of_current_model = order_of_current_model + 1
    model_name = models_names[order_of_current_model]
    method = model$method
    tune_grid = model$tune_grid
    weights = model$weights 
    pre_process = model$pre_process
    model_start = Sys.time()
    model_result = train(   interpolated_acetic~. ,
                            data = data,
                            method = method,
                            trControl = train_control,
                            tuneGrid = tune_grid,
                            tuneLength = tune_length,
                            preProcess = pre_process,
                            weights = weights
    )
    if ( is.null( test ) ) 
    {
      collectors = c( collectors , list(model = model_result) )
      names(collectors)[order_of_current_model] = model_name
      next
    }
    else
    {
      Predict = predict( model_result, newdata = test )
      y = test$interpolated_acetic
      rmse =  mean( (y - Predict )**2 )**0.5
      mean_per_error = sum( abs(y - Predict ) ) / sum( y )
      collectors = c( collectors , 
                      list(  model = list( model = model_result , rmse = rmse , mean_per_error = mean_per_error)  )
      )
      names(collectors)[order_of_current_model] = model_name
    }
    model_end = Sys.time()
    cat("model:",method,"costs:", (model_end - model_start),"\n")
    cat("\n")
  }
  total_end = Sys.time()
  cat("jobs finished with total time:", (total_end - total_start) , "\n")
  cat("\n")
  if ( cpu != 0  )
  {
    stopCluster( cl )
    removeParallel()
  }
  return (collectors)
}

# t is the input of EFT of online data; t_k is a sequence of EFT within offline data 

get_single_weight = function( t , t_ofline) 
{
  weight = sum(  abs( t - t_ofline ) )
  return (weight)
}
# input a sequence of t
get_data_weights = function( t , t_ofline , scale = TRUE)
{
  weights = c()
  for (i in 1:length(t) )
  {
    weight = 1 / get_single_weight( t[i] , t_ofline )
    weights = c( weights , weight)
  }
  if (scale) 
  {
    weights = weights / sum( weights )
  }
  return (weights)
}

cat("data loaded successfully \n")

