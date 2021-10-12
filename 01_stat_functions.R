

### HELP FUNCTIONS
## Model performance test - Root Mean Square error
#' @param Xa - observed vector
#' @param Xb - predicted vector
#' @return RMSE 

RMSE  <-  function(Xa,Xb){
  N <- length(Xa)
  rmse <- (sum((Xa-Xb)^2)/N)^0.5
  return(rmse)
}
MAE  <-  function(Xa,Xb){
  N <- length(Xa)
  me <- (sum((Xa-Xb)^1)/N)^1
  return(me)
}

# model efficiency
ME  <-  function(Xa,Xb){
  me <- 1- (sum((Xa-Xb)^2)/sum((Xa-mean(Xa))^2))
  return(me)
}
