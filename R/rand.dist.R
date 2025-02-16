rand.dist <- function(exactt.object, beta0, variable.name = NULL){
  
  if(!methods::is(exactt.object, "exactt")){
   stop('"exactt.object" must be of class "exactt".')
  }
  
  rand.dist.list <- vector("list")
  
  for(name in names(exactt.object$geometry)){
    if(name %in% variable.name){
      
      num.temp <- sapply(beta0,
                         function(x){
                           exactt.object$geometry[[name]]$line.data.num$b + 
                             exactt.object$geometry[[name]]$line.data.num$m * x
                         })
      
      if(!is.null(exactt.object$geometry[[name]]$line.data.denom)){
        denom.temp <- sapply(exactt.object$geometry[[name]]$line.data.denom,
                             function(x){
                               sqrt(predict(x, beta0))
                             }) |>
          t()
        
        rand.dist.temp <- num.temp/denom.temp
      } else{
        rand.dist.temp <- num.temp
      }
      
      if(is.null(exactt.object$call$side)){
        rand.dist.temp <- abs(rand.dist.temp)
      }
      
      colnames(rand.dist.temp) <- paste0("beta.null.", beta0)
      
      rand.dist.list[[name]] <- data.frame(rand.dist.temp)
      
    }
  }
  
  return(rand.dist.list)
}
