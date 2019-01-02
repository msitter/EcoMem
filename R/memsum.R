#####################################################################
#####################################################################
# Summarize ecomem output
#####################################################################
#####################################################################

memcoef = function(x,cred.int=0.95,verbose=TRUE,...){
  
  n.chains = length(x$post.samps)
  params = names(x$post.samps[[1]])
  quantiles = (1-cred.int)/2 + c(0,cred.int)
  
  summ = function(y){
    if (is.null(dim(y))){
      tab = t(c(mean(y),median(y),quantile(y,quantiles)))
      } else {
      tab = t(apply(y,2,function(z)c(mean(z),median(z),quantile(z,quantiles))))
      }
    tab = as.data.frame(tab)
    names(tab) = 
      c("mean","median",paste("q",quantiles[1],sep=""),paste("q",quantiles[2],sep=""))
    return(tab)
  }
  
  # summarize beta
  coef.tab = cbind(data.frame(var=x$pred.vars,
                              summ(do.call("rbind",
                                           lapply(x$post.samps,
                                                  function(y)y$beta)))))
  
  # summarize sig.y
  if ("sig.y" %in% params){
    
    coef.tab = rbind(coef.tab,
                     data.frame(var="sig.y",summ(unlist(lapply(x$post.samps,
                                                           function(y)y$sig.y)))))
    
  }
  
  # summarize memory functions
  if ("w" %in% params){
    
    for (i in 1:length(mem.vars)){
      coef.tab = 
        rbind(coef.tab,data.frame(var=paste(
          mem.vars[i],":w",1:ncol(x$post.samps[[1]]$w[[i]])-1,
          sep=""),summ(
            do.call("rbind",lapply(x$post.samps,
                                   function(y)y$w[[i]])))))
    }
    
    
  }
  
  if (verbose == TRUE){
    print(coef.tab,digits=3,row.names=FALSE)
  }
  
  return(coef.tab)
  
}
