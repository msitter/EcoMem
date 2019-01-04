#####################################################################
#####################################################################
# Converts MCMC samples from ecomem functions to MCMC list for use
# with coda
#####################################################################
#####################################################################

mem2mcmc = function(x,...){

  params = names(x$post.samps[[1]])

  mcmc.samps = coda::mcmc.list(lapply(x$post.samps,function(y){
    samps = y$beta
    param.names = x$pred.vars
    if ("sig.y" %in% params){
      samps = cbind(samps,y$sig.y)
      param.names = c(param.names,"sig.y")
    }
    if ("w" %in% params){
      for (i in 1:length(x$mem.vars)){
        samps = cbind(samps,y$w[[i]])
        param.names = c(param.names,paste(
          x$mem.vars[i],":w",1:ncol(y$w[[i]])-1,sep=""))
      }
    }
    if ("tau.sq" %in% params){
      samps = cbind(samps,y$tau.sq)
      param.names = c(param.names,paste(
        x$mem.vars,":tau.sq",sep=""))
    }
    colnames(samps) = param.names
    return(coda::mcmc(samps))
  }))

  return(mcmc.samps)

}
