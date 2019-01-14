#####################################################################
#####################################################################
# Plot memory functions
#####################################################################
#####################################################################

plotmem = function(x,cred.int=0.95,...){

  n.chains = length(x$post.samps)
  params = names(x$post.samps[[1]])
  quantiles = (1-cred.int)/2 + c(0,cred.int)

  summ = function(y){
    tab = as.data.frame(t(apply(y,2,function(z)c(mean(z),median(z),quantile(z,quantiles)))))
    names(tab) =
      c("mean","median","lwr","upr")
    return(tab)
  }

  # summarize memory functions
  w.post = do.call("rbind",lapply(1:length(x$mem.vars),function(i){
    cbind(data.frame(var=x$mem.vars[i],
                     lag=1:ncol(x$post.samps[[1]]$w[[i]])-1,
          stringsAsFactors=FALSE),
          summ(do.call("rbind",lapply(x$post.samps,function(y)y$w[[i]]))))
  }))

  # plot memory functions
  p = ggplot2::ggplot(aes(x=lag,y=mean,ymin=lwr,ymax=upr),data=w.post) +
    geom_ribbon(aes(fill=paste(100*cred.int,"% Cred. Int.",sep=""))) +
    scale_fill_manual(values=c("lightgrey")) +
    # geom_line(aes(color="Post. mean"),size=0.3) +
    geom_point(aes(color="Post. mean")) +
    scale_color_manual(values="black") +
    xlab("Lag") + ylab("Weight") +
    facet_wrap(~var,scale="free") +
    theme_bw() +
    theme(legend.title=element_blank(),
          legend.spacing.y=unit(0.1,"pt"))

  print(p)

}
