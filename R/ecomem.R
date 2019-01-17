#####################################################################
#####################################################################
# Ecological memory function - ecomem(...)
#####################################################################
#####################################################################

ecomem = function(formula,data,mem.vars,
                  L,timeID,groupID=NA,starting=NULL,
                  smooth=NULL,n.post=1000,thin=10,burn.in=5000,
                  n.step=5,n.chains=3,parallel=TRUE,max.cpu=NULL,
                  inputs.only=FALSE,...){

  options(stringsAsFactors=FALSE)
  snow::setDefaultClusterOptions(type="SOCK")

  ######################################################################
  #### Check for unused arguments ######################################
  ######################################################################

  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
    if(! i %in% formal.args)
      warning("'",i, "' is not an argument")
  }

  ######################################################################

  ######################################################################
  #### Parse formula inputs ############################################
  ######################################################################

  if(missing(formula)){stop("model formula must be specified")}

  if(class(formula) == "formula"){
    resp = as.character(formula[[2]])
    main = attr(terms(formula),"term.labels")[attr(terms(formula),"order")==1]
    mod.order = max(attr(terms(formula),"order"))
    # Check model order
    if (mod.order>2){stop("only 2nd order interactions allowed")}
    # Check for interaction terms
    if (mod.order==2){
      inter.terms = attr(terms(formula),"term.labels")[attr(terms(formula),"order")==2]
      inter.vars = strsplit(inter.terms,":")
      inter = TRUE
    } else {
      inter = FALSE
      inter.terms = NULL
      inter.vars = NULL
    }
  } else {
    stop("model formula is misspecified")
  }

  ######################################################################

  ######################################################################
  #### Define lagged covariates ########################################
  ######################################################################

  # Check inputs
  if(missing(data)){stop("no model data frame specified")}
  if(missing(mem.vars)){
    mem.vars = NULL
    warning("no memory variables specified")
  } else {
    if (!isTRUE(is.character(mem.vars))){stop("mem.vars must be a character")}
  }
  if(missing(timeID)){stop("time variable must be specified")}
  if (!isTRUE(is.character(timeID))){stop("timeID must be a character")}
  if(!isTRUE(all(data[,timeID]==floor(data[,timeID])))){stop("time variable must only contain integer values")}
  if(!is.na(groupID)){
    if (!isTRUE(is.character(groupID))){stop("groupID must be a character")}
  }

  p.mem = length(mem.vars)

  # Check L input
  if(missing(L)){stop("L must be specified")}
  check.L = length(L)
  if (check.L==1){
    L = rep(L,p.mem)
  } else {
    if (check.L!=p.mem){stop("L is missing for one or more memory covariates")}
  }
  names(L) = mem.vars

  # Check if fixed values are provided for smoothing parameters
  if (!is.null(smooth)){
    if (length(smooth)==1){
      smooth = rep(smooth,p.mem)
    } else {
      if (length(smooth)!=p.mem){
        stop("smoothing parameter is missing for one or more memory covariates")
      }
    }
    names(smooth) = mem.vars
  }

  # Define groups
  if (is.na(groupID)){
    data$groupID = rep(1,nrow(data))
    groupID = "groupID"
  }
  group = data[,groupID]
  group.idx = sort(unique(group))
  n.group = length(group.idx)

  # # Check var.type
  # if (!is.null(var.type)){
  #   if (!all(var.type%in%c("C","D"))){stop("one or more memory covariates are mis-categorized")}
  #   check.var.type = length(var.type)
  #   if (check.var.type==1){
  #     var.type = rep(var.type,p.mem)
  #   } else {
  #     if (check.var.type!=p.mem){stop("var.type is missing for one or more memory covariates")}
  #   }
  # } else {
  #   var.type = rep("C",p.mem)
  # }

  # Assign variable types to covariates
  aux.vars = main[!main%in%mem.vars]
  aux.vars.C = aux.vars[which(apply(as.matrix(data[,aux.vars]),
                                    2,class)%in%c("numeric","integer"))]
  aux.vars.D = aux.vars[!aux.vars%in%aux.vars.C]
  if (!all(apply(as.matrix(data[,mem.vars]),2,class)%in%c("numeric","integer"))){
    stop("memory variables must be numeric or integer")
  }
  mem.vars.D = mem.vars[which(apply(as.matrix(data[,mem.vars]),2,function(x)
    all(na.omit(x)%in%c(0,1))))]
  mem.vars.C = mem.vars[!mem.vars%in%mem.vars.D]

  # aux.vars = main[!main%in%mem.vars]
  # aux.vars.C = aux.vars[which(apply(as.matrix(data[,aux.vars]),
  #                                   2,class)%in%c("numeric","integer"))]
  # aux.vars.D = aux.vars[!aux.vars%in%aux.vars.C]
  # mem.vars.C = mem.vars[var.type=="C"]
  # nC = length(mem.vars.C)
  # mem.vars.D = mem.vars[var.type=="D"]
  # nD = length(mem.vars.D)

  # Determine interaction terms including memory variables
  if (inter==TRUE){
    inter.terms = unlist(lapply(inter.vars,function(x){
      if (!any(x%in%mem.vars)){
        d = NULL
      } else {
        if (any(x%in%aux.vars.D)){
          mem.var.idx = which(mem.vars%in%x)
          aux.var.idx = which(aux.vars.D%in%x)
          var = aux.vars.D[aux.var.idx]
          aux.var.pos = which(x==var)
          if (!is.factor(data[,var])){
            tmp = factor(data[,var])
          } else {
            tmp = data[,var]
          }
          if (aux.var.pos==1){
            d = paste(paste(var,levels(tmp)[-1],sep=""),mem.vars[mem.var.idx],sep=":")
          } else {
            d = paste(mem.vars[mem.var.idx],paste(var,levels(tmp)[-1],sep=""),sep=":")
          }
        } else {
          d = paste(x,collapse=":")
        }
      }
      return(d)
    }))

    inter.vars = lapply(1:length(inter.terms),function(i){
      unlist(strsplit(inter.terms[i],":"))
    })

  }

  # # Check discrete vars and calculate counts
  # if (nD > 0){
  #   if (!all(data[,mem.vars.D]%in%c(NA,0,1))){stop("Non-binary discrete memory covariates specified")}
  #   D = array(NA,dim=c(n.group,nD))
  #   for (i in 1:n.group){
  #     for (j in 1:nD){
  #       tmp = data[data[,groupID]==group.idx[i],mem.vars.D[j]]
  #       D[i,j] = sum(tmp,na.rm=TRUE)
  #     }
  #   }
  #   colnames(D) = mem.vars.D
  # }

  data = data[order(data[,groupID],data[,timeID]),]
  scaled.X = scale(data[,c(mem.vars.C,aux.vars.C)])
  # scale.factors = attr(scaled.X,"scaled:scale")
  if (dim(scaled.X)[2]==1){
    data[,c(mem.vars.C,aux.vars.C)] = as.numeric(scaled.X)
    # names(scale.factors) = c(mem.vars.C,aux.vars.C)
  } else {
    data[,c(mem.vars.C,aux.vars.C)] = scaled.X
  }
  mod.data = data

  # # Define function to calculate time since disturbance
  # tsD = function(t,v,max,n.col){
  #   tmp = outer(t,t[v==1&!is.na(v)],"-")
  #   tmp[tmp<0|tmp>max] = 9999
  #   tmp = cbind(tmp,matrix(9999,nrow(tmp),n.col-ncol(tmp)))
  #   return(tmp)
  # }

  # # Define function to weight discrete data
  # wtD = function(x,w){
  #   tmp = matrix(w[x+1],dim(x))
  #   tmp[is.na(tmp)] = 0
  #   wtd.val = apply(tmp,1,sum)
  #   return(wtd.val)
  # }

  #################################################
  #### Fix time lag NA for discrete covariates ####
  #################################################
  x.mem.all.obs = lapply(1:length(L),function(i){
    do.call("rbind",lapply(1:n.group,function(j){
      tmp.dat = data[data[,groupID]==group.idx[j],]
      # if (var.type[i]=="C"){
        x.lag.mat = t(sapply(1:nrow(tmp.dat),function(k){
          if (!is.na(tmp.dat[k,resp])){
            v = (0:(-L[i])) + tmp.dat[k,timeID]
            if (all(v%in%tmp.dat[,timeID])){
              x.vals = tmp.dat[match(v,tmp.dat[,timeID]),which(names(tmp.dat)==mem.vars[i])]
            } else {
              x.vals = rep(NA,L[i]+1)
            }
          } else {
            x.vals = rep(NA,L[i]+1)
          }
          return(x.vals)
        }))
        storage.mode(x.lag.mat) = "double"
      # } else {
      #   x.lag.mat = tsD(tmp.dat[,timeID],tmp.dat[,mem.vars[i]],L[i],
      #                   max(D[,mem.vars[i]],na.rm=TRUE))
      #   for (k in 1:nrow(tmp.dat)){
      #     if (!is.na(tmp.dat[k,resp])){
      #       v = (0:(-L[i])) + tmp.dat[k,timeID]
      #       if (!all(v%in%tmp.dat[,timeID])){
      #         x.lag.mat[k,] = rep(NA,ncol(x.lag.mat))
      #       }
      #     } else {
      #       x.lag.mat[k,] = rep(NA,ncol(x.lag.mat))
      #     }
      #   }
      #   storage.mode(x.lag.mat) = "integer"
      # }
      return(x.lag.mat)
    }))
  })

  drop.idx = sort(unique(
    c(unlist(lapply(x.mem.all.obs,function(x){
      which(apply(x,1,function(y)any(is.na(y)))==T)
    })),which(apply(as.matrix(data[,aux.vars]),1,
                    function(z)any(is.na(z)))))))

  mod.data[drop.idx,resp] = NA
  data = data[-drop.idx,]
  n = nrow(data)
  x.mem = lapply(x.mem.all.obs,function(x){
    x[-drop.idx,]
  })
  names(x.mem) = mem.vars
  group = group[-drop.idx]
  if (any(as.numeric(table(group))==0)){
    warning("no data for one or more groups")
  }
  timeframe = min(data[,timeID]):max(data[,timeID])

  ######################################################################

  ######################################################################
  #### Create model inputs #############################################
  ######################################################################

  #### Create data inputs ##############################################

  ### Form design matrix ###
  X = model.matrix(formula,data)
  p = ncol(X)
  pred.vars = colnames(X)
  storage.mode(X) = "double"
  ### Memory function inputs ###
  # Define basis functions
  bf = list()
  for (j in 1:length(L)){
    t.s = (0:L[j])/L[j]
    time = data.frame(t=0:L[j],t.s=t.s)
    n.knots = L[j] + 1
    CRbasis = mgcv::smoothCon(mgcv::s(t.s,k=n.knots,bs="cr"),data=time,knots=NULL,absorb.cons=TRUE,
                              scale.penalty=TRUE)
    RE = diag(ncol(CRbasis[[1]]$S[[1]]))
    bf[[j]] = list(S=CRbasis[[1]]$S[[1]]+(1E-07)*RE,
                   H=CRbasis[[1]]$X,
                   k=ncol(CRbasis[[1]]$X),
                   U=chol(CRbasis[[1]]$S[[1]]+(1E-07)*RE))
  }
  # Define smoothing parameters (if provided)
  if(!is.null(smooth)){
    inputs = list(y=as.double(data[,resp]),X.fix=X,n=as.integer(n),
                  p=as.integer(p),mem.vars=as.character(mem.vars),
                  # var.type=as.character(var.type),
                  p.mem=as.integer(p.mem),L=as.integer(L),x.lag=x.mem,
                  bf=bf,tau.sq=as.double(smooth),update.smooth=FALSE,
                  inter=inter,inter.terms=inter.terms,inter.vars=inter.vars,
                  n.post=as.integer(n.post),thin=as.integer(thin),
                  burn.in=as.integer(burn.in),n.chains=as.integer(n.chains),
                  n.step=as.integer(n.step))
  } else {
    inputs = list(y=as.double(data[,resp]),X.fix=X,n=as.integer(n),
                  p=as.integer(p),mem.vars=as.character(mem.vars),
                  # var.type=as.character(var.type),
                  p.mem=as.integer(p.mem),L=as.integer(L),x.lag=x.mem,
                  bf=bf,update.smooth=TRUE,inter=inter,inter.terms=inter.terms,
                  inter.vars=inter.vars,n.post=as.integer(n.post),thin=as.integer(thin),
                  burn.in=as.integer(burn.in),n.chains=as.integer(n.chains),
                  n.step=as.integer(n.step))
  }

  ######################################################################

  ######################################################################
  #### Define priors ###################################################
  ######################################################################

  nu = 4
  A = 0.1
  a.y = 1e-07
  b.y = 1e07
  sig2.0 = 1e07

  priors = list(nu=nu,A=A,a.y=a.y,b.y=b.y,sig2.0=sig2.0)

  ######################################################################

  ######################################################################
  #### Generate starting values ########################################
  ######################################################################

  # Initial starting value check
  if(!is.null(starting) & length(starting)!=n.chains){
    starting = NA
    warning("Starting values not provided for each MCMC chain,\nreverting to default starting values")
  }

  if (!is.null(starting)){
    # Check starting value inputs
    if(!is.list(starting)){stop("starting values must be specified as a list")}
    for (i in 1:n.chains){
      # beta
      if (! "beta" %in% names(starting[[i]])){
        stop(paste("starting values not provided for beta in chain",i,sep=" "))
      }
      if (length(starting[[i]]$beta)!=inputs[[1]]$p){stop("beta starting values are wrong length")}
      storage.mode(starting[[i]]$beta) = "double"
      # tau.sq
      if (is.null(smooth)){
        if (! "tau.sq" %in% names(starting[[i]])){
          stop(paste("starting values not provided for tau.sq in chain",i,sep=" "))
        }
        if (length(starting[[i]]$tau.sq)!=inputs[[1]]$p.mem){
          stop("tau.sq starting values are wrong length")
        }
      }
      # sig.y
      if (! "sig.y" %in% names(starting[[i]])){
        stop(paste("starting values not provided for sig.y in chain",i,sep=" "))
      }
      if (length(starting[[i]]$sig.y)!=1){stop("sig.y starting value must be a scalar")}
      storage.mode(starting[[i]]$sig.y) = "double"
      # mem
      if (! "mem" %in% names(starting[[i]])){
        stop(paste("starting values not provided for memory functions in chain",i,sep=" "))
      }
      if (!is.list(starting[[i]]$mem)){stop("memory function starting values must be a list")}
      if (length(starting[[i]]$mem)!=inputs[[1]]$p.mem){
        stop("starting values missing for one or more memory variables")
      }
      names(starting[[i]]$mem) = inputs[[1]]$mem.vars
      for (j in 1:inputs[[1]]$p.mem){
        if (! "eta" %in% names(starting[[i]]$mem[[j]])){
          stop(paste("eta starting values missing for",inputs[[1]]$mem.vars[j],
                     "in chain",i,sep=" "))
        }
        if (length(starting[[i]]$mem[[j]]$eta)!=inputs[[1]]$bf[[j]]$k){
          stop(paste("eta starting values for",inputs[[1]]$mem.vars[j],
                     "in chain",i,"incorrect length",sep=" "))
        }
        storage.mode(starting[[i]]$mem[[j]]$eta) = "double"
        if (! "w" %in% names(starting[[i]]$mem[[j]])){
          stop(paste("weight starting values missing for",inputs[[1]]$mem.vars[j],
                     "in chain",i,sep=" "))
        }
        if (length(starting[[i]]$mem[[j]]$w)!=(inputs[[1]]$L[j]+1)){
          stop(paste("weight starting values for",inputs[[1]]$mem.vars[j],
                     "in chain",i,"incorrect length",sep=" "))
        }
        storage.mode(starting[[i]]$mem[[j]]$w) = "double"
      }
      if (is.null(smooth)){
        starting[[i]]$tau = sqrt(starting[[i]]$tau.sq)
        starting[[i]]$tau.sq = NULL
        starting[[i]]$tau.tune = rep(0.1,inputs[[1]]$p.mem)
      }
      X=inputs[[1]]$X.fix
      X[,inputs[[1]]$mem.vars] = sapply(1:inputs[[1]]$p.mem,function(j){
        # if (inputs[[1]]$var.type[j]=="C"){
          inputs[[1]]$x.lag[[j]]%*%starting[[i]]$mem[[j]]$w
        # } else {
        #   wtD(inputs[[1]]$x.lag[[j]],starting[[i]]$mem[[j]]$w)
        # }
      })
      if (inter==TRUE){
        X[,inter.terms] = sapply(1:length(inter.vars),function(j){
          apply(X[,inter.vars[[j]]],1,prod)
        })
      }
      starting[[i]]$X = X
      starting[[i]]$sig.tune = 0.1
    }
  } else {
    X = inputs$X.fix
    start.mod = lm(inputs$y~X-1)
    starting = lapply(1:n.chains,function(i){
      beta = rnorm(inputs$p,coef(start.mod),0.1)
      if(is.null(smooth)){
        tau = 10^(-sample(2:4,inputs$p.mem,replace=T)/2)
        tau.tune = rep(0.1,inputs$p.mem)
      }
      sig.y = exp(rnorm(1,log(sd(start.mod$residuals)),0.1))
      sig.tune = 0.1
      X.start = X
      if (inputs$p.mem > 0){
        mem = lapply(1:inputs$p.mem,function(j){
          uni.wts = rep(1/(inputs$L[j]+1),inputs$L[j]+1)
          eta = rnorm(inputs$bf[[j]]$k,coef(lm(uni.wts~inputs$bf[[j]]$H-1)),0.1)
          tmp = exp(inputs$bf[[j]]$H%*%eta)
          w = tmp/sum(tmp)
          return(list(eta=eta,w=w))
        })
        names(mem) = inputs$mem.vars
        X.start[,inputs$mem.vars] = sapply(1:inputs$p.mem,function(j){
          # if (inputs$var.type[j]=="C"){
            inputs$x.lag[[j]]%*%mem[[j]]$w
          # } else {
          #   wtD(inputs$x.lag[[j]],mem[[j]]$w)
          # }
        })
        if (inter==TRUE){
          X.start[,inter.terms] = sapply(1:length(inter.vars),function(j){
            apply(X.start[,inter.vars[[j]]],1,prod)
          })
        }
      } else {
        mem = NULL
      }
      if (is.null(smooth)){
        return(list(beta=beta,tau=tau,tau.tune=tau.tune,sig.y=sig.y,
                    sig.tune=sig.tune,mem=mem,X=X.start))
      } else {
        return(list(beta=beta,sig.y=sig.y,sig.tune=sig.tune,mem=mem,X=X.start))
      }
    })
  }

  ######################################################################

  ######################################################################
  #### Define MCMC inputs ##############################################
  ######################################################################

  mcmc.inputs = lapply(1:n.chains,function(i){
    list(inputs=inputs,priors=priors,starting=starting[[i]],chain=i)
  })
  names(mcmc.inputs) = paste("chain",1:n.chains,sep="")

  ######################################################################

  ######################################################################
  #### Summarize inputs ################################################
  ######################################################################

  if(!isTRUE(inputs.only)){
    cat("\n","---Starting MCMC---","\n",
        "Number of chains:",n.chains,"\n",
        "Running in parallel:",parallel,"\n",
        "Number of posterior samples:",n.post,"\n",
        "Thin interval:",thin,"\n",
        "Burn-in length:",burn.in,"\n",
        "Memory variables:",mem.vars,"\n",
        "Number of groups:",n.group,"\n",
        "Study period time points:",timeframe[1],
        "to",timeframe[length(timeframe)],"\n")
  } else {
    cat("\n","---Generating MCMC inputs---","\n",
        "Running in parallel:",parallel,"\n",
        "Number of posterior samples:",n.post,"\n",
        "Thin interval:",thin,"\n",
        "Burn-in length:",burn.in,"\n",
        "Memory variables:",mem.vars,"\n",
        "Number of groups:",n.group,"\n",
        "Study period time points:",timeframe[1],
        "to",timeframe[length(timeframe)],"\n")
  }

  ######################################################################

  ######################################################################
  #### Run MCMC ########################################################
  ######################################################################

  if (!isTRUE(inputs.only)){
    if (isTRUE(parallel)){
      cat("\n","track progress using 'track-ecomem.txt'","\n")
      if(!is.null(max.cpu)){
        snowfall::sfInit(parallel=T,cpus=max.cpu,slaveOutfile="track-ecomem.txt")
        snowfall::sfClusterSetupRNG()
        mod.out = snowfall::sfClusterApply(mcmc.inputs,ecomemMCMC)
        snowfall::sfStop()
        names(mod.out) = paste("chain",1:n.chains,sep="")
      } else {
        snowfall::sfInit(parallel=T,cpus=n.chains,slaveOutfile="track-ecomem.txt")
        snowfall::sfClusterSetupRNG()
        mod.out = snowfall::sfClusterApply(mcmc.inputs,ecomemMCMC)
        snowfall::sfStop()
        names(mod.out) = paste("chain",1:n.chains,sep="")
      }
    } else {
      mod.out = lapply(mcmc.inputs,function(x){
        ecomemMCMC(x)
      })
    }
  }

  ######################################################################

  ######################################################################
  #### Return output ###################################################
  ######################################################################

  if (isTRUE(inputs.only)){
    out = list(inputs=mcmc.inputs,data=mod.data,n=n,
               pred.vars=pred.vars,
               mem.vars=mem.vars)
  } else {
    if (n.chains>1){
      out = list(post.samps=mod.out,data=mod.data,n=n,
                 pred.vars=pred.vars,
                 mem.vars=mem.vars)
    } else {
      out = list(post.samps=mod.out[[1]],data=mod.data,n=n,
                 pred.vars=pred.vars,
                 mem.vars=mem.vars)
    }
  }

  class(out) = "ecomem"

  ######################################################################

  return(out)

} # End function
