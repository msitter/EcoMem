# Ecological Memory Function
# Gaussian Process Weights
# April 13, 2018

ecomemMCMC = function(x){

  #### Define general functions ####

  # logit transform
  logit = function(x,a,b){
    log((x-a)/(b-x))
  }

  # logit backtransformation
  logit.b = function(y,a,b){
    b - ((b-a)/(1+exp(y)))
  }

  # log-determinant of Jacobian matrix
  # for logit transform
  logit.J = function(x,a,b){
    log(x-a) + log(b-x)
  }

  # MVN function
  Rmvn = function(mu,S){
    x = rnorm(length(mu))
    y = mu + crossprod(chol(S),x)
    return(as.numeric(y))
  }

  # # weight discrete data
  # wtD = function(x,w){
  #   tmp = matrix(w[x+1],dim(x))
  #   tmp[is.na(tmp)] = 0
  #   wtd.val = apply(tmp,1,sum)
  #   return(wtd.val)
  # }

  #### Define log-likelihood functions ####

  # Log-likelihood data
  ll.y = function(Z){
    l = -(1/(sig.y^2))*crossprod(y-Z%*%beta)
    return(as.numeric(l))
  }

  # Log-target density tau
  ll.tau = function(v,k,z,S){
    l = -(k-1)*log(v) - (1/(2*(v^2)))*crossprod(z,S%*%z) -
      ((nu+1)/2)*log(1+(1/nu)*((v/A)^2))
    return(as.numeric(l))
  }

  # Log-target density sig.y
  ll.sig = function(v){
    l = -n*log(v) - (1/(2*(v^2)))*crossprod(y-X%*%beta) +
      logit.J(v,a.y,b.y)
    return(as.numeric(l))
  }

  #### Load MCMC inputs ####

  chain = x$chain

  # Load inputs
  for (i in 1:length(x$inputs)){
    tmp = x$inputs[[i]]
    assign(names(x$inputs)[i],tmp)
  }

  # Load priors
  for (i in 1:length(x$priors)){
    tmp = x$priors[[i]]
    assign(names(x$priors)[i],tmp)
  }

  # Load starting values
  for (i in 1:length(x$starting)){
    tmp = x$starting[[i]]
    assign(names(x$starting)[i],tmp)
  }

  #### Setup MCMC ####

  # MCMC tracking/adaptive parameters
  n.iter = n.post*thin + burn.in
  n.save = (n.iter-burn.in)/thin
  iter2save = seq(burn.in+1,n.iter,thin)
  n.block = floor(0.05*n.iter)
  track = seq(n.block,n.iter,n.block)
  n.batch = 50
  adapt.step = seq(n.batch,n.iter,n.batch)

  # Define data arrays
  beta.sim = array(NA,dim=c(n.save,p))
  sig.y.sim = rep(NA,n.save)
  tau.track = array(0,dim=c(n.iter,p.mem))
  sig.y.track = rep(0,n.iter)
  if (p.mem > 0){
    eta.sim = wts.sim = wts.track = list()
    for (i in 1:p.mem){
      eta.sim[[i]] = array(NA,dim=c(n.save,bf[[i]]$k))
      wts.sim[[i]] = array(NA,dim=c(n.save,L[i]+1))
    }
    names(eta.sim) = names(wts.sim) = mem.vars
    X.sim = array(NA,dim=c(n.save,n*p))
    tau.sq.sim = array(NA,dim=c(n.save,p.mem))
    theta.sim = array(NA,dim=c(n.save,p.mem))
  }

  theta = rep(NA,p.mem)

  #### Run MCMC sampler ####

  for (iter in 1:n.iter){

    ###########################
    #### Update parameters ####
    ###########################

    #### Update memory functions ####

    if (p.mem > 0){

      #### Update antecedent weights ####

      # Loop over p.mem
      for (i in 1:p.mem){
        eta.curr = mem[[i]]$eta
        w.curr = mem[[i]]$w
        X.curr = X
        for (j in 1:n.step){
          ll.thresh = ll.y(X.curr) + log(runif(1))
          xi = backsolve((1/tau[i])*bf[[i]]$U,rnorm(bf[[i]]$k))
          theta[i] = runif(1)*(2*pi)
          I = c(theta[i]-(2*pi),theta[i])
          eta.star = (cos(theta[i]))*eta.curr + (sin(theta[i]))*xi
          tmp = exp(bf[[i]]$H%*%eta.star)
          w.star = as.numeric(tmp/sum(tmp))
          X.star = X
          # if (var.type[i]=="C"){
            X.star[,mem.vars[i]] = x.lag[[i]]%*%w.star
          # } else {
          #   X.star[,mem.vars[i]] = wtD(x.lag[[i]],w.star)
          # }
          if (inter==TRUE){
            for (l in 1:length(inter.vars)){
              X.star[,inter.terms[l]] = apply(X.star[,inter.vars[[l]]],1,prod)
            }
          }
          if (any(is.nan(w.star))){
            ll.star = "skip"
          } else {
            ll.star = ll.y(X.star)
          }
          if (ll.star!="skip" & ll.star > ll.thresh){
            eta.curr = eta.star
            w.curr = w.star
            X.curr = X.star
          } else {
            if (theta[i] < 0){
              I[1] = theta[i]
            } else if (theta[i] > 0){
              I[2] = theta[i]
            } else {
              warning("antecedent weights not updated",immediate.=TRUE)
              break()
            }
            val = 1
            while (val>0){
              theta[i] = runif(1,I[1],I[2])
              eta.star = (cos(theta[i]))*eta.curr + (sin(theta[i]))*xi
              tmp = exp(bf[[i]]$H%*%eta.star)
              w.star = as.numeric(tmp/sum(tmp))
              X.star = X
              # if (var.type[i]=="C"){
                X.star[,mem.vars[i]] = x.lag[[i]]%*%w.star
              # } else {
              #   X.star[,mem.vars[i]] = wtD(x.lag[[i]],w.star)
              # }
              if (inter==TRUE){
                for (l in 1:length(inter.vars)){
                  X.star[,inter.terms[l]] = apply(X.star[,inter.vars[[l]]],1,prod)
                }
              }
              if (any(is.nan(w.star))){
                ll.star = "skip"
              } else {
                ll.star = ll.y(X.star)
              }
              if (ll.star!="skip" & ll.star > ll.thresh){
                eta.curr = eta.star
                w.curr = w.star
                X.curr = X.star
                val = 0
              } else if (theta[i] < 0){
                I[1] = theta[i]
              } else if (theta[i] > 0){
                I[2] = theta[i]
              } else {
                warning("antecedent weights not updated",immediate.=TRUE)
                break()
              }
            }
          }
        }
        mem[[i]]$eta = eta.curr
        mem[[i]]$w = w.curr
        X = X.curr
      }

      #### Update smoothing parameters ####

      if (update.smooth==TRUE){
        for (i in 1:p.mem){
          tau.star = exp(rnorm(1,log(tau[i]),tau.tune[i]))
          r = exp(ll.tau(tau.star,bf[[i]]$k,mem[[i]]$eta,bf[[i]]$S)-
                    ll.tau(tau[i],bf[[i]]$k,mem[[i]]$eta,bf[[i]]$S))
          if (r > runif(1)){
            tau[i] = tau.star
            tau.track[iter,i] = 1
          }
        }
      }

    }

    #### Update regression coefficients ####

    S = chol2inv(chol(crossprod(X)/(sig.y^2) + diag(p)/sig2.0))
    s = crossprod(X,y)/(sig.y^2)
    beta = Rmvn(S%*%s,S)

    #### Update residual std. dev. ####

    sig.y.star = logit.b(rnorm(1,logit(sig.y,a.y,b.y),sig.tune),a.y,b.y)
    r = exp(ll.sig(sig.y.star)-ll.sig(sig.y))
    if (r > runif(1)){
      sig.y = sig.y.star
      sig.y.track[iter] = 1
    }

    #########################
    #### Adaptation Step ####
    #########################

    if (iter %in% adapt.step){
      delta.n = min(0.1,1/sqrt(which(adapt.step %in% iter)))

      #### Update tuning variance for tau ####
      if (update.smooth==TRUE){
        tau.rate = apply(as.matrix(tau.track[(iter-(n.batch-1)):iter,]),2,mean)
        for (i in 1:p.mem){
          if (tau.rate[i]>0.44){
            tau.tune[i] = exp(log(tau.tune[i])+delta.n)
          } else if (tau.rate[i]<0.44){
            tau.tune[i] = exp(log(tau.tune[i])-delta.n)
          }
        }
      }

      #### Update tuning variance for sig.y ####
      sig.y.rate = mean(sig.y.track[(iter-(n.batch-1)):iter])
      if (sig.y.rate>0.44){
        sig.tune = exp(log(sig.tune) + delta.n)
      } else if (sig.y.rate<0.44){
        sig.tune = exp(log(sig.tune) - delta.n)
      }

    }

    ######################
    #### Save & Track ####
    ######################

    #### Save Samples ####

    if (iter %in% iter2save){
      idx = which(iter2save %in% iter)
      beta.sim[idx,] = beta
      sig.y.sim[idx] = sig.y
      if (p.mem > 0){
        for (i in 1:p.mem){
          eta.sim[[i]][idx,] = mem[[i]]$eta
          wts.sim[[i]][idx,] = mem[[i]]$w
        }
        if (update.smooth==TRUE){
          tau.sq.sim[idx,] = tau^2
        }
        X.sim[idx,] = as.numeric(t(X))
      }
    }

    #### Track progress ####

    if (iter %in% track){
      cat("\n","Chain",chain,"of",n.chains,
          "\n",paste(iter/n.iter*100,"%",sep=""),"complete","\n")
    }

  } # end MCMC loop

  if (p.mem > 0){
    if (update.smooth==TRUE){
      return(list(beta=beta.sim,sig.y=sig.y.sim,eta=eta.sim,w=wts.sim,
                  X=X.sim,tau.sq=tau.sq.sim))
    } else {
      return(list(beta=beta.sim,sig.y=sig.y.sim,eta=eta.sim,w=wts.sim,
                  X=X.sim))
    }
  } else {
    return(list(beta=beta.sim,sig.y=sig.y.sim))
  }

} # end function
