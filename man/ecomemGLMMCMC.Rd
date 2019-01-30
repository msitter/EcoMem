\name{ecomemGLMMCMC}
\alias{ecomemGLMMCMC}
\title{Markov chain Monte Carlo algorithm to estimate ecological memory within a univariate
generalized linear model}

\description{Implements Markov chain Monte Carlo (MCMC) simulation to fit a Bayesian hierarchical
generalized linear regression model estimating ecological memory.}

\usage{
ecomemGLMMCMC(x)
}

\arguments{
  \item{x}{an \code{ecomem} class list containing inputs, prior parameters, and starting
  values for a univariate generalized linear regression model to estimate ecological
  memory.}
}

\value{
  A named list containing posterior samples for model parameters including \code{beta}.
  If memory variables are specified, the list includes samples for \code{eta},
  \code{w}, \code{X}, and \code{tau.sq} (only if \code{smooth = NULL}). \code{eta} and \code{w}
  are named lists with basis function coefficients and weights for each memory variable
  specified. The \eqn{i^{th}} row of \code{X} contains covariate values weighted based
  on the \eqn{i^{th}} posterior sample of ecological memory function weights formed by converting
  the weighted design matrix to a vector by rows.
}

\author{
  Malcolm S. Itter \email{malcolm.itter@helsinki.fi}
}

\seealso{
  \code{\link{ecomemGLM}}
}


\examples{
\dontrun{

## Simulate some time series data

####################################
#### Binomial response #############
####################################

## Define logit backtransformation
inv.logit = function(x){
  1/(1+exp(-x))
}

## Simulate some time series data

set.seed(1)

# Maximum lag
L = 8

# Sample size
n = 200

# Total number of time points
TT = n + L

# Generate weights
w = exp(-0.7*(0:L))/sum(exp(-0.7*(0:L)))

# Simulate covariate values
x = scale(rnorm(TT,(1:TT)/5000+0.02*cos(2*pi*(1:TT)/TT),2.4),center=FALSE)

# Form matrix containing lagged covariate values
x.lag = t(sapply(1:n,function(i){
  x[(L+1):1 + i - 1]
}))

# Calculate weighted covariate values
x.tilde = x.lag\%*\%w

# Form design matrix
X = cbind(rep(1,n),x.tilde)

# Set model parameter values
beta = c(-0.8,0.75)

# Simulate from ecological memory model
n.trials = sample(1:5,n,replace=T)
y = c(rep(NA,L),rbinom(n,n.trials,inv.logit(X\%*\%beta)))

# Form model data frame
data = data.frame(time=1:TT,x=x,y=y,n.trials=c(rep(NA,L),n.trials))

## Generate ecological memory model inputs

n.chains = 2

mod = ecomemGLM(y~x,family="binomial",data=data,mem.vars="x",
                L=L,timeID="time",offset="n.trials",
                n.chains=n.chains,inputs.only=TRUE)

## Run ecomemMCMC

mod.fit = lapply(1:n.chains, function(i){
  ecomemGLMMCMC(mod$inputs[[i]])
})

## Assess model convergence

mod$post.samps = mod.fit
plot(mem2mcmc(mod),ask=T)

## Plot memory function

plotmem(mod)

}
}
