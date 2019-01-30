\name{ecomem}
\alias{ecomem}
\title{Function for quantifying ecological memory within a univariate regression model}

\description{Estimates ecological memory functions for a specified subset of covariates and a continuous response within a linear model framework. The function is a wrapper to
\code{\link{ecomemMCMC}} which implements Markov chain Monte Carlo (MCMC) simulation to fit a Bayesian hierarchical model estimating ecological memory.}

\usage{
ecomem(formula, data, mem.vars, L, timeID, groupID = NA,
       starting = NULL, smooth = NULL, n.post = 1000, thin = 10,
       burn.in = 5000, n.step = 5, n.chains = 3, parallel = TRUE,
       max.cpu = NULL, inputs.only = FALSE, ...)
}

\arguments{
  \item{formula}{an object of class \code{\link[stats]{formula}} providing
  a symbolic representation of the linear regression model to be fit.}
  \item{data}{a data frame including the response, explanatory, group, and time index
  variables.}
  \item{mem.vars}{an optional character vector indicating the subset of explanatory
  variables for which ecological memory functions should be estimated (memory
  variables).}
  \item{L}{an integer or an integer vector specifying the maximum lag for memory
  variables. If a vector is specified, its dimension must match the number of memory
  variables specified in \code{mem.vars}.}
  \item{timeID}{a character value indicating the time index variable within
  \code{data}.}
  \item{groupID}{an optional character value indicating a group index variable
  within \code{data}. Only one \code{groupID} variable is allowed.}
  \item{starting}{optional list of starting values for model parameters. If
  specified, the length of the list must match \code{n.chains}. Each element
  of the \code{starting} list must be a named list with names and dimensions
  corresponding to model parameters. Valid names include \code{beta}, \code{sig.y},
  \code{mem}, and \code{tau.sq} (if \code{smooth = NULL}). The \code{mem} element,
  must be a list with length equal to the number of memory variables with elements
  \code{eta} and \code{w} providing starting values for spline basis function
  coefficients and memory function weights, respectively.}
  \item{smooth}{an optional vector with length equal to the number of memory
  variables specifying penalty parameters for penalized spline regression used
  to estimate ecological memory functions. If included, \code{tau.sq} is fixed
  within the model.}
  \item{n.post}{an optional scalar indicating the number of desired posterior
  MCMC samples post \code{burn.in}. Total number of samples is given by
  \code{burn.in} + \code{thin}*\code{n.post}. The default is set to 1000.}
  \item{thin}{an optional scalar indicating the thinning interval for posterior
  samples. The default is set to 10.}
  \item{burn.in}{an optional scalar indicating the number of burn-in samples for
  MCMC chains. The default is set to 5000.}
  \item{n.step}{an optional scalar indicating the number of basis function
  coefficient updates taken per iteration of the MCMC sampler. The default is set
  to 5.}
  \item{n.chains}{an optional scalar indicating the number of MCMC chains. The
  default is set to 3.}
  \item{parallel}{an optional \code{logical} indicating whether MCMC chains
  should be run in parallel. If \code{TRUE} (default), \code{n.chains} are
  run in parallel using \code{\link{snowfall}}.}
  \item{max.cpu}{an optional scalar indicating the maximum number of CPUs to
  utilize. Only evaluated if \code{parallel = TRUE}. The default CPU number is
  equal to \code{n.chains}.}
  \item{inputs.only}{an optional \code{logical} indicating whether model should be
  fit using MCMC within the \code{ecomem} function. If \code{TRUE} inputs for
  \code{\link{ecomemMCMC}} are returned, but no posterior samples. The default is
  set to \code{FALSE}.}
  \item{...}{currently no additional arguments.}
}

\details{
  Unequally-spaced time points are not supported. The \code{timeID}
  variable should include equally-spaced, sequential time points within each group.
  Modeled time periods for different groups may be non-overlapping, but the spacing
  of time points \emph{must} be consistent among groups. Time points with \code{NA}
  values for the response and/or explanatory variables are dropped from the analysis.
  \code{NA} values are not supported in the \code{timeID} and \code{groupID}
  index variables.

  Explanatory variable observations for \code{L} time points prior to the
  response are needed to fit the regression model when estimating ecological
  memory functions for a subset of covariates. If explanatory variable observations
  are not available for \code{L} time points prior to the response, the observation
  is not evaluated as part of the regression model. This is particularly important
  during model initialization. If explanatory variable observations
  are available prior to the initial response (first time point), the first
  \code{\link[base]{max(L)}} responses for each group (if \code{groupID} is
  non-\code{NULL}) within \code{data} should be set to \code{NA}. Otherwise,
  the first \code{\link[base]{max(L)}} response observations within each
  group will not be evaluated as part of the regression model.

\value{
  An object of class \code{ecomem}, which is a list with the following elements:
  \item{post.samps}{a list of posterior samples with length equal to \code{n.chains}.
  Each element is a named list with names corresponding to model parameters
  (\code{beta}, \code{sig.y}, \code{eta}, \code{w}, \code{X.tilde}, \code{tau.sq}).
  Returned only if \code{inputs.only = FALSE} (default).}
  \item{inputs}{a list of MCMC inputs with length equal to \code{n.chains}. Each
  element is a named list containing \code{inputs}, \code{priors}, and
  \code{starting} to be passed to \code{\link{ecomemMCMC}}. Returned only if
  \code{inputs.only = TRUE}.}
  \item{data}{a data frame containing data used to fit the regression model including
  a response variable with the first \code{\link[base]{max(L)}} within each group set
  to \code{NA} and scaled values for all continuous explanatory variables.}
  \item{n}{an integer value indicating the number of reponse observations used to fit
  the regression model across all groups.}
}

\author{
  Malcolm S. Itter \email{malcolm.itter@helsinki.fi}
}

\seealso{
  \code{\link{ecomemGLM}}
}

\examples{
\dontrun{
####################################
#### Continuous covariate ##########
####################################

## Simulate some time series data

set.seed(1)

# Maximum lag
L = 10

# Sample size
n = 1000

# Total number of time points
TT = n + L

# Generate weights
w = exp(-0.5*(0:L))/sum(exp(-0.5*(0:L)))

# Simulate a single continuous covariate
x = scale(rnorm(TT,(1:TT)/1000+0.05*cos(
    2*pi*10*(1:TT)/TT),2.4),center=FALSE)

# Form matrix containing lagged covariate values
x.lag = t(sapply(1:n,function(i){
  x[(L+1):1 + i - 1]
}))

# Calculate weighted covariate values
x.tilde = x.lag\%*\%w

# Form design matrix
X = cbind(rep(1,n),x.tilde)

# Set model parameter values
beta = c(1.3,-0.4)
sig.y = 0.1

# Simulate from ecological memory model
y = c(rep(NA,L),rnorm(n,X\%*\%beta,sig.y))

# Form model data frame
data = data.frame(time=1:TT,x=x,y=y)

## Run ecomem function

mod = ecomem(y~x,data=data,mem.vars="x",L=L,timeID="time")

## Assess model convergence

post.samps = mem2mcmc(mod)
plot(post.samps,ask=T)

## Process results

coef.summ = memsum(mod)

p = plotmem(mod)
print(p)

####################################
#### Discrete covariate ############
####################################

## Define function to calculate time since events
tsD = function(t,v,max,n.col){
  tmp = outer(t,t[v==1],"-")
  tmp[tmp<0|tmp>max] = 9999
  tmp = cbind(tmp,matrix(9999,nrow(tmp),n.col-ncol(tmp)))
  return(tmp)
}

## Define function to weight event data
wtD = function(x,w){
  tmp = matrix(w[x+1],dim(x))
  tmp[is.na(tmp)] = 0
  wtd.val = apply(tmp,1,sum)
  return(wtd.val)
}

## Simulate some time series data

set.seed(1)

# Maximum lag
L = 5

# Sample size
n = 500

# Total number of time points
TT = n + L

# Generate weights
w = exp(-0.95*(0:L))/sum(exp(-0.95*(0:L)))

# Simulate a single binary covariate
x = rbinom(TT,1,0.05)

# Calculate event counts
nx = sum(x)

# Form matrix with time since event values
tsd.vals = tsD(1:TT,x,L,nx)[(L+1):TT,]

# Calculate weighted covariate values
x.tilde = wtD(tsd.vals,w)

# Form design matrix
X = cbind(rep(1,n),x.tilde)

# Define model parameters
beta = c(2.5,3.4)
sig.y = 0.2

# Simulate from ecological memory model
y = c(rep(NA,L),rnorm(n,X\%*\%beta,sig.y))

# Form model data frame
data = data.frame(time=1:TT,x=x,y=y)

## Run ecomem function

mod = ecomem(y~x,data=data,mem.vars="x",var.type="D",
      L=L,timeID="time")

## Assess model convergence

post.samps = mem2mcmc(mod)
plot(post.samps,ask=T)

## Process results

coef.summ = memsum(mod)

p = plotmem(mod)
print(p)
}
}
