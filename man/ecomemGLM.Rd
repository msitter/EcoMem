\name{ecomemGLM}
\alias{ecomemGLM}
\title{Function for quantifying ecological memory within a univariate generalized linear model}

\description{Estimates ecological memory functions for a specified subset of covariates and a discrete response within a generalized linear model framework. The function is a wrapper to
\code{\link{ecomemGLMMCMC}} which implements Markov chain Monte Carlo (MCMC) simulation to fit a Bayesian hierarchical model estimating ecological memory.}

\usage{
ecomem(formula, family = "binomial", data, mem.vars, var.type = NULL,
       L, timeID, groupID = NA, offset, starting = NULL,
       smooth = NULL, n.post = 1000, thin = 10, burn.in = 5000,
       n.step = 5, n.chains = 3, parallel = TRUE, max.cpu = NULL,
       inputs.only = FALSE, ...)
}

\arguments{
  \item{formula}{an object of class \code{\link[stats]{formula}} providing
  a symbolic representation of the linear regression model to be fit.}
  \item{family}{exponential family of the response. Currently supports only
  \code{binomial} (default) and \code{poisson}, respectively applying logit and
  log link functions.}
  \item{data}{a data frame including the response, explanatory, group, and time index
  variables.}
  \item{mem.vars}{an optional character vector indicating the subset of explanatory
  variables for which ecological memory functions should be estimated (memory
  variables).}
  \item{var.type}{an optional character vector specifying the structure of memory
  variables. If non-\code{NULL}, values should be "\code{C}" for continuous or
  "\code{D}" for discrete (i.e., binary). The default setting is all memory
  variables are continuous.}
  \item{L}{an integer or an integer vector specifying the maximum lag for memory
  variables. If a vector is specified, its dimension must match the number of memory
  variables specified in \code{mem.vars}.}
  \item{timeID}{a character value indicating the time index variable within
  \code{data}.}
  \item{groupID}{an optional character value indicating a group index variable
  within \code{data}. Only one \code{groupID} variable is allowed.}
  \item{offset}{an optional character value indicating an offset variable within
  \code{data}. For \code{binomial} data, the offset indicates the number of trials (or
  total size for proportional data).}
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

If \code{var.type = "D"}, the explanatory variable must consist of binary (0/1)
integer values with 1 indicating the occurrence of an event (e.g.,
disturbance/stress event) and 0 indicating no event occurrence.
}
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
\item{scale.factors}{a named numeric vector containing the factors used to scale
explanatory variables. Names correspond to model variables.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{

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
n = 2000

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

## Run ecomem function

mod = ecomemGLM(y~x,family="binomial",data=data,mem.vars="x",L=L,
                timeID="time",offset="n.trials")

## Assess model convergence

post.samps = mem2mcmc(mod)
plot(post.samps,ask=T)

## Process results

coef.summ = memsum(mod)

p = plotmem(mod)
print(p)

####################################
#### Poisson response ##############
####################################

## Simulate some time series data

set.seed(1)

# Maximum lag
L = 8

# Sample size
n = 2000

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
beta = c(-1.3,0.95)

# Simulate from ecological memory model
effort = runif(n,0,10)
y = c(rep(NA,L),rpois(n,exp(X\%*\%beta)*effort))

# Form model data frame
data = data.frame(time=1:TT,x=x,y=y,effort=c(rep(NA,L),effort))

## Run ecomem function

mod = ecomemGLM(y~x,family="poisson",data=data,mem.vars="x",L=L,
                timeID="time",offset="effort",n.step=12,thin=15)

## Assess model convergence

post.samps = mem2mcmc(mod)
plot(post.samps,ask=T)

## Process results

coef.summ = memsum(mod)

p = plotmem(mod)
print(p)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
