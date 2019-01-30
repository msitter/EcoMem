\name{ecomemMCMC}
\alias{ecomemMCMC}
\title{Markov chain Monte Carlo algorithm to estimate ecological memory within a univariate
regression model}

\description{Implements Markov chain Monte Carlo (MCMC) simulation to fit a Bayesian hierarchical
linear regression model estimating ecological memory.}

\usage{
  ecomemMCMC(x)
}

\arguments{
  \item{x}{an \code{ecomem} class list containing inputs, prior parameters, and starting
  values for a univariate regression model to estimate ecological memory.}
}

\value{
  A named list containing posterior samples for each model parameter including: \code{beta} and
  \code{sig.y}. If memory variables are specified, the list includes samples for \code{eta},
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
  \code{\link{ecomem}}
}

\examples{
\dontrun{

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

## Generate ecological memory model inputs

n.chains = 2

mod = ecomem(y~x,data=data,mem.vars="x",L=L,timeID="time",
                n.chains=n.chains, inputs.only=TRUE)

## Run ecomemMCMC

mod.fit = lapply(1:n.chains, function(i){
  ecomemMCMC(mod$inputs[[i]])
})

## Assess model convergence

mod$post.samps = mod.fit
plot(mem2mcmc(mod),ask=T)

## Plot memory function

plotmem(mod)

}
}
