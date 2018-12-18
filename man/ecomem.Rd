\name{ecomem}
\alias{ecomem}
\title{Ecological memory function}

\description{Estimates ecological memory functions for a specified subset of covariates based
on a continuous response using a linear model framework. The function is a wrapper to
\code{\link{ecomemMCMC}} which implements Markov chain Monte Carlo (MCMC) simulation to fit the
ecological memory model.}

\usage{
ecomem(formula, data, mem.vars, var.type = NULL, L, timeID, groupID = NA,
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
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
