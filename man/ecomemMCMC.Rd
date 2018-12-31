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
