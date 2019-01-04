\name{mem2mcmc}
\alias{mem2mcmc}
\title{Converts ecomem posterior samples to MCMC list}

\description{Converts ecomem posterior samples to named
\code{\link[coda]{mcmc.list}} object for use in the
\code{\link[coda]{coda}} package.}

\usage{
  memsum(x, ...)
}

\arguments{
  \item{x}{an object of \code{\link[ecomem]{ecomem}} providing posterior
  samples for all ecological memory model parameters.}
  \item{...}{currently no additional arguments.}
}

\value{
\code{\link[coda]{mcmc.list}} object containing all ecological memory
model parameters for use in the \code{\link[coda]{coda}} package.
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
