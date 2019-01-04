\name{memsum}
\alias{memsum}
\title{Summarize ecological memory function parameters}

\description{Summarizes marginal posterior distribution of ecological
  memory model parameters based on \code{ecomem} object.}

\usage{
  memsum(x, cred.int = 0.95, verbose = TRUE, ...)
}

\arguments{
  \item{x}{an object of \code{\link[ecomem]{ecomem}} providing posterior
  samples for all ecological memory model parameters.}
  \item{cred.int}{numeric value between 0 and 1 indicating the credible
  interval to return in the posterior summary.}
  \item{verbose}{logical indicating whether the marginal posterior summary
  should be printed to the console.}
  \item{...}{currently no additional arguments.}
}

\value{
A data frame including marginal posterior distribution summary statistics for
each ecological memory function model parameter. Parameters are listed within
the \code{var} variable. Summary statistics include the mean, median, and
lower and upper quantiles defining the specified credible interval.
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
