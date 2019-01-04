\name{plotmem}
\alias{plotmem}
\title{Plot ecological memory functions}

\description{Plots ecological memory function weights for all memory
variables included in ecological memory model.}

\usage{
  memsum(x, cred.int = 0.95, ...)
}

\arguments{
  \item{x}{an object of \code{\link[ecomem]{ecomem}} providing posterior
  samples for all ecological memory model parameters.}
  \item{cred.int}{numeric value between 0 and 1 indicating the credible
  interval to be displayed within plot.}
  \item{...}{currently no additional arguments.}
}

\value{
Summary plot of ecological memory function weights. The memory variable
name is displayed as the title to each subplot. Time lags correspond to
\code{L} passed to either \code{\link[ecomem]{ecomem}} or
\code{\link[ecomem]{ecomemGLM}}.
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
