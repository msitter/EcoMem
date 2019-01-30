\name{memsum}
\alias{memsum}
\title{Summarize ecological memory function parameters}

\description{Summarizes marginal posterior distribution of ecological
  memory model parameters based on \code{ecomem} object.}

\usage{
  memsum(x, cred.int = 0.95, verbose = TRUE, ...)
}

\arguments{
  \item{x}{an object of \code{ecomem} providing posterior
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

\author{
  Malcolm S. Itter \email{malcolm.itter@helsinki.fi}
}

\seealso{
  \code{\link{ecomem}}, \code{\link{ecomemGLM}}
}
