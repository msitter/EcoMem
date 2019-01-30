\name{mem2mcmc}
\alias{mem2mcmc}
\title{Converts ecomem posterior samples to MCMC list}

\description{Converts ecomem posterior samples to named
\code{\link[coda]{mcmc.list}} object for use in the
\code{\link[coda]{coda}} package.}

\usage{
  mem2mcmc(x, ...)
}

\arguments{
  \item{x}{an object of \code{ecomem} providing posterior
  samples for all ecological memory model parameters.}
  \item{...}{currently no additional arguments.}
}

\value{
  \code{\link[coda]{mcmc.list}} object containing all ecological memory
  model parameters for use in the \code{\link[coda]{coda}} package.
}

\author{
  Malcolm S. Itter \email{malcolm.itter@helsinki.fi}
}

\seealso{
  \code{\link{ecomem}}, \code{\link{ecomemGLM}}
}
