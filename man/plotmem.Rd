\name{plotmem}
\alias{plotmem}
\title{Plot ecological memory functions}

\description{Plots ecological memory function weights for all memory
variables included in ecological memory model.}

\usage{
  plotmem(x, cred.int = 0.95, ...)
}

\arguments{
  \item{x}{an object of \code{ecomem} providing posterior
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

\author{
  Malcolm S. Itter \email{malcolm.itter@helsinki.fi}
}

\seealso{
  \code{\link{ecomem}}, \code{\link{ecomemGLM}}
}
