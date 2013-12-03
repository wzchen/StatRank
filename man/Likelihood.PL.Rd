\name{Likelihood.PL}
\alias{Likelihood.PL}
\title{A faster Likelihood for PL}
\usage{
  Likelihood.PL(Data, parameter)
}
\arguments{
  \item{Data}{ranking data}

  \item{parameter}{Mean of Exponential Distribution}
}
\value{
  log likelihood
}
\description{
  A faster Likelihood for PL
}
\examples{
data(Data.Test)
parameter = Generate.RUM.Parameters(5, "exponential")
Likelihood.PL(Data.Test, parameter)
}

