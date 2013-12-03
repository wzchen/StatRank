\name{Likelihood.RUM}
\alias{Likelihood.RUM}
\title{Likelihood for general RUM}
\usage{
  Likelihood.RUM(Data, parameter, dist = "exp", range = 10,
    res = 0.1)
}
\arguments{
  \item{Data}{ranking data}

  \item{parameter}{Mean of Exponential Distribution}

  \item{dist}{norm or norm.fixedvariance}

  \item{range}{range}

  \item{res}{res}
}
\value{
  log likelihood
}
\description{
  Likelihood for general RUM
}
\examples{
data(Data.Test)
parameter = Generate.RUM.Parameters(5, "normal")
Likelihood.RUM(Data.Test,parameter, "norm")
parameter = Generate.RUM.Parameters(5, "exponential")
Likelihood.RUM(Data.Test,parameter, "exp")
}

