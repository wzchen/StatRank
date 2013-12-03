\name{Generate.RUM.Ranks}
\alias{Generate.RUM.Ranks}
\title{Generate observation of ranks given parameters}
\usage{
  Generate.RUM.Ranks(parameter, m, n, distribution)
}
\arguments{
  \item{parameter}{list of parameters from the
  Generate.RUM.Parameters function}

  \item{m}{number of alternatives}

  \item{n}{number of agents}

  \item{distribution}{can be either 'normal' or
  'exponential'}
}
\value{
  a matrix of observed rankings
}
\description{
  Given a list of parameters (generated via the Generate
  RUM Parameters function), generate random utilities from
  these models and then return their ranks
}
\examples{
parameter= Generate.RUM.Parameters(10, "normal")
Generate.RUM.Ranks(parameter,m=10,n=5,"normal")
parameter= Generate.RUM.Parameters(10, "exponential")
Generate.RUM.Ranks(parameter,m=10,n=5,"exponential")
}

