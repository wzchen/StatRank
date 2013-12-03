\name{Generate.Pairwise.Probabilities}
\alias{Generate.Pairwise.Probabilities}
\title{Creates pairwise matrices to compare inference results with the empirical pairwise probabilities}
\usage{
  Generate.Pairwise.Probabilities(Data.pairs, Parameters,
    get.pairwise.prob, name.of.method)
}
\arguments{
  \item{Data.pairs}{datas broken into pairs}

  \item{Parameters}{The Parameter element of a result from
  an Estimation function}

  \item{get.pairwise.prob}{function that we use to generate
  the pairwise probability of beating}

  \item{name.of.method}{names of the alternatives}
}
\value{
  none
}
\description{
  Creates pairwise matrices to compare inference results
  with the empirical pairwise probabilities
}
\examples{
library(ggplot2)
library(gridExtra)
data(Data.Test)
Data.Test.pairs <- Breaking(Data.Test, "full")
Parameters <- Estimation.PL.GMM(Data.Test.pairs, 5)$Parameters
PL.Pairwise.Prob <- function(a, b) a$Mean / (a$Mean + b$Mean)
Generate.Pairwise.Probabilities(Data.Test.pairs, Parameters, PL.Pairwise.Prob, "PL on Test Data")
}

