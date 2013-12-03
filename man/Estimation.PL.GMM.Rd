\name{Estimation.PL.GMM}
\alias{Estimation.PL.GMM}
\title{GMM Method for PL Estimation}
\usage{
  Estimation.PL.GMM(Data.pairs, m, prior = 0)
}
\arguments{
  \item{Data.pairs}{data broken up into pairs}

  \item{m}{number of alternatives}

  \item{prior}{magnitude of fake observations input into
  the model}
}
\value{
  Estimated mean parameters for distribution of underlying
  exponential
}
\description{
  GMM Method for PL Estimation
}
\examples{
data(Data.Test)
Data.Test.pairs <- Breaking(Data.Test, "full")
Estimation.PL.GMM(Data.Test.pairs, 5)
}

