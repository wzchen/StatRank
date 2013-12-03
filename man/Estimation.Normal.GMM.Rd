\name{Estimation.Normal.GMM}
\alias{Estimation.Normal.GMM}
\title{GMM Method for Normal RUM model Estimation}
\usage{
  Estimation.Normal.GMM(Data.pairs, m, prior = 0)
}
\arguments{
  \item{Data.pairs}{data broken up into pairs}

  \item{m}{number of alternatives}

  \item{prior}{magnitude of fake observations input into
  the model}
}
\value{
  Estimated mean parameters for distribution of underlying
  normal (variance is fixed at 1)
}
\description{
  GMM Method for Normal RUM model Estimation
}
\examples{
data(Data.Test)
Data.Test.pairs <- Breaking(Data.Test, "full")
Estimation.Normal.GMM(Data.Test.pairs, 5)
}

