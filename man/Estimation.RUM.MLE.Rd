\name{Estimation.RUM.MLE}
\alias{Estimation.RUM.MLE}
\title{Performs parameter estimation for a Random Utility Model}
\usage{
  Estimation.RUM.MLE(Data, iter = 10, dist)
}
\arguments{
  \item{Data}{data in either partial or full rankings}

  \item{iter}{number of EM iterations to run}

  \item{dist}{underlying distribution. Can be "normal",
  "normal.fixedvariance", "exp"}
}
\value{
  parameters of the latent RUM distributions
}
\description{
  This function supports RUMs 1) Normal 2) Normal with
  fixed variance (fixed at 1) 3) Exponential
}
\examples{
data(Data.Test)
Estimation.RUM.MLE(Data.Test, iter = 2, dist="norm")
}

