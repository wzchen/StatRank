\name{Estimation.PL.MLE}
\alias{Estimation.PL.MLE}
\title{Performs parameter estimation for a Random Utility Model using an MM algorithm}
\usage{
  Estimation.PL.MLE(Data, iter = 10)
}
\arguments{
  \item{Data}{data in either partial or full rankings}

  \item{iter}{number of MM iterations to run}
}
\value{
  list of estimated means (Gamma) and the log likelihoods
}
\description{
  Performs parameter estimation for a Random Utility Model
  using an MM algorithm
}
\examples{
data(Data.Test)
Estimation.PL.MLE(Data.Test)
}

