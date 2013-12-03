\name{generateC}
\alias{generateC}
\title{Generate a matrix of pairwise wins}
\usage{
  generateC(Data.pairs, m, Zemel = FALSE, prior = 0)
}
\arguments{
  \item{Data.pairs}{the data broken up into pairs}

  \item{m}{the total number of alternatives}

  \item{Zemel}{whether or not this function should use the
  rank differences instead of indicators for when i beats
  j}

  \item{prior}{the intial "fake data" that you want to
  include in C. A prior of 1 would mean that you initially
  "observe" that all alternatives beat all other
  alternatives exactly once.}
}
\value{
  a Count matrix of how many times alternative i has beat
  alternative j
}
\description{
  This function takes in data that has been broken up into
  pair format. The user is given a matrix C, where element
  C[i, j] represents exactly how many times alternative i
  has beaten alternative j
}
\examples{
data(Data.Test)
Data.Test.pairs <- Breaking(Data.Test, "full")
generateC(Data.Test.pairs, 5, Zemel = FALSE, prior = 0)
}

