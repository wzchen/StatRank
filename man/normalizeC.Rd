\name{normalizeC}
\alias{normalizeC}
\title{Converts a pairwise count matrix into a probability matrix}
\usage{
  normalizeC(C)
}
\arguments{
  \item{C}{original pairwise count matriix}
}
\value{
  a pairwise probability matrix
}
\description{
  Converts a pairwise count matrix into a probability
  matrix
}
\examples{
C= matrix(c(2,4,3,5),2,2)
normalizeC(C)
}

