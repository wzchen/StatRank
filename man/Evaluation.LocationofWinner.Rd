\name{Evaluation.LocationofWinner}
\alias{Evaluation.LocationofWinner}
\title{Calculates the location of the True winner in the estimated ranking}
\usage{
  Evaluation.LocationofWinner(EstimatedRank, TrueRank)
}
\arguments{
  \item{EstimatedRank}{estimated ranking}

  \item{TrueRank}{true ranking}
}
\value{
  The location of the true best in the estimated rank
}
\description{
  Calculates the location of the True winner in the
  estimated ranking
}
\examples{
rank1 <- scramble(1:10)
rank2 <- scramble(1:10)
Evaluation.LocationofWinner(rank1, rank2)
}

