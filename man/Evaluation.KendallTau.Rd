\name{Evaluation.KendallTau}
\alias{Evaluation.KendallTau}
\title{Calculates the Kendall Tau correlation between two ranks}
\usage{
  Evaluation.KendallTau(rank1, rank2)
}
\arguments{
  \item{rank1}{two rankings. Order does not matter}

  \item{rank2}{two rankings. Order does not matter}
}
\value{
  The Kendall Tau correlation
}
\description{
  Calculates the Kendall Tau correlation between two ranks
}
\examples{
rank1 <- scramble(1:10)
rank2 <- scramble(1:10)
Evaluation.KendallTau(rank1, rank2)
}

