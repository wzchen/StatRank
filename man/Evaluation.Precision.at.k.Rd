\name{Evaluation.Precision.at.k}
\alias{Evaluation.Precision.at.k}
\title{Calculates the Average Precision at k}
\usage{
  Evaluation.Precision.at.k(EstimatedRank, RelevanceLevel,
    k)
}
\arguments{
  \item{EstimatedRank}{estimated ranking}

  \item{RelevanceLevel}{score for the document}

  \item{k}{positive that we want to run this algorithm for}
}
\value{
  The AP at k for this estimation and relevance level
}
\description{
  Calculates the Average Precision at k
}
\examples{
EstimatedRank <- scramble(1:10)
RelevanceLevel <- runif(10)
Evaluation.Precision.at.k(EstimatedRank, RelevanceLevel, 5)
}

