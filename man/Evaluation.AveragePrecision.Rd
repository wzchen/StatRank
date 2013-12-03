\name{Evaluation.AveragePrecision}
\alias{Evaluation.AveragePrecision}
\title{Calculates the Average Precision}
\usage{
  Evaluation.AveragePrecision(EstimatedRank,
    RelevanceLevel)
}
\arguments{
  \item{EstimatedRank}{estimated ranking}

  \item{RelevanceLevel}{score for the document}
}
\value{
  The AP for this estimation and relevance level
}
\description{
  Calculates the Average Precision
}
\examples{
EstimatedRank <- scramble(1:10)
RelevanceLevel <- runif(10)
Evaluation.AveragePrecision(EstimatedRank, RelevanceLevel)
}

