\name{Evaluation.NDCG}
\alias{Evaluation.NDCG}
\title{Calculates the Normalized Discounted Cumluative Gain}
\usage{
  Evaluation.NDCG(EstimatedRank, RelevanceLevel)
}
\arguments{
  \item{EstimatedRank}{estimated ranking}

  \item{RelevanceLevel}{score for the document}
}
\value{
  The NDCG for this estimation and relevance level
}
\description{
  Calculates the Normalized Discounted Cumluative Gain
}
\examples{
EstimatedRank <- scramble(1:10)
RelevanceLevel <- runif(10)
Evaluation.NDCG(EstimatedRank, RelevanceLevel)
}

