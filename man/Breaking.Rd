\name{Breaking}
\alias{Breaking}
\title{Breaks full or partial orderings into pairwise comparisons}
\usage{
  Breaking(Data, method, k = NULL)
}
\arguments{
  \item{Data}{data in either full or partial ranking
  format}

  \item{method}{- can be full, adjacent, top or
  top.patrixl}

  \item{k}{This applies to the top method, choose which top
  k to focus on}
}
\value{
  Pairwise breakings, where the three columns are winner,
  loser and rank distance (latter used for Zemel)
}
\description{
  Given full or partial orderings, this function will
  generate pairwise comparise Options 1. full - All
  available pairwise comparisons. This is used for partial
  rank data where the ranked objects are a random subset of
  all objects 2. adjacent - Only adjacent pairwise
  breakings 3. top - also takes in k, will break within top
  k and will also generate pairwise comparisons comparing
  the top k with the rest of the data 4. top.partial - This
  is used for partial rank data where the ranked
  alternatives are preferred over the non-ranked
  alternatives
}
\examples{
data(Data.Test)
Data.Test.pairs <- Breaking(Data.Test, "full")
}

