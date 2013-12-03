\name{scores.to.order}
\alias{scores.to.order}
\title{Converts scores to a ranking}
\usage{
  scores.to.order(scores)
}
\arguments{
  \item{scores}{the scores (e.g. means) of a set of
  alternatives}
}
\value{
  an ordering of the index of the winner, second place,
  etc.
}
\description{
  takes in vector of scores (with the largest score being
  the one most prefered) and returns back a vector of
  WINNER, SECOND PLACE, ... LAST PLACE
}
\examples{
scores <- Generate.RUM.Parameters(10, "exponential")$Mean
scores.to.order(scores)
}

