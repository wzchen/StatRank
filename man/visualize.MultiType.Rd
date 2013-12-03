\name{visualize.MultiType}
\alias{visualize.MultiType}
\title{Multitype Random Utility visualizer}
\usage{
  visualize.MultiType(multitype.output, names)
}
\arguments{
  \item{multitype.output}{output from a multitype fitter}

  \item{names}{names of alternatives}
}
\value{
  none
}
\description{
  Multitype Random Utility visualizer
}
\examples{
library(ggplot2)
library(grid)
Data.Tiny <- matrix(c(1, 2, 3, 3, 2, 1, 1, 2, 3), ncol = 3, byrow = TRUE)
multitype.output <- Estimation.RUM.MultiType.MLE(Data.Tiny, iter = 2, dist = "norm", ratio = .5)
names <- 1:5
plots <- visualize.MultiType(multitype.output, names)
# the list of plots is appropriate for passing into the multiplot function at
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# as input to the plotlist argument
}

