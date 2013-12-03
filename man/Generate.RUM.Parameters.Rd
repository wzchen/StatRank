\name{Generate.RUM.Parameters}
\alias{Generate.RUM.Parameters}
\title{Parameter Generation for a RUM model}
\usage{
  Generate.RUM.Parameters(m, distribution)
}
\arguments{
  \item{m}{number of sets of parameters to be drawn}

  \item{distribution}{either 'normal' or 'exponential'}
}
\value{
  a list of RUM parameters
}
\description{
  Exponential models mean parameters are drawn from a
  uniform distribution Normal models, mean and standard
  devaition parameters are drawn from a standard unifrom
}
\examples{
Generate.RUM.Parameters(10, "normal")
Generate.RUM.Parameters(10, "exponential")
}

