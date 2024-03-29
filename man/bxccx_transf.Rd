% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TBoCo.R
\name{bxccx_transf}
\alias{bxccx_transf}
\title{Title}
\usage{
bxccx_transf(bas)
}
\arguments{
\item{bas}{matrix to transform ( each columns is transform using a boxCox transformation)
and then it is scaled}
}
\value{
the transformed matrix ( after a scaling step)
}
\description{
Title
}
\examples{
x <- rnorm(100)
y <- rpois(100,3)
z <- rnorm(100,10,2)
z[z<0] <- 0
bas <- cbind(x = exp(x), y = exp(y), z = z)
bxccx_transf(bas)
}
