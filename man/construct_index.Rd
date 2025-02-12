% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AdjustmentFunctions.R
\name{construct_index}
\alias{construct_index}
\title{Frencuency and proportions of constructs -- construct_index()}
\usage{
construct_index(wimp)
}
\arguments{
\item{wimp}{Subject's WimpGrid. It must be a "wimp" S3 object
imported by the \code{\link{importwimp}} function.}
}
\value{
A matrix with the frequency and proportion of congruents, discrepants
       , dilemmatics and undefined constructs.
}
\description{
This function calculates frequency and proportion of
             congruents, discrepants, dilemmatics and undefined constructs.
}
\examples{

construct_index(example.wimp)

}
\author{
Alejandro Sanfeliciano
}
