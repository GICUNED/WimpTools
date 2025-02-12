% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SystemDynamicsFunctions.R
\name{pcsd_summary}
\alias{pcsd_summary}
\title{PCSD summary -- pcsd_summary()}
\usage{
pcsd_summary(scn)
}
\arguments{
\item{scn}{A "scn" S3 class object, the result of the \code{\link{scenariomatrix}} function.
This object contains the scenario matrix and relevant metadata.}
}
\value{
Returns a matrix with the PCSD summary.
}
\description{
This function returns a summary of the PCSD. It informs us the
initial and final value of each construct and the difference between them.
}
\examples{
# Example Scenario Matrix (scn)
example_scn <- scenariomatrix(
  wimp = example.wimp,
  infer = "self dynamics",
  thr = "saturation",
  act.vector = c(0.1, -0.2, 0.3, -0.1, 0)
)

# Calculate PCSD summary for each construct
pcsd_summary(example_scn)

}
