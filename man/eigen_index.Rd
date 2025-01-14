% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CentralityFunctions.R
\name{eigen_index}
\alias{eigen_index}
\title{Eigenvalue Centrality index -- eigen_index()}
\usage{
eigen_index(wimp, matrix = "weights", num.vectors = 2)
}
\arguments{
\item{wimp}{wimp An object of class 'wimp' (weighted implications grid)}

\item{matrix}{A character string specifying which matrix to use for the centrality analysis. Accepted values are
'direct', 'weights', or 'implications'. Default is 'implications'.}

\item{num.vectors}{An integer specifying the number of eigenvectors to use for computing centrality scores.}
}
\value{
A dataframe containing the constructs' names and their respective centrality scores.
}
\description{
This function calculates centrality scores for constructs within a `wimp` object,
             based on the eigenvalue decomposition of a specified adjacency matrix from the WIMP scores.
             It supports analyzing centrality using the 'direct', 'weights', or 'implications' matrices.
             The centrality calculation is performed over the specified number of eigenvectors.
}
\examples{

eigen_index(example.wimp)

}
\author{
Carlos Hurtado
}
