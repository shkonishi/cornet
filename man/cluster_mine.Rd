% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_mine.R
\name{cluster_mine}
\alias{cluster_mine}
\title{The reshape of mine output}
\usage{
cluster_mine(cl_dat)
}
\arguments{
\item{cl_dat}{A list of dataframes, result of 'cornet::cluster_mat'}
}
\value{
A list of dataframes, result of minerva::mine soreted by 'TIC'.
}
\description{
The result of mine is multi matrix of all node pairs, and that is too many memory usage using big data.
   this program retuns dataframe as edge list and atrributes of these edge.
}
\examples{
# sample data, result of 'cluster_mat'
# data("cluster_dat")
#
}
