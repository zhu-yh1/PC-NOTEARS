% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{PC_nt}
\alias{PC_nt}
\title{PC+NOTEARS}
\usage{
PC_nt(data, lambda1 = 0.01, alpha = 0.05)
}
\arguments{
\item{lambda1}{lambda for l1 penalty, default 0.01.}

\item{alpha}{alpha for PC algorithm, default 0.05.}

\item{X}{a numeric data matrix with samples as rows and variables as columns.}

\item{losses}{"logistic" for binary variables and "l2" for continuous variables.}

\item{no_parent}{constraint that certain variables do not have parent nodes.}

\item{no_edge}{a 0,1 matrix where 1 species no_edge.}

\item{max.iter}{default 100.}

\item{h.tol}{default 1e-6.}

\item{rho.max}{default 1e+6.}

\item{w.threshold}{default 0.01.}
}
\value{
Adjacency matrix.
}
\description{
PC+NOTEARS
}
