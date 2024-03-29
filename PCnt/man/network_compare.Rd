% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{network_compare}
\alias{network_compare}
\title{Compare two networks with same nodes embedding}
\usage{
network_compare(
  output_dag,
  truth,
  gvarType,
  gvarShape,
  manual_colors = NULL,
  topE = 200,
  seed = 1
)
}
\arguments{
\item{output_dag}{adjacency matrix.}

\item{truth}{number of top edges selected. default 150.}

\item{gvarType}{description.}

\item{gvarShape}{description.}

\item{manual_colors}{node colors.}

\item{topE}{default 200.}

\item{seed}{default 1.}
}
\description{
Compare two networks with same nodes embedding
}
