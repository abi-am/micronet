% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{construct_network}
\alias{construct_network}
\title{construct_network}
\usage{
construct_network(
  ps,
  level.columns = NULL,
  tax.lvl = "Species",
  subset.columns = list(),
  rabun.thr = 0,
  prev.thr = 0,
  rho.pos.thr = 0.8,
  rho.neg.thr = NULL,
  p.adj.thr = 0.01,
  n.bootstrap.batches = 1000
)
}
\arguments{
\item{ps}{phyloseq object}

\item{level.columns}{Vector of column names categories of which should be treated separately, for example c('Treatment'), Default: NULL}

\item{tax.lvl}{Taxonomic level, Default: 'Species'}

\item{subset.columns}{named list for subsetting the phyloseq object based on sample data, where the names are column names and values are vectors of values to keep, Default: list()}

\item{rabun.thr}{threshold which is used for filtering out abundances with lower relative abundance, Default: 0}

\item{prev.thr}{threshold which is used for filtering out abundances with lower prevalence, Default: 0}

\item{rho.pos.thr}{threshold which is used for filtering out pairs with lower positive correlation, Default: 0.8}

\item{rho.neg.thr}{threshold which is used for filtering out pairs with higher negative correlation, if not set -rho.pos.thr is taken, Default: NULL}

\item{p.adj.thr}{threshold which is used for filtering out pairs with larger p adjusted value, Default: 0.01}

\item{n.bootstrap.batches}{number of times to bootstrap, Default: 1000}
}
\value{
named list, containing corr.table, network.object, aggregated.genera.species which are the correlation table, igraph oject and aggregated genera if tax.lvl is 'Species' respectively.
}
\description{
Generates co-occurrence network.
}
