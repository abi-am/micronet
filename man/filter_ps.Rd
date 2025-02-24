% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{filter_ps}
\alias{filter_ps}
\title{filter_ps}
\usage{
filter_ps(
  ps,
  tax.lvl = "Species",
  relative.abundance = F,
  rabun.thr = 0,
  prev.thr = 0,
  subset.columns = list(),
  level.columns = NULL,
  species.names.modify = T
)
}
\arguments{
\item{ps}{phyloseq object}

\item{tax.lvl}{tax level, Default: 'Species'}

\item{relative.abundance}{make the abundance relative if set TRUE, Default: F}

\item{rabun.thr}{threshold which is used for filtering out abundances with lower relative abundance, Default: 0}

\item{prev.thr}{threshold which is used for filtering out abundances with lower prevalence, Default: 0}

\item{subset.columns}{named list for subsetting the phyloseq object based on sample data, where the names are column names and values are vectors of values to keep, Default: list()}

\item{level.columns}{columns based on which to separately the phyloseq object and separately do prevalence filtering, Default: NULL}

\item{species.names.modify}{if TRUE, changes the taxa_names names to consider the aggregation when tax.lvl is 'Species', Default: T}
}
\value{
filtered phyloseq object
}
\description{
Filter the phyloseq object
}
