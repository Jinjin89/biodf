% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_enrich.R
\name{fun_enrich_patch}
\alias{fun_enrich_patch}
\title{enrichment analysis for input_dta using enriched,}
\usage{
fun_enrich_patch(
  input_dta,
  input_methods = "ora",
  outfile = NULL,
  pval = 0.2,
  logFC = "logFC",
  gene = "ID",
  term2gene = c("term", "gene"),
  enrich_db = list(hallmark = biodata::db_hallmark, kegg = biodata::db_kegg_v2, go =
    biodata::db_go, reactome = biodata::db_reactome, reactome_metabolism =
    biodata::db_reactome_metabolism, biocarta = biodata::db_biocarta),
  return_list = list()
)
}
\arguments{
\item{input_dta}{genes vector for ora, data.frame for gsea}

\item{input_methods}{ora or gsea}

\item{outfile}{the outfile position,default is NULL, not saving the results}

\item{pval}{pval cutoff for pvalue and qvalue}

\item{logFC}{logFC}

\item{gene}{gene}

\item{term2gene}{term2gene}

\item{enrich_db}{db list, must be names list}

\item{return_list}{return_list}
}
\value{
list of different enrichment obj
}
\description{
enrichment analysis for input_dta using enriched,
}
