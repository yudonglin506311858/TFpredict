\name{get_promoter_sequences}
\alias{get_promoter_sequences}
\title{Get promoter sequences for target genes}
\usage{
get_promoter_sequences(
  target_gene,
  species = c("human", "mouse"),
  promoter_start = 2000,
  promoter_end = 200
)
}
\arguments{
\item{target_gene}{A character vector of gene symbols}

\item{species}{A character string specifying species ("human" or "mouse")}

\item{promoter_start}{Distance upstream of TSS (default: 2000)}

\item{promoter_end}{Distance downstream of TSS (default: 200)}
}
\value{
A DNAStringSet object containing promoter sequences
}
\description{
This function retrieves promoter sequences for given target genes from the specified
genome build. The promoter region is defined relative to the transcription start site (TSS).
}
\examples{
promoter_seq <- get_promoter_sequences("GATA1", "human", 2000, 200)
}
