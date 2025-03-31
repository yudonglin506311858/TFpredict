\name{TFpredict}
\alias{TFpredict}
\title{Predict transcription factors binding to promoter regions}
\usage{
TFpredict(
  target_gene,
  species = c("human", "mouse"),
  promoter_start = 2000,
  promoter_end = 200,
  min_score = "85\%"
)
}
\arguments{
\item{target_gene}{A character vector of gene symbols}

\item{species}{A character string specifying species ("human" or "mouse")}

\item{promoter_start}{Distance upstream of TSS (default: 2000)}

\item{promoter_end}{Distance downstream of TSS (default: 200)}

\item{min_score}{Minimum score threshold for PWM matching (default: "85\%")}
}
\value{
A data.frame with predicted transcription factors and their binding sites
}
\description{
This function predicts transcription factors that may bind to the promoter regions
of target genes by scanning the promoter sequences against the JASPAR database
using position weight matrices (PWMs).
}
\examples{
TF <- TFpredict(target_gene = "GATA1", species = "human",
                promoter_start = 2000, promoter_end = 200)
}
