\name{preGramm}
\alias{preGramm}
\title{
Preproccess the metabolome data and microbiome data
}
\description{
Preproccess the input data of metabolome and microbiome. Missing values may be imputed and filled (KNN method). Metabolome data and microbiomedata may be normalized and transformed by logarithm transformation and centered log-ratio (CLR) algorithm.
}
\usage{
preGramm(A,B,metaNor = TRUE,rarefaction = FALSE)
}
\arguments{
  \item{A}{
The metabolome data under pretreatment (SummarizedExperiment object).
}
  \item{B}{
 The microbiome data under pretreatment (SummarizedExperiment object).
}
  \item{metaNor}{
  Should metabolome data normalized? Using normalization when your metabolites are qualitative; and no normalization when the metabolites are quantitative. Default:TRUE.
  }
    \item{rarefaction}{
Resample an OTU table such that all samples have the same library size. Here refers to a repeated sampling procedure to assess species richness, first proposed in 1968 by Howard Sanders.(see wikipedia for more detail.)
Default:FALSE.
}
}
\value{
\item{x}{Metabolome preprocessed data}
\item{y}{Microbiome preprocessed data}
}
\references{
Gloor, G. B., Macklaim, J. M., Pawlowsky-Glahn, V., Egozcue, J. J., Microbiome Datasets Are Compositional: And This Is Not Optional. Front. Microbiol. 2017, 8 (2224).
}
\author{
Mengci Li, Dandan Liang, Tianlu Chen and Wei Jia
}
\seealso{
\link{naiveGramm} for naive correlation method;
\link{nlfitGramm} for nonlinear fitting;
\link{Gramm}: the whole strategy of this method.
}
\examples{
data("metabolites")
data("microbes")
preGramm(metabolites,microbes)
}

