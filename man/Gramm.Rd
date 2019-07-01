\name{Gramm}
\alias{Gramm}
\title{
Get the association between metabolites and microbes
}
\description{
The entire strategy to get the association between metabolites and microbes, using linear and nonlinear methods, and plot the regression figures.
}
\usage{
Gramm(A,B,C,metaNor,r,alpha)
}
\arguments{
   \item{A}{
A SummarizedExperiment object contains data of metabolome, where rows represent features of metabolites and columns represent samples.
}
  \item{B}{
A SummarizedExperiment object contains data of microbiome, where rows represent features of microbes and columns represent samples.
}
  \item{C}{
A SummarizedExperiment object contains data of covariates, where rows represent features of covariates and columns represent samples.
}
  \item{metaNor}{
  Should metabolome data normalized? Using normalization when your metabolites are qualitative; and no normalization when the metabolites are quantitative. Default:TRUE.
  }
  \item{r}{
The linear regression coefficients threshold for using nonlinear method. Default: 0.5.
}
 \item{alpha}{
  The linear regression p-value threshold for using nonlinear method.Default: 0.05.
}
}
\value{
  \item{pretreatment}{The result of pretreatment}
  \item{correlation}{The result of correlation, see \link{naiveGramm} for detail }
A file named  "R value top 10 pairs.pdf" will be created automatically (corrlation coefficient top 10 pairs) .
}
\references{
Gloor, G. B., Macklaim, J. M., Pawlowsky-Glahn, V., Egozcue, J. J., Microbiome Datasets Are Compositional: And This Is Not Optional. Front. Microbiol. 2017, 8 (2224).
Chambers, J. M. (1992) Linear models. Chapter 4 of Statistical Models in S eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.
D. Reshef, Y. Reshef, H. Finucane, S. Grossman, G. McVean, P. Turnbaugh, E. Lander, M. Mitzenmacher, P. Sabeti. (2011) Detecting novel associations in large datasets. Science 334, 6062.
D. Albanese, M. Filosi, R. Visintainer, S. Riccadonna, G. Jurman, C. Furlanello. minerva and minepy: a C engine for the MINE suite and its R, Python and MATLAB wrappers. Bioinformatics (2013) 29(3): 407-408.
}
\author{
Mengci Li, Dandan Liang, Tianlu Chen and Wei Jia
}
\seealso{
\link{preGramm} for pretreatment;\link{nlfitGramm} for nonlinear fitting;\link{naiveGramm} for naive correlation method.
}
\examples{
data("metabolites");data("microbes");data("covariates")
Gramm(metabolites,microbes,covariates)
}
