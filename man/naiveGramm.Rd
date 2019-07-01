\name{naiveGramm}
\alias{naiveGramm}
\title{
Corrlation methods of Gramm
}
\description{
A function to test the association among metabolites and microbes. Using linear (linear regression) or nonlinear (Maximal Information Coefficient) methods.
}
\usage{
naiveGramm(x,y,z,r,alpha)
}
\arguments{
  \item{x}{
A SummarizedExperiment object contains data of preproccessed data, where rows represent features of metabolites and columns represent samples.
}
  \item{y}{
A SummarizedExperiment object contains data of microbiome preproccessed data , where rows represent features of microbes and columns represent samples.
}
  \item{z}{
An optional SummarizedExperiment object containing the confounders.
}
  \item{r}{
The linear regression coefficients threshold for using nonlinear method. Default: 0.5.
}
 \item{alpha}{
  The linear regression p-value threshold for using nonlinear method.Default: 0.05.
}
}
\value{
  \item{r}{Correlation coefficients.}
  \item{p}{Correlation p values.}
  \item{type}{Correlation methods for computing (linear or nonlinear).}

}
\references{
Chambers, J. M. (1992) Linear models. Chapter 4 of Statistical Models in S eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.
D. Reshef, Y. Reshef, H. Finucane, S. Grossman, G. McVean, P. Turnbaugh, E. Lander, M. Mitzenmacher, P. Sabeti. (2011) Detecting novel associations in large datasets. Science 334, 6062.
D. Albanese, M. Filosi, R. Visintainer, S. Riccadonna, G. Jurman, C. Furlanello. minerva and minepy: a C engine for the MINE suite and its R, Python and MATLAB wrappers. Bioinformatics (2013) 29(3): 407-408.
}
\author{
Mengci Li, Dandan Liang, Tianlu Chen and Wei Jia
}
\seealso{
\link{preGramm} for pretreatment;\link{nlfitGramm} for nonlinear fitting;\link{Gramm}: the whole strategy of this method.
}
\examples{
data("metabolites")
data("microbes")
data("covariates")
naiveGramm(metabolites,microbes,covariates)
}
