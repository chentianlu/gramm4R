\name{nlfitGramm}
\alias{nlfitGramm}
\title{
Plot the nonlinear fitting of input data
}
\description{
Plot, draw regression line and confidence interval, and show regression equation, R-square and P-value. The function includes the following models in the latest version: "line2P" (formula as: y=a*x+b), "line3P" (y=a*x^2+b*x+c), "log2P" (y=a*ln(x)+b), "exp2P" (y=a*exp(b*x)),"exp3P" (y=a*exp(b*x)+c), "power2P" (y=a*x^b), "power3P" (y=a*x^b+c), and "S" (y=a/(1+exp((x-b)/c))). Of each correlation pairs the figure who has the highest R-square among these models will be saved in a "pdf" file at working dirctory.
}
\usage{
nlfitGramm(X,Y)
}
\arguments{

  \item{X}{
Metabolome preproccessed data  (SummarizedExperiment object).
}
  \item{Y}{
Microbiome preproccessed data  (SummarizedExperiment object).
}
}
\details{
nothing

}
\value{
A file named "curve fit.pdf" will be created automatically.
}
\references{
Bates, D. M., and Watts, D. G. (2007) Nonlinear Regression Analysis and its Applications. Wiley.

Greenwell B. M., and Schubert-Kabban, C. M. (2014) investr: An R Package for Inverse Estimation. The R Journal, 6(1), 90-100.
}
\author{
Mengci Li, Dandan Liang, Tianlu Chen and Wei Jia
}

\seealso{
\link{naiveGramm} for naive correlation method;
\link{preGramm} for pretreatment;
\link{Gramm}: the whole strategy of this method.
}
\examples{
data("metabolites")
data("microbes")
nlfitGramm(metabolites,microbes)
}
\keyword{Nonlinear Fitting }

