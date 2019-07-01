\name{covariates}
\alias{covariates}
\docType{data}
\title{
Test data of covariates
}
\description{
The test data for examples at function\link{naiveGramm},  \link{Gramm},\link{nlfitGramm}
}
\usage{data("covariates")}
\format{
  The format is:
Formal class 'SummarizedExperiment' [package "SummarizedExperiment"] with 5 slots
  ..@ colData        :Formal class 'DataFrame' [package "S4Vectors"] with 6 slots
  .. .. ..@ rownames       : chr [1:42] "B1" "B2" "B3" "B4" ...
  .. .. ..@ nrows          : int 42
  .. .. ..@ listData       : Named list()
  .. .. ..@ elementType    : chr "ANY"
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ assays         :Formal class 'SimpleAssays' [package "SummarizedExperiment"] with 1 slot
  .. .. ..@ data:Formal class 'SimpleList' [package "S4Vectors"] with 4 slots
  .. .. .. .. ..@ listData       :List of 1
  .. .. .. .. .. ..$ counts: int [1:3, 1:42] 25 1 90 26 1 70 27 1 68 28 ...
  .. .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. .. .. .. ..$ : NULL
  .. .. .. .. .. .. .. ..$ : chr [1:42] "B1" "B2" "B3" "B4" ...
  .. .. .. .. ..@ elementType    : chr "ANY"
  .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. ..@ metadata       : list()
  ..@ NAMES          : NULL
  ..@ elementMetadata:Formal class 'DataFrame' [package "S4Vectors"] with 6 slots
  .. .. ..@ rownames       : NULL
  .. .. ..@ nrows          : int 3
  .. .. ..@ listData       :List of 1
  .. .. .. ..$ X: Factor w/ 3 levels "age","gender",..: 1 2 3
  .. .. ..@ elementType    : chr "ANY"
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ metadata       : list()
}

\examples{
data(covariates)
## maybe str(covariates) ; plot(covariates) ...
}
\keyword{datasets}
