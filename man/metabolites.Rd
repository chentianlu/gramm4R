\name{metabolites}
\alias{metabolites}
\docType{data}
\title{
Test data of metabolome
}
\description{
The test data for examples at function \link{naiveGramm},  \link{Gramm},\link{nlfitGramm}
}
\usage{data("metabolites")}
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
  .. .. .. .. .. ..$ counts: num [1:40, 1:42] 8.75 3.31 5.09 1.49 6.2 8.81 2.37 1.76 4.47 2.18 ...
  .. .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. .. .. .. ..$ : NULL
  .. .. .. .. .. .. .. ..$ : chr [1:42] "B1" "B2" "B3" "B4" ...
  .. .. .. .. ..@ elementType    : chr "ANY"
  .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. ..@ metadata       : list()
  ..@ NAMES          : NULL
  ..@ elementMetadata:Formal class 'DataFrame' [package "S4Vectors"] with 6 slots
  .. .. ..@ rownames       : NULL
  .. .. ..@ nrows          : int 40
  .. .. ..@ listData       :List of 1
  .. .. .. ..$ X: Factor w/ 40 levels "Decadienylcarnitine",..: 16 35 31 25 3 4 1 8 9 6 ...
  .. .. ..@ elementType    : chr "ANY"
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ metadata       : list()
}
\examples{
data(metabolites)
## maybe str(metabolites) ; plot(metabolites) ...
}
\keyword{datasets}
