\name{microbes}
\alias{microbes}
\docType{data}
\title{
Test data of microbiome
}
\description{
The test data for examples at function\link{naiveGramm},  \link{Gramm},\link{nlfitGramm}
}
\usage{data("microbes")}
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
  .. .. .. .. .. ..$ counts: num [1:85, 1:42] 1 1 13 2 4 8 11 24 6 1 ...
  .. .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. .. .. .. ..$ : NULL
  .. .. .. .. .. .. .. ..$ : chr [1:42] "B1" "B2" "B3" "B4" ...
  .. .. .. .. ..@ elementType    : chr "ANY"
  .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. ..@ metadata       : list()
  ..@ NAMES          : NULL
  ..@ elementMetadata:Formal class 'DataFrame' [package "S4Vectors"] with 6 slots
  .. .. ..@ rownames       : NULL
  .. .. ..@ nrows          : int 85
  .. .. ..@ listData       :List of 1
  .. .. .. ..$ X: Factor w/ 84 levels "[Ruminococcus] spp",..: 74 24 1 12 20 21 29 70 65 20 ...
  .. .. ..@ elementType    : chr "ANY"
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ metadata       : list()
}
\examples{
data(microbes)
## maybe str(microbes) ; plot(microbes) ...
}
\keyword{datasets}
