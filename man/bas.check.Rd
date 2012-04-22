\name{bas.check}
\alias{bas.check}
\title{ Utility Function to Detect Erroneously Included Basis Functions }
\description{
  The purpose of this function is to detect repeated (redundant) basis functions across libraries, as well as zero columns and constant columns.
}
\usage{
bas.check(baseslist)
}
\arguments{
  \item{baseslist}{ List of basis function matrices. Each matrix of which is a library. }
}
\value{
  \item{selindex}{ List of tag vectors. Each vector corresponds to a library indicating which columns are excluded (indexed by "1"s). }
}
\keyword{ smooth }
