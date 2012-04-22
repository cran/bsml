\name{summary.bsml}
\alias{summary.bsml}
\title{ Summary Function for BSML Objects }
\description{
  This function provides a quick summary for objects from BSML, HAS, BSML-C and BSML-S.
}
\usage{
\method{summary}{bsml}(object, ...)
}
\arguments{
  \item{object}{ Output object from BSML, HAS, BSML-C OR BSML-S. }
  \item{\dots}{ Not used. }
}
\value{
  \item{call}{Function call.}
  \item{nb}{Number of bases in the final model.}
  \item{dof}{Degree of Freedom for the model.}
  \item{sigma}{Estimated standard error.}
}
\keyword{ smooth }
