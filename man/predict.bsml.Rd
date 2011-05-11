\name{predict.bsml}
\alias{predict.bsml}
\title{ Predict Function for BSML Objects }
\description{
  This function can be used to make predictions using the output objects from HAS, BSML-C and BSML-S. Confidence intervals can also be obtained using this function.
}
\usage{
\method{predict}{bsml}(object, bases.include=NULL, new.baseslist=NULL, confint=F, alpha=0.05, alg="a", addbas="cv", ...)
}
\arguments{
  \item{object}{ Output object from HAS, BSML-C OR BSML-S. }
  \item{bases.include}{ Of length (number of null bases) + (number of non-null libraries). Vector indicating the null bases and non-null libraries to include in the prediction. }
  \item{new.baseslist}{ Newly generated list of libraries for prediction. The original list when fitting can be used. For libraries not supplied, NULL must be in place of each them in the list. }
  \item{confint}{ True or False. If True, a confidence band of significance level alpha will be calculated. It will take extra computation time. The default is False.}
  \item{alpha}{ Level of significance. The default value is 0.05. }
  \item{alg}{ Algorithm to use. There are two options: "a" or "b", both of which use bootstrap. Algorithm "a" (default) uses models with the same number of basis functions as that for the final model for the perturbed responses in the Monte Carlo simulation. Algorithm "b" finds the optimal model size using MSE for the perturbed responses and starts with a "big" model that can be adjusted by the "addbas" argument also described here. }
  \item{addbas}{ Used with alg="b" only. Number of basis to be added to the initial fit to creat a "big" model. It can be a number or "cv" (default), when the CV approach is used to identify a "big" model. }
  \item{\dots}{ Not used. }
}
\value{
  \item{fit}{ Predicted values.}
  \item{lower}{ Lower bound of confidence band.}
  \item{upper}{ Upper bound of confidence band.}
  \item{allfits}{ Fits from which the percentile-based confidence band is derived from.}
  \item{alist}{ Lists of bases functions from which fits are calculated.}
  \item{afits}{ All fits based on the original design points. They can be used to derive a confidence band if original design points are used.}
}
\examples{
##############################
# Heavisine Function Example #
##############################
n=128
x=seq(0,1,length=n)
pt=(rep(1,n)\%o\%x)[,-c(n)]
f=2.2*(4*sin(4*pi*x)-sign(x-.3)-sign(.72-x))
set.seed(123)
sigma=3
y=f+sigma*rnorm(n)
lib1=stdz(cbind(rep(1,length=n)))
lib2=stdz(1*cbind((x-pt)>0))
baseslist=list(lib1,lib2)
bsmlc.obj=bsmlc(y,baseslist)
# To make predictions using only lib2
bsmlc.lib2.fit=predict(bsmlc.obj,c(0,1),baseslist)
}
\keyword{ smooth }
