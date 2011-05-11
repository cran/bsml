\name{has}
\alias{has}
\title{ Hybrid Adaptive Splines }
\description{
  This is a re-implementation of the HAS methodology with added flexibility of using more than one library and using estimated GDF values instead of a fixed IDF value in model selection.
}
\usage{
has(y,baseslist,maxbas=30,backward=NULL,nobascheck=F,control=NULL)
}
\arguments{
  \item{y}{ Response vector. }
  \item{baseslist}{ List of libraries. }
  \item{maxbas}{ Maximum number of bases to be selected. }
  \item{backward}{ Backward elimination criterion. Available options are "AIC", "BIC", "Cp" and "GCV". Default is NULL, when no backward elimination is carried out. }
  \item{nobascheck}{ TRUE or FALSE. Default is FALSE. When TRUE, no safeguarding will be used to detect erroneously inputted basis functions. This is to save time when a loop is run. Please use bas.check() to help clean up the libraries before using this option. }
  \item{control}{ Pass control options to control.bsml(). Must be entered as a list. There are 3 control option groups: has.control, gdf.control and gcvpk.control, each of which must be entered as a list. When an argument in an option group needs to be changed, only this argument has to be included in the option group list to be passed along. Other arguments will take default values.}
}
\value{
  \item{coefficients}{Coefficients}
  \item{residuals}{Residuals}
  \item{fitted.values}{Fitted Values}
  \item{chosen.bases.full}{Full list of indices of chosen bases, excluding null bases.}
  \item{chosen.bases.trim}{List of indices of chosen bases in the final model, excluding null bases.}
  \item{chosen.bases.matrix}{Matrix of chosen bases in final model, excluding null bases.}
  \item{null.bases.matrix}{Null bases used in fitting the model.}
  \item{lib.size}{Size of each library used in fitting the model.}
  \item{score}{Scores used to determine optimal model complexity. The first score is for the null model.}
  \item{idf}{Inflated Degrees of Freedom for each chosen basis}
  \item{gdf}{Cumulated Generalized Degrees of Freedom for the corresponding number of chosen bases in the model. Not applicable for HAS}
  \item{nb}{Total number of bases in the model, including both null and chosen bases.}
  \item{sigma_sq}{Estimated variance.}
  \item{dof}{Degrees of Freedom for the model.}
  \item{index.bas}{Indices of chosen bases in one-dimensional designations. For internal use.}
  \item{y}{ Response vector. }
  \item{baseslist}{ List of libraries. }
  \item{method}{ It should have the value "has".}
  \item{maxbas}{ Maximum number of bases to be selected. }
  \item{backward}{ Backward elimination criterion. Possible options are "AIC", "BIC", "Cp" and "GCV". Default is NULL, when no backward elimination is carried out. }
  \item{nobascheck}{ TRUE or FALSE. Default is FALSE. When TRUE, no safeguarding will be used to detect erroneously inputted basis functions. This is to save time when a loop is run. }
  \item{control}{ Control options passed to control.bsml().}
}
\references{ 
Zhen Luo, Grace Wahba. Hybrid Adaptive Splines. Journal of the American Statistical Association, Vol. 92, No. 437. (Mar., 1997), pp. 107-116.
}
\author{ Junqing Wu \email{wjqpu@yahoo.com} , Jeff Sklar, Wendy Meiring and Yuedong Wang \email{yuedong@pstat.ucsb.edu} }
\seealso{
\code{\link{stdz}},\code{\link{bsmlc}},\code{\link{bsmls}},\code{\link{control.bsml}}
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
has.obj=has(y,baseslist)
plot(x,y)
lines(x,has.obj$fit,col="red")
}
\keyword{ smooth }
