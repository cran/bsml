\name{bsmlc}
\alias{bsmlc}
\title{ Basis Selection from Multiple Libraries Combined }
\description{
  This is the variant of the BSML methodology that combines all libraries into one and selects one basis only at a time during forward selection. Model complexity is assessed using either Generalized Degree of Freedom or Covariance Penalty.
}
\usage{
bsmlc(y,baseslist,maxbas=30,sub.maxbas=c(),backward=NULL,nobascheck=F,approach="gdf",control=NULL)
}
\arguments{
  \item{y}{ Response vector. }
  \item{baseslist}{ List of libraries. }
  \item{maxbas}{ Maximum number of bases to be selected. }
  \item{sub.maxbas}{ Maximum number of bases to be pre-selected for each library. }
  \item{backward}{ Backward elimination criterion. Available options are "AIC", "BIC", "Cp" and "GCV". Default is NULL, when no backward elimination is carried out. }
  \item{nobascheck}{ TRUE or FALSE. Default is FALSE. When TRUE, no safeguarding will be used to detect erroneously inputted basis functions. This is to save time when a loop is run. Please use bas.check() to help clean up the libraries before using this option. }
  \item{approach}{ GDF approach or Covariance Penalty approach. }
  \item{control}{ Pass control options to control.bsml(). Must be entered as a list. There are 4 control option groups: covpen.control, gdf.control, gcvpk.control and cv.control, each of which must be entered as a list. When an argument in an option group needs to be changed, only this argument has to be included in the option group list to be passed along. Other arguments will take default values. }
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
  \item{gdf}{Cumulated Generalized Degrees of Freedom for the corresponding number of chosen bases in the model. See Ref. 2 below.}
  \item{nb}{Total number of bases in the model, including both null and chosen bases.}
  \item{sigma_sq}{Estimated variance.}
  \item{dof}{Degrees of Freedom for the model.}
  \item{allisb}{ Lists of all bases selected during the Monte Carlo simulation.}
  \item{ptb}{ All perturbations used during the Monte Carlo simulation.}
  \item{allfits}{ All fits calculated during the Monte Carlo simulation.}
  \item{y}{ Response vector. }
  \item{baseslist}{ List of libraries. }
  \item{method}{ It has the value "bsmlc".}
  \item{maxbas}{ Maximum number of bases to be selected. }
  \item{sub.maxbas}{ Maximum number of bases to be pre-selected for each library. Specified using a vector. The first value is for the first non-null library.}
  \item{backward}{ Backward elimination criterion. Possible options are "AIC", "BIC", "Cp" and "GCV". Default is NULL, when no backward elimination is carried out. }
  \item{nobascheck}{ TRUE or FALSE. Default is FALSE. When TRUE, no safeguarding will be used to detect erroneously inputted basis functions. This is to save time when a loop is run. }
  \item{approach}{ GDF approach or Covariance Penalty approach. }
  \item{control}{ Control options passed to control.bsml().}
}
\references{ 
1) Zhen Luo, Grace Wahba. Hybrid Adaptive Splines. Journal of the American Statistical Association, Vol. 92, No. 437. (Mar., 1997), pp. 107-116.

2) Jianming Ye. On Measuring and Correcting the Effects of Data Mining and Model Selection. Journal of the American Statistical Association, Vol. 93, No. 441. (Mar., 1998), pp. 120-131.

3) Efron. The Estimation of Prediction Error: Covariance Penalties and Cross-Validation. Journal of the American Statistical Association, Vol. 99, No. 467. (Sep., 2004), pp. 619-632.
}
\author{ Junqing Wu \email{wjqpu@yahoo.com} , Jeff Sklar, Wendy Meiring and Yuedong Wang \email{yuedong@pstat.ucsb.edu} }
\seealso{
\code{\link{stdz}},\code{\link{has}},\code{\link{bsmls}},\code{\link{control.bsml}}
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
plot(x,y)
lines(x,bsmlc.obj$fit,col="green")
}
\keyword{ smooth }
