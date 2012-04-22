\name{bsml}
\alias{bsml}
\title{ Basis Selection from Multiple Libraries Ver. 1.5-1 }
\description{
  This is the implementation of the namesake methodology. To fit a nonparametric or semiparametric regression model, multiple libraries may be used to approximate the mean function. The program will adaptively select basis functions from each library and find the optimal model complexity by using the Generalized Degrees of Freedom or Covariance Penalty.
}
\usage{
bsml(y,baseslist,method="bsmlc",maxbas=30,sub.maxbas=c(),backward=NULL,nobascheck=F,bsml.approach="gdf",control=NULL)
}
\arguments{
  \item{y}{ Response vector. }
  \item{baseslist}{ List of libraries. }
  \item{method}{ Three methods available: HAS (Hybrid Adaptive Splines, see Ref. 1 below), BSML-C and BSML-S. Use "has", "bsmlc" and "bsmls" to specify.}
  \item{maxbas}{ Maximum number of bases to be selected. }
  \item{sub.maxbas}{ Maximum number of bases to be pre-selected for each library. Specify using a vector. Do not include a zero for the null library. }
  \item{backward}{ Backward elimination criterion. Available options are "AIC", "BIC", "Cp" and "GCV". Default is NULL, when no backward elimination is carried out. }
  \item{nobascheck}{ TRUE or FALSE. Default is FALSE. When TRUE, no safeguarding will be used to detect erroneously inputted basis functions. This is to save time when a loop is run. Please use bas.check() to help clean up the libraries before using this option. }
  \item{bsml.approach}{ GDF approach or Covariance Penalty approach. For BSML-C or BSML-S only. }
  \item{control}{ Pass control options to control.bsml(). Must be entered as a list. There are 5 control option groups, each of which must be entered as a list. When an argument in an option group needs to be changed, only this argument has to be included in the option group list to be passed along. Other arguments will take default values. For HAS, the available control option groups are: has.control, gdf.control and gcvpk.control. For BSML-C, the available control option groups are: covpen.control, gdf.control, gcvpk.control and cv.control. For BSML-S, the available control option groups are the same as those for BSML-C except an additional one: dpc.control, which controls the Doubly-Penalized Criterion used in the forward selection part of the procedure. }
}
\value{
  \item{coefficients}{ Coefficients}
  \item{residuals}{ Residuals}
  \item{fitted.values}{ Fitted Values}
  \item{chosen.bases.full}{ Full list of indices of chosen bases, excluding null bases.}
  \item{chosen.bases.trim}{ List of indices of chosen bases in the final model, excluding null bases.}
  \item{chosen.bases.matrix}{ Matrix of chosen bases in final model, excluding null bases.}
  \item{null.bases.matrix}{ Null bases used in fitting the model.}
  \item{lib.size}{ Size of each library used in fitting the model.}
  \item{score}{ Scores used to determine optimal model complexity. The first score is for the null model.}
  \item{idf}{ Inflated Degrees of Freedom for each chosen basis}
  \item{gdf}{ Cumulated Generalized Degrees of Freedom for the corresponding number of chosen bases in the model. See Ref. 2 below.}
  \item{nb}{ Total number of bases in the model, including both null and chosen bases.}
  \item{sigma_sq}{ Estimated variance.}
  \item{dof}{ Degrees of Freedom for the model.}
  \item{index.bas}{ HAS only. Indices of chosen bases in one-dimensional designations. For internal use.}
  \item{allisb}{ BSML only. Lists of all bases selected during the Monte Carlo simulation.}
  \item{ptb}{ BSML only. All perturbations used during the Monte Carlo simulation.}
  \item{allfits}{ BSML only. All fits calculated during the Monte Carlo simulation.}
  \item{y}{ Response vector. }
  \item{baseslist}{ List of libraries. }
  \item{method}{ Three methods possible: HAS (Hybrid Adaptive Splines, see Ref. 1 below), BSML-C and BSML-S. Represented by "has", "bsmlc" and "bsmls", respectively.}
  \item{maxbas}{ Maximum number of bases to be selected. }
  \item{sub.maxbas}{ Maximum number of bases to be pre-selected for each library. Specified using a vector. The first value is for the first non-null library.}
  \item{backward}{ Backward elimination criterion. Possible options are "AIC", "BIC", "Cp" and "GCV". Default is NULL, when no backward elimination is carried out. }
  \item{nobascheck}{ TRUE or FALSE. Default is FALSE. When TRUE, no safeguarding will be used to detect erroneously inputted basis functions. This is to save time when a loop is run. }
  \item{approach}{ For BSML-C or BSML-S only. GDF approach or Covariance Penalty approach. }
  \item{control}{ Control options passed to control.bsml().}
}
\references{ 
1) Zhen Luo, Grace Wahba. Hybrid Adaptive Splines. Journal of the American Statistical Association, Vol. 92, No. 437. (Mar., 1997), pp. 107-116.

2) Jianming Ye. On Measuring and Correcting the Effects of Data Mining and Model Selection. Journal of the American Statistical Association, Vol. 93, No. 441. (Mar., 1998), pp. 120-131.

3) Efron. The Estimation of Prediction Error: Covariance Penalties and Cross-Validation. Journal of the American Statistical Association, Vol. 99, No. 467. (Sep., 2004), pp. 619-632.
}
\author{ Junqing Wu \email{wjqpu@yahoo.com} , Jeff Sklar, Wendy Meiring and Yuedong Wang \email{yuedong@pstat.ucsb.edu} }
\seealso{
\code{\link{stdz}},\code{\link{has}},\code{\link{bsmlc}},\code{\link{bsmls}},\code{\link{control.bsml}}
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
has.obj=bsml(y,baseslist,method="has")
bsmlc.obj=bsml(y,baseslist,method="bsmlc")
bsmls.obj=bsml(y,baseslist,method="bsmls")
plot(x,y)
lines(x,has.obj$fit,col="red")
lines(x,bsmlc.obj$fit,col="green")
lines(x,bsmls.obj$fit,col="blue")
}
\keyword{ smooth }
