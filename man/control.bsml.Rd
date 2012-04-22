\name{control.bsml}
\alias{control.bsml}
\title{ Function to control all the options for HAS, BSML-C and BSML-S }
\description{
  The function generate all the options for the HAS and BSML methodologies.
}
\usage{
control.bsml(has.control=list(crit="GCV",idf=1.2,ridge=F),
			cv.control=list(cvfold=10),
			covpen.control=list(nboot=100,tauhat="Rice",crit="CIC",ct=2),
			gdf.control=list(nrep=100,tauhat="Rice",crit="GCV",ct=2),
			gcvpk.control=list(delta=.1,job=1005,lamlim=c(-6,2),tau=10,ntbl=100,maxtbl=100),
			dpc.control=list(c1=NA,c2=NA)
			)
}
\arguments{
  \item{has.control}{ Includes three controlling arguments: "crit" for basis selection criterion for HAS; Available options are "AIC", "BIC", "Cp" and "GCV". "idf" for IDF value for HAS; Default is 1.2; If value NA is supplied, it will be estimated. "ridge" for use ridge regression or not. TRUE or FALSE. }
  \item{cv.control}{ Includes only "cvfold", which controls the number of folds for Cross-Validation. The default value is 10.}
  \item{covpen.control}{ Includes 3 items: "nboot" controls the number of bootstrap samples; "tauhat" controls the value or the method of estimation for standard deviation for generating perturbations for bootstrap; "ct" controls the coefficient in AIC}
  \item{gdf.control}{ Controlling parameters for GDF estimation: "nrep" is for number of repetitions of the Monte Carlo simulation; "tauhat" controls the value or the method of estimation for standard deviation for generating perturbations for calculating the GDF. }
  \item{gcvpk.control}{ Controlling parameters for HAS. }
  \item{dpc.control}{ Coefficients used in Doubly Penalized Criterion (DPC) used in BSML-S. Default values are c1=log(n) and c2=2. }
}
\value{
  \item{has.control}{ Includes three controlling arguments: "crit" for basis selection criterion for HAS; Available options are "AIC", "BIC", "Cp" and "GCV". "idf" for IDF value for HAS; Default is 1.2; If value NA is supplied, it will be estimated. "ridge" for use ridge regression or not. TRUE or FALSE. }
  \item{cv.control}{ Includes only "cvfold", which controls the number of folds for Cross-Validation. The default value is 10.}
  \item{covpen.control}{ Includes 4 items: "nboot" controls the number of bootstrap samples; "tauhat" controls the value or the method of estimation for standard deviation for generating perturbations for bootstrap; "crit" controls the model selection criterion to use, available options are "GCV" and "CIC" (default); "ct" controls the coefficient in CIC. }
  \item{gdf.control}{ Controlling parameters for GDF estimation: "nrep" is for number of repetitions of the Monte Carlo simulation; "tauhat" controls the value or the method of estimation for standard deviation for generating perturbations for calculating the GDF; "crit" controls the model selection criterion to use, available options are "GCV" (default) and "CIC"; "ct" controls the coefficient in CIC. }
  \item{gcvpk.control}{ Controlling parameters for HAS. }
  \item{dpc.control}{ Coefficients used in Doubly Penalized Criterion (DPC) used in BSML-S. Default values are c1=log(n) and c2=2. }
}
\keyword{ smooth }
