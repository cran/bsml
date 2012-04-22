###########################################
# Basis Selection from Multiple Libraries #
# Updated Apr 22 2012					  #
# Version 1.5-1                           #
###########################################

bsml <- function(y,baseslist,method="bsmlc",maxbas=30,
				 sub.maxbas=c(),backward=NULL,nobascheck=F,
                 bsml.approach="gdf",control=NULL
                 ) 
{

   # Let idf=NA if you want idf to be estimated from the data
   #### Select Spatially Adaptive Regression Method ####
   # has = Hybrid Adaptive Splines (with option to estimate IDF)
   # bsmlc = BSML Libraries Combined
   # bsmls = BSML Libraries Separate

   if (method=="bsmls")  { 
       est <- bsmls(y,baseslist,maxbas=maxbas,
					  sub.maxbas=sub.maxbas,backward=backward,nobascheck=nobascheck,
					  approach=bsml.approach,control=control)
   }


   if (method=="bsmlc")  {
       est <- bsmlc(y,baseslist,maxbas=maxbas,
					  sub.maxbas=sub.maxbas,backward=backward,nobascheck=nobascheck,
					  approach=bsml.approach,control=control)
   }


   if (method=="has")  {
       est <- has(y,baseslist,maxbas=maxbas,
					  backward=backward,nobascheck=nobascheck,
					  control=control)
   }

   est$call=match.call()
   return(est)

}


control.bsml <- function(has.control=list(crit="GCV",idf=1.2,ridge=F),
						cv.control=list(cvfold=10),
						covpen.control=list(nboot=100,tauhat="Rice",crit="CIC",ct=2),
						gdf.control=list(nrep=100,tauhat="Rice",crit="GCV",ct=2),
						gcvpk.control=list(delta=.1,job=1005,lamlim=c(-6,2),tau=10,ntbl=100,maxtbl=100),
						dpc.control=list(c1=NA,c2=NA)
						) {
						list(has.control=has.control,cv.control=cv.control,covpen.control=covpen.control,gdf.control=gdf.control,gcvpk.control=gcvpk.control,dpc.control=dpc.control)
						}


assign.control <- function(control) {		
ctrl.vals <- control.bsml()
if (!is.null(control)) {
	if (!is.list(control)) stop("Control options must be entered as a list. \n")
	for (i in 1:length(control)) {
		if (!is.list(control[[i]])) stop("Each control option must be entered as a list. \n")
		for (j in names(control[[i]])) {
		flag=T
			for (k in 1:6) {
			if (names(control)[i]==names(ctrl.vals)[k]) {
				for (h in names(ctrl.vals[[k]])) {
				if (j==h) {
					ctrl.vals[[k]][[sum((h==names(ctrl.vals[[k]]))*(1:length(ctrl.vals[[k]])))]]=control[[i]][[sum((h==names(control[[i]]))*(1:length(control[[i]])))]]
					flag=F
					}
				}
			}
			}
		if (flag) warning(paste("Option ",j," in group ",names(control)[i]," is invalid.","\n"))
		}
	}
}
return(ctrl.vals)
}


print.bsml <- function(x, ...)
{
Call=x$call
Total_Number_of_Bases=x$nb
Chosen_Bases_List="None. Only null bases are kept in the final model.\n"
if (!is.null(x$chosen.bases.trim)) Chosen_Bases_List=as.matrix(x$chosen.bases.trim)
Generalized_Degree_of_Freedom=x$dof
Estimated_Variance=x$sigma_sq
cat("Call:\n")
print(Call)
cat("\nTotal Number of Bases in Model: (including null bases)\n")
print(Total_Number_of_Bases)
cat("\nChosen Bases List: (excluding null bases)\n")
print(Chosen_Bases_List)
cat("\nGeneralized Degree of Freedom:\n")
print(Generalized_Degree_of_Freedom)
cat("\nEstimated Variance:\n")
print(Estimated_Variance)
}


print.summary.bsml <- function(x, ...)
{
Call=x$call
Total_Number_of_Bases=x$nb
Generalized_Degree_of_Freedom=x$dof
Standard_Error=x$sigma
cat("Call:\n")
print(Call)
cat("\nTotal Number of Bases in Model:\n")
print(Total_Number_of_Bases)
cat("\nGeneralized Degree of Freedom:\n")
print(Generalized_Degree_of_Freedom)
cat("\nStandard Error:\n")
print(Standard_Error)
}


summary.bsml <- function(object, ...)
{
bsuml <- list(call=object$call, nb=object$nb, dof=object$dof, sigma=sqrt(object$sigma_sq))
class(bsuml) <- "summary.bsml"
return(bsuml)
}


predict.bsml <- function(object, bases.include=NULL, new.baseslist=NULL, confint=F, alpha=0.05, alg="a", addbas="cv", bootrep=200, ...)
{
if (is.null(object$chosen.bases.matrix)) warning("Null model is the best model.","\n")
baseslist=new.baseslist
if (is.null(baseslist)) baseslist=object$baseslist
if (!is.list(baseslist)) stop("Basis functions must be entered as a list.","\n")   
if (length(baseslist)!=length(object$lib.size)) stop("Number of libraries different from that used in the original fit.","\n")
for (i in 1:length(baseslist)) {
  if (!is.null(baseslist[[i]])) {
    if (ncol(as.matrix(baseslist[[i]]))!=object$lib.size[i]) warning("Size of new library ", i, " different from original.", "\n")
  }
}

null.lib=object$null.bases.matrix
vec.len=length(object$lib.size)-1
coef.chosen=object$coef
coef.chosen[is.na(coef.chosen)]=0
coef.null=0
if (is.null(bases.include)) {
bases.include=rep(0,length(object$lib.size))
bases.include[unique(object$chosen.bases.trim[,1])]=1
bases.include=bases.include[-1]
if (!is.null(object$null.bases.matrix)) bases.include=c(rep(1,object$lib.size[1]),bases.include)
}
include.chosen=bases.include
include.null=0
if (!is.null(null.lib)) {
nnull=ncol(null.lib)
vec.len=nnull+vec.len
coef.null=coef.chosen[1:nnull]
coef.chosen=coef.chosen[(nnull+1):length(coef.chosen)]
include.null=bases.include[1:nnull]
include.chosen=include.chosen[(nnull+1):length(include.chosen)]
}
if (length(bases.include)!=vec.len) stop("bases.include has the wrong length.","\n")

fit.null<-0
if (sum(include.null)>0) {
nu.co=coef.null[include.null==1]
if (is.null(baseslist[[1]])) stop("Null library is empty.", "\n")
nu.mat=baseslist[[1]][,include.null==1]
fit.null=nu.mat%*%cbind(nu.co)
}
if (is.null(object$chosen.bases.matrix)) return(fit.null)

sep.coef<-sep.bases<-sep.matrix<-list()
chosen=object$chosen.bases.trim
uni.lib=unique(chosen[,1])
if (ncol(chosen)==1) {
uni.lib=chosen[1]
sep.bases[[uni.lib]]=chosen[2]
if (!is.null(baseslist[[uni.lib]])) sep.matrix[[uni.lib]]=baseslist[[uni.lib]][,sep.bases[[uni.lib]]]
sep.coef[[uni.lib]]=coef.chosen
}
if (ncol(chosen)>1) {
for (i in uni.lib) {
sep.bases[[i]]=chosen[chosen[,1]==i,2]
if (!is.null(baseslist[[i]])) sep.matrix[[i]]=baseslist[[i]][,sep.bases[[i]]]
sep.coef[[i]]=coef.chosen[chosen[,1]==i]
}
}

fit<-fit.chosen<-flag<-0
if (sum(include.chosen)>0) {
lib.pick=2:length(object$lib.size)
lib.pick=c(lib.pick[include.chosen==1])
if (is.null(baseslist[[(lib.pick[1])]])) stop("Library ", lib.pick[1], " is empty.", "\n")
fit.chosen=rep(0,nrow(baseslist[[(lib.pick[1])]]))
for (i in lib.pick) {
	flag=1
	if (length(sep.matrix)<i || is.null(sep.matrix[[i]])) {
	warning("Library ", i, " is empty or its bases not chosen.", "\n")
	flag=0
	}
	if (flag==1) fit.chosen=fit.chosen+sep.matrix[[i]]%*%cbind(sep.coef[[i]])
	}
}


fit=fit.null+fit.chosen

if (confint==F) return(fit)

unbiasedsigma=sqrt(object$sigma_sq)
if (object$approach=="gdf") control=list(gdf.control=list(tauhat=0.5*unbiasedsigma,nrep=bootrep))
else if (object$approach=="covpen") control=list(covpen.control=list(tauhat=0.5*unbiasedsigma,nboot=bootrep))
else control=NULL
ctrl.vals <- object$control
if (!is.null(control)) {
	if (!is.list(control)) stop("Control options must be entered as a list. \n")
	for (i in 1:length(control)) {
		if (!is.list(control[[i]])) stop("Each control option must be entered as a list. \n")
		for (j in names(control[[i]])) {
		flag=T
			for (k in 1:6) {
			if (names(control)[i]==names(ctrl.vals)[k]) {
				for (h in names(ctrl.vals[[k]])) {
				if (j==h) {
					ctrl.vals[[k]][[sum((h==names(ctrl.vals[[k]]))*(1:length(ctrl.vals[[k]])))]]=control[[i]][[sum((h==names(control[[i]]))*(1:length(control[[i]])))]]
					flag=F
					}
				}
			}
			}
		if (flag) warning(paste("Option ",j," in group ",names(control)[i]," is invalid.","\n"))
		}
	}
}

if (alg=="a") {
unbiasedfit=object$fit
initfit=bsml(y=unbiasedfit,baseslist=object$baseslist,maxbas=object$maxbas,sub.maxbas=object$sub.maxbas,backward=object$backward,nobascheck=object$nobascheck,bsml.approach=object$approach,control=ctrl.vals,method=object$method)
afits=initfit$allfits
alist<-allfits<-NULL
if (initfit$method=="bsmls") {
	for (i in 1:dim(initfit$allisb)[3]) {
	modelsize=nrow(object$chosen.bases.trim)
	tofit=rbind(initfit$chosen.bases.trim[-modelsize,],initfit$allisb[,,i][modelsize,])
	alist[[i]]=tofit
	onefit=pdt(initfit,bases.include,new.baseslist,tofit,(unbiasedfit+initfit$ptb[,i]))
	allfits=cbind(allfits,onefit)
	}
	lower=apply(allfits,1,quantile,prob=alpha/2,na.rm=T)
	upper=apply(allfits,1,quantile,prob=1-alpha/2,na.rm=T)
} else {
	for (i in 1:dim(initfit$allisb)[3]) {
	modelsize=nrow(object$chosen.bases.trim)
	tofit=initfit$allisb[,,i][1:modelsize,]
	alist[[i]]=tofit
	onefit=pdt(initfit,bases.include,new.baseslist,tofit,(unbiasedfit+initfit$ptb[,i]))
	allfits=cbind(allfits,onefit)
	}
	lower=apply(allfits,1,quantile,prob=alpha/2,na.rm=T)
	upper=apply(allfits,1,quantile,prob=1-alpha/2,na.rm=T)
}
} else if (alg=="b") {
if (is.numeric(addbas)) {
bigmodelsize=min(object$maxbas-object$lib.size[1],nrow(object$chosen.bases.trim)+addbas)
unbiasedfit=pdt(object,NULL,NULL,object$chosen.bases.full[1:bigmodelsize,],object$y)
}
else if (addbas=="cv") unbiasedfit=bsml(y=object$y,baseslist=object$baseslist,maxbas=object$maxbas,sub.maxbas=object$sub.maxbas,backward=object$backward,nobascheck=object$nobascheck,bsml.approach="cv",control=ctrl.vals,method=object$method)$fit
initfit=bsml(y=unbiasedfit,baseslist=object$baseslist,maxbas=object$maxbas,sub.maxbas=object$sub.maxbas,backward=object$backward,nobascheck=object$nobascheck,bsml.approach=object$approach,control=ctrl.vals,method=object$method)
afits=initfit$allfits
alist<-allfits<-NULL
if (initfit$method=="bsmls") {
	for (i in 1:dim(initfit$allisb)[3]) {
	fitall=initfit$allfits[,,i]
	alldiff=NULL
	for (j in ncol(fitall)) alldiff=cbind(alldiff,fitall[,j]-unbiasedfit)
	scores=apply(alldiff^2,2,sum)/length(object$y)
	scores=scores[(object$lib.size[1]+1):length(scores)]
	modelsize=order(scores)[1]
	tofit=rbind(initfit$chosen.bases.full[1:(modelsize-1),],initfit$allisb[,,i][modelsize,])
	alist[[i]]=tofit
	onefit=pdt(initfit,bases.include,new.baseslist,tofit,(unbiasedfit+initfit$ptb[,i]))
	allfits=cbind(allfits,onefit)
	}
	lower=apply(allfits,1,quantile,prob=alpha/2,na.rm=T)
	upper=apply(allfits,1,quantile,prob=1-alpha/2,na.rm=T)
} else {
	for (i in 1:dim(initfit$allisb)[3]) {
	fitall=initfit$allfits[,,i]
	alldiff=NULL
	for (j in ncol(fitall)) alldiff=cbind(alldiff,fitall[,j]-unbiasedfit)
	scores=apply(alldiff^2,2,sum)/length(object$y)
	scores=scores[(object$lib.size[1]+1):length(scores)]
	modelsize=order(scores)[1]
	tofit=initfit$allisb[,,i][1:modelsize,]
	alist[[i]]=tofit
	onefit=pdt(initfit,bases.include,new.baseslist,tofit,(unbiasedfit+initfit$ptb[,i]))
	allfits=cbind(allfits,onefit)
	}
	lower=apply(allfits,1,quantile,prob=alpha/2,na.rm=T)
	upper=apply(allfits,1,quantile,prob=1-alpha/2,na.rm=T)
}
} else {
warning(paste("Algorithm ",alg," does NOT exist!","\n"))
lower<-upper<-allfits<-alist<-afits<-NA
}

return(list(fit=fit,lower=lower,upper=upper,allfits=allfits,alist=alist,afits=afits))
}


pdt <- function(object, bases.include=NULL, new.baseslist=NULL, bases.chosen,y)
{
baseslist=new.baseslist
if (is.null(baseslist)) baseslist=object$baseslist
options(warn=-1)
bc=bases.chosen
fm=object$null.bases.matrix
if (ncol(bc)>1) {
  for (i in 1:nrow(bc)) {
    fm=cbind(fm,object$baseslist[[bc[i,1]]][,bc[i,2]])
  }
}
if (ncol(bc)==1) {
  fm=cbind(fm,object$baseslist[[bc[1]]][,bc[2]])
}
new.obj=glm(y~fm-1)

null.lib=object$null.bases.matrix
vec.len=length(object$lib.size)-1
coef.chosen=new.obj$coef
coef.chosen[is.na(coef.chosen)]=0
coef.null=0
if (is.null(bases.include)) {
bases.include=rep(0,length(object$lib.size))
bases.include[unique(bases.chosen[,1])]=1
bases.include=bases.include[-1]
if (!is.null(object$null.bases.matrix)) bases.include=c(rep(1,object$lib.size[1]),bases.include)
}
include.chosen=bases.include
include.null=0
if (!is.null(null.lib)) {
nnull=ncol(null.lib)
vec.len=nnull+vec.len
coef.null=coef.chosen[1:nnull]
coef.chosen=coef.chosen[(nnull+1):length(coef.chosen)]
include.null=bases.include[1:nnull]
include.chosen=include.chosen[(nnull+1):length(include.chosen)]
}
if (length(bases.include)!=vec.len) stop("bases.include has the wrong length.","\n")

fit.null<-0
if (sum(include.null)>0) {
nu.co=coef.null[include.null==1]
if (is.null(baseslist[[1]])) stop("Null library is empty.", "\n")
nu.mat=baseslist[[1]][,include.null==1]
fit.null=nu.mat%*%cbind(nu.co)
}
if (is.null(object$chosen.bases.matrix)) return(fit.null)

sep.coef<-sep.bases<-sep.matrix<-list()
chosen=bases.chosen
uni.lib=unique(chosen[,1])
if (ncol(chosen)==1) {
uni.lib=chosen[1]
sep.bases[[uni.lib]]=chosen[2]
if (!is.null(baseslist[[uni.lib]])) sep.matrix[[uni.lib]]=baseslist[[uni.lib]][,sep.bases[[uni.lib]]]
sep.coef[[uni.lib]]=coef.chosen
}
if (ncol(chosen)>1) {
for (i in uni.lib) {
sep.bases[[i]]=chosen[chosen[,1]==i,2]
if (!is.null(baseslist[[i]])) sep.matrix[[i]]=baseslist[[i]][,sep.bases[[i]]]
sep.coef[[i]]=coef.chosen[chosen[,1]==i]
}
}

fit<-fit.chosen<-flag<-0
if (sum(include.chosen)>0) {
lib.pick=2:length(object$lib.size)
lib.pick=c(lib.pick[include.chosen==1])
if (is.null(baseslist[[(lib.pick[1])]])) stop("Library ", lib.pick[1], " is empty.", "\n")
fit.chosen=rep(0,nrow(baseslist[[(lib.pick[1])]]))
for (i in lib.pick) {
	flag=1
	if (length(sep.matrix)<i || is.null(sep.matrix[[i]])) {
	warning("Library ", i, " is empty or its bases not chosen.", "\n")
	flag=0
	}
	if (flag==1) fit.chosen=fit.chosen+sep.matrix[[i]]%*%cbind(sep.coef[[i]])
	}
}

fit=fit.null+fit.chosen
options(warn=1)
return(fit)
}


stdz <- function(x) {
x <- as.matrix(x)
check.zero <- 1:ncol(x)
for (j in 1:ncol(x)) {
  check.zero[j] <- sum(abs(x[,j]))       
#  if (check.zero[j]==0) cat("Column", j, "contains all zeros.","\n")
}
if (min(check.zero)==0) warning("At least one basis function contains all zeros.","\n")     
n=nrow(x)
numcol=ncol(x)
stdz.proc <- .Fortran("STDZ",
	X=as.double(x),
	N=as.integer(n),
	NUMCOL=as.integer(numcol)
	)[c("X")]
x=matrix(stdz.proc$X,n,numcol)
return(x)
}


bas.check <- function(baseslist) {
	basel <- liblib <- basbas <- selex <- NULL
	selindex <- list()
    starter <- ifelse(is.null(baseslist[[1]]),2,1)
    nnull <- ifelse(is.null(baseslist[[1]]),0,ncol(as.matrix(baseslist[[1]])))
	lib.size <- rep(0,length(baseslist))
    for (i in starter:length(baseslist)) {
		lib.size[i] <- ncol(as.matrix(baseslist[[i]]))
	}
    for (i in starter:length(baseslist)) {
	basel <- cbind(basel,baseslist[[i]])
	selindex[[i]]<-rep(0,lib.size[i])
    liblib <- c(liblib,rep(i,ncol(as.matrix(baseslist[[i]]))))
    basbas <- c(basbas,1:ncol(as.matrix(baseslist[[i]])))
	}
	basel <- stdz(basel)
	n <- nrow(as.matrix(basel))
	numcol <- ncol(as.matrix(basel))
	selex=rep(0,numcol)
	for (j in 1:numcol) {
	 if ((var(basel[,j])==0) && (j>nnull)) {
	  selex[j]=1
	 }
	}
	bascheck.proc <- .Fortran("BASCHECK",
							  BASEL=as.double(basel),
							  SELEX=as.integer(selex),
							  N=as.integer(n),
							  NUMCOL=as.integer(numcol)
							 )[c("SELEX")]
	selex=bascheck.proc$SELEX
	if (sum(selex)>0) {
	flg=which(selex>0)
		for (k in 1:length(flg)) {
			selindex[[liblib[flg[k]]]][basbas[flg[k]]]=1
		}
	}
	for (j in 1:numcol) {
	 if (selindex[[liblib[j]]][basbas[j]]==1) {
	  warning("Library ", liblib[j], " Column ", basbas[j], " is zero or constant or repeated.","\n")
	 }
	}
	return(selindex)
}


###########################
# Hybrid Adaptive Splines #
###########################


has <- function(y,baseslist,maxbas=30,
				backward=NULL,nobascheck=F,
                control=NULL) 
{
	ctrl.vals <- assign.control(control)
	has.control=ctrl.vals$has.control
	gdf.control=ctrl.vals$gdf.control
	gcvpk.control=ctrl.vals$gcvpk.control

	#### Check if the libraries are entered as a list ####
    if (!is.list(baseslist)) 
      stop("Basis functions must be entered as a list.","\n")   

	nobs <- length(y)
    starter <- ifelse(is.null(baseslist[[1]]),2,1)
    for (i in starter:length(baseslist)) baseslist[[i]] <- as.matrix(baseslist[[i]])
	
	# Check if the number of rows of all libraries are the same
	for (i in starter:length(baseslist)) {
	lib.dim=nrow(as.matrix(baseslist[[i]]))
	if (lib.dim!=nobs) stop("Library ", i, " has the wrong number of rows.", "\n")
	}

	crit=has.control$crit
	idf=has.control$idf
	ridge=has.control$ridge
    coefficients <- coef <- residuals <- fitted.values <- index.bas <- isb <- nb <- NULL
    fits <- y
    nrep <- gdf.control$nrep
    tauhat <- gdf.control$tauhat
    vmu <- 1
    control <- gcvpk.control
    nnull <- ifelse(is.null(baseslist[[1]]),0,ncol(baseslist[[1]]))
    bas <- basel <- liblib <- basbas <- iwk <- NULL
	selindex <- list()
    for (i in 1:length(baseslist)) {
    bas <- cbind(bas,baseslist[[i]])
    }
	nbas <- ncol(bas)
	
    # Double-checking if total number of bases is greater than maxbas
    totalbases <- 0
	lib.size <- rep(0,length(baseslist))
    for (i in starter:length(baseslist)) {
		lib.size[i] <- ncol(baseslist[[i]])
		totalbases <- totalbases+ncol(baseslist[[i]])
	}
    if (totalbases<maxbas) stop("Total number of basis functions in all libraries is less than MAXBAS.","\n")  
	
    # Check if any two columns (bases) among the libraries are the same
    for (i in starter:length(baseslist)) selindex[[i]]<-rep(0,lib.size[i])
	if (nobascheck==F) {
	options(warn=-1)
	selindex=bas.check(baseslist)
	options(warn=1)
	}
	for (i in 1:length(baseslist)) iwk = c(iwk,selindex[[i]])
	if (sum(iwk)>0) warning("At least on basis is zero or constant or repeated. Use bas.check() to find out which ones.","\n")
	if (abs(sum(iwk-1))<maxbas) stop("Total number of UNIQUE basis functions is less than MAXBAS.","\n")
	
    if (is.na(idf)) {
      if (tauhat == "Rice") tauhat <- sqrt(mean(diff(y)^2)/2)
	else if (tauhat == "Robust") {
	  gaps=abs(diff(y))
	  gaps.lb=quantile(gaps,0.025)
	  gaps.ub=quantile(gaps,0.975)
	  gaps=gaps*(1*(gaps>=gaps.lb & gaps<=gaps.ub))
	  tauhat <- sqrt(mean((gaps[gaps>0])^2)/2)
	}
      idf.all <- bf.cost2(y,bas,iwk,nnull=nnull,maxbas=maxbas,wkamax=maxbas,nrep=nrep,tauhat=tauhat)$idfhat
	idf <- mean(idf.all)
    }

    nf <- sys.nframe()
    ctl <- eval(control,sys.frame(which=nf))
    if (crit=="GCV") tmpcrit <- 1
    if (crit=="AIC") tmpcrit <- 2
    if (crit=="BIC") tmpcrit <- 3
    if (crit=="Cp")  tmpcrit <- 4
    cond <- 1
    if (is.null(backward) & ridge==T) cond <- 2

    isb <- coef <- score <- rss <- seq(0, maxbas)
    yy <- hasfit <- adiag <- svals <- prefits <- y 
    hasdev <- hasdf <- nb <- info <- 0
    iwork <- rep(1,2*maxbas-nnull)
    baswk <- bas 
    wka <- des <- predes <- matrix(0,nobs,maxbas) 
    baswork <- matrix(0,nobs,ncol(bas))
    sigma <- matrix(0,maxbas,maxbas) 
    tbl <- matrix(0,ctl$maxtbl,3)
    lwa <- (maxbas-nnull)*(maxbas-2*nnull+2+nobs)+maxbas+nobs 
    work <- rep(0,lwa)
    fits <- u <- h <- z <- zk <- rep(0,nobs)
    fitswk <- matrix(0,nobs,maxbas)

    hasfit <- .Fortran("has1",
            NOBS = as.integer(nobs),
            Y = as.double(y),
            BAS = as.double(bas),
            NNULL = as.integer(nnull),
            NBAS = as.integer(nbas),
            MAXBAS = as.integer(maxbas),
            CRIT = as.integer(tmpcrit),
            IDF = as.double(idf),
            DELTA = as.double(ctl$delta), 
            COND = as.integer(cond),
            NB = as.integer(nb),
            ISB = as.integer(isb),
            COEF = as.double(coef),
            HASFIT = as.double(hasfit),
            HASDEV = as.double(hasdev),
            HASDF = as.double(hasdf),
            SCORE = as.double(score),
            RSS = as.double(rss),
            INFO = as.integer(info),
            LAMLIM = as.double(ctl$lamlim),
            MAXTBL = as.integer(ctl$maxtbl),
            NTBL = as.integer(ctl$ntbl),
            JOB = as.integer(ctl$job),
            TAU = as.double(ctl$tau),
            IWK = as.integer(iwk), 
            IWORK = as.integer(iwork), 
            YY = as.double(yy), 
            BASWK = as.double(baswk), 
            WKA = as.double(wka), 
            DES = as.double(des), 
            SIGMA = as.double(sigma), 
            ADIAG = as.double(adiag),
            SVALS = as.double(svals),
            TBL = as.double(tbl),
            WORK = as.double(work),
            VMU = as.integer(vmu),
            U = as.double(u),
            H = as.double(h),
            FITSWK = as.double(fitswk),
            Z = as.double(z),
            ZK = as.double(zk))[c("NB", "ISB", "COEF", 
            "HASFIT", "HASDEV", "HASDF", "SCORE", "RSS", "INFO","ADIAG","IDF")]

    # if (hasfit$INFO > 0) stop("GCVPACK/dsnsm error")

    idf <- hasfit$IDF
    isb <- hasfit$ISB
    isb <- isb[1:(length(isb)-1)]
    nb <- hasfit$NB
    selbas <- bas[,isb[1:nb]]
    A.diag <- hasfit$ADIAG	
	scr <- hasfit$SCORE
	scr <- scr[1:(length(scr)-1)]
	
    if (is.null(backward)) {
      if (ridge == T & nb > nnull) {
        coef <- hasfit$COEF
        fits <- hasfit$HASFIT
        sd <- hasfit$HASDEV
        df <- hasfit$HASDF
        A.diag <- hasfit$ADIAG  
      }
      if (ridge == F | nb == nnull) {
        fit1 <- glm(y~selbas-1)
        fits <- fit1$fitted.values
        coef <- fit1$coef
        sd <- sqrt(fit1$deviance/fit1$df.residual)
        df <- nb
        A.diag <- hat(selbas)
      }     
    }

    if (!is.null(backward) & nb == nnull) {
      fit1 <- glm(y~selbas-1)
      fits <- fit1$fitted.values
      coef <- fit1$coef
      sd <- sqrt(fit1$deviance/fit1$df.residual)
      df <- nb
      A.diag <- hat(selbas)	
    }
	                        
    if (!is.null(backward) & nb > nnull) {

      tmpbas <- bas[,isb[1:nb]]
      cstep <- regsubsets(tmpbas,y,force.in=c(1:nnull),method="backward",intercept=F)

      if (nb==maxbas) isb <- isb[cstep$vorder]
      if (nb<maxbas) isb <- isb[c(cstep$vorder,(nb+1):maxbas)]
      tmp <- cstep$rss[1:nb]
      if (backward == "GCV") scr <- tmp/nobs/(1-(1:nb)/nobs)**2
      if (backward == "AIC") scr <- nobs*log(tmp/nobs) + 2*(1:nb)
      if (backward == "BIC") scr <- nobs*log(tmp/nobs) + log(nobs)*(1:nb)
      if (backward == "Cp") scr <- tmp + 2*(1:nb)*tmp[nb]/(nobs-nb)
      tmp1 <- scr==min(scr)
      cum <- 0         
      cumsum <- 0
      while (cumsum==0) {
        cum <- cum+1
        cumsum <- sum(tmp1[1:cum])
      }
      nb <- cum
      selbas <- bas[,isb[1:nb]]
      if (ridge == F | nb == nnull) {
        fit1 <- glm(y~selbas-1)
        fits <- fit1$fitted.values
        coef <- fit1$coef
        sd <- sqrt(fit1$deviance/fit1$df.residual)
        df <- nb
      }
      if (ridge == T & nb > nnull) {
        pen.mat <- bdiag(list(matrix(0,nnull,nnull),selbas[1:(nb-nnull),(nnull+1):nb]+ctl$delta*diag(rep(1,nb-nnull))))
        fit1 <- snsm1(selbas, y, S=pen.mat)
        fits <- fit1$fitted.values
        pen.info <- fit1$penalty.info
        sd <- sqrt(pen.info[3]/(nobs-(nobs-pen.info[4]-nnull)*idf-nnull))
        df <- nobs - pen.info[4]
        coef <- fit1$coef
	  A.diag <- fit1$adiag
      }
    }

    # Identify the library index and basis index of the chosen bases
    chosen.bases.index <- chosen.bases.matrix <- chosen.bases.full <- chosen.bases.trim <- NULL
    if (nnull==0) index.sel <- isb else index.sel <- isb[-(1:nnull)]-nnull
    L <- length(baseslist)
    liblen <- libnum <- NULL
    for (j in 2:L) {
      liblen <- c(liblen,ncol(as.matrix(baseslist[[j]])))
      libnum <- c(libnum,1:ncol(as.matrix(baseslist[[j]])))
    }
    libindex <- rep(2:L,liblen)
    lib.num <- libindex[index.sel]
    bas.num <- libnum[index.sel]  
    chosen.bases.full <- cbind(lib.num,bas.num)
    if (nb>nnull) {
	chosen.bases.trim <- cbind(chosen.bases.full[1:(nb-nnull),])
    chosen.bases.matrix <- selbas[,(nnull+1):ncol(selbas)]
    }
	if (!is.null(backward)) score <- scr[nnull:nb] else score <- scr[nnull:length(scr)]
	
	# Summary information
	residuals <- y-fits
	sigma_sq <- sum(residuals^2)/(nobs-idf*nb)
	dof <- idf*nb
	
	if (!is.null(chosen.bases.full)) chosen.bases.full=as.matrix(chosen.bases.full)
	if (!is.null(chosen.bases.trim)) chosen.bases.trim=as.matrix(chosen.bases.trim)
	if (!is.null(chosen.bases.matrix)) chosen.bases.matrix=as.matrix(chosen.bases.matrix)
	est <- list(coefficients=coef,residuals=residuals,fitted.values=fits,chosen.bases.full=chosen.bases.full,chosen.bases.trim=chosen.bases.trim,chosen.bases.matrix=chosen.bases.matrix,null.bases.matrix=baseslist[[1]],lib.size=lib.size,score=score,gdf=NA,idf=idf,nb=nb,sigma_sq=sigma_sq,dof=dof,index.bas=isb,y=y,baseslist=baseslist,maxbas=maxbas,backward=backward,nobascheck=nobascheck,control=ctrl.vals)
	est$method="has"
	est$call <- match.call()
	class(est) <- "bsml"
    return(est)

}

    # Function for ridge regression
    snsm1 <- function(X, Y, method = 1, S = diag(ncol(X))) {
    # op=0 using GML; op= 1 using GCV; op=2 using UBR
    call <- match.call()
    n <- dim(S)[1]
    storage.mode(S) <- "double"
    rk.penalty.matrix <- .Fortran("dchdc", x = S , as.integer(n), as.integer(n), double(n), piv = as.integer(rep(0,n)), as.integer(1), rank = integer(1))$rank

    #print(rk.penalty.matrix)
        npar <- ncol(X)
        nobs <- nrow(X)
        tau <- 1
        lamlim <- c(-20, 20)
        ldtbl <- ntbl <- 1000
        tbl <- matrix(rep(0, ldtbl * 3), ncol = 3)
        auxtbl <- matrix(rep(0, 9), ncol = 3)
        nnull <- npar - rk.penalty.matrix
        lwa <- (npar - nnull) * (npar - 2 * nnull + 2 + nobs) + npar + nobs
        liwa <- 2 * npar - nnull
        job <- 1005    # tbl GCV data,
        z <- .Fortran("dsnsm1",
                X = as.double(X),
                LdX = as.integer(nrow(X)),
                Y = as.double(Y),
                S = as.double(S),
                LdS = as.integer(nrow(S)),
                nobs = as.integer(nobs),
                npar = as.integer(npar),
                nnull = as.integer(nnull),
                adiag = double(nobs),
                tau = as.double(tau),
                lamlim = as.double(lamlim),
                ntbl = as.integer(ntbl),
                penalty.info = double(5),
                iout = integer(3),
                coef = double(npar),
                svals = double(npar - nnull),
                gcv.info = as.double(tbl),
                Ldtbl = as.integer(ldtbl),
                auxtbl = as.double(auxtbl),
                iwork = integer(liwa),
                liwa = as.integer(liwa),
                work = double(lwa),
                lwa = as.integer(lwa),
                job = as.integer(job),
                info = integer(1),
                vm = as.integer(method))

        return(list(coefficients = z$coef, residuals = (Y - z$Y),
              fitted.values = z$Y, penalty.matrix = S,
              penalty.info = z$penalty.info, gcv.info = z$gcv.info,
              aux.info = z$aux.info, call = call,
              info = z$info, adiag = z$adiag))
    }

######################################################
# Basis Selection from Multiple Libraries - Combined #
######################################################


bsmlc <- function(y,baseslist,maxbas=30,
				  sub.maxbas=c(),backward=NULL,nobascheck=F,
                  approach="gdf",control=NULL) 
{
	ctrl.vals <- assign.control(control)
	covpen.control=ctrl.vals$covpen.control
	gdf.control=ctrl.vals$gdf.control
	gcvpk.control=ctrl.vals$gcvpk.control
	cv.control=ctrl.vals$cv.control
	cvfold=cv.control$cvfold

	#### Check if the libraries are entered as a list ####
    if (!is.list(baseslist)) 
      stop("Basis functions must be entered as a list.","\n")   
	ct <- ifelse(approach=="gdf",gdf.control$ct,covpen.control$ct)
	n <- nobs <- length(y)
    starter <- ifelse(is.null(baseslist[[1]]),2,1)
    for (i in starter:length(baseslist)) baseslist[[i]] <- as.matrix(baseslist[[i]])
	
	# Check if the number of rows of all libraries are the same
	for (i in starter:length(baseslist)) {
	lib.dim=nrow(as.matrix(baseslist[[i]]))
	if (lib.dim!=nobs) stop("Library ", i, " has the wrong number of rows.", "\n")
	}

	allisb <- allisbraw <- coefficients <- coef <- residuals <- res <- fitted.values <- fitted <- score <- score.null <- gdf <- gdf.new <- idf <- idf.new <- nb <- NULL
	chosen.bases.full <- chosen.bases.trim <- chosen.bases.matrix <- NULL
	score <- rep(0,maxbas)
	repp <- ifelse(approach=="gdf",gdf.control$nrep,covpen.control$nboot)
	tauhat <- ifelse(approach=="gdf",gdf.control$tauhat,covpen.control$tauhat)
	crit <- ifelse(approach=="gdf",gdf.control$crit,covpen.control$crit)
	nnull <- ifelse(is.null(baseslist[[1]]),0,ncol(baseslist[[1]]))
    for (i in 2:length(baseslist)) colnames(baseslist[[i]]) <- as.character(1:ncol(baseslist[[i]]))
    allfits <- array(rep(0,nobs*maxbas*gdf.control$nrep),dim=c(nobs,maxbas,gdf.control$nrep))
	ptb <- matrix(rep(0,nobs*gdf.control$nrep),nobs,gdf.control$nrep)

    # If a subset of bases of length sub.maxbas is requested for 
    # the libraries, then HAS is used to identify the best subsets of knots.

    if (!is.null(sub.maxbas)) {
	  orig.index <- matrix(length(sub.maxbas)*max(sub.maxbas),max(sub.maxbas),length(sub.maxbas))
      for (i in 1:length(sub.maxbas)) {
        if (sub.maxbas[i] < ncol(baseslist[[i+1]])) {
          has.obj <- has(y,list(baseslist[[1]],baseslist[[i+1]]),maxbas=nnull+sub.maxbas[i],control=list(gdf.control=gdf.control,gcvpk.control=gcvpk.control))
          baseslist[[i+1]] <- baseslist[[i+1]][,(has.obj$index.bas[(nnull+1):length(has.obj$index.bas)]-nnull)]
	   orig.index[,i]=as.numeric(colnames(baseslist[[i+1]]))
        }
      }
    }

    bas <- iwk <- NULL
	selindex <- list()
    for (i in 1:length(baseslist)) bas <- cbind(bas,baseslist[[i]])

	# Double-checking if total number of subknots chosen is greater than maxbas
    totalbases <- 0
	lib.size <- rep(0,length(baseslist))
    for (i in starter:length(baseslist)) {
		lib.size[i] <- ncol(baseslist[[i]])
		totalbases <- totalbases+ncol(baseslist[[i]])
	}
    if (totalbases<maxbas) stop("Total number of basis functions in all libraries is less than MAXBAS.","\n")  

    # Check if any two columns (bases) among the libraries are the same
    for (i in starter:length(baseslist)) selindex[[i]]<-rep(0,lib.size[i])
	if (nobascheck==F) {
	options(warn=-1)
	selindex=bas.check(baseslist)
	options(warn=1)
	}
	for (i in 1:length(baseslist)) iwk = c(iwk,selindex[[i]])
	if (sum(iwk)>0) warning("At least on basis is zero or constant or repeated. Use bas.check() to find out which ones.","\n")
	if (abs(sum(iwk-1))<maxbas) stop("Total number of UNIQUE basis functions is less than MAXBAS.","\n")
	
    # Compute Rice Estimator of standard deviation

    if (tauhat == "Rice") {
	  tauhat <- sqrt(mean(diff(y)^2)/2)
	  est.var=tauhat^2
	}
	else if (tauhat == "Robust") {
	  gaps=abs(diff(y))
	  gaps.lb=quantile(gaps,0.025)
	  gaps.ub=quantile(gaps,0.975)
	  gaps=gaps*(1*(gaps>=gaps.lb & gaps<=gaps.ub))
	  tauhat <- sqrt(mean((gaps[gaps>0])^2)/2)
	  est.var=tauhat^2
	}
    else if (tauhat == "MAD") {
      tauhat <- median(abs(diff(y)))
	  est.var=tauhat^2
    }
	else if (is.numeric(tauhat)) {
	  est.var=tauhat^2
	}
#	if (approach=="covpen") tautilde=2*tauhat
#	else tautilde=tauhat

    startpt <- ifelse(is.null(baseslist[[1]]),1,nnull)
    tmp <- bf.cost2(y,bas,iwk,nnull=nnull,maxbas=maxbas,wkamax=maxbas,nrep=repp,tauhat=tauhat,approach=approach)
	allfits <- tmp$fits
	ptb <- tmp$ptb
    gdf <- tmp$gdfhat
	if (approach=="gdf" || approach=="covpen") allisbraw <- as.matrix(tmp$icb)
    if (approach=="gdf") gdf <- gdf[gdf<n]
    gdf.new <- gdf[(nnull+1):length(gdf)]
    idf.new <- tmp$idfhat[(nnull+1):length(gdf)]
    index.sel <- tmp$isb[(nnull+1):length(gdf)]
    rss <- tmp$rss[startpt:length(gdf)]
    gdf <- gdf[startpt:length(gdf)]
    if (approach=="gdf") {
		if (crit=="GCV") score <- rss/(n-gdf)^2
		else if (crit=="CIC") score <- rss+ct*gdf*est.var
		else stop(paste("Criterion ",crit," is not supported!","\n",sep=""))
	}
    else if (approach=="covpen") {
		if (crit=="CIC") score <- rss+ct*gdf
		else if (crit=="GCV") score <- rss/(n-(gdf/est.var))^2
		else stop(paste("Criterion ",crit," is not supported!","\n",sep=""))
	}
	else if (approach=="cv") {
		cvfold=cvfold
		sel.mat0<-score<-NULL
		sel.mat0=bas[,tmp$isb[1:maxbas]]
		for (j in startpt:ncol(sel.mat0)) {
			rand <- sample(n)%%cvfold + 1
			sel.mat=as.matrix(sel.mat0[,1:j])
			mse.pred<-0
			for (i in 1:length(unique(rand))) {
				coef=lm(y[rand!=i]~sel.mat[rand!=i,]-1)$coef
                usemat=sel.mat[rand==i,!is.na(coef)]  # 100802 YW correction (WM) for leave-1-out CV
				usecoef=coef[!is.na(coef)]
                if (length(usecoef)==1) predd=as.vector(usemat*usecoef)   # 100802 YW (WM)
                if (length(usecoef)>1) predd=as.vector(usemat%*%usecoef)  # 100802 YW (WM)

### OR: (YW alternative 100802) # predd=ifelse(rep(length(usecoef)==1,nrow(usemat)),usemat*usecoef,usemat%*%usecoef)   
###     (alternative not tested in sims)
                         
				mse.pred <- mse.pred+sum((predd-y[rand==i])^2)
			}
			score=c(score,mse.pred/n)
		}
	}
    score[is.na(score)] <- max(score[!is.na(score)])
    kk <- min((startpt:(length(score)+startpt-1))[score==min(score)])
    selbas <- bas[,tmp$isb[1:kk]]
    fit1 <- glm(y~selbas-1)
    fitted <- fit1$fitted
    coef <- fit1$coef
    res <- fit1$res

    # Identify the library index and basis index of the chosen bases
    begin <- ifelse(is.null(baseslist[[1]]),2,1)
    L <- length(baseslist)
    liblen <- libnum <- NULL
    for (j in begin:L) {
      liblen <- c(liblen,ncol(baseslist[[j]]))
      libnum <- c(libnum,1:ncol(baseslist[[j]]))
    }
    libindex <- rep(begin:L,liblen)
    lib.num <- libindex[index.sel]
    bas.num <- libnum[index.sel]  
    chosen.bases.full <- cbind(lib.num,bas.num)
	if (kk>nnull) {
    chosen.bases.trim <- cbind(chosen.bases.full[1:(kk-nnull),1],chosen.bases.full[1:(kk-nnull),2])
    chosen.bases.matrix <- selbas[,(nnull+1):kk]
	}
    if (!is.null(sub.maxbas)) {
      for (i in 1:length(sub.maxbas)) {
        libnum.full=chosen.bases.full[,1]
        temp.full=chosen.bases.full[libnum.full==(i+1),2]
        chosen.bases.full[libnum.full==(i+1),2]=orig.index[temp.full,i]
		if (kk>nnull) {
        libnum.trim=chosen.bases.trim[,1]
        temp.trim=chosen.bases.trim[libnum.trim==(i+1),2]
        chosen.bases.trim[libnum.trim==(i+1),2]=orig.index[temp.trim,i]
		}
      }
    }
    nb=kk
#	score=crit

	# Convert all index sequences into two-column array format
	if (approach=="gdf" || approach=="covpen") {
	allisb=array(0,dim=c(nrow(allisbraw)-nnull,2,ncol(allisbraw)))
	for (i in 1:ncol(allisbraw)) {
	begin <- ifelse(is.null(baseslist[[1]]),2,1)
    L <- length(baseslist)
    liblen <- libnum <- NULL
    for (j in begin:L) {
      liblen <- c(liblen,ncol(baseslist[[j]]))
      libnum <- c(libnum,1:ncol(baseslist[[j]]))
    }
    libindex <- rep(begin:L,liblen)
    lib.num <- libindex[allisbraw[,i][(nnull+1):nrow(allisbraw)]]
    bas.num <- libnum[allisbraw[,i][(nnull+1):nrow(allisbraw)]]
	allisb[,,i]=cbind(lib.num,bas.num)
	}
	}
	
	
    # Backward Elimination (output of GDFs and IDF will not be affected)
    if (!is.null(backward) && !is.null(chosen.bases.full)) {
    if (nnull==0) importance=regsubsets(selbas,y,method="backward",intercept=F)
    else importance=regsubsets(selbas,y,force.in=c(1:nnull),method="backward",intercept=F)

    imp.rss <- importance$rss[nnull:kk]
    gdf.trim <- gdf[1:kk]
    gdf.trim <- gdf.trim[importance$vorder[1:kk]]
	gdf.trim <- gdf.trim[nnull:kk]
  
    if (backward == "GCV") scr <- imp.rss/(nobs-(gdf.trim))^2
    if (backward == "AIC") scr <- nobs*log(imp.rss/nobs) + 2*(gdf.trim)
    if (backward == "BIC") scr <- nobs*log(imp.rss/nobs) + log(nobs)*(gdf.trim)
    if (backward == "Cp")  scr <- imp.rss + 2*(gdf.trim)*imp.rss[kk-nnull]/(nobs-(gdf.trim))

    trim.again <- order(scr)[1]-1
	chosen.bases.trim.neworder <- chosen.bases.trim[(importance$vorder[(nnull+1):kk]-nnull),]
    chosen.bases.trim <- chosen.bases.trim[(importance$vorder[(nnull+1):kk]-nnull),][1:trim.again,]
	# The following line controls whether to preserve the original bases order
	chosen.bases.full <- rbind(chosen.bases.trim.neworder,chosen.bases.full[-(1:(nb-nnull)),])
    chosen.bases.matrix <- chosen.bases.matrix[,(importance$vorder[(nnull+1):kk]-nnull)][,1:trim.again]
    final.bases.backward <- cbind(baseslist[[1]],chosen.bases.matrix)
    fit <- glm(y~final.bases.backward-1)
    fitted <- fit$fitted
    coef <- fit$coef
    res <- fit$res
    nb=trim.again+nnull
	score=scr
    }
	
	# Summary information
	if (approach=="gdf") {
	gen_deg_fre <- ifelse(is.null(baseslist[[1]]),gdf.new[nb],gdf.new[(nb-nnull)])
	sigma_sq <- sum(res^2)/(nobs-gen_deg_fre)
	}
	else {
	gen_deg_fre <- ifelse(is.null(baseslist[[1]]),gdf.new[nb],gdf.new[(nb-nnull)])/est.var
	sigma_sq <- sum(res^2)/(nobs-gen_deg_fre)
	}
	if (!is.null(chosen.bases.full)) chosen.bases.full=as.matrix(chosen.bases.full)
	if (!is.null(chosen.bases.trim)) chosen.bases.trim=as.matrix(chosen.bases.trim)
	if (!is.null(chosen.bases.matrix)) chosen.bases.matrix=as.matrix(chosen.bases.matrix)
	est <- list(coefficients=coef,residuals=res,fitted.values=fitted,chosen.bases.full=chosen.bases.full,chosen.bases.trim=chosen.bases.trim,chosen.bases.matrix=chosen.bases.matrix,null.bases.matrix=baseslist[[1]],lib.size=lib.size,score=score,gdf=gdf.new,idf=idf.new,nb=nb,sigma_sq=sigma_sq,dof=gen_deg_fre,allisb=allisb,ptb=ptb,allfits=allfits,y=y,baseslist=baseslist,maxbas=maxbas,sub.maxbas=sub.maxbas,backward=backward,nobascheck=nobascheck,approach=approach,control=ctrl.vals)
	est$method <- "bsmlc"
	est$call <- match.call()
	class(est) <- "bsml"
    return(est)
}

######################################################
# Basis Selection from Multiple Libraries - Separate #
######################################################


bsmls <- function(y,baseslist,maxbas=30,
				  sub.maxbas=c(),backward=NULL,nobascheck=F,
				  approach="gdf",control=NULL) 
{
	ctrl.vals <- assign.control(control)
	covpen.control=ctrl.vals$covpen.control
	gdf.control=ctrl.vals$gdf.control
	gcvpk.control=ctrl.vals$gcvpk.control
	cv.control=ctrl.vals$cv.control
	cvfold=cv.control$cvfold
	dpc.control=ctrl.vals$dpc.control

	#### Check if the libraries are entered as a list ####
  if (!is.list(baseslist)) 
    stop("Basis functions must be entered as a list.","\n")   
  
  c1=dpc.control$c1
  if (is.na(c1)) c1=log(length(y))
  c2=dpc.control$c2
  if (is.na(c2)) c2=2
  ct <- ifelse(approach=="gdf",gdf.control$ct,covpen.control$ct)
  nobs=length(y)
  starter <- ifelse(is.null(baseslist[[1]]),2,1)
  for (i in starter:length(baseslist)) baseslist[[i]] <- as.matrix(baseslist[[i]])
  
  # Check if the number of rows of all libraries are the same
  for (i in starter:length(baseslist)) {
    lib.dim=nrow(as.matrix(baseslist[[i]]))
    if (lib.dim!=nobs) stop("Library ", i, " has the wrong number of rows.", "\n")
  }

  coefficients <- coef <- residuals <- res <- fitted.values <- fitted <- score <- score.null <- gdf <- idf <- nb <- NULL
  chosen.bases.full <- chosen.bases.trim <- chosen.bases.matrix <- NULL
  bas <- basel <- liblib <- basbas <- iwk <- NULL
  selindex <- list()
  n=length(y)
  repp <- ifelse(approach=="gdf",gdf.control$nrep,covpen.control$nboot)
  tauhat <- ifelse(approach=="gdf",gdf.control$tauhat,covpen.control$tauhat)
  crit <- ifelse(approach=="gdf",gdf.control$crit,covpen.control$crit)
  nnull=ifelse(is.null(baseslist[[1]]),0,ncol(baseslist[[1]]))
  updatevar=0
  for (i in 2:length(baseslist)) colnames(baseslist[[i]]) <- as.character(1:ncol(baseslist[[i]]))
  allfits <- array(rep(0,nobs*maxbas*gdf.control$nrep),dim=c(nobs,maxbas,gdf.control$nrep))

  # If a subset of bases of length sub.maxbas is requested for 
  # the libraries, then HAS is used to identify the best subsets of knots.

  if (!is.null(sub.maxbas)) {
	orig.index <- matrix(length(sub.maxbas)*max(sub.maxbas),max(sub.maxbas),length(sub.maxbas))
    for (i in 1:length(sub.maxbas)) {
      if (sub.maxbas[i] < ncol(baseslist[[i+1]])) {
        has.obj <- has(y,list(baseslist[[1]],baseslist[[i+1]]),maxbas=nnull+sub.maxbas[i],control=list(gdf.control=gdf.control,gcvpk.control=gcvpk.control))
        baseslist[[i+1]] <- baseslist[[i+1]][,(has.obj$index.bas[(nnull+1):length(has.obj$index.bas)]-nnull)]
	 orig.index[,i]=as.numeric(colnames(baseslist[[i+1]]))
      }
    }
  }


  for (i in 1:length(baseslist)) bas <- cbind(bas,baseslist[[i]])
  # Double-checking if total number of subknots chosen is greater than maxbas
  totalbases <- 0
  size.of.lib <- rep(0,length(baseslist))
  for (i in starter:length(baseslist)) {
	size.of.lib[i] <- ncol(baseslist[[i]])
	totalbases <- totalbases+ncol(baseslist[[i]])
  }
  if (totalbases<maxbas) stop("Total number of basis functions in all libraries is less than MAXBAS.","\n")  

    # Check if any two columns (bases) among the libraries are the same
    for (i in starter:length(baseslist)) selindex[[i]]<-rep(0,size.of.lib[i])
	if (nobascheck==F) {
	options(warn=-1)
	selindex=bas.check(baseslist)
	options(warn=1)
	}
	for (i in 1:length(baseslist)) iwk = c(iwk,selindex[[i]])
	if (sum(iwk)>0) warning("At least on basis is zero or constant or repeated. Use bas.check() to find out which ones.","\n")
	if (abs(sum(iwk-1))<maxbas) stop("Total number of UNIQUE basis functions is less than MAXBAS.","\n")
	
  #### Rice Estimator of the variance if sigma^2 not supplied ####
  if (tauhat == "Rice") {
  est.var=mean(diff(y)^2)/2
  tauhat=sqrt(est.var)
  }
  if (tauhat == "Robust") {
    gaps=abs(diff(y))
    gaps.lb=quantile(gaps,0.025)
    gaps.ub=quantile(gaps,0.975)
    gaps=gaps*(1*(gaps>=gaps.lb & gaps<=gaps.ub))
    tauhat <- sqrt(mean((gaps[gaps>0])^2)/2)
    est.var=tauhat^2
  }
  if (tauhat == "MAD") {
    tauhat=median(abs(diff(y)))
    est.var=tauhat^2
  }
  if (is.numeric(tauhat)) {
#  s.hat=tauhat
  est.var=tauhat^2
  }
#  if (approach=="covpen") tautilde=s.hat
#  else tautilde=0.5*s.hat

  rss.sel<-bas.gdf<-rep(0,maxbas)
  lib=baseslist
  new.sel=matrix(0,length(y),length(baseslist))
  lib.cost<-cum.bas.cost<-rss<-bas.cost<-lib.size<-dpc<-index.new<-rep(0,length(baseslist))
  lib.index=list()
  for (j in starter:length(baseslist)) {
  lib.size[j]=abs(sum(selindex[[j]]-1))
  lib.index[[j]]=1:ncol(as.matrix(baseslist[[j]]))
  }
  lib.var<-new.source<-sel.bases<-NULL
  lib.num<-bas.num<-rep(0,maxbas)
  chosen.bases.index=cbind(lib.num,bas.num)
  df.var<-gdf<-0
  cp.null.score<-score.null<-NA
  chosen.bases.full=NULL
  allisb <- allisbraw <- libisb <- wka.tmp <- wkall.tmp <- yy.tmp <- yyptb.tmp <- NULL

  if (approach=="gdf") {
  # Pre-evaluation of library cost for each library for up to "maxbas" number of bases selected
  gdf.pre<-idf.pre<-matrix(0,maxbas,length(baseslist)-1)
  for (w in 2:length(baseslist)) {
    iwk=c(selindex[[1]],selindex[[w]])
    pre.eval=bf.cost2(y,cbind(baseslist[[1]],baseslist[[w]]),iwk,nnull=nnull,maxbas=min(maxbas,abs(sum(iwk-1))),wkamax=maxbas,nrep=repp,tauhat=tauhat)
    gdf.pre[,w-1][1:length(pre.eval$gdf)]=pre.eval$gdf
    idf.pre[,w-1][1:length(pre.eval$idf)]=pre.eval$idf
  }
  }

  if (nnull!=0) {
  iwk=selindex[[1]]
  init.eval=bf.cost2(y,baseslist[[1]],iwk,nnull=nnull,maxbas=nnull,wkamax=maxbas,nrep=repp,tauhat=tauhat,approach=approach)
  allfits[,1:nnull,]=init.eval$fits
  wka=init.eval$wka
  wkall=init.eval$wkall
  ptb=init.eval$ptb
  yy=init.eval$yy
  yyptb=init.eval$yyptb
  if (crit=="CIC") cp.null.score=init.eval$rss[nnull]+ct*init.eval$gdf[nnull]
  else if (crit=="GCV") cp.null.score=init.eval$rss[nnull]/(n-init.eval$gdf[nnull]/est.var)^2
  }

  if (nnull==0) {
  iwk=NULL
  for (i in 1:length(baseslist)) iwk = c(iwk,selindex[[i]])
  init.eval=bf.cost2(y,bas,iwk,nnull=nnull,maxbas=1,wkamax=maxbas,nrep=repp,tauhat=tauhat,approach=approach)
  allfits[,1,]=init.eval$fits
  wka=init.eval$wka
  wkall=init.eval$wkall
  ptb=init.eval$ptb
  yy=init.eval$yy
  yyptb=init.eval$yyptb
  begin <- 2
  L <- length(baseslist)
  liblen <- libnum <- NULL
  for (j in begin:L) {
    liblen <- c(liblen,ncol(baseslist[[j]]))
    libnum <- c(libnum,1:ncol(baseslist[[j]]))
  }
  libindex <- rep(begin:L,liblen)
  index.sel <- init.eval$isb
  pick <- lib.num <- libindex[index.sel]
  index.new[pick] <- bas.num <- libnum[index.sel]  
  sel.bases=cbind(sel.bases,as.matrix(lib[[pick]])[,index.new[pick]])
  chosen.bases.index[1,]=c(pick,lib.index[[pick]][index.new[pick]])
  bas.gdf[1]=init.eval$gdf
  rss.sel[1]=init.eval$rss
  lib.size[pick]=0
  if (ncol(lib[[pick]])>1) {
  lib[[pick]]=as.matrix(lib[[pick]])[,-index.new[pick]]
  lib.index[[pick]]=lib.index[[pick]][-index.new[pick]]
  selindex[[pick]]=selindex[[pick]][-index.new[pick]]
  lib.size[pick]=abs(sum(selindex[[pick]]-1))
  }
  }

  # Choose a candidate basis from each library and compute the corresponding basis cost and library cost
  i=ifelse(is.null(baseslist[[1]]),2,1)
  if (approach=="gdf") cond=((i <= (maxbas-nnull)) && (gdf < length(y)))
  else cond=(i <= (maxbas-nnull))
  
  while (cond) {
  # Calculate costs for bases from different libraries
  for (j in 2:length(baseslist)) {
#    dpc[j]=Inf
	rss[j]<-cum.bas.cost[j]<-lib.cost[j]<-Inf	
    if (lib.size[j] >= 1) {
	new.source=cbind(baseslist[[1]],sel.bases,lib[[j]])
	iwk=c(rep(0,ncol(as.matrix(cbind(baseslist[[1]],sel.bases)))),selindex[[j]])
	bas.eval=pick.bas(yy,yyptb,new.source,iwk,wka=wka,wkall=wkall,wkamax=maxbas,ptb=ptb,nnull=nnull+i-1,nrep=repp,tauhat=tauhat,approach=approach)
	libisb[[j]]=bas.eval$icb-(nnull+i-1)
	wka.tmp[[j]]=bas.eval$wka
	wkall.tmp[[j]]=bas.eval$wkall
	yy.tmp[[j]]=bas.eval$yy
	yyptb.tmp[[j]]=bas.eval$yyptb
	cum.bas.cost[j]=bas.eval$gdfhat
	if (approach=="gdf") {
	bas.cost[j]=max(cum.bas.cost[j]-(nnull+i-1),1)
	lib.cost[j]=idf.pre[nnull+ncol(baseslist[[j]])-ncol(as.matrix(lib[[j]]))+1,j-1]
	}
	rss[j]=bas.eval$rss
	new.sel[,j]=new.source[,bas.eval$isb]
    index.new[j]=bas.eval$isb-(nnull+i-1)
    }
  }

  # Update variance estimate (not used)
  if (updatevar==1) {
  lib.var=cbind(baseslist[[1]],sel.bases,new.sel[,-1])
  df.var=nnull+i-1+sum(bas.cost)
  if (df.var>=length(y)) warning("Effective degree of freedom exceeds sample size.","\n")
  lm.var=lm(y~lib.var)
  est.var=sum((y-lm.var$fit)^2/(length(y)-df.var))
#  tauhat=sqrt(est.var)
  }

  if (approach=="gdf") {
  # Compute the Doubly Penalized Criterion (DPC)
  for (j in 2:length(baseslist)) {
    dpc[j]=rss[j]+c1*est.var*(cum.bas.cost[j])+c2*est.var*lib.cost[j]
#	print(c(rss[j],est.var*cum.bas.cost[j],est.var*lib.cost[j]))
  }
  dpc2=dpc[2:length(dpc)]
  pick=order(dpc2)[1]+1
  }
  else if (approach=="cv") {
  # Compute Cross Validation Mean Square Errors
  cvfold=cvfold
  sel.mat<-cv<-NULL
  rand <- sample(n)%%cvfold + 1
  for (j in 2:length(baseslist)) {
	sel.mat=cbind(baseslist[[1]],sel.bases,new.sel[,j])
	mse.pred<-0
	for (k in 1:length(unique(rand))) {
	  coef=lm(y[rand!=k]~sel.mat[rand!=k,]-1)$coef
      usemat=sel.mat[rand==k,!is.na(coef)]  # 100802 YW correction (WM) for leave-1-out CV
	  usecoef=coef[!is.na(coef)]
      if (length(usecoef)==1) predd=as.vector(usemat*usecoef)   # 100802 YW (WM)
      if (length(usecoef)>1) predd=as.vector(usemat%*%usecoef)  # 100802 YW (WM)

### OR: (YW alternative 100802) # predd=ifelse(rep(length(usecoef)==1,nrow(usemat)),usemat*usecoef,usemat%*%usecoef)
###     (alternative not tested in sims)

	  mse.pred <- mse.pred+sum((predd-y[rand==k])^2)
	}
	  cv=c(cv,mse.pred/n)
  }
  pick=order(cv)[1]+1  
  }
  else if (approach=="covpen") {
  cp.score=NULL
  for (j in 2:length(baseslist)){
    if (crit=="CIC") cp.score=c(cp.score,rss[j]+ct*cum.bas.cost[j])
	else if (crit=="GCV") cp.score=c(cp.score,rss[j]/(n-cum.bas.cost[j]/est.var)^2)
	else stop(paste("Criterion ",crit," is not supported!","\n",sep=""))
  }
  pick=order(cp.score)[1]+1  
  }
  sel.bases=cbind(sel.bases,cbind(lib[[pick]])[,index.new[pick]])

  iwk=rep(0,ncol(as.matrix(cbind(baseslist[[1]],sel.bases))))
  pick.eval=pick.bas(yy,yyptb,cbind(baseslist[[1]],sel.bases),iwk,wka=wka,wkall=wkall,wkamax=maxbas,ptb=ptb,nnull=nnull+i-1,nrep=repp,tauhat=tauhat,approach=approach)
  allfits[,i+nnull,]=pick.eval$fits

#  allisbraw=rbind(allisbraw,c(pick,lib.index[[pick]][libisb[[pick]]]))
  if (approach=="gdf" || approach=="covpen") allisbraw=rbind(allisbraw,c(pick,lib.index[[pick]][pick.eval$icb-(nnull+i-1)]))
#  difdif=c(pick,lib.index[[pick]][libisb[[pick]]])-c(pick,lib.index[[pick]][pick.eval$icb-(nnull+i-1)])
#  print(difdif)
  wka=wka.tmp[[pick]]
  wkall=pick.eval$wkall
  yy=yy.tmp[[pick]]
  yyptb=pick.eval$yyptb
  chosen.bases.index[i,]=c(pick,lib.index[[pick]][index.new[pick]])
#  print(chosen.bases.index[i,])
  rss.sel[i]=rss[pick]

  if (approach=="gdf") {
  # GDF calculation #
  bas.gdf[i]=ifelse(i==1,(bas.cost[pick]+nnull),bas.cost[pick])
  if (i>1) bas.gdf[i]=bas.gdf[i]+bas.gdf[i-1]
  gdf=bas.gdf[i]  
  }
  else if (approach=="cv") {
  bas.gdf[i]=min(cv)
  }
  else if (approach=="covpen") {
  bas.gdf[i]=min(cp.score)
  }

  lib.size[pick]=0
  if (ncol(as.matrix(lib[[pick]]))>1) {
  lib[[pick]]=lib[[pick]][,-index.new[pick]]
  lib.index[[pick]]=lib.index[[pick]][-index.new[pick]]
  selindex[[pick]]=selindex[[pick]][-index.new[pick]]
  lib.size[pick]=abs(sum(selindex[[pick]]-1))
  }

  i=i+1
  if (approach=="gdf") cond=((i <= (maxbas-nnull)) && (gdf < length(y)))
  else cond=(i <= (maxbas-nnull))
  }
  maxbas.sel=i-1
#  maxbas.sel=length(bas.gdf[1:(i-1)])

  # Trimming back
  if (!is.null(lib[[1]])) {
    fit.null=glm(y~lib[[1]]-1)
    chosen.bases.full=cbind(chosen.bases.index[1:maxbas.sel,])
    if (approach=="gdf") {
	  if (crit=="GCV") score.null=(fit.null$dev)/(nobs-ncol(lib[[1]]))^2
	  else if (crit=="CIC") score.null=fit.null$dev + ct*ncol(lib[[1]])*est.var
	}
	else if (approach=="cv") {
	  cvfold=cvfold
	  sel.mat<-NULL
	  rand <- sample(n)%%cvfold + 1
	  sel.mat=lib[[1]]
	  mse.pred<-0
	  for (k in 1:length(unique(rand))) {
	  coef=lm(y[rand!=k]~sel.mat[rand!=k,]-1)$coef
      usemat=sel.mat[rand==k,!is.na(coef)]  # 100802 YW correction (WM) for leave-1-out CV
	  usecoef=coef[!is.na(coef)]
      if (length(usecoef)==1) predd=as.vector(usemat*usecoef)   # 100802 YW (WM)
      if (length(usecoef)>1) predd=as.vector(usemat%*%usecoef)  # 100802 YW (WM)

### OR: (YW alternative 100802) # predd=ifelse(rep(length(usecoef)==1,nrow(usemat)),usemat*usecoef,usemat%*%usecoef)
###     (alternative not tested in sims)

	  mse.pred <- mse.pred+sum((predd-y[rand==k])^2)
	  }
	  score.null=mse.pred/n
	}
	else if (approach=="covpen") score.null=cp.null.score
    fitted=fit.null$fitted
    coef=fit.null$coef
    res=fit.null$res
    nb=nnull
  }

  score=rep(0,maxbas.sel)
  if (approach=="gdf") {
  for (i in 1:maxbas.sel) {
	if (crit=="GCV") score[i]=rss.sel[i]/(length(y)-bas.gdf[i])^2
	else if (crit=="CIC") score[i]=rss.sel[i] + ct*bas.gdf[i]*est.var
	else stop(paste("Criterion ",crit," is not supported!","\n",sep=""))
	}
  }
  else if (approach=="cv" || approach=="covpen") score=bas.gdf[1:maxbas.sel]
  trim=order(score)[1]
  if (is.na(score.null) || (!is.na(score.null)&&(score.null>score[trim]))) {
  final.bases=cbind(baseslist[[1]],sel.bases[,1:trim])
  fit=glm(y~final.bases-1)
  fitted=fit$fitted
  coef=fit$coef
  res=fit$res
  gdf=bas.gdf[1:maxbas.sel]
  idf=(bas.gdf-nnull)/(maxbas.sel-nnull)
  idf=idf[1:maxbas.sel]
  idf[idf<0]=0
  chosen.bases.full=cbind(chosen.bases.index[1:maxbas.sel,])
  chosen.bases.trim=cbind(chosen.bases.index[1:trim,])
  if (!is.null(sub.maxbas)) {
    for (i in 1:length(sub.maxbas)) {
      libnum.full=chosen.bases.full[,1]
      temp.full=chosen.bases.full[libnum.full==(i+1),2]
      chosen.bases.full[libnum.full==(i+1),2]=orig.index[temp.full,i]
      libnum.trim=chosen.bases.trim[,1]
      temp.trim=chosen.bases.trim[libnum.trim==(i+1),2]
      chosen.bases.trim[libnum.trim==(i+1),2]=orig.index[temp.trim,i]
    }
  }
  chosen.bases.matrix=sel.bases[,1:trim]
  nb=trim+nnull
  }
  score=c(score.null,score)

  # Construct two-column array output for all index sequences
  if (approach=="gdf" || approach=="covpen") {
  allisb=array(0,dim=c(nrow(allisbraw),2,ncol(allisbraw)-1))
  for (i in 2:ncol(allisbraw)) allisb[,,i-1]=cbind(allisbraw[,1],allisbraw[,i])
  }
  
  # Backward Elimination (output of GDFs and IDFs will not be affected)
  if (!is.null(backward) && !is.null(chosen.bases.full)) {
    if (nnull==0) importance=regsubsets(final.bases,y,method="backward",intercept=F)
    else importance=regsubsets(final.bases,y,force.in=c(1:nnull),method="backward",intercept=F)

    imp.rss=importance$rss[nnull:(nnull+trim)]
    gdf.trim=bas.gdf[1:(nnull+trim)]
    gdf.trim=gdf.trim[importance$vorder[1:(nnull+trim)]]
	gdf.trim=gdf.trim[nnull:(nnull+trim)]
  
    if (backward == "GCV") scr <- imp.rss/(nobs-(gdf.trim))^2
    if (backward == "AIC") scr <- nobs*log(imp.rss/nobs) + 2*(gdf.trim)
    if (backward == "BIC") scr <- nobs*log(imp.rss/nobs) + log(nobs)*(gdf.trim)
    if (backward == "Cp")  scr <- imp.rss + 2*(gdf.trim)*imp.rss[trim]/(nobs-(gdf.trim))

    trim.again=order(scr)[1]-1
	chosen.bases.trim.neworder=chosen.bases.trim[(importance$vorder[(nnull+1):(nnull+trim)]-nnull),]
    chosen.bases.trim=chosen.bases.trim[(importance$vorder[(nnull+1):(nnull+trim)]-nnull),][1:trim.again,]
	# The following line controls whether to preserve the original bases order
	chosen.bases.full <- rbind(chosen.bases.trim.neworder,chosen.bases.full[-(1:(nb-nnull)),])
    chosen.bases.matrix=chosen.bases.matrix[,(importance$vorder[(nnull+1):(nnull+trim)]-nnull)][,1:trim.again]
    final.bases.backward=cbind(baseslist[[1]],chosen.bases.matrix)
    fit=glm(y~final.bases.backward-1)
    fitted=fit$fitted
    coef=fit$coef
	res=fit$res
    nb=trim.again+nnull
	score=scr
  }

  # Summary information
  if (approach=="gdf") {
  gen_deg_fre <- ifelse(is.null(baseslist[[1]]),gdf[nb],gdf[(nb-nnull)])
  sigma_sq <- sum(res^2)/(n-gen_deg_fre)
  }
  else {
  gen_deg_fre <- ifelse(is.null(baseslist[[1]]),gdf[nb],gdf[(nb-nnull)])/est.var
  sigma_sq <- sum(res^2)/(n-gen_deg_fre)
  }
  if (!is.null(chosen.bases.full)) chosen.bases.full=as.matrix(chosen.bases.full)
  if (!is.null(chosen.bases.trim)) chosen.bases.trim=as.matrix(chosen.bases.trim)
  if (!is.null(chosen.bases.matrix)) chosen.bases.matrix=as.matrix(chosen.bases.matrix)
  est <- list(coefficients=coef,residuals=res,fitted.values=fitted,chosen.bases.full=chosen.bases.full,chosen.bases.trim=chosen.bases.trim,chosen.bases.matrix=chosen.bases.matrix,null.bases.matrix=baseslist[[1]],lib.size=size.of.lib,score=score,gdf=gdf,idf=idf,nb=nb,sigma_sq=sigma_sq,dof=gen_deg_fre,allisb=allisb,ptb=ptb,allfits=allfits,y=y,baseslist=baseslist,maxbas=maxbas,sub.maxbas=sub.maxbas,backward=backward,nobascheck=nobascheck,approach=approach,control=ctrl.vals)
  est$method="bsmls"
  est$call <- match.call()
  class(est) <- "bsml"
  return(est)
}

#########################
# Basis Function Cost 2 #
# Updated Oct 13, 2009  #
#########################

# Selects up to maxbas number of basis functions and calculates the corresponding
# GDF's and IDF's. Two approaches are available: the one described in Ye's paper
# or the covariance penalty approach.
# Calls either LIBCOST2 (Ye) or LIBCOST3 (Covariance penalty)


bf.cost2 <- function(y,bas,iwk=NA,nnull=NA,maxbas=NA,wkamax=NA,nrep=100,tauhat=NA,approach="gdf")
{
    nobs <- length(y)
    nbas <- ncol(as.matrix(bas))
 
    idfhat <- gdfhat <- rep(0,maxbas)
    if (nbas < maxbas) stop("nbas < maxbas")	
    isb <- rss <- seq(0, maxbas)
    yy <- y
	if (ncol(as.matrix(bas))!=length(iwk)) stop("Number of columns in library is not equal to length of flags vector!","\n")
    baswk <- bas 
    wka <- matrix(0,nobs,wkamax)
    wkall <- array(rep(0,nobs*wkamax*nrep),dim=c(nobs,wkamax,nrep))
    icb <- matrix(0,maxbas,nrep)
    ptb <- matrix(0,nobs,nrep) 
    u <- h <- z <- zk <- rep(0,nobs)
    fits <- array(rep(0,nobs*maxbas*nrep),dim=c(nobs,maxbas,nrep))
    fitswk <- matrix(0,nobs,maxbas)
    yyptb <- matrix(0,nobs,nrep)
    if (approach=="gdf") routine="libcost2"
    else if (approach=="covpen") {
	routine="libcost3"
	}
	else if (approach=="cv") {
		routine="libcost2"
		nrep=1
	}

    hasfit <- .Fortran(routine,
            NOBS = as.integer(nobs),
            Y = as.double(y),
            BAS = as.double(bas),
            NNULL = as.integer(nnull),
            NBAS = as.integer(nbas),
            MAXBAS = as.integer(maxbas),
            IDFHAT = as.double(idfhat),
	      nrep = as.integer(nrep),
	      GDFHAT = as.double(gdfhat),
            TAUHAT = as.double(tauhat),
            ISB = as.integer(isb),
            RSS = as.double(rss),
            IWK = as.integer(iwk), 
            YY = as.double(yy), 
            BASWK = as.double(baswk), 
            WKA = as.double(wka), 
		WKALL = as.double(wkall),
		WKAMAX = as.integer(wkamax),
		ICB = as.double(icb),
		PTB = as.double(ptb),
            U = as.double(u),
            H = as.double(h),
            FITS = as.double(fits),
            FITSWK = as.double(fitswk),
            Z = as.double(z),
            ZK = as.double(zk),
		YYPTB = as.double(yyptb)
	      )[c("IDFHAT","GDFHAT","ISB","IWK","RSS","ICB","PTB","FITS","U","H","Z","ZK","WKA","WKALL","BASWK","YY","FITSWK","YYPTB")]
    rss <- hasfit$RSS[1:maxbas]
    isb <- hasfit$ISB[1:maxbas]
    iwk <- hasfit$IWK
    idfhat <- hasfit$IDFHAT 
    gdfhat <- hasfit$GDFHAT
    icb <- matrix(as.integer(hasfit$ICB),ncol=nrep,byrow=F)
    ptb <- matrix(hasfit$PTB,ncol=nrep,byrow=F)
    fits <- array(hasfit$FITS,dim=c(nobs,maxbas,nrep))
    u=hasfit$U
    h=hasfit$H
    z=hasfit$Z
    zk=hasfit$ZK
    wka=matrix(hasfit$WKA,nrow=nobs)
    wkall=array(hasfit$WKALL,dim=c(nobs,wkamax,nrep))
    basnew=matrix(hasfit$BASWK,nrow=nobs)
    yy=hasfit$YY
    yyptb=matrix(hasfit$YYPTB,nrow=nobs)
    last.fit=matrix(hasfit$FITSWK,nobs)

    return(list(isb=isb,icb=icb,iwk=iwk,rss=rss,gdfhat=gdfhat,icb=icb,ptb=ptb,idfhat=idfhat,fits=fits,u=u,h=h,z=z,zk=zk,wka=wka,wkall=wkall,basnew=basnew,yy=yy,yyptb=yyptb,last.fit=last.fit))
}

################
# Pick Basis   #
# Oct 13, 2009 #
################

# Determines cost of adding a basis from a library given that a set of bases
# have already been selected for the model.
# Calls Fortran subroutine PICKBAS.


pick.bas <- function(y,yyptb,bas,iwk,wka,wkall,wkamax,ptb,nnull,nrep,tauhat,approach="gdf")
{
    nobs <- length(y)
    nbas <- ncol(bas)
    idfhat <- gdfhat <- 0
    isb <- rss <- 0
    yy <- y
    baswk <- bas 
    u <- h <- z <- zk <- rep(0,nobs)
    fits <- matrix(0,nobs,nrep)
    fitswk <- rep(0,nobs)
	icb <- rep(0,nrep)
	if (ncol(as.matrix(bas))!=length(iwk)) stop("Number of columns in library is not equal to length of flags vector!","\n")
	if (approach=="gdf") routine="PICKBAS"
	else if (approach=="covpen") routine="PICKBAS2"
	else if (approach=="cv") {
	routine="PICKBAS"
	nrep=1
	}

    hasfit <- .Fortran(routine,
            NOBS = as.integer(nobs),
            Y = as.double(y),
            BAS = as.double(bas),
            NNULL = as.integer(nnull),
            NBAS = as.integer(nbas),
            IDFHAT = as.double(idfhat),
	      nrep = as.integer(nrep),
	      GDFHAT = as.double(gdfhat),
            TAUHAT = as.double(tauhat),
            ISB = as.integer(isb),
			ICB = as.integer(icb),
            RSS = as.double(rss),
            IWK = as.integer(iwk), 
            YY = as.double(yy), 
            BASWK = as.double(baswk), 
            WKA = as.double(wka), 
		WKALL = as.double(wkall),
		WKAMAX = as.integer(wkamax),
		PTB = as.double(ptb),
            U = as.double(u),
            H = as.double(h),
            FITS = as.double(fits),
            FITSWK = as.double(fitswk),
            Z = as.double(z),
            ZK = as.double(zk),
		YYPTB = as.double(yyptb)
	      )[c("IDFHAT","GDFHAT","ISB","ICB","IWK","RSS","FITS","U","H","Z","ZK","WKA","WKALL","BASWK","YY","YYPTB","FITSWK")]

    rss <- hasfit$RSS
    isb <- hasfit$ISB
	icb <- hasfit$ICB
    iwk <- hasfit$IWK
    idfhat <- hasfit$IDFHAT 
    gdfhat <- hasfit$GDFHAT
    fits <- matrix(hasfit$FITS,nobs)
    u=hasfit$U
    h=hasfit$H
    z=hasfit$Z
    zk=hasfit$ZK
    wka=matrix(hasfit$WKA,nobs)
    wkall=array(hasfit$WKALL,dim=c(nobs,wkamax,nrep))
    basnew=matrix(hasfit$BASWK,nobs)
    yy=hasfit$YY
    yyptb=matrix(hasfit$YYPTB,nobs)
    last.fit=hasfit$FITSWK

    return(list(isb=isb,icb=icb,iwk=iwk,rss=rss,gdfhat=gdfhat,idfhat=idfhat,fits=fits,u=u,h=h,z=z,zk=zk,wka=wka,wkall=wkall,wkamax=wkamax,basnew=basnew,yy=yy,yyptb=yyptb,last.fit=last.fit))
}
