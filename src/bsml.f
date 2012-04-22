c gcvpack.new
c
c Modified by jeff sklar 8-24-05
c
c Collection of fortran routines from gcvpack, and additional routines
c that allow ridge regression with GML and UBR options for smoothing
c parameter selection.
c Additional routines include: dgcvnew.f dgmin1.f dsnsm1.f dumin1.f dvl1.f 
c dvl21.f dvl31.f dvlop1.f dvmin1.f 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c modified by jeff sklar 1-18-2003
      subroutine dgcvnew (x,ldx,y,sigma,ldsigm,nobs,npar,nnull,sgpvt,
     * dcaux,ldcaux,adiag,lamlim,ntbl,dout,coef,tbl,ldtbl,auxtbl,work,
     * lwa,job,info,vm)
      integer ldx,ldsigm,nobs,npar,nnull,sgpvt(npar),ldcaux,ntbl,ldtbl,
     * lwa,info,job, vm
      double precision x(ldx,npar),y(nobs),sigma(ldsigm,npar),
     * dcaux(ldcaux),adiag(nobs),lamlim(2),dout(4),coef(npar),
     * tbl(ldtbl,3),auxtbl(3,3),work(lwa)
c
c Purpose: determine the generalized cross validation estimate of the 
c	smoothing parameter and fit model parameters for a penalized
c 	least squares problem with a semi-norm smoothing matrix.
c
c On Entry:
c   x(ldx,npar)		design matrix as returned by ddcom
c   ldx			leading dimension of x as declared in the 
c			calling program, must be at least max(nobs,npar)
c   y(nobs)		response vector
c   sigma(ldsigm,npar)	Cholesky factor of sigma as returned by ddcom
c   ldsigm		leading dimension of sigma as
c			declared in the calling program
c   nobs		number of observations
c   npar		number of parameters
c   nnull		dimension of the null space of sigma
c   sgpvt(npar)		permuted indices from the pivoted Cholesky
c			decomposition of the semi-norm matrix
c   dcaux(ldcaux)	auxiliary vector from ddcom
c   ldcaux		length of ldcaux. Must be at least
c			(npar-nnull)**2+2*npar-nnull
c   adiag(nobs)	 	"true" y values on entry if predictive mse if 
c			requested 
c   lamlim(2)		limits on lambda hat search (in log10(nobs*
c			lambda)	scale) if user input limits are 
c			requested. if lamlim(1) = lamlim(2) then lamhat
c			is set to (10**lamlim(1))/nobs
c   ntbl		number of evenly spaced values for 
c			log10(nobs*lambda) to be used in the initial 
c			grid search for lambda hat
c			if ntbl = 0 only a golden ratio search will be 
c			done and tbl is not referenced, if ntbl > 0
c			there will be ntbl rows returned in tbl
c   ldtbl		leading dimension of tbl as declared in the 
c			calling program	
c   job			integer with decimal expansion abc
c			if a is nonzero then predictive mse is computed
c			   using adiag as true y
c			if b is nonzero then user input limits on search
c			   for lambda hat are used
c			if c is nonzero then diagonal of the hat matrix
c			   is calculated
c
c On Exit:
c   y(nobs)		predicted values
c   adiag(nobs)		diagonal elements of the hat matrix if requested
c   lamlim(2)		limits on lambda hat search 
c			(in log10(nobs*lambda) scale)
c   dout(4)		contains:
c  			1  lamhat   generalized cross validation
c				    estimate of the smoothing parameter
c			2  penlty   smoothing penalty
c			3  rss	    residual sum of squares
c			4  tr(I-A)  trace of I - A
c   coef(npar)		coefficient estimates
c   tbl(ldtbl,3)	column	contains
c			  1 	grid of log10(nobs*lambda) 
c			  2  	V(lambda)
c			  3     R(lambda) if requested
c   auxtbl(3,3)		auxiliary table
c			1st row contains:
c			    log10(nobs*lamhat), V(lamhat) and  
c			    R(lamhat) if requested
c			    where lamhat is the gcv estimate of lambda
c			2nd row contains:
c			    0, V(0) and  R(0) if requested
c			3rd row contains:
c			    0, V(infinity) and R(infinity) if requested
c   info		error indicator
c			  0 : successful completion
c			 -1 : log10(nobs*lamhat) <= lamlim(1) 
c			      (not fatal)
c			 -2 : log10(nobs*lamhat) >= lamlim(2) 
c			      (not fatal)
c			  1 : dimension error
c			  2 : error in ntbl
c			  3 : ldcaux (length of dcaux) is too small
c			  4 : lwa (length of work) is too small
c			  5 : lamlim(1) > lamlim(2)
c			 10 < info < 20 : 10 + nonzero info returned
c					  from dvlop
c			 20 < info < 30 : 20 +nonzero info returned
c					  from dcfcr
c
c Work Arrays:
c   work(lwa) 		double precision work vector
c   lwa			length of work as declared in the calling 
c			program must be at least (npar-nnull)+nobs
c
c Subprograms Called Directly:
c	Gcvpack - drsap dvlop dpmse dcfcr dpdcr ddiag
c
c Subprograms Called Indirectly:
c	Gcvpack - dvmin dvl
c	Linpack - dqrsl dtrsl
c	Blas    - dcopy ddot dgemv
c	Other   - dprmut dset
c
c $Header: dgcv.f,v 2.100.1.1 86/10/07 12:49:16 lindstrom Exp $
c
      double precision addend,nlamht,ssqw2
      integer mp1,mnh,pmh,minnp,pp1,pmhp1,wsize,npsing,i,jpmse,jlaml,
     * jadiag,sinfo
c
      jpmse = job/100
      jlaml = mod(job,100)/10
      jadiag = mod(job,10)
c			check dimensions
      sinfo = 0
      info = 0
      pmh = npar - nnull
      if ((nobs .le. 0) .or. (npar .le. 0) .or. (nnull .lt. 0) 
     * .or. (pmh .le. 0)) then
         info = 1
         return
      endif
      if ((ntbl .lt. 0) .or. (ntbl .gt. ldtbl)) then
         info = 2
         return
      endif
      wsize = pmh**2 + 2*npar - nnull
      if (ldcaux .lt. wsize) then
         info = 3
         return
      endif
      wsize = pmh + nobs
      if (lwa .lt. wsize) then
         info = 4
         return
      endif
      if (jlaml .ne. 0 .and. (lamlim(1) .gt. lamlim(2))) then
         info = 5
         return
      endif
      pmhp1 = pmh + 1
      pp1 = npar + 1
c			calculate npsing
      minnp = min(nobs-nnull,pmh)
      do 30 i=1,minnp
         if (dcaux(npar+i)**2 .gt. 0) npsing = i
   30 continue
c			apply rotations to y
      mp1 = nnull + 1
      mnh = nobs - nnull
      call drsap (x(1,pmhp1),ldx,nobs,nnull,dcaux(pmhp1),x(mp1,1),ldx,
     * mnh,npsing,y,ssqw2,addend,work)
c			minimize V(lambda)
c     vm added 4/27/01 11:26PM
      call dvlop1 (y(mp1),dcaux(pp1),nobs,nnull,npsing,addend,ssqw2,
     * lamlim,ntbl,nlamht,tbl,ldtbl,auxtbl,dout(3),jlaml,info,vm)
      if (info .gt. 0) then
         info = info + 10
         return
      endif
      if (info .lt. 0) sinfo = info
c			calculate predictive mse
      if (jpmse .ne. 0) then
         call dpmse(x(1,pmhp1),ldx,nobs,nobs,nnull,dcaux(pmhp1),
     *    dcaux(pp1),npsing,x(mp1,1),ldx,y,y(mp1),ntbl,adiag,tbl,ldtbl,
     *    auxtbl,work)
      endif
c			calculate coefficients
      call dcfcr (x(1,pmhp1),ldx,nnull,sigma,ldsigm,npar,
     * pmh,dcaux,sgpvt,x,ldx,dcaux(pp1),npsing,dcaux(npar+pmhp1),pmh, 
     * nlamht,y,y(mp1),coef,dout(2),work,info)
      if (info .gt. 0) then
         info = info + 20
         return
      endif
c			calculate predicted values 
      call dpdcr (x(1,pmhp1),ldx,nobs,nnull,dcaux(pmhp1),dcaux(pp1),
     * npsing,x(mp1,1),ldx,nlamht,y,y(mp1),y,work)
      if (jadiag .ne. 0) then
          call ddiag (x(1,pmhp1),ldx,nobs,nnull,dcaux(pmhp1),
     * 	  dcaux(pp1),npsing,x(mp1,1),ldx,nlamht,adiag,work)
      endif
      if (sinfo .lt. 0) info = sinfo
      return
      end
      double precision function dgmin1(lower,upper,svals,z,npsing,
     * ntbl,tbl,ldtbl,vlamht,info)
      integer npsing,ntbl,ldtbl,info
      double precision lower,upper,svals(npsing),z(npsing),
     * tbl(ldtbl,3),vlamht
c
c Purpose: evaluate V(lambda) for a grid of ln(nobs*lambda) values 
c	between	lower and upper, store these in the array tbl, and find
c	minimizer of v.
c
c On Entry:
c   lower		lower bound of interval (in nobs*ln(lambda)
c			scale) over which V(lambda) is to be minimized
c   upper		upper bound of interval (in nobs*ln(lambda)
c			scale) over which V(lambda) is to be minimized
c   svals(npsing)	singular values 
c   z(npsing)		data vector in canonical coordinates
c   npsing		number of positive elements of svals 
c   ntbl		number of evenly spaced values for 
c			ln(nobs*lambda)	to be used in the initial grid 
c			search for lambda hat
c			if ntbl = 0 only a golden ratio search will be 
c			done and tbl is not referenced, if ntbl > 0
c			there will be ntbl rows returned in tbl
c   ldtbl		leading dimension of tbl as declared in the 
c			calling program	
c   
c On Exit:
c   tbl(ldtbl,3)	column	contains
c			  1 	grid of ln(nobs*lambda) 
c			  2  	V(lambda)
c   vlamht		V(lambda hat)
c   dvmin		ln(nobs*lambda hat)
c   info		error indicator
c			  0 : successful completion
c			 -1 : dvmin <= lower (not fatal)
c			 -2 : dvmin >= upper (not fatal)
c
c Subprograms Called Directly:
c	Gcvpack - dvl
c
c $Header: dmin.f,v 2.100.1.1 86/10/07 12:59:51 lindstrom Exp $
c
c last modified by Yu-Chieh Yang Jul 17 16:01:27 PDT 2001
c argument setting is similar to dvmin.f except for dvl2

      double precision a,b,c,d,vc,vd,del,k1,k2,x,v
      integer j,jmin,k
      double precision dvl21
c				null interval
      if (lower .eq. upper) then
	 dvmin1 = lower
	 info = -1
	 vlamht = dvl21(lower,svals,z,npsing)
	 do 10 j = 1, ntbl
	    tbl(j,1) = lower
	    tbl(j,2) = vlamht
   10    continue
	 return
      end if
c				non-null interval
      info = 0
      a = lower
      b = upper
      if (ntbl .eq. 1) then
	 x = (a + b)/2
	 tbl(1,1) = x
	 tbl(1,2) = dvl21(x,svals,z,npsing)
      else if (ntbl .ge. 2) then
c			do grid search
	 v=dvl21(lower,svals,z,npsing)*2.0d0
	 del=(upper-lower)/(ntbl-1)
	 do 20 j = 1, ntbl
	    tbl(j,1) = lower + (j - 1) * del
	    tbl(j,2) = dvl21(tbl(j,1),svals,z,npsing)
	    if (tbl(j,2) .le. v) then
	       jmin = j
	       v = tbl(j,2)
	    endif
   20    continue	      
	 a=tbl(jmin,1)-del
	 b=tbl(jmin,1)+del
      end if
c			do golden ratio search			
      k1=(3.0d0-dsqrt(5.0d0))/2.0d0
      k2=(dsqrt(5.0d0)-1)/2.0d0
      c = a + k1*(b - a)
      d = a + k2*(b - a)
      vc = dvl21(c,svals,z,npsing)
      vd = dvl21(d,svals,z,npsing)
      do 30 k=1,50
	 if (vd .lt. vc) then
	    a = c
	    c = d
	    d = a + k2*(b - a)
	    vc = vd
	    vd = dvl21(d,svals,z,npsing)
	 else
	    b = d
	    d = c
	    c = a + k1*(b - a)
	    vd = vc
	    vc = dvl21(c,svals,z,npsing)
	 end if
   30 continue
      x=(a+b)/2
      if (x .le. lower) info = -1
      if (x .ge. upper) info = -2
      vlamht=dvl21(x,svals,z,npsing)
      dgmin1 = x
      return
      end
c modified by jeff sklar 1-19-2003
      subroutine dsnsm1 (x,ldx,y,sigma,ldsigm,nobs,npar,nnull,adiag,
     * tau,lamlim,ntbl,dout,iout,coef,svals,tbl,ldtbl,auxtbl,
     * iwork,liwa,work,lwa,job,info,vm)
      integer ldx,ldsigm,nobs,npar,nnull,ntbl,iout(3),ldtbl,liwa,
     * iwork(liwa),lwa,job,info,vm
      double precision x(ldx,npar),y(nobs),sigma(ldsigm,npar),
     * adiag(nobs),tau,lamlim(2),dout(5),coef(npar),svals(*),
     * tbl(ldtbl,3),auxtbl(3,3),work(lwa)
c
c Purpose: determine the generalized cross validation estimate of the 
c	smoothing parameter and fit model parameters for a penalized 
c	least squares problem with a semi-norm smoothing matrix.
c
c On Entry:
c   x(ldx,npar)		design matrix
c   ldx			leading dimension of x as declared in the 
c			calling program, must be at least max(nobs,npar)
c   y(nobs)		response vector
c   sigma(ldsigm,npar)	symmetric matrix that defines the semi-norm
c   ldsigm		leading dimension of sigma as declared
c			in the calling program
c   nobs		number of observations
c   npar		number of parameters
c   nnull		dimension of the null space of sigma
c   adiag(nobs)		"true" y values on entry if computation of 
c			predictive mse is requested
c   lamlim(2)		limits on lambda hat search (in log10(nobs*
c			lambda) scale) if user input limits are 
c			requested if lamlim(1) = lamlim(2) then lamhat
c			is set to (10**lamlim(1))/nobs
c   tau			multiplier controlling the amount of truncation
c			if truncation is requested (try tau = 1 
c			to start then try 10 and 100)
c   ntbl		number of evenly spaced values for 
c			log10(nobs*lambda) to be used in the initial 
c			grid search for lambda hat
c			if ntbl = 0 only a golden ratio search will be 
c			done and tbl is not referenced, if ntbl > 0
c			there will be ntbl rows returned in tbl
c   ldtbl		leading dimension of tbl as declared in the 
c			calling program	
c   job			integer with decimal expansion abcd
c			if a is nonzero then truncation is used
c			if b is nonzero then predictive mse is computed
c			   using adiag as true y
c			if c is nonzero then user input limits on search
c			   for lambda hat are used
c			if d is nonzero then the diagonal of the hat 
c			   matrix is calculated
c
c On Exit:
c   x(ldx,npar)		overwritten with many intermediate results
c   y(nobs)		predicted values
c   sigma(ldsigm,npar)	overwritten with the QR decomposition of the 
c			Cholesky factor of sigma
c   adiag(nobs)		diagonal elements of the hat matrix if requested
c   lamlim(2)		limits on lambda hat search 
c			(in log10(nobs*lambda) scale)
c   dout(5)		contains:
c  			1 lamhat   generalized cross validation 
c				   estimate of the smoothing parameter
c			2 penlty   smoothing penalty
c			3 rss	   residual sum of squares
c			4 tr(I-A)  trace of I - A
c			5 truncation ratio = 1/(1+(normk/(nobs*lamhat)))
c				   where normk = norm(R - R sub k)**2
c   iout(3)		contains:
c			1  npsing   number of positive singular
c				    values
c				    if info indicates nonzero info in 
c				    dsvdc then iout(1) contains info as
c				    it was returned from dsvdc
c			2  npar	    number of parameters
c			3  nnull    size of the null space of sigma
c   coef(npar)		coefficient estimates
c   svals(npar-nnull)	first npsing entries contain singular values 
c			of the matrix j2 
c			if info indicates nonzero info in dsvdc then 
c			svals is as it was returned from dsvdc
c   tbl(ldtbl,3)	column	contains
c			  1 	grid of log10(nobs*lambda) 
c			  2  	V(lambda)
c			  3     R(lambda) if requested
c   auxtbl(3,3)		auxiliary table
c			1st row contains:
c			    log10(nobs*lamhat), V(lamhat) and  
c			    R(lamhat) if requested
c			    where lamhat is the gcv estimate of lambda
c			2nd row contains:
c			    0, V(0) and  R(0) if requested
c			3rd row contains:
c			    0, V(infinity) and R(infinity) if requested
c   info		error indicator
c			   0 : successful completion
c			  -3 : nnull is too small (not fatal)
c			  -2 : log10(nobs*lamhat) >= lamlim(2) 
c			       (not fatal)
c			  -1 : log10(nobs*lamhat) <= lamlim(1)
c			       (not fatal)
c			   1 : dimension error 	
c			   2 : lwa (length of work) is too small
c			   3 : liwa (length of iwork) is too small
c			   4 : error in ntbl or tau
c			  100< info <200 : 100 + nonzero info returned
c					   from ddcom
c			  200< info <300 : 200 + nonzero info returned
c					   from dgcv
c
c Work Arrays:
c   work(lwa)		double precision work vector
c   lwa			length of work as declared in the calling 
c			program 
c			must be at least 
c			(npar-nnull)*(npar-2*nnull+2+nobs)+npar+nobs
c   iwork(liwa)		integer work vector
c   liwa		length of iwork as declared in the calling 
c			program
c			must be at least 2*npar - nnull 
c
c Subprograms Called Directly:
c	Gcvpack - ddcom dgcv
c
c Subprograms Called Indirectly:
c	Gcvpack - dcrtz dsgdc dcfcr drsap dvlop dtsvdc
c		  dpmse dvmin dvl dzdc dpdcr ddiag
c	Linpack - dchdc dqrdc dqrsl dtrsl dsvdc dtrco
c	Blas    - dcopy ddot dgemv dswap
c	Other 	- dcpmut dprmut dset 
c
c $Header: dsnsm.f,v 2.100.1.2 86/11/19 09:24:39 lindstrom Exp $
c
      integer bigp,pmh,pp1,lwa2,wsize,ldcaux,job1,job2,npsing,sinfo
      sinfo = 0
      info = 0
      iout(2) = npar
c			check dimensions
      pmh = npar - nnull
      if ((nobs .le. 0) .or. (npar .le. 0) .or. (nnull .lt. 0) 
     * .or.(pmh .le. 0)) then
         info = 1
         return
      endif
      wsize = pmh**2 +2*npar - nnull + pmh*(nobs - nnull) + pmh + nobs
      if (lwa .lt. wsize) then
         info = 2
         return
      endif
      wsize = npar + pmh
      if (liwa .lt. wsize) then
         info = 3
         return
      endif
      if ((ntbl .lt. 0) .or. (ntbl .gt. ldtbl) .or. (tau .lt. 0)) then
         info = 4
         return
      endif
c 			first p positions of iwork contain sgpvt
c 			first bigp-1  positions of work contain dcaux 
c			iwork vector starts at pp1
      pp1 = npar + 1
c 	work vector starts at bigp
      bigp = pmh**2 + 2*npar - nnull + 1
      ldcaux = pmh**2 + 2*npar - nnull
      lwa2 = lwa - ldcaux
      job1=job/1000
      job2=mod(job,1000)
c			decompose sigma and design matrix
      call ddcom (x,ldx,sigma,ldsigm,nobs,npar,nnull,tau,npsing,
     * svals,iwork,work,ldcaux,dout(5),work(bigp),lwa2,iwork(pp1),pmh,
     * job1,info)
      iout(1) = npsing
      iout(3) = nnull
      if (info .gt. 0) then
         info = info + 100
         return
      endif
      if (info .lt. 0) sinfo = info
c			compute lambda hat and other parameters
c	vm added by Yu-Chieh Yang 4/27/01 11:27PM
      call dgcvnew(x,ldx,y,sigma,ldsigm,nobs,npar,nnull,iwork,work,
     * ldcaux,adiag, lamlim,ntbl,dout,coef,tbl,ldtbl,auxtbl,work(bigp),
     * lwa2,job2,info,vm)
c      call dgcv(x,ldx,y,sigma,ldsigm,nobs,npar,nnull,iwork,work,ldcaux,
c     * adiag, lamlim,ntbl,dout,coef,tbl,ldtbl,auxtbl,work(bigp),lwa2,
c     * job2,info)

      if (info .gt. 0) then
         info = info + 200
         return
      endif
      dout(5) = 1.0d0/(1.0d0+dout(5)/nobs/dout(1))
      if (sinfo .lt. 0) info = sinfo
      return
      end
c modified by jeff sklar 1-18-2003
      double precision function dumin1(lower,upper,svals,z,npsing,
     * ntbl,tbl,ldtbl,vlamht,info)
      integer npsing,ntbl,ldtbl,info
      double precision lower,upper,svals(npsing),z(npsing),
     * tbl(ldtbl,3),vlamht
c
c Purpose: evaluate V(lambda) for a grid of ln(nobs*lambda) values 
c	between	lower and upper, store these in the array tbl, and find
c	minimizer of v.
c
c On Entry:
c   lower		lower bound of interval (in nobs*ln(lambda)
c			scale) over which V(lambda) is to be minimized
c   upper		upper bound of interval (in nobs*ln(lambda)
c			scale) over which V(lambda) is to be minimized
c   svals(npsing)	singular values 
c   z(npsing)		data vector in canonical coordinates
c   npsing		number of positive elements of svals 
c   ntbl		number of evenly spaced values for 
c			ln(nobs*lambda)	to be used in the initial grid 
c			search for lambda hat
c			if ntbl = 0 only a golden ratio search will be 
c			done and tbl is not referenced, if ntbl > 0
c			there will be ntbl rows returned in tbl
c   ldtbl		leading dimension of tbl as declared in the 
c			calling program	
c   
c On Exit:
c   tbl(ldtbl,3)	column	contains
c			  1 	grid of ln(nobs*lambda) 
c			  2  	V(lambda)
c   vlamht		V(lambda hat)
c   dvmin		ln(nobs*lambda hat)
c   info		error indicator
c			  0 : successful completion
c			 -1 : dvmin <= lower (not fatal)
c			 -2 : dvmin >= upper (not fatal)
c
c Subprograms Called Directly:
c	Gcvpack - dvl3
c
c $Header: dvmin.f,v 2.100.1.1 86/10/07 12:59:51 lindstrom Exp $
c

      double precision a,b,c,d,vc,vd,del,k1,k2,x,v
      integer j,jmin,k
      double precision dvl31
c				null interval
      if (lower .eq. upper) then
	 dvmin1 = lower
	 info = -1
	 vlamht = dvl31(lower,svals,z,npsing)
	 do 10 j = 1, ntbl
	    tbl(j,1) = lower
	    tbl(j,2) = vlamht
   10    continue
	 return
      end if
c				non-null interval
      info = 0
      a = lower
      b = upper
      if (ntbl .eq. 1) then
	 x = (a + b)/2
	 tbl(1,1) = x
	 tbl(1,2) = dvl31(x,svals,z,npsing)
      else if (ntbl .ge. 2) then
c			do grid search
	 v=dvl31(lower,svals,z,npsing)*2.0d0
	 del=(upper-lower)/(ntbl-1)
	 do 20 j = 1, ntbl
	    tbl(j,1) = lower + (j - 1) * del
	    tbl(j,2) = dvl31(tbl(j,1),svals,z,npsing)
	    if (tbl(j,2) .le. v) then
	       jmin = j
	       v = tbl(j,2)
	    endif
   20    continue	      
	 a=tbl(jmin,1)-del
	 b=tbl(jmin,1)+del
      end if
c			do golden ratio search			
      k1=(3.0d0-dsqrt(5.0d0))/2.0d0
      k2=(dsqrt(5.0d0)-1)/2.0d0
      c = a + k1*(b - a)
      d = a + k2*(b - a)
      vc = dvl31(c,svals,z,npsing)
      vd = dvl31(d,svals,z,npsing)
      do 30 k=1,50
	 if (vd .lt. vc) then
	    a = c
	    c = d
	    d = a + k2*(b - a)
	    vc = vd
	    vd = dvl31(d,svals,z,npsing)
	 else
	    b = d
	    d = c
	    c = a + k1*(b - a)
	    vd = vc
	    vc = dvl31(c,svals,z,npsing)
	 end if
   30 continue
      x=(a+b)/2
      if (x .le. lower) info = -1
      if (x .ge. upper) info = -2
      vlamht=dvl31(x,svals,z,npsing)
      dumin1 = x
      return
      end
      double precision function dvl1(lgnlam,svals,z,npsing)
      integer npsing
      double precision lgnlam,svals(npsing),z(npsing)
c
c Purpose: evaluate the cross-validation function with a semi-norm.
c
c On Entry:
c   lgnlam		log10(nobs*lambda) where lambda is the value of
c			lambda for which V is evaluated
c   svals(npsing)	singular values 
c   z(npsing)		data vector in canonical coordinates
c   npsing		number of positive svals 
c
c On Exit:
c   dvl			V(lambda)
c
c $Header: dvl.f,v 2.100.1.1 86/10/07 12:59:22 lindstrom Exp $
c
      integer j
      double precision nlam,numrtr,denom,factor
c
      common / gcvcom / addend,rss,tria,n,h
      integer n,h
      double precision rss,tria,addend
c     			see dvlop1 for definition of common block 
c			variables
c
      nlam = 10**lgnlam
      numrtr = addend
      denom = dble(n - h - npsing)
      do 10 j = 1,npsing 
         factor = 1.0d0/(1.0d0 + (svals(j)**2)/nlam)
         numrtr = numrtr + (factor*z(j))**2
         denom = denom + factor
   10 continue
      rss=numrtr
      tria=denom
      dvl1=dble(n)*numrtr/denom**2
      return
      end
      double precision function dvl21(lgnlam,svals,z,npsing)
      integer npsing
      double precision lgnlam,svals(npsing),z(npsing)
c
c Purpose: evaluate the cross-validation function with a semi-norm.
c dvl.f in GMLpack has been modified slightly 00/04/22
c
c On Entry:
c   lgnlam		log10(nobs*lambda) where lambda is the value of
c			lambda for which V is evaluated
c   svals(npsing)	singular values 
c   z(npsing)		data vector in canonical coordinates
c   npsing		number of positive svals 
c
c On Exit:
c   dvl			V(lambda)
c
c $Header: dvl.f,v 2.100.1.1 86/10/07 12:59:22 lindstrom Exp $
c
      integer j
c     double precision nlam,numrtr,denom,factor
c last modified by Yu-Chieh Yang, Jul 17 16:04:14 PDT 2001
      double precision nlam,numrtr,denom,factor,numr,denr
c
      common / gcvcom / addend,rss,tria,n,h
      integer n,h
      double precision rss,tria,addend
c     	see dvlop for definition of common block variables
c
c definition for GML as follows
c
c      y^2-z^2 + sum_i( ((n*lambda)/(di^2+n*lambda)) *zi^2) 
c  M = --------------------------------------------------------
c        n*( prod_i( (n*lambda)/(di^2+n*lambda)  )  )^(1/n)
c
c for easy maintainece, new arguments added : numr denr for GML 
c
      nlam = 10**lgnlam
      numrtr = addend
      numr = addend
      denr = 1.0
      denom = dble(n - h - npsing)
      do 10 j = 1,npsing 
         factor = 1.0d0/(1.0d0 + (svals(j)**2)/nlam)         
         numrtr = numrtr + (factor*z(j))**2
         denom = denom + factor
         numr = numr + factor*z(j)**2
         denr = denr*factor  
   10 continue
      rss=numrtr 
      tria=denom
      dvl21 = numr/( dble(n)*denr**(1.0/dble(n))) 
      return
      end
c modified by Jeff Sklar 1-18-2003
      double precision function dvl31(lgnlam,svals,z,npsing)
      integer npsing
      double precision lgnlam,svals(npsing),z(npsing)
c
c Purpose: evaluate the cross-validation function with a semi-norm.
c
c On Entry:
c   lgnlam		log10(nobs*lambda) where lambda is the value of
c			lambda for which V is evaluated
c   svals(npsing)	singular values 
c   z(npsing)		data vector in canonical coordinates
c   npsing		number of positive svals 
c
c On Exit:
c   dvl			V(lambda)
c
c $Header: dvl.f,v 2.100.1.1 86/10/07 12:59:22 lindstrom Exp $
c
      integer j
      double precision nlam,numrtr,denom,factor
c
      common / gcvcom / addend,rss,tria,n,h
      integer n,h
      double precision rss,tria,addend
c     			see dvlop for definition of common block 
c			variables
c
      nlam = 10**lgnlam
      numrtr = addend
      denom = dble(n - h - npsing)
      do 10 j = 1,npsing 
         factor = 1.0d0/(1.0d0 + (svals(j)**2)/nlam)
         numrtr = numrtr + (factor*z(j))**2
         denom = denom + factor
   10 continue
c
c .modified by Yu-Chieh Yang, Tue Jul  3 23:38:55 PDT 2001
c .keep above block as in dvl.f w/o any change, using above symbols
c            1           2   numrtr
c  UBR(AIC)= -*numrtr +  -*(--------)*[n-denom]
c            n           n   denom
c
c .then only modify command dvl3 below
      rss=numrtr
      tria=denom
c     dvl=dble(n)*numrtr/denom**2
      dvl31=(1.0/dble(n))*(numrtr + 2*(numrtr/denom)*(dble(n)-denom)) 
      return
      end
c modified by jeff sklar 1-18-2003
      subroutine dvlop1(z,svals,nobs,nnull,npsing,inadd,ssqw2,
     * lamlim,ntbl,nlamht,tbl,ldtbl,auxtbl,dout,job,info,vm)
      integer nobs,nnull,npsing,ntbl,ldtbl,job,info,vm
      double precision z(npsing),svals(npsing),inadd,ssqw2,lamlim(2),
     * nlamht,tbl(ldtbl,3),auxtbl(3,3),dout(2)
c
c Purpose: determine the optimal lambda for the generalized cross 
c	validation function given singular values and the data vector 
c	in canonical coordinates.
c
c On Entry:
c   z(npsing)		data vector in canonical coordinates
c   svals(npsing)	singular values 
c   nobs		number of observations
c   nnull		dimension of the null space of sigma
c   npsing		number of positive elements of svals 
c   inadd		constant term in expression for V
c   ssqw2		squared length of w2
c   lamlim(2)		limits on lambda hat search (in log10(nobs*
c			lambda) scale) if user input limits are 
c			requested. if lamlim(1) = lamlim(2) then nlamht
c			is set to 10**lamlim(1)
c   ntbl		number of evenly spaced values for 
c			log10(nobs*lambda) to be used in the initial 
c			grid search for lambda hat
c			if ntbl = 0 only a golden ratio search will be 
c			done and tbl is not referenced, if ntbl > 0
c			there will be ntbl rows returned in tbl
c   ldtbl		leading dimension of tbl as declared in the 
c			calling program	
c   job      		if job is nonzero then user input limits on 
c			lambda hat search are used
c
c On Exit:
c   lamlim(2)		limits on lambda hat search 
c			(in log10(nobs*lambda) scale)
c   nlamht		nobs*(lambda hat) where lambda hat is the gcv 
c			estimate of lambda
c   tbl(ldtbl,3)	column	contains
c			  1 	grid of log10(nobs*lambda) 
c			  2  	V(lambda)
c   auxtbl(3,3)		auxiliary table
c			1st row contains:
c			    log10(nobs*lambda hat), V(lambda hat)
c			2nd row contains:
c			    0, V(0) 
c			3rd row contains:
c			    0, V(infinity) 
c   dout(2)		contains:
c			1  rss
c			2  tr(I-A)
c   info		error indicator
c			  0 : successful completion
c			 -1 : log10(nlamht) <= lamlim(1) (not fatal)
c			 -2 : log10(nlamht) >= lamlim(2) (not fatal)
c			  1 : svals(1) = 0.0d0
c			  2 : npsing is incorrect
c			  3 : lamlim(1) > lamlim(2)
c
c Subprograms Called Directly:
c	Gcvpack - dvmin 
c
c Subprograms Called Indirectly:
c	Gcvpack - dvl
c
c $Header: dvlop.f,v 2.100.1.1 86/10/07 12:59:38 lindstrom Exp $
c
      integer i,k
      double precision vlamht,w
c	dgmin, dumin added by Yu-Chieh Yang 4/27/01 11:13PM
      double precision dvmin1, dgmin1, dumin1
c
      common / gcvcom / addend,rss,tria,n,h
      integer n,h
      double precision addend,rss,tria,machpr,one
c
      info = 0
      one = 1.0d0
      machpr = 1.0d0
   10 machpr = machpr/2.0d0
      if (one .lt. 1.0d0 + machpr) goto 10
      machpr = machpr*2.0d0
c
      n=nobs
      h=nnull
      addend = inadd
      if (svals(1) .eq. 0.0d0) then
         info = 1
         return
      endif
      k = 0
      do 20 i = 1,npsing
         if (svals(i) .gt. 0) then
            k = i
         endif
   20 continue
      if (k .ne. npsing) then
	 info = 2
    	 return
      endif
      if (job .ne. 0 .and. (lamlim(1) .gt. lamlim(2))) then
         info = 3
         return
      endif
      if (job .eq. 0) then
         lamlim(2) = 2.0d0*dlog10(svals(1))+2.0d0
         lamlim(1) = 2.0d0*dlog10(svals(npsing))-2.0d0
      endif
c	modification added 4/27/01 11:18PM
c       revisit Wed Jul  4 00:00:49 PDT 2001
      if ( vm .eq. 0 ) then
      nlamht = dgmin1 (lamlim(1),lamlim(2),svals,z,npsing,ntbl,tbl,
     * ldtbl,vlamht,info)
      endif
      if ( vm .eq. 1 ) then  
      nlamht = dvmin1 (lamlim(1),lamlim(2),svals,z,npsing,ntbl,tbl,
     * ldtbl,vlamht,info)
      endif
      if ( vm .eq. 2 ) then
      nlamht = dumin1 (lamlim(1),lamlim(2),svals,z,npsing,ntbl,tbl,
     * ldtbl,vlamht,info)
      endif
c	end modified block // 4/27/01 11:19PM
      dout(1) = rss
      dout(2) = tria
c			compute auxtbl
      auxtbl(1,1)=nlamht
      auxtbl(1,2)=vlamht
c			lambda = 0
      auxtbl(2,1)=0.0d0
      auxtbl(2,2)=0.0d0
      if ((nobs-nnull) .ne. npsing) then
         auxtbl(2,2)=inadd*(nobs)/(nobs-nnull-npsing)**2
      endif
      if ((nobs-nnull) .eq. npsing) then
         w=0.0d0
         do 30 i=npsing,1,-1
c           w=w+(z(i)*svals(npsing)**2/(svals(i)**2))**2
            w=w+(z(i)*(svals(npsing)/svals(i))**2)**2
   30    continue
         auxtbl(2,2)=nobs*w
         w=0.0d0
         do 40 i=npsing,1,-1
            w=w+(svals(npsing)/svals(i))**2
   40    continue
         auxtbl(2,2)=auxtbl(2,2)/(w**2)
      endif
c			lambda = infinity
      auxtbl(3,1)=0.0d0
      auxtbl(3,2)=ssqw2/(nobs - nnull)
      nlamht = 10**nlamht
      return
      end
      double precision function dvmin1(lower,upper,svals,z,npsing,
     * ntbl,tbl,ldtbl,vlamht,info)
      integer npsing,ntbl,ldtbl,info
      double precision lower,upper,svals(npsing),z(npsing),
     * tbl(ldtbl,3),vlamht
c
c Purpose: evaluate V(lambda) for a grid of ln(nobs*lambda) values 
c	between	lower and upper, store these in the array tbl, and find
c	minimizer of v.
c
c On Entry:
c   lower		lower bound of interval (in nobs*ln(lambda)
c			scale) over which V(lambda) is to be minimized
c   upper		upper bound of interval (in nobs*ln(lambda)
c			scale) over which V(lambda) is to be minimized
c   svals(npsing)	singular values 
c   z(npsing)		data vector in canonical coordinates
c   npsing		number of positive elements of svals 
c   ntbl		number of evenly spaced values for 
c			ln(nobs*lambda)	to be used in the initial grid 
c			search for lambda hat
c			if ntbl = 0 only a golden ratio search will be 
c			done and tbl is not referenced, if ntbl > 0
c			there will be ntbl rows returned in tbl
c   ldtbl		leading dimension of tbl as declared in the 
c			calling program	
c   
c On Exit:
c   tbl(ldtbl,3)	column	contains
c			  1 	grid of ln(nobs*lambda) 
c			  2  	V(lambda)
c   vlamht		V(lambda hat)
c   dvmin		ln(nobs*lambda hat)
c   info		error indicator
c			  0 : successful completion
c			 -1 : dvmin <= lower (not fatal)
c			 -2 : dvmin >= upper (not fatal)
c
c Subprograms Called Directly:
c	Gcvpack - dvl
c
c $Header: dvmin.f,v 2.100.1.1 86/10/07 12:59:51 lindstrom Exp $
c

      double precision a,b,c,d,vc,vd,del,k1,k2,x,v
      integer j,jmin,k
      double precision dvl1
c				null interval
      if (lower .eq. upper) then
	 dvmin1 = lower
	 info = -1
	 vlamht = dvl1(lower,svals,z,npsing)
	 do 10 j = 1, ntbl
	    tbl(j,1) = lower
	    tbl(j,2) = vlamht
   10    continue
	 return
      end if
c				non-null interval
      info = 0
      a = lower
      b = upper
      if (ntbl .eq. 1) then
	 x = (a + b)/2
	 tbl(1,1) = x
	 tbl(1,2) = dvl1(x,svals,z,npsing)
      else if (ntbl .ge. 2) then
c			do grid search
	 v=dvl1(lower,svals,z,npsing)*2.0d0
	 del=(upper-lower)/(ntbl-1)
	 do 20 j = 1, ntbl
	    tbl(j,1) = lower + (j - 1) * del
	    tbl(j,2) = dvl1(tbl(j,1),svals,z,npsing)
	    if (tbl(j,2) .le. v) then
	       jmin = j
	       v = tbl(j,2)
	    endif
   20    continue	      
	 a=tbl(jmin,1)-del
	 b=tbl(jmin,1)+del
      end if
c			do golden ratio search			
      k1=(3.0d0-dsqrt(5.0d0))/2.0d0
      k2=(dsqrt(5.0d0)-1)/2.0d0
      c = a + k1*(b - a)
      d = a + k2*(b - a)
      vc = dvl1(c,svals,z,npsing)
      vd = dvl1(d,svals,z,npsing)
      do 30 k=1,50
	 if (vd .lt. vc) then
	    a = c
	    c = d
	    d = a + k2*(b - a)
	    vc = vd
	    vd = dvl1(d,svals,z,npsing)
	 else
	    b = d
	    d = c
	    c = a + k1*(b - a)
	    vd = vc
	    vc = dvl1(c,svals,z,npsing)
	 end if
   30 continue
      x=(a+b)/2
      if (x .le. lower) info = -1
      if (x .ge. upper) info = -2
      vlamht=dvl1(x,svals,z,npsing)
      dvmin1 = x
      return
      end
c gcvpack.f, (Aug 10 1995, 17:20)
c this file contains all the subroutines in GCVPACK (release 2).
c ----------------------------------------------------------------------
      subroutine dcfcr (fg,ldfg,nnull,qr,ldqr,npar,pmh,qraux,
     * sgpvt,capz,ldcapz,svals,npsing,v,ldv,nlamht,w1,z,coef,penlty,
     * work,info)
      integer ldfg,nnull,ldqr,npar,pmh,sgpvt(npar),ldcapz,npsing,ldv,
     * info
      double precision fg(ldfg,nnull),qr(ldqr,pmh),
     * qraux(pmh),capz(ldcapz,*),svals(npsing),v(ldv,npsing),nlamht,
     * w1(nnull),z(npsing),coef(npar),penlty,work(nnull)
c
c Purpose: determine the coefficients for a given value of nlamht
c	and vectors z and w1.
c
c On Entry:
c   fg(ldfg,nnull)	information on the Householder transformations
c   			that define f and g
c   ldfg		leading dimension of fg as declared in the 
c			calling	program
c   nnull		number of columns in g
c   qr(ldqr,pmh)	information on the Householder transformations 
c			that define q and r
c   ldqr		leading dimension of qr as declared
c			in the calling program
c   npar		number of rows in q
c   pmh			number of columns in r
c   qraux(pmh)		auxiliary information on the
c			qr Householder transformations
c   sgpvt(npar)		permuted indices from the pivoted Cholesky 
c			decomposition of the matrix which defines the 
c			semi-norm 
c   capz(ldcapz,pmh)	first part of the rotated design matrix
c   ldcapz		leading dimension of capz as declared
c			in the calling program
c   svals(npsing)	singular values 
c   npsing		number of positive singular values
c   v(ldv,npsing)	right singular vectors corresponding to svals
c   ldv			leading dimension of v as declared
c			in the calling program
c   nlamht		nobs*lambda hat
c   w1(nnull)		leading part of rotated response vector
c   z(npsing)		u'w2
c
c On Exit:
c   z(npsing)		g = [(D**2 +nlamht)**-1]Dz
c   coef(npar)		estimated coefficients
c   penlty		smoothness penalty which equals	gamma'gamma
c   info		error indicator
c			  0 : successful completion
c			  1 : error in dtrco, g is singular
c			  2 : error in dtrsl, r is singular
c Work arrays
c   work(nnull)		double precision work vector
c
c   Subprograms Called Directly:
c	Linpack - dqrsl dtrsl dtrco
c	Blas    - ddot dgemv
c	Other   - dprmut 
c
c $Header: dcfcr.f,v 2.100.1.1 86/10/07 12:47:09 lindstrom Exp $
c
      integer i
      double precision dummy(1),machpr,one,rcond
      double precision ddot
c
      info = 0
      one = 1.0d0
      machpr = 1.0d0
   10 machpr = machpr/2.0d0
      if (one .lt. 1.0d0 + machpr) goto 10
      machpr = machpr*2.0d0
c			form gamma-hat and penlty
      do 20 i = 1,npsing
         z(i) = z(i)*svals(i)/(svals(i)**2 + nlamht)
   20 continue
      call dgemv('N',pmh,npsing,1.0d0,v,ldv,z,1,0.0d0,coef,1)
      penlty = ddot(pmh,coef,1,coef,1)
c			form delta
      call dcopy(nnull,w1,1,coef(pmh+1),1)
      call dgemv('N',nnull,pmh,-1.0d0,capz,ldcapz,coef,1,1.0d0,
     *   coef(pmh+1),1)
c			check condition number of g
      if (nnull .ne. 0) then
          call dtrco(fg,ldfg,nnull,rcond,work,1)
          if (rcond .le. machpr*100) then
             info = 1
             return
          endif
          call dtrsl (fg,ldfg,nnull,coef(pmh+1),01,info)
      endif
c			form theta from gamma and delta
      call dtrsl (qr,ldqr,pmh,coef,11,info)
      if (info .ne. 0) then
         info = 2
         return
      endif
      call dqrsl (qr,ldqr,npar,pmh,qraux,coef,coef,dummy,dummy,dummy,
     * dummy,10000,info)
      call dprmut (coef,npar,sgpvt,1)
      return
      end
      subroutine dcfcr1(fg,ldfg,ncts1,fgaux,u,ldu,f1kf2,ldfkf,nuobs,
     * svals,npsing,nlamht,w1,z,coef,penlty,work,info)
      integer ldfg,ldfkf,ncts1,ldu,nuobs,npsing,info
      double precision fg(ldfg,ncts1),fgaux(ncts1),u(ldu,*),
     * f1kf2(ldfkf,*),svals(npsing),nlamht,w1(ncts1),z(npsing),coef(*),
     * penlty,work(*)
c
c Purpose: determine the coefficients for a given value of nlamht
c	and vectors z and w1.
c
c On Entry:
c   fg(ldfg,ncts1)	information on the Householder transformations
c   			that define f and g
c   ldfg		leading dimension of fg as declared in the 
c			calling	program
c   ncts1		number of columns in g
c   fgaux(ncts1)	auxiliary information on the fg Householder
c			transformations
c   u(ldu,npsing)	left singular vectors corresponding to svals
c   ldu	    		leading dimension of u as declared in the 
c			calling	program
c   f1kf2(ldfkf,nuobs-ncts1) f1 k f2
c   ldfkf		leading dimension of f1kf2 as declared
c			in the calling program
c   nuobs		number of rows in fg
c   svals(npsing)	singular values of f2'k f2
c   npsing		number of positive singular
c   nlamht		nobs*(lambda hat)
c   w1(ncts1)		leading part of rotated response vector
c   z(npsing)		u'w2
c
c On Exit:
c   z(npsing)		g = [ (D**2 +nlamht)**-1 ] D z
c   coef(nuobs+ncts1)	estimated coefficients
c   penlty		smoothness penalty which equals	gamma'gamma
c   info		error indicator
c			  0 : successful completion
c			  1 : error in dtrco, g is singular
c
c Work Arrays:
c   work(nuobs-ncts1)  	double precision work vector
c
c   Subprograms Used:
c      Linpack - dqrsl dtrsl 
c      Blas    - ddot dcopy dgemv
c
c $Header: dcfcr1.f,v 2.100.1.1 86/10/07 12:47:31 lindstrom Exp $
c
      integer i,j,nmnct,nctp1,locinf
      double precision dummy(1),machpr,one,rcond
      double precision ddot
c
      info = 0
      one = 1.0d0
      machpr = 1.0d0
   10 machpr = machpr/2.0d0
      if (one .lt. 1.0d0 + machpr) goto 10
      machpr = machpr*2.0d0
c
      nmnct = nuobs - ncts1
      nctp1 = ncts1 + 1
c			form g and penalty
      do 20 j = 1,npsing
         z(j) = z(j)*svals(j)/(svals(j)**2 + nlamht)
   20 continue
      penlty = ddot(npsing,z,1,z,1)
c			z now contains g
c			compute xi 
      do 30 i = 1,npsing
         coef(i) = z(i)/svals(i)
   30 continue
      call dgemv('N',nmnct,npsing,1.0d0,u,ldu,coef,1,0.0d0,
     *  work,1)
      do 40 j = 1,ncts1
         coef(ncts1+j) = 0.0d0
   40 continue
      call dcopy(nmnct,work,1,coef(2*ncts1+1),1)
      call dqrsl(fg,ldfg,nuobs,ncts1,fgaux,coef(nctp1),coef(nctp1),
     * dummy,dummy,dummy,dummy,10000,locinf)
c			compute beta
      call dcopy(ncts1,w1,1,coef,1)
      call dgemv('N',ncts1,nmnct,-1.0d0,f1kf2,ldfkf,work,1,1.0d0,
     *  coef,1)
c			check condition number of g
      call dtrco(fg,ldfg,ncts1,rcond,work,1)
      if (rcond .le. machpr*100) then
         info = 1
         return
      endif
      call dtrsl (fg,ldfg,ncts1,coef,01,info)
      return
      end
      subroutine dcpmut (x,ldx,nobs,npar,jpvt,job)
      integer ldx,nobs,npar,jpvt(npar),job
      double precision x(ldx,npar)
c
c Purpose: permute the columns of the matrix x according to the index 
c	vector jpvt (either forward or backward permutation).
c
c On Entry:
c   x(ldx,npar)		matrix whose columns are to be permuted
c   ldx			leading dimension of x as declared
c			in the calling program
c   nobs		number of rows of x used
c   npar		number of columns of x
c   jpvt(npar)		indices of the permutation
c   job			indicator of forward or backward permutation
c			if job = 0 forward permutation  
c				x(jpvt(i)) moved to x(i)
c			else backward permutation
c				x(i) moved to x(jpvt(i))
c On Exit:
c   x(ldx,npar)		matrix with columns permuted
c
c Subprograms Called Directly
c     Blas	- dswap
c
c  Written:	Yin Ling	U. of Maryland, August,1978
c
c $Header: dcpmut.f,v 2.100.1.1 86/10/07 12:47:45 lindstrom Exp $
c
      integer i,j,k
c
      if (npar .le. 1) then
         return
      endif
      do 10 j = 1,npar
         jpvt(j) = -jpvt(j)
   10 continue
      if (job .eq. 0) then
c		forward permutation
         do 30 i = 1,npar 
            if (jpvt(i) .gt. 0) then
               goto 30
            endif
            j = i
            jpvt(j) = -jpvt(j)
            k = jpvt(j)
c           while
   20       if (jpvt(k) .lt. 0) then
               call dswap (nobs,x(1,j),1,x(1,k),1)
               jpvt(k) = -jpvt(k)
               j = k
               k = jpvt(k)
               goto 20
c           endwhile
            endif
   30    continue
      endif
      if (job .ne. 0) then
c		backward permutation
         do 50 i = 1,npar 
            if (jpvt(i) .gt. 0) then
               goto 50
            endif
            jpvt(i) = -jpvt(i)
            j = jpvt(i)
c           while
   40       if (j .ne. i) then
               call dswap (nobs,x(1,i),1,x(1,j),1)
               jpvt(j) = -jpvt(j)
               j = jpvt(j)
               goto 40
c           endwhile
            endif
   50    continue
      endif
      return
      end
      subroutine dcrtz (x,ldx,nobs,npar,pmh,qr,ldqr,qraux,sgpvt,work,
     * info)
      integer ldx,nobs,npar,pmh,ldqr,sgpvt(npar),info
      double precision x(ldx,npar),qr(ldqr,pmh),qraux(pmh),work(npar)
c
c Purpose: create z from x using permutations and householder 
c	   transformation from sigma
c
c On Entry:
c   x(ldx,npar)		design matrix
c   ldx		     	leading dimension of x as declared in the 
c			calling	program
c   nobs   		number of observations
c   npar      		number of parameters
c   pmh    		npar minus the dimension of the null space of 
c			the semi-norm
c   qr(ldqr,pmh)	information on the qr decomposition of the 
c			Cholesky factor of the semi-norm
c   ldqr    		leading dimension of qr as declared in the 
c			calling	program
c   qraux(pmh)		auxiliary information on the qr decomposition 
c			of the factor of the semi-norm
c   sgpvt(npar) 	permuted indices from the pivoted Cholesky
c			decomposition of the semi-norm matrix
c
c On Exit:
c   x(ldx,npar)		z (x after householder transformation)
c   info    		error indicator
c			  0 : successful completion
c			  1 : error in dtrsl, r is singular
c
c Work Arrays:
c   work(npar)		double precision work vector
c
c Subprograms Called Directly:
c	Linpack - dqrsl dtrsl 
c	Blas    - dcopy 
c	Other   - dcpmut 
c
c Subprograms Called Indirectly:
c	Blas    - dswap
c
c $Header: dcrtz.f,v 2.100.1.1 86/10/07 12:47:59 lindstrom Exp $
c
      double precision dummy(1)
      integer i,locinf
c
      info = 0
c			permute columns of x according to sgpvt
      call dcpmut (x,ldx,nobs,npar,sgpvt,0)
c			apply Householder transformations to rows of
c			permuted x to generate z1 and z2
      do 20 i = 1,nobs 
         call dcopy (npar,x(i,1),ldx,work,1)
         call dqrsl (qr,ldqr,npar,pmh,qraux,work,dummy,work,dummy,dummy,
     *    dummy,01000,locinf)
         call dtrsl (qr,ldqr,pmh,work,01,info)
         if (info .ne. 0) then
            info = 1
            return
         endif
         call dcopy (npar,work,1,x(i,1),ldx)
   20 continue
      return
      end
      subroutine dctsx(ts1,ldts1,ncts1,s2,lds2,ncov2,nobs,tbsb1,ldtbsb,
     * nb,sigma,ldsigm,x,ldx,npar,fgaux,work)
      integer ldts1,ncts1,lds2,ncov2,nobs,ldtbsb,nb,ldsigm,ldx,npar
      double precision ts1(ldts1,ncts1),s2(lds2,*),tbsb1(ldtbsb,ncts1),
     * sigma(ldsigm,*),x(ldx,*),fgaux(ncts1),work(nb)
c
c Purpose: compute x and sigma from t,s and k.
c
c On Entry:
c   ts1(nobs,ncts1)	[t:s1]
c   ldts1		leading dimension of ts1 as declared in the
c			calling	program 
c   ncts1	 	number of columns in [t:s1]
c   s2(lds2,ncov2) 	columns contain the covariates which do not 
c			duplicate the replication pattern of des
c   lds2		leading dimension of s2 as declared in the
c			calling	program 
c   ncov2		number of columns in s2
c   nobs		number of observations
c   tbsb1(ldtbsb,ncts1)	unique rows of t and s1 (or [t:s1] created from
c			basis functions)
c   ldtbsb		leading dimension of tbsb1 as declared in 
c			calling program
c   nb			number of basis functions
c   sigma		ku (ku must start in the ncov2+1 row and
c			the ncov2+1 column)
c   ldsigm		leading dimension of sigma as
c			declared in the calling program
c   x(ldx,npar)		k (k must start in the 1st row and column)
c   ldx			leading dimension of x as
c			declared in the calling program
c
c On Exit:
c   tbsb1(ldtbsb,ncts1)	qr decomposition of tbsb1
c   sigma(ldsigm,npar)  symmetric matrix that defines the semi-norm
c				[0:   0  ]
c				[0:f2'ku f2]
c   x(ldx,npar)		design matrix [t:s1:s2:kf2]
c   npar		number of parameters (nb+ncov2)
c   fgaux(ncts1)	auxiliary vector for qr decomposition of tbsb1
c
c Work Arrays:
c   work(nb)		double precision work vector
c
c Subprograms Called Directly
c	Linpack - dqrdc dqrsl
c	Blas    - dcopy
c	Other   - dset dftkf
c
c Subprograms Called Indirectly
c	Linpack - dqrsl
c	Blas    - dcopy
c
c $Header: dctsx.f,v 2.100.1.1 86/10/07 12:48:06 lindstrom Exp $
c
      double precision dummy(1)
      integer i,locinf,idummy(1)
c
      npar = nb + ncov2
      call dqrdc(tbsb1,ldtbsb,nb,ncts1,fgaux,idummy,dummy,0)
c			calculate k f2  put in last nb -ncts1 columns 
c			of x
      do 10 i=1,nobs 
         call dcopy(nb,x(i,1),ldx,work,1)
         call dqrsl(tbsb1,ldtbsb,nb,ncts1,fgaux,work,dummy,work,dummy,
     *    dummy,dummy,01000,locinf)
         call dcopy(nb-ncts1,work(ncts1+1),1,x(i,ncov2+ncts1+1),ldx)
   10 continue
c			copy [t:s1] into first ncts1 columns of x
      do 20 i=1,ncts1
         call dcopy(nobs,ts1(1,i),1,x(1,i),1)
   20 continue
c			copy s2 into next ncov2 columns of x
      do 30 i=1,ncov2
         call dcopy(nobs,s2(1,i),1,x(1,ncts1+i),1)
   30 continue
c			calculate sigma
      call dftkf(tbsb1,ldtbsb,nb,ncts1,fgaux,sigma(ncov2+1,ncov2+1),
     * ldsigm,work)
      do 40 i = 1,ncts1+ncov2
         call dset(npar,0.0d0,sigma(1,i),1)
   40 continue
      do 50 i = ncts1+ncov2+1,npar
         call dset(ncts1+ncov2,0.0d0,sigma(1,i),1)
   50 continue
      return
      end
      subroutine ddcom (x,ldx,sigma,ldsigm,nobs,npar,nnull,tau,
     * npsing,svals,sgpvt,dcaux,ldcaux,normk,work,lwa,iwork,liwa,
     * job,info)
      integer ldx,ldsigm,nobs,npar,nnull,npsing,sgpvt(npar),
     * ldcaux,lwa,liwa,iwork(liwa),job,info
      double precision x(ldx,npar),sigma(ldsigm,npar),tau,svals(*),
     * dcaux(ldcaux),normk,work(lwa)
c
c Purpose: decompose sigma and the design matrix.
c
c On Entry:
c   x(ldx,npar)		design matrix
c   ldx			leading dimension of x as declared in the 
c			calling program. Must be at least max(npar,nobs)
c   sigma(ldsigm,npar)	symmetric matrix that defines the semi-norm
c   ldsigm		leading dimension of sigma as
c			declared in the calling program
c   nobs		number of observations
c   npar		number of parameters
c   nnull		dimension of the null space of sigma
c   tau			multiplier controlling the amount of truncation
c			if truncation is requested
c   ldcaux		length of auxiliary vector as declared in the 
c			calling program. 
c			Must be at least (npar-nnull)**2+2*npar-nnull
c   job 		job nonzero use truncation 
c			job = 0 no truncation
c
c On Exit:
c   x(ldx,npar)		overwritten with many intermediate results
c			must be	passed to dgcv intact
c   sigma(ldsigm,npar)	overwritten with the qr decomposition of the 
c			Cholesky factor of sigma
c			must be passed intact to dgcv
c   nnull		if input nnull is too small it is replaced by
c			larger value such that sigma has rank npar-nnull
c   npsing       	number of singular values calculated if info = 0
c			otherwise npsing contains nonzero info from 
c			dsvdc.
c   svals(npar-nnull)	singular values of the matrix j2 if info = 0
c			if info from dsvdc is nonzero then svals 
c			is as it was returned from dsvdc.
c   sgpvt(npar)		permuted indices from the Cholesky decomposition
c			with pivots of sigma
c			must be passed intact to dgcv 
c   dcaux(ldcaux)	auxiliary vector, must be passed intact to dgcv 
c   normk		frobenius norm of the k by k lower right sub 
c			matrix of j2
c   info		error indicator
c			   0 : successful completion
c			  -3 : nnull is too small (not fatal)
c			   1 : dimension error or tau < 0
c			   2 : ldcaux (length of dcaux) is too small
c			   3 : lwa (length of work) is too small
c			   4 : liwa (length of iwork) is too small
c			   10 < info < 20 : 10 + nonzero info returned
c					    from dsgdc
c			   20 < info < 30 : 20 + nonzero info returned 
c					    from dcrtz
c			   30 < info < 40 : 30 + nonzero info returned 
c					    from dzdc
c
c Work Arrays:
c   work(lwa)		double precision work vector
c   lwa			length of work as declared in the calling 
c			program Must be at least 
c			(npar-nnull)*(nobs-nnull+1)+nobs
c   iwork(liwa)		integer work vector
c   liwa		length of iwork as declared in the calling
c			program 
c			must be at least npar-nnull
c
c Subprograms Called Directly:
c	Gcvpack - dsgdc dcrtz dzdc
c	Blas	- dcopy
c
c Subprograms Called Indirectly:
c	Gcvpack - dtsvdc
c	Linpack - dtrco dchdc dqrdc dqrsl dtrsl dsvdc
c	Blas    - dcopy dswap ddot
c	Other   - dprmut dcpmut dset
c
c $Header: ddcom.f,v 2.100.1.1 86/10/07 12:48:26 lindstrom Exp $
c
      integer pmh,pp1,pmhp1,wsize,nmh,sinfo
c
      sinfo = 0
      info = 0
c			check dimensions
      nmh = nobs - nnull
      pmh = npar - nnull
      pmhp1 = pmh + 1
      pp1 = npar + 1
      if ((nobs .le. 0) .or. (npar .le. 0) .or. (nnull .lt. 0) 
     * .or. (pmh .le. 0) .or. (ldx .lt. nobs).or. (ldx .lt. npar)
     * .or. (tau .lt. 0)) then
         info = 1
         return
      endif
      wsize = pmh**2 + 2*npar -nnull
      if (ldcaux .lt. wsize) then
         info = 2
         return
      endif
      wsize = pmh*nmh + pmh + nobs
      if (lwa .lt. wsize) then
         info = 3
         return
      endif
      if (liwa .lt. pmh) then
         info = 4
         return
      endif
c			decompose sigma
      call dsgdc (sigma,ldsigm,npar,nnull,dcaux,sgpvt,info)
      if (info.gt. 0) then
         info = info + 10
         return
      endif
      if (info .lt. 0) then
	  sinfo = info
	  pmh = npar - nnull
          pmhp1 = pmh + 1
      endif
c			create z
      call dcrtz(x,ldx,nobs,npar,pmh,sigma,ldsigm,dcaux,sgpvt,work,info)
      if (info .gt. 0) then
         info = info + 20
         return
      endif
c			create j and decompose j2
      call dzdc(x,ldx,nobs,npar,pmh,tau,dcaux(pmhp1),dcaux(pp1),
     * npsing,dcaux(npar+pmhp1),pmh,normk,work,lwa,iwork,liwa,job,info)
c			copy svals
      call dcopy(pmh,dcaux(pp1),1,svals,1)
      if (info .gt. 0) then
         info = info + 30
         return
      endif
      if (sinfo .lt. 0) info = sinfo
      return
      end
      subroutine ddiag(fg,ldfg,nobs,nnull,fgaux,svals,npsing,u,ldu,
     * nlamht,adiag,work)
      integer ldfg,nobs,nnull,npsing,ldu
      double precision fg(ldfg,nnull),fgaux(nnull),svals(npsing),
     * u(ldu,npsing),nlamht,adiag(nobs),work(*)
c
c Purpose: determine the diagonal of hat matrix for nobs*lamhat
c
c On Entry:
c   fg(ldfg,nnull)	information on the Householder transformations
c			that define f and g
c   ldfg		leading dimension of fg as declared in the 
c			calling	program
c   nobs		number of rows in f
c   nnull		number of columns in g
c   fgaux(nnull)	auxiliary information on the fg Householder
c			transformations
c   svals(npsing)	singular values 
c   npsing		number of positive singular values 
c   u(ldu,npsing)	left singular vectors corresponding to svals
c   ldu	    		leading dimension of u as declared in the 
c			calling	program
c   nlamht		nobs*lambda hat
c
c On Exit:
c   adiag(nobs)		diagonal elements of the hat matrix if requested
c
c Work Arrays:
c   work(nobs+npsing)	double precision work vector
c
c Subprograms Called Directly:
c	Linpack - dqrsl 
c	Blas    - ddot dgemv
c	Other   - dset
c
c $Header: ddiag.f,v 2.100.1.1 86/10/07 13:28:57 lindstrom Exp $
c
      integer i,j,hp1,locinf,nmh,np1
      double precision dummy(1)
      double precision ddot
c
      np1 = nobs + 1
      hp1 = nnull + 1
      nmh = nobs - nnull
c			form adiag
      do 20 i = 1,nobs 
         call dset(nobs,0.0d0,work,1)
         work(i)=1.0d0
         call dqrsl(fg,ldfg,nobs,nnull,fgaux,work,dummy,work,dummy, 
     *    dummy,dummy,01000,locinf)
         adiag(i)=ddot(nnull,work,1,work,1)
	 call dgemv('T',nmh,npsing,1.0d0,u,ldu,work(hp1),1,0.0d0,
     *    work(np1),1)
         do 10 j=1,npsing
            work(nobs+j)=work(nobs+j)*svals(j)/dsqrt(svals(j)**2+nlamht)
   10    continue
         adiag(i)=adiag(i) + ddot(npsing,work(np1),1,work(np1),1)
   20 continue
      return
      end
      subroutine dftkf(fg,ldfg,nrf,ncg,fgaux,kk,ldkk,work)
      integer ldfg,nrf,ncg,ldkk
      double precision fg(ldfg,ncg),fgaux(ncg),kk(ldkk,nrf),work(nrf)
c
c Purpose: create f'k f.
c
c On Entry:
c   fg(ldfg,ncg)	qr decomposition of [t:s1]
c   ldfg		leading dimension of fg as declared in the
c   			calling program 
c   nrf 		number of rows in f
c   ncg			number of columns in g
c   fgaux(ncg)		auxiliary information on the qr decomposition
c			of [t:s1]
c   kk(ldkk,nrf) 	k
c   ldkk		leading dimension of kk as declared in the
c   			calling program 
c
c On Exit:
c   kk(ldkk,nrf)  	f'k f
c
c Work Array:
c   work(nrf)		double precision work vector
c
c Subprograms Called Directly:
c	Linpack - dqrsl
c	Blas    - dcopy
c
c $Header: dftkf.f,v 2.100.1.1 86/10/07 12:48:38 lindstrom Exp $
c
      double precision dummy(1)
      integer i,locinf
c	  		calculate k f, store in kk
      do 10 i=1,nrf 
         call dcopy(nrf,kk(i,1),ldkk,work,1)
         call dqrsl(fg,ldfg,nrf,ncg,fgaux,work,dummy,work,dummy,dummy,
     *    dummy,01000,locinf)
         call dcopy(nrf,work,1,kk(i,1),ldkk)
   10 continue
c	  		calculate f'k f
      do 20 i=1,nrf 
         call dqrsl(fg,ldfg,nrf,ncg,fgaux,kk(1,i),dummy,kk(1,i),dummy,
     *    dummy,dummy,01000,locinf)
   20 continue
      return
      end
      subroutine dgcv(x,ldx,y,sigma,ldsigm,nobs,npar,nnull,sgpvt,dcaux,
     * ldcaux,adiag,lamlim,ntbl,dout,coef,tbl,ldtbl,auxtbl,work,lwa,
     * job,info)
      integer ldx,ldsigm,nobs,npar,nnull,sgpvt(npar),ldcaux,ntbl,ldtbl,
     * lwa,info,job
      double precision x(ldx,npar),y(nobs),sigma(ldsigm,npar),
     * dcaux(ldcaux),adiag(nobs),lamlim(2),dout(4),coef(npar),
     * tbl(ldtbl,3),auxtbl(3,3),work(lwa)
c
c Purpose: determine the generalized cross validation estimate of the 
c	smoothing parameter and fit model parameters for a penalized
c 	least squares problem with a semi-norm smoothing matrix.
c
c On Entry:
c   x(ldx,npar)		design matrix as returned by ddcom
c   ldx			leading dimension of x as declared in the 
c			calling program, must be at least max(nobs,npar)
c   y(nobs)		response vector
c   sigma(ldsigm,npar)	Cholesky factor of sigma as returned by ddcom
c   ldsigm		leading dimension of sigma as
c			declared in the calling program
c   nobs		number of observations
c   npar		number of parameters
c   nnull		dimension of the null space of sigma
c   sgpvt(npar)		permuted indices from the pivoted Cholesky
c			decomposition of the semi-norm matrix
c   dcaux(ldcaux)	auxiliary vector from ddcom
c   ldcaux		length of ldcaux. Must be at least
c			(npar-nnull)**2+2*npar-nnull
c   adiag(nobs)	 	"true" y values on entry if predictive mse if 
c			requested 
c   lamlim(2)		limits on lambda hat search (in log10(nobs*
c			lambda)	scale) if user input limits are 
c			requested. if lamlim(1) = lamlim(2) then lamhat
c			is set to (10**lamlim(1))/nobs
c   ntbl		number of evenly spaced values for 
c			log10(nobs*lambda) to be used in the initial 
c			grid search for lambda hat
c			if ntbl = 0 only a golden ratio search will be 
c			done and tbl is not referenced, if ntbl > 0
c			there will be ntbl rows returned in tbl
c   ldtbl		leading dimension of tbl as declared in the 
c			calling program	
c   job			integer with decimal expansion abc
c			if a is nonzero then predictive mse is computed
c			   using adiag as true y
c			if b is nonzero then user input limits on search
c			   for lambda hat are used
c			if c is nonzero then diagonal of the hat matrix
c			   is calculated
c
c On Exit:
c   y(nobs)		predicted values
c   adiag(nobs)		diagonal elements of the hat matrix if requested
c   lamlim(2)		limits on lambda hat search 
c			(in log10(nobs*lambda) scale)
c   dout(4)		contains:
c  			1  lamhat   generalized cross validation
c				    estimate of the smoothing parameter
c			2  penlty   smoothing penalty
c			3  rss	    residual sum of squares
c			4  tr(I-A)  trace of I - A
c   coef(npar)		coefficient estimates
c   tbl(ldtbl,3)	column	contains
c			  1 	grid of log10(nobs*lambda) 
c			  2  	V(lambda)
c			  3     R(lambda) if requested
c   auxtbl(3,3)		auxiliary table
c			1st row contains:
c			    log10(nobs*lamhat), V(lamhat) and  
c			    R(lamhat) if requested
c			    where lamhat is the gcv estimate of lambda
c			2nd row contains:
c			    0, V(0) and  R(0) if requested
c			3rd row contains:
c			    0, V(infinity) and R(infinity) if requested
c   info		error indicator
c			  0 : successful completion
c			 -1 : log10(nobs*lamhat) <= lamlim(1) 
c			      (not fatal)
c			 -2 : log10(nobs*lamhat) >= lamlim(2) 
c			      (not fatal)
c			  1 : dimension error
c			  2 : error in ntbl
c			  3 : ldcaux (length of dcaux) is too small
c			  4 : lwa (length of work) is too small
c			  5 : lamlim(1) > lamlim(2)
c			 10 < info < 20 : 10 + nonzero info returned
c					  from dvlop
c			 20 < info < 30 : 20 +nonzero info returned
c					  from dcfcr
c
c Work Arrays:
c   work(lwa) 		double precision work vector
c   lwa			length of work as declared in the calling 
c			program must be at least (npar-nnull)+nobs
c
c Subprograms Called Directly:
c	Gcvpack - drsap dvlop dpmse dcfcr dpdcr ddiag
c
c Subprograms Called Indirectly:
c	Gcvpack - dvmin dvl
c	Linpack - dqrsl dtrsl
c	Blas    - dcopy ddot dgemv
c	Other   - dprmut dset
c
c $Header: dgcv.f,v 2.100.1.1 86/10/07 12:49:16 lindstrom Exp $
c
      double precision addend,nlamht,ssqw2
      integer mp1,mnh,pmh,minnp,pp1,pmhp1,wsize,npsing,i,jpmse,jlaml,
     * jadiag,sinfo
c
      jpmse = job/100
      jlaml = mod(job,100)/10
      jadiag = mod(job,10)
c			check dimensions
      sinfo = 0
      info = 0
      pmh = npar - nnull
      if ((nobs .le. 0) .or. (npar .le. 0) .or. (nnull .lt. 0) 
     * .or. (pmh .le. 0)) then
         info = 1
         return
      endif
      if ((ntbl .lt. 0) .or. (ntbl .gt. ldtbl)) then
         info = 2
         return
      endif
      wsize = pmh**2 + 2*npar - nnull
      if (ldcaux .lt. wsize) then
         info = 3
         return
      endif
      wsize = pmh + nobs
      if (lwa .lt. wsize) then
         info = 4
         return
      endif
      if (jlaml .ne. 0 .and. (lamlim(1) .gt. lamlim(2))) then
         info = 5
         return
      endif
      pmhp1 = pmh + 1
      pp1 = npar + 1
c			calculate npsing
      minnp = min(nobs-nnull,pmh)
      do 30 i=1,minnp
         if (dcaux(npar+i)**2 .gt. 0) npsing = i
   30 continue
c			apply rotations to y
      mp1 = nnull + 1
      mnh = nobs - nnull
      call drsap (x(1,pmhp1),ldx,nobs,nnull,dcaux(pmhp1),x(mp1,1),ldx,
     * mnh,npsing,y,ssqw2,addend,work)
c			minimize V(lambda)
      call dvlop (y(mp1),dcaux(pp1),nobs,nnull,npsing,addend,ssqw2,
     * lamlim,ntbl,nlamht,tbl,ldtbl,auxtbl,dout(3),jlaml,info)
      dout(1)=nlamht/nobs
      if (info .gt. 0) then
         info = info + 10
         return
      endif
      if (info .lt. 0) sinfo = info
c			calculate predictive mse
      if (jpmse .ne. 0) then
         call dpmse(x(1,pmhp1),ldx,nobs,nobs,nnull,dcaux(pmhp1),
     *    dcaux(pp1),npsing,x(mp1,1),ldx,y,y(mp1),ntbl,adiag,tbl,ldtbl,
     *    auxtbl,work)
      endif
c			calculate coefficients
      call dcfcr (x(1,pmhp1),ldx,nnull,sigma,ldsigm,npar,
     * pmh,dcaux,sgpvt,x,ldx,dcaux(pp1),npsing,dcaux(npar+pmhp1),pmh, 
     * nlamht,y,y(mp1),coef,dout(2),work,info)
      if (info .gt. 0) then
         info = info + 20
         return
      endif
c			calculate predicted values 
      call dpdcr (x(1,pmhp1),ldx,nobs,nnull,dcaux(pmhp1),dcaux(pp1),
     * npsing,x(mp1,1),ldx,nlamht,y,y(mp1),y,work)
      if (jadiag .ne. 0) then
          call ddiag (x(1,pmhp1),ldx,nobs,nnull,dcaux(pmhp1),
     * 	  dcaux(pp1),npsing,x(mp1,1),ldx,nlamht,adiag,work)
      endif
      if (sinfo .lt. 0) info = sinfo
      return
      end
      subroutine dgcv1(fkf,ldfkf,y,nuobs,nobs,fg,ldfg,ncts1,fgaux,svals,
     * adiag,lamlim,ssqrep,ntbl,dout,coef,tbl,ldtbl,auxtbl,work,lwa,
     * job,info)
      integer ldfkf,nuobs,nobs,ldfg,ncts1,ntbl,ldtbl,lwa,job,info
      double precision fkf(ldfkf,nuobs),y(nuobs),fg(ldfg,ncts1),
     * fgaux(ncts1),svals(*),adiag(nuobs),lamlim(2),ssqrep,dout(4),
     * coef(*),tbl(ldtbl,3),auxtbl(3,3),work(lwa)
c
c Purpose: determine the generalized cross validation estimate of the 
c	smoothing parameter and fit model parameters for a semi-norm
c	thin plate spline model.
c
c On Entry:
c   fkf(ldfkf,nuobs) 	intermediate results as created by dsgdc1
c   ldfkf		leading dimension of fkf as declared in the 
c			calling program
c   y(nuobs)		B1'y 
c   nuobs		number of rows in fg
c   nobs		number of observations
c   fg(ldfg,ncts1)	qr decomposition of [t:s1]	
c   ldfg		leading dimension of fg as
c			declared in the calling program
c   ncts1		number of columns in [t:s1]	
c   fgaux(ncts1)	auxiliary information on the decomposition of 
c			[t:s1]
c   svals(nuobs-ncts1) 	positive singular values of the Cholesky factor
c			of f2'k f2
c   adiag(nuobs)	B1'(true y) if predictive mse is requested
c   lamlim(2)		limits on lambda hat search (in log10(nobs*
c			lambda) scale) if user input limits are 
c			requested. if lamlim(1) = lamlim(2) then lamhat
c			is set to (10**lamlim(1))/nobs
c   ssqrep		sum of squares for replication
c   ntbl		number of evenly spaced values for 
c			log10(nobs*lambda) to be used in the initial 
c			grid search for lambda hat
c			if ntbl = 0 only a golden ratio search will be 
c			done and tbl is not referenced, if ntbl > 0
c			there will be ntbl rows returned in tbl
c   ldtbl		leading dimension of tbl as declared in the 
c			calling program	
c   job			integer with decimal expansion abc
c			if a is nonzero then predictive mse is computed
c			   using adiag as true y
c			if b is nonzero then user input limits on search
c			   for lambda hat are used
c			if c is nonzero then diagonal of the hat matrix
c			   is calculated
c
c On Exit:
c   y(nuobs)		B1'(predicted values)
c   adiag(nuobs)	diagonal elements of the hat matrix if requested
c   lamlim(2)		limits on lambda hat search 
c			(in log10(nobs*lambda) scale)
c   dout(4)		contains:
c  			1  lamhat   generalized cross validation 
c				    estimate of the smoothing parameter
c			2  penlty   smoothing penalty
c			3  rss	    residual sum of squares
c			4  tr(I-A)  trace of I - A
c   coef(nuobs+ncts1) 	estimated coefficients
c   tbl(ldtbl,3)	column	contains
c			  1 	grid of log10(nobs*lambda) 
c			  2  	V(lambda)
c			  3     R(lambda) if requested
c   auxtbl(3,3)		auxiliary table
c			1st row contains:
c			    log10(nobs*lamhat), V(lamhat) and  
c			    R(lamhat) if requested
c			    where lamhat is the gcv estimate of lambda
c			2nd row contains:
c			    0, V(0) and  R(0) if requested
c			3rd row contains:
c			    0, V(infinity) and R(infinity) if requested
c   info		error indicator
c			  0 : successful completion
c			 -1 : log10(nobs*lamhat) <= lamlim(1) 
c			      (not fatal)
c			 -2 : log10(nobs*lamhat) >= lamlim(2)
c			      (not fatal)
c			  1 : dimension error
c			  2 : error in ntbl
c			  3 : lwa (length of work) is too small
c			  4 : lamlim(1) > lamlim(2)
c			 10 < info < 20 : 10 + nonzero info returned 
c					  from dvlop
c			 20 < info < 30 : 20 + nonzero info returned
c					  from dcfcr1
c
c Working Storage:
c   work(lwa) 		double precision work vector
c   lwa			length of work as declared in the calling 
c			program
c			must be at least nuobs-ncts1+nobs
c
c Subprograms Called Directly:
c       Gcvlib - drsap dvlop dpmse dcfcr1 dpdcr ddiag
c
c Subprograms Called Indirectly:
c	Gcvlib  - dvl vmin
c	Linpack - dqrsl dtrsl
c	Blas    - ddot dcopy dgemv
c
c $Header: dgcv1.f,v 2.100.1.1 86/10/07 12:49:43 lindstrom Exp $
c
      double precision addend,nlamht,ssqw2
      integer npsing,i,jpmse,jlaml,jadiag,nmnct,nctp1,sinfo
c
c
      sinfo = 0
      info = 0
      nmnct = nuobs - ncts1
      nctp1=ncts1+1
      jpmse = job/100
      jlaml = mod(job,100)/10
      jadiag = mod(job,10)
c			check dimensions
      if ((nuobs.le.0).or.(ncts1 .le. 0).or.(nmnct .le. 0)) then
         info = 1
         return
      endif
      if ((ntbl .lt. 0) .or. (ntbl .gt. ldtbl)) then
         info = 2
         return
      endif
      if (lwa .lt. nobs+nuobs - nmnct) then
         info = 3
         return
      endif
      if (jlaml .ne. 0 .and. (lamlim(1) .gt. lamlim(2))) then
         info = 4
         return
      endif
c			calculate npsing
      do 30 i=1,nmnct
         if (svals(i)**2 .gt. 0.0d0) npsing = i
   30 continue
c			apply rotations to y
      call drsap(fg,ldfg,nuobs,ncts1,fgaux,fkf(nctp1,nctp1),ldfkf,nmnct,
     * npsing,y,ssqw2,addend,work)
      addend = addend + ssqrep
      ssqw2 = ssqw2 + ssqrep
c			minimize V(lambda)
      call dvlop (y(nctp1),svals,nobs,ncts1,npsing,addend,ssqw2,
     * lamlim,ntbl,nlamht,tbl,ldtbl,auxtbl,dout(3),jlaml,info)
      dout(1)=nlamht/nobs
      if (info .gt. 0) then
         info = info + 10
         return
      endif
      if (info .lt. 0)  sinfo = info
c			calculate predictive mse
      if (jpmse .ne. 0) then
         call dpmse(fg,ldfg,nuobs,nobs,ncts1,fgaux,svals,npsing,
     *    fkf(nctp1,nctp1),ldfkf,y,y(nctp1),ntbl,adiag,tbl,ldtbl,
     *    auxtbl,work)
      endif
c			calculate coefficients
      call dcfcr1(fg,ldfg,ncts1,fgaux,fkf(nctp1,nctp1),ldfkf,
     * fkf(1,nctp1),ldfkf,nuobs,svals,npsing,nlamht,y,y(nctp1),coef,
     * dout(2),work,info)
      if (info .gt. 0) then
         info = info + 20
         return
      endif
      call dpdcr (fg,ldfg,nuobs,ncts1,fgaux,svals,npsing,
     * fkf(nctp1,nctp1),ldfkf,nlamht,y,y(nctp1),y,work)
      if (jadiag .ne. 0) then
          call ddiag(fg,ldfg,nuobs,ncts1,fgaux,svals,npsing,
     *     fkf(nctp1,nctp1),ldfkf,nlamht,adiag,work)
      endif
      if (sinfo .lt. 0) info = sinfo
      return
      end
Caveat receptor.  (Jack) dongarra@anl-mcs, (Eric Grosse) research!ehg
Compliments of netlib   Sun Jul  6 09:34:18 CDT 1986
C
C***********************************************************************
C
C     File of the  DOUBLE PRECISION  Level 2 BLAS routines:
C
C      DGEMV, DGBMV, DSYMV, DSBMV, DSPMV, DTRMV, DTBMV, DTPMV,
C      DGER , DSYR , DSPR ,
C      DSYR2, DSPR2,
C      DTRSV, DTBSV, DTPSV.
C
C     See:
C
C        Dongarra J. J., Du Croz J. J., Hammarling S. and Hanson R. J..
C        A proposal for an extended set of Fortran Basic Linear Algebra
C        Subprograms. Technical Memorandum No.41 (revision 1),
C        Mathematics and Computer Science Division, Argone National
C        Laboratory, 9700 South Cass Avenue, Argonne, Illinois 60439,
C        USA or NAG Technical Report TR4/85, Nuemrical Algorithms Group
C        Inc., 1101 31st Street, Suite 100, Downers Grove, Illinois
C        60515-1263, USA.
C
C***********************************************************************
C
      SUBROUTINE DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
      CHARACTER*1        TRANS
      INTEGER            M, N, LDA, INCX, INCY
      DOUBLE PRECISION   ALPHA, A( LDA, * ), X( * ), BETA, Y( * )
*
*  Purpose
*  =======
*
*  DGEMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*
*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
*
*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y
*.
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the leading dimension of A as
*           declared in the calling (sub) program. LDA must be at least
*           m.
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*  BETA   - REAL            .
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry with BETA non-zero, the incremented array Y
*           must contain the vector y. On exit, Y is overwritten by the
*           updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y.
*           Unchanged on exit.
*
*
*  Note that TRANS, M, N and LDA must be such that the value of the
*  LOGICAL variable OK in the following statement is true.
*
*     OK = ( ( TRANS.EQ.'N' ).OR.( TRANS.EQ.'n' ).OR.
*    $       ( TRANS.EQ.'T' ).OR.( TRANS.EQ.'t' ).OR.
*    $       ( TRANS.EQ.'C' ).OR.( TRANS.EQ.'c' )     )
*    $     .AND.
*    $     ( M.GE.0 )
*    $     .AND.
*    $     ( N.GE.0 )
*    $     .AND.
*    $     ( LDA.GE.M )
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 30-August-1985.
*     Sven Hammarling, Nag Central Office.
*
      INTEGER            I     , IX    , IY    , J     , JX    , JY
      INTEGER            KX    , KY    , LENX  , LENY
      DOUBLE PRECISION   ONE   ,         ZERO
      PARAMETER        ( ONE   = 1.0D+0, ZERO  = 0.0D+0 )
      DOUBLE PRECISION   TEMP
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.
     $    ( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     Set LENX and LENY, the lengths of the vectors x and y.
*
      IF( ( TRANS.EQ.'N' ).OR.( TRANS.EQ.'n' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y  and set up the start points in X and Y if
*     the increments are not both unity.
*
      IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
         IF( BETA.NE.ONE )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         END IF
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( LENX - 1 )*INCX
         END IF
         IF( INCY.GT.0 )THEN
            KY = 1
         ELSE
            KY = 1 - ( LENY - 1 )*INCY
         END IF
         IF( BETA.NE.ONE )THEN
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( ( TRANS.EQ.'N' ).OR.( TRANS.EQ.'n' ) )THEN
*
*        Form  y := alpha*A*x + y.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 60, J = 1, N
               IF( X( J ).NE.ZERO )THEN
                  TEMP = ALPHA*X( J )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
   60       CONTINUE
         ELSE
            JX = KX
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
*
*        Form  y := alpha*A'*x + y.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 100, J = 1, N
               TEMP = ZERO
               DO 90, I = 1, M
                  TEMP  = TEMP + A( I, J )*X( I )
   90          CONTINUE
               Y( J ) = Y( J ) + ALPHA*TEMP
  100       CONTINUE
         ELSE
            JY = KY
            DO 120, J = 1, N
               TEMP = ZERO
               IX   = KX
               DO 110, I = 1, M
                  TEMP  = TEMP + A( I, J )*X( IX )
                  IX    = IX   + INCX
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
      RETURN
*
*     End of DGEMV .
*
      END
      subroutine dmakek(m,n,dim,des,lddes,nb,desb,lddesb,kk,ldkk)
      integer m,n,dim,lddes,nb,lddesb,ldkk
      double precision des(lddes,dim),desb(lddesb,dim),kk(ldkk,nb)
c
c Purpose: create the k matrix.
c   
c On Entry:
c   m			order of the derivatives in the penalty
c   n			number of rows in des
c   dim			dimension of the space to be splined 
c   des(lddes,dim)	variables to be splined		
c   lddes		leading dimension of des as declared in the 
c			calling	program
c   nb			number of rows in desb 
c   desb(lddesb,dim)	positions of unique design points or basis 
c			functions
c   lddesb		leading dimension of desb as declared in the 
c			calling	program
c   ldkk		leading dimension of kk as declared in the
c			calling	program
c On Exit:
c   kk(ldkk,nb)		k matrix
c
c Subprograms Called:
c	Other   - fact
c
c $Header: dmakek.f,v 2.100.1.1 86/10/07 12:50:38 lindstrom Exp $
c
      integer i,j,k,fact
      double precision tauij,expo,theta,t,pi
c
c			t to be used in computation of theta
      pi = 4.0d0*atan(1.0d0)
      t = 2 ** (2*m) * pi**(dim/2.0d0) * fact(m-1)
c
c 			exponent for tauij
      expo = m - (dim / 2.0d0)
      if (dim .eq. 2*(dim/2)) then
c			1+dim odd
         theta = 1.0 / (0.5 * t * fact (m-dim/2))
         if ((2*m+dim) .eq. 4*((2*m+dim)/4)) theta = -theta
         do 30 i=1,n 
            do 20 j=1,nb 
               tauij = 0
               do 10 k=1,dim
                  tauij = tauij + (des(i,k)-desb(j,k))**2
   10          continue
               if (tauij .eq. 0.0d0) then
                  kk(i,j) = 0.0d0
               else
                  kk(i,j) = theta*tauij**expo * 0.5 * log(tauij)
               endif
   20       continue
   30    continue
      else
c			1+dim even
c			compute theta
c			compute gamma(dim/2 - m)
         j = (1 - (dim-2*m)) / 2
         theta = sqrt(pi)
         do 40 i=1,j
	      theta = -theta / (i - 0.5d0)
   40    continue
	 theta = theta / t

         do 70 i=1,n 
            do 60 j=1,nb 
               tauij = 0
               do 50 k=1,dim
                  tauij = tauij + (des(i,k)-desb(j,k))**2
   50          continue
               if (tauij .eq. 0.0d0) then
                  kk(i,j) = 0.0d0
               else
                  kk(i,j) = theta*tauij**expo
               endif
   60       continue
   70    continue
      endif
      end
      subroutine dmaket(m,n,dim,des,lddes,s1,lds1,ncov1,npoly,t,ldt,
     * wptr,info)
      integer m,n,dim,lddes,lds1,ncov1,npoly,ldt,wptr(dim),info
      double precision des(lddes,dim),s1(lds1,*),t(ldt,*)
c
c Purpose: create t matrix and append s1 to it.
c
c On Entry:
c   m			order of the derivatives in the penalty
c   n			number of rows in des
c   dim			number of columns in des
c   des(lddes,dim)	variables to be splined
c   lddes		leading dimension of des as declared in the 
c			calling program
c   s1(lds1,ncov1)	covariates which duplicate the replication 
c			structure of des
c   lds1		leading dimension of s1 as declared in the 
c			calling program
c   ncov1		number of columns in s1
c   ldt			leading dimension of t as declared in the 
c			calling program
c
c On Exit:
c   npoly		dimension of polynomial part of spline
c   t(ldt,npoly+ncov1)	[t:s1]
c   info 		error indication
c   			   0 : successful completion
c		 	   1 : error in creation of t
c Work Arrays:
c   wptr(dim)		integer work vector
c
c Subprograms Called Directly:
c	Blas  - dcopy
c	Other - mkpoly
c
c $Header: dmaket.f,v 2.100.1.3 86/11/21 11:32:19 lindstrom Exp $
c
      integer i,j,k,tt,nt,bptr,eptr
      integer mkpoly
c
      info = 0
      npoly = mkpoly(m,dim)
      call dset(n,1.0d0,t(1,1),1)
      nt = 1
      if (npoly .gt. 1) then
          do 10 j=1,dim 
             nt = j + 1
             wptr(j) = nt
             call dcopy(n,des(1,j),1,t(1,nt),1)
   10     continue
c
c     get cross products of x's in null space for m>2
c
c     WARNING: do NOT change next do loop unless you fully understand:
c              This first gets x1*x1, x1*x2, x1*x3, then
c              x2*x2, x2*x3, and finally x3*x3 for dim=3,n=3
c              wptr(1) is always at the beginning of the current
c	       level of cross products, hence the end of the
c	       previous level which is used for the next.
c	       wptr(j) is at the start of xj * (previous level)
c
          do 50 k=2,m-1 
             do 40 j=1,dim 
                bptr = wptr(j)
                wptr(j) = nt + 1
                eptr = wptr(1) - 1
                do 30 tt=bptr,eptr 
                   nt = nt + 1
                   do 20 i=1,n
                      t(i,nt) = des(i,j) * t(i,tt)
   20              continue
   30           continue
   40        continue
   50     continue
          if (nt .ne. npoly) then
	      info = 1
	      return
          endif
      endif
c			append s1 to t
      do 60 i = 1,ncov1
         call dcopy(n,s1(1,i),1,t(1,nt+i),1)
   60 continue
      end
      subroutine dpdcr(fg,ldfg,nobs,nnull,fgaux,svals,npsing,u,ldu, 
     * nlamht,w1,g,pred,work)
      integer ldfg,nobs,nnull,npsing,ldu
      double precision fg(ldfg,nnull),fgaux(nnull),svals(npsing),
     * u(ldu,npsing),nlamht,w1(nnull),g(npsing),pred(nobs),
     * work(*)
c
c Purpose: determine the predicted responses for a given value of 
c	nobs*lamhat and vectors g and w1.
c
c On Entry:
c   fg(ldfg,nnull)	information on the Householder transformations
c			that define f and g
c   ldfg		leading dimension of fg as declared in the 
c			calling	program
c   nobs		number of rows in f
c   nnull		number of columns in g
c   fgaux(nnull)	auxiliary information on the fg Householder
c			transformations
c   svals(npsing)	singular values 
c   npsing		number of positive singular values 
c   u(ldu,npsing)	left singular vectors corresponding to svals
c   ldu	    		leading dimension of u as declared in the 
c			calling	program
c   nlamht		nobs*lambda hat
c   w1(nnull)		leading part of rotated response vector
c   g(npsing)		(D**2 + nlamht*I)*-1 Dz
c
c On Exit:
c   pred(nobs)		predicted responses
c
c Work Arrays:
c   work(nobs+npsing)	double precision work vector
c
c Subprograms Called Directly:
c	Linpack - dqrsl 
c	Blas    - dcopy dgemv
c
c $Header: dpdcr.f,v 2.100.1.1 86/10/07 12:51:16 lindstrom Exp $
c
      integer i,locinf,nmh,np1
      double precision dummy(1)
c
c
      np1 = nobs + 1
      nmh = nobs - nnull
c			form the response vector
      call dcopy (nnull,w1,1,pred,1)
      call dcopy (npsing,g,1,work(np1),1)
      do 10 i = 1,npsing
         work(nobs+i) = work(nobs+i)*svals(i)
   10 continue
      call dgemv('N',nmh,npsing,1.0d0,u,ldu,work(np1),1,0.0d0,
     *  pred(nnull+1),1)
      call dqrsl (fg,ldfg,nobs,nnull,fgaux,pred,pred,dummy,dummy,dummy,
     * dummy,10000,locinf)
      return
      end
      subroutine dpmse(fg,ldfg,nuobs,nobs,nnull,fgaux,svals,npsing,u,
     * ldu,w1,z,ntbl,adiag,tbl,ldtbl,auxtbl,work)
      integer ldfg,nuobs,nobs,nnull,npsing,ldu,ntbl,ldtbl
      double precision fg(ldfg,nnull),fgaux(nnull),svals(npsing),
     * u(ldu,npsing),w1(nnull),z(npsing),adiag(nuobs),tbl(ldtbl,3),
     * auxtbl(3,3),work(npsing)
c
c Purpose: determine the predictive mean squared error for each lambda 
c	value in tbl.
c
c On Entry:
c   fg(ldfg,nnull)	information on the Householder transformations 
c			that define f and g
c   ldfg		leading dimension of fg as declared
c			in the calling program
c   nuobs		number of rows in f
c   nnull		number of columns in g
c   fgaux(nnull)	auxiliary information on the fg Householder 
c			transformations 
c   svals(npsing)	singular values 
c   npsing		number of singular values
c   u(ldu,npsing)	left singular vectors corresponding to svals
c   ldu			leading dimension of u as declared in the 
c			calling program
c   w1(nnull)		leading part of rotated response vector
c   z(npsing)		u'w2
c   ntbl		number of rows in tbl
c   adiag(nuobs)	"true" y values 
c   tbl(ldtbl,3)	column	contains
c			  1 	grid of log10(nobs*lambda) 
c   ldtbl		leading dimension of tbl as declared in the 
c			calling program	
c   auxtbl(3,3)		auxiliary table
c			auxtbl(1,1) contains log10(nobs*lamhat) where
c			lamhat is the gcv estimate of lambda
c
c On Exit:
c   tbl(ldtbl,3)	column	contains
c			  1 	grid of log10(nobs*lambda) 
c			  3     R(lambda) 
c   auxtbl(3,3)		auxiliary table
c			3rd column contains:
c			    [R(lamhat) , R(0), R(infinity)]'
c
c Work Arrays:
c   work(npsing)	double precision work vector
c
c Subprograms Called Directly:
c      Linpack - dqrsl 
c      Blas    - ddot dgemv
c
c $Header: dpmse.f,v 2.100.1.1 86/10/07 12:51:29 lindstrom Exp $
c
      integer i,nmh,k,locinf
      double precision dummy(1),nlam,wrk1,addtru,wrk
      double precision ddot
c
      nmh = nuobs - nnull
      addtru = 0.0d0
      call dqrsl(fg,ldfg,nuobs,nnull,fgaux,adiag,dummy,adiag,dummy, 
     * dummy,dummy,01000,locinf)
c			the first nnull positions of adiag now contain 
c			w1 true the last nuobs-nnull positions contain 
c			w2 true
      do 10 i = 1,nnull
         addtru=addtru + (w1(i)-adiag(i))**2
   10 continue
      addtru = addtru + ddot(nmh,adiag(nnull+1),1,adiag(nnull+1),1)
      call dgemv('T',nmh,npsing,1.0d0,u,ldu,adiag(nnull+1),1,0.0d0,
     *  work,1)
      addtru = addtru - ddot(npsing,work,1,work,1)
c			addtru contains ||w1 - (w1 true)||**2 +
c			||w2 true||**2 - ||z true||**2
c			work contains z true
c
c			compute predictive mse for each lambda in tbl
      do 30 k = 1,ntbl 
         nlam = 10**tbl(k,1)
         wrk=0.0d0
         do 20 i=1,npsing 
            wrk1=(svals(i)**2)/(svals(i)**2+nlam)
            wrk = wrk + (work(i)-z(i)*wrk1)**2
   20    continue
         tbl(k,3)=(addtru+wrk)/nobs
   30 continue
c			add pred. mse for lambda hat to auxtbl
      wrk=0.0d0
      nlam=10**auxtbl(1,1)
      do 40 i=1,npsing 
         wrk1=(svals(i)**2)/(svals(i)**2+nlam)
         wrk = wrk + (work(i)-z(i)*wrk1)**2
   40 continue
      auxtbl(1,3)=(addtru+wrk)/nobs
c			add pmse for lambda = 0
      wrk=0.0d0
      do 50 i=1,npsing 
         wrk = wrk + (work(i)-z(i))**2
   50 continue
      auxtbl(2,3)=(addtru+wrk)/nobs
c			add pmse for lambda = infinity
      auxtbl(3,3)=(addtru+ddot(npsing,work,1,work,1))/nobs
      return
      end
      subroutine dpred(pdes,ldpdes,npred,dim,m,desb,lddesb,ndesb,ps,
     * ldps,ncov1,ncov2,coef,npar,pred,work,lwa,iwork,info) 
      integer ldpdes,npred,dim,m,lddesb,ndesb,ldps,ncov1,ncov2,npar,lwa,
     * iwork(dim),info
      double precision pdes(ldpdes,dim),desb(lddesb,dim),ps(ldps,*),
     * coef(npar),pred(npred),work(lwa)
c
c   Purpose: determine predicted values at the locations in pdes and ps
c
c  On Entry:
c   pdes(ldpdes,dim) 	prediction design for splined variables
c   ldpdes		leading dimension of pdes as declared in the 
c			calling	program 
c   npred		number of rows in pdes
c   desb(lddesb,dim) 	locations for the basis functions
c			(returned from dtpss and dptpss in the 
c			variable des)
c   lddesb		leading dimension of desb as declared in the
c			calling	program 
c   ndesb		number of rows in desb
c   dim			number of columns in desb
c   m			order of the derivatives in the penalty
c   ps(ldps,ncov1+ncov2) prediction covariates corresponding to pdes
c   ldps		leading dimension of ps as declared in the
c			calling	program 
c   ncov1		number of covariates which duplicate the 
c			replication structure of pdes
c   ncov2		number of covariates which do not duplicate the 
c			replication structure of pdes
c   coef(npar)		coefficient estimates  [delta':xi']'
c   npar		ndesb + (m+dim-1 choose dim) + ncov1 + ncov2
c
c On Exit:
c   pred(npred)		predicted values
c   info		error indicator
c			  0 : successful completion
c			  1 : dimension error
c			  2 : error in npar,ncov1,ncov2,m or dim
c			  3 : lwa too small
c			  4 : error in dmaket
c			
c
c Working Storage:
c   work(lwa)		double precision work vector
c   lwa			length of work vector
c			must be at least npred*(nct+ndesb)
c			where nct = (m+dim-1 choose dim)
c   iwork(dim)		integer work vector
c
c Subprograms Called Directly:
c    Gcvpack - dmaket dmakek
c    Blas    - dgemv
c
c Subprograms Called Indirectly:
c    Blas    - dcopy
c    Other   - fact mkpoly
c
c $Header: dpred.f,v 2.100.1.1 86/10/07 12:51:40 lindstrom Exp $
c
      double precision dummy
      integer nct,p1,p1p1,npoly
      integer mkpoly 
c
      nct = mkpoly(m,dim) 
      if ((ndesb .le. 0) .or. (nct .le. 0) .or. (m .le. 0) .or. 
     * (dim .le. 0) .or. 2*m - dim .le. 0) then
	 info = 1
	 return
      endif
      if (npar .ne. ndesb + nct + ncov1 + ncov2) then
         info = 2
         return
      endif
      if (lwa .lt. npred*(nct+ndesb)) then
	 info = 3
	 return
      endif
c			first npred*nct positions of work contain t
      p1 = npred*nct
c			next npred*ndesb positions of work contain k
      p1p1 = p1 + 1
c changed desb from dummy 7/28/03 by Jeff Sklar
      call dmaket(m,npred,dim,pdes,ldpdes,desb,1,0,npoly,work(1),npred,
     * iwork,info)
      if (info .ne. 0) then
         info = 4
         return
      endif
      call dmakek(m,npred,dim,pdes,ldpdes,ndesb,desb,lddesb,work(p1p1),
     * npred)
c			compute predicted values
      call dgemv('N',npred,nct,1.0d0,work,npred,coef,1,0.0d0,
     * pred,1)
      call dgemv('N',npred,ncov1+ncov2,1.0d0,ps,ldps,coef(nct+1),1,
     * 1.0d0,pred,1)
      call dgemv('N',npred,ndesb,1.0d0,work(p1+1),npred,
     * coef(nct+ncov1+ncov2+1),1,1.0d0,pred,1)
      return
      end
      subroutine dprmut (x,npar,jpvt,job)
      integer npar,jpvt(npar),job
      double precision x(npar)
c
c Purpose: permute the elements of the array x according to the index 
c	vector jpvt (either forward or backward permutation).
c
c On Entry:
c   x(npar)		array to be permuted
c   npar		size of x (and jpvt)
c   jpvt		indices of the permutation
c   job			indicator of forward or backward permutation
c			if job = 0 forward permutation  
c				x(jpvt(i)) moved to x(i)
c			if job is nonzero backward permutation 
c				x(i) moved to x(jpvt(i))
c On Exit:
c   x(npar)		array with permuted entries
c
c   Written:	Yin Ling	U. of Maryland, August,1978
c
c $Header: dprmut.f,v 2.100.1.1 86/10/07 12:51:58 lindstrom Exp $
c
      integer i,j,k
      double precision t
c
      if (npar .le. 1) then
         return
      endif
      do 10 j = 1,npar
         jpvt(j) = -jpvt(j)
   10 continue
      if (job .eq. 0) then
c		forward permutation
         do 30 i = 1,npar 
            if (jpvt(i) .gt. 0) then
               goto 30
            endif
            j = i
            jpvt(j) = -jpvt(j)
            k = jpvt(j)
c           while
   20       if (jpvt(k) .lt. 0) then
               t = x(j)
               x(j) = x(k)
               x(k) = t
               jpvt(k) = -jpvt(k)
               j = k
               k = jpvt(k)
               goto 20
c           endwhile
            endif
   30    continue
      endif
      if (job .ne. 0 ) then
c			backward permutation
         do 50 i = 1,npar 
            if (jpvt(i) .gt. 0) then
               goto 50
            endif
            jpvt(i) = -jpvt(i)
            j = jpvt(i)
c           while
   40       if (j .ne. i) then
               t = x(i)
               x(i) = x(j)
               x(j) = t
               jpvt(j) = -jpvt(j)
               j = jpvt(j)
               goto 40
c           endwhile
            endif
   50    continue
      endif
      return
      end
      subroutine dptpss(des,lddes,nobs,dim,m,s,lds,ncov1,ncov2,y,ntbl,
     * adiag,lamlim,dout,iout,coef,svals,tbl,ldtbl,auxtbl,work,
     * lwa,iwork,liwa,job,info)
      integer lddes,nobs,dim,m,lds,ncov1,ncov2,ntbl,iout(4),ldtbl,lwa,
     * liwa,iwork(liwa),job,info
      double precision des(lddes,dim),s(lds,*),y(nobs),adiag(nobs),
     * lamlim(2),dout(4),coef(*),svals(*),tbl(ldtbl,3),
     * auxtbl(3,3),work(lwa)
c
c Purpose: determine the generalized cross validation estimate of the 
c	smoothing parameter and fit model parameters for a partial thin
c	plate spline model.
c
c On Entry:
c   des(lddes,dim) 	design for the variables to be splined
c   lddes		leading dimension of des as declared in the
c   			calling program 
c   nobs		number of observations
c   dim			number of columns in des
c   m			order of the derivatives in the penalty
c   s(lds,ncov1+ncov2)	design for the covariates
c			first ncov1 columns contain covariates which
c			  duplicate the replication structure of des
c			next ncov2 columns contain covariates which
c			  do not duplicate the replication structure of
c			  des
c   lds			leading dimension of s as declared in the
c   			calling program 
c   ncov1		number of covariates which duplicate the 
c			replication structure of des
c   ncov2		number of covariates which do not duplicate the 
c			replication structure of des
c   y(nobs)		response vector
c   ntbl		number of evenly spaced values for 
c			log10(nobs*lambda) to be used in the initial 
c			grid search for lambda hat
c			if ntbl = 0 only a golden ratio search will be 
c			done and tbl is not referenced, if ntbl > 0
c			there will be ntbl rows returned in tbl
c   adiag(nobs)	 	"true" y values on entry if predictive mse is 
c			requested
c   lamlim(2)		limits on lambda hat search (in log10(nobs*
c			lambda)	scale) if user input limits are 
c			requested. if lamlim(1) = lamlim(2) then lamhat 
c			is set to (10**lamlim(1))/nobs
c   ldtbl		leading dimension of tbl as declared in the 
c			calling program	
c   job			integer with decimal expansion abc
c			if a is nonzero then predictive mse is computed
c			   using adiag as true y
c			if b is nonzero then user input limits on search
c			   for lambda hat are used
c			if c is nonzero then adiag will be calculated
c On Exit:
c   des(lddes,dim)	unique rows of des
c   y(nobs)		predicted values
c   adiag(nobs)		diagonal elements of the hat matrix if requested
c   lamlim(2)		limits on lambda hat search 
c			(in log10(nobs*lambda) scale)
c   dout(4)		contains:
c  			1  lamhat   generalized cross validation 
c				    estimate of the smoothing parameter
c			2  penlty   smoothing penalty
c			3  rss      residual sum of squares
c			4  tr(I-A)  trace of I - A
c   iout(4)		contains:
c			1  npsing   number of positive singular values
c				    if info indicates nonzero info 
c				    from dsvdc then npsing contains
c				    info as it was returned from dsvdc
c			2  npar	    number of parameters 
c				    (npar = nuobs + nnull)
c			3  nnull    size of the null space of sigma
c				    (m+dim-1 choose dim)+ncov1+ncov2
c			4  nuobs    number of unique rows in des
c   coef(npar)		coefficient estimates [beta':alpha':delta']'
c			coef must have a dimension of at least 
c			nuobs+nnull
c   svals(npsing)	singular values, svals must have a dimension,
c			of at least nuobs-nnull.
c			if info indicates nonzero info in dsvdc then 
c			svals is as returned from dsvdc.
c   tbl(ldtbl,3)	column	contains
c			  1 	grid of log10(nobs*lambda) 
c			  2  	V(lambda)
c			  3     R(lambda) if requested
c   auxtbl(3,3)		auxiliary table
c			1st row contains:
c			    log10(nobs*lamhat), V(lamhat) and  
c			    R(lamhat) if requested
c			    where lamhat is the gcv estimate of lambda
c			2nd row contains:
c			    0, V(0) and  R(0) if requested
c			3rd row contains:
c			    0, V(infinity) and R(infinity) if requested
c   info		error indicator
c			   0 : successful completion
c			  -1 : log10(nobs*lamhat) <= lamlim(1) 
c			       (not fatal)
c			  -2 : log10(nobs*lamhat) >= lamlim(2) 
c			       (not fatal)
c			   1 : dimension error 	
c			   2 : error in dreps, the first ncov1 columns
c			       of s do not duplicate the replication 
c			       structure of des
c			   3 : lwa (length of work) is too small
c			   4 : liwa (length of iwork) is too small
c			   5 : error in dmaket
c			   6 : sigma is rank deficient
c			  1000< info : 1000 + nonzero info returned from
c				       dsnsm
c
c Working Storage:
c   work(lwa)		double precision work vector
c   lwa			length of work as declared in the calling 
c			program 
c			must be at least lwa1 + lwa2 where 
c			  lwa1 = (nnull-ncov2)*(nobs+nuobs+1)
c				 +npar*(nobs+npar)
c			  lwa2 = (npar-nnull)*(npar-2*nnull+2+nobs)
c			 	 +npar+nobs
c			
c			
c   iwork(liwa)		integer work vector
c   liwa		length of the iwork as declared in the calling 
c			program
c			must be at least 3*nobs - (nnull - ncov2)
c
c Subprograms Called Directly:
c 	Gcvpack - dreps dmaket duni dmakek dctsx dsnsm
c	Linpack - dqrdc dqrsl 
c	Blas    - dcopy
c	Other 	- dprmut dset prmut mkpoly
c
c Subprograms Called Indirectly:
c	Gcvpack - dcrtz ddcom dgcv dsgdc dtsvdc drsap ddiag
c		  dvlop dvlop dpmse dcfcr dpdcr dvmin dvl dzdc
c	Linpack - dchdc dqrdc dqrsl dtrsl dsvdc dtrco
c	Blas    - dcopy ddot dgemv dswap
c	Other 	- dcpmut dprmut dset dftkf fact mkpoly
c
c $Header: dptpss.f,v 2.100.1.2 86/11/21 12:22:18 lindstrom Exp $
c
      double precision dummy(1),dout2(5)
      integer i,ncts1,jadiag,p1,p1p1,p2,p2p1,p3,p3p1,p4,p4p1,p5,p5p1,
     * ip1,ip1p1,ip2,ip2p1,nuobs,lwa2,liwa2,npoly,nnull,snnpar,q1,
     * wsize,lwa1,sinfo,ldx,iout2(4)
      integer mkpoly
c
c
      info = 0
      sinfo = 0
      ncts1 = mkpoly(m,dim) + ncov1
      nnull = ncts1 + ncov2
      iout(3) = nnull
      jadiag = mod(job,10)
c			check dimensions
      if((nobs .le. 0) .or. (m .le. 0) .or. (dim .le. 0) .or.
     * (ntbl .lt. 0) .or. (ntbl .gt. ldtbl) .or. 
     * (2*m - dim .le. 0)) then
         info = 1
         return
      endif
c			set up pointers for iwork
c  			first nobs positions of iwork contain order
      ip1 = nobs
c			next nobs positions of iwork contain xrep
      ip1p1 = ip1 + 1
      ip2 = ip1 + nobs
c			rest of iwork is an integer work vector
      ip2p1 = ip2 + 1
      liwa2 = liwa - 2*nobs
      call dreps(des,lddes,nobs,dim,s,lds,ncov1,ncov2,dummy,iwork,nuobs,
     * iwork(ip1p1),0,info)
      iout(2) = nuobs+ncts1+ncov2
      iout(4) = nuobs
      if (info .ne. 0) then
         info = 2
         return
      endif
c			undo permutation of des, s and xrep
      do 10 i = 1,dim
          call dprmut(des(1,i),nobs,iwork,1)
   10 continue
      do 20 i = 1,ncov1+ncov2
          call dprmut(s(1,i),nobs,iwork,1)
   20 continue
      call prmut(iwork(ip1p1),nobs,iwork,1)
c
c			check size of work vectors
      snnpar = nuobs + ncov2
      lwa1 = ncts1*(nobs+nuobs+1)+snnpar*(nobs+snnpar)
      lwa2 = (snnpar-nnull)*(snnpar-2*nnull+2+nobs)+snnpar+nobs
      wsize = lwa1+lwa2
      if (lwa .lt. wsize) then
         info = 3
         return
      endif
      wsize = 3*nobs - ncts1
      if (liwa .lt. wsize) then
         info = 4
         return
      endif
c			set up pointers for work
c
c	     name     runs from  to		
c	     ----     ---------  --
c	     [t:s1]      1	 p1	  p1 = nobs*ncts1	
c	     [tu:s1u]    p1p1    p2       p2 = p1 + nuobs*ncts1
c	      x	         p2p1    p3       p3 = p2 + max(nobs,snnpar)
c								*snnpar
c	      sigma      p3p1    p4	  p4 = p3 + snnpar**2
c	      fgaux      p4p1    p5	  p5 = p4 + ncts1
c	      working    p5p1    p5+(snnpar-nnull)*(snnpar-2+nnull+
c					2+nobs)+snnpar+nobs
c
c		after the call to duni ku runs from q1 to p4
c
      p1 = nobs*ncts1
      p1p1 = p1 + 1
      p2 = p1 + nuobs*ncts1
      p2p1 = p2 + 1
      ldx=max(nobs,snnpar)
      p3 = p2 + ldx*snnpar
      p3p1 = p3 + 1
      q1 = p3 + ncov2*snnpar+ncov2+1
      p4 = p3 + snnpar**2
      p4p1 = p4 + 1
      p5 = p4 + ncts1
      p5p1 = p5 + 1
      lwa2 = lwa - ( ncts1*(nobs+nuobs+1) + snnpar*(nobs+snnpar))
c			make [t:s1] and k
      if (m .ne. 1) then
          call dmaket(m,nobs,dim,des,lddes,s,lds,ncov1,npoly,work,nobs,
     *     iwork(ip2p1),info)
          if (info .ne. 0) then
             info = 5
             return
          endif
c			put unique rows of des into des
          call duni(des,lddes,nobs,dim,iwork(ip1p1),des,lddes)
c			put k into x
          call dmakek(m,nobs,dim,work(1+nobs),nobs,nuobs,des,lddes,
     *     work(p2p1),ldx)
      else
c			put unique rows of des in 1st nuobs 
c			positions of work
          call duni(des,lddes,nobs,dim,iwork(ip1p1),work,1)
c			put k into x
          call dmakek(m,nobs,dim,work,nobs,nuobs,des,lddes,
     *     work(p2p1),ldx)
c			copy unique rows of des from work to des 
          call dcopy(nuobs,work,1,des,1)
c			set t = 1
	  call dset(nobs,1.0d0,work,1)
      endif
c			put unique rows of k into (ncov2+1,ncov2+1)
c			position of sigma
      call duni(work(p2p1),ldx,nobs,nuobs,iwork(ip1p1),work(q1),snnpar)
c			compute unique rows of [t:s1] 
      call duni(work(1),nobs,nobs,ncts1,iwork(ip1p1),work(p1p1),nuobs)
c			make x and sigma
      call dctsx(work(1),nobs,ncts1,s(1,ncov1+1),lds,ncov2,nobs,
     * work(p1p1),nuobs, nuobs,work(p3p1),snnpar,work(p2p1),ldx,snnpar,
     * work(p4p1),work(p5p1))
c
      call dsnsm(work(p2p1),ldx,y,work(p3p1),snnpar,nobs,snnpar,nnull,
     * adiag,1.0d0,lamlim,ntbl,dout2,iout2,coef,svals,tbl,ldtbl,auxtbl,
     * iwork(ip2p1),liwa2,work(p5p1),lwa2,job,info)
      call dcopy(4,dout2,1,dout,1)
      iout(3) = nnull
      iout(1) = iout2(1)
      if (info .eq. -3) then
	  info = 6
	  return
      endif
      if (info .gt. 0) then
	 info = info + 1000
	 return
      endif
      if (info .lt. 0) sinfo = info
c			form psi = f2 zeta
      call dset(ncts1,0.0d0,work(p5p1),1)
      call dcopy(nuobs-ncts1,coef(ncts1+ncov2+1),1,work(p5+ncts1+1),1)
      call dqrsl(work(p1p1),nuobs,nuobs,ncts1,work(p4p1),work(p5p1),
     * coef(ncts1+ncov2+1),dummy,dummy,dummy,dummy,10000,info)
      if (sinfo .lt. 0) info = sinfo
      return
      end
      subroutine dreps(des,lddes,nobs,dim,s,lds,ncov1,ncov2,c1,order,
     * nuobs,xrep,job,info)
      integer lddes,nobs,dim,lds,ncov1,ncov2,order(nobs),nuobs,
     * xrep(nobs),info,job
      double precision des(lddes,dim),s(lds,*),c1(*)
c
c Purpose: sort des and s, compute degrees of freedom for replication
c	and c1.
c
c On Entry:
c   des(lddes,dim) 	design for the variables to be splined
c			should be entered in lexicographical order 
c			(smallest to largest) if possible for efficient
c			computing
c   lddes		leading dimension of des as declared in the
c   			calling program 
c   nobs		number of observations
c   dim			number of columns in des
c   s(lds,ncov1+ncov2) 	design for the covariates
c   lds			leading dimension of s as declared in the
c   			calling program 
c   ncov1		number of covariates which duplicate the 
c			replication structure of des
c   ncov2		number of covariates which do not duplicate
c			replication structure of des
c   job			if job is nonzero then c1 is computed
c   			if job = 0 then c1 is not referenced
c
c On Exit:
c   des(lddes,dim) 	des sorted lexicographically
c   s(lds,ncov1+ncov2) 	s, sorted to correspond to des
c   c1(nuobs)		if job is nonzero then c1(i) the square root of
c			the number of replicates of the ith sorted 
c			design point 
c   order(nobs)		order of the sorted des
c   nuobs		number of unique rows in des
c   xrep(nobs)		xrep(i) = 1 if the ith sorted design point is a
c			replicate, 0 if not
c   info		error indicator
c			   0 : successful completion
c			   1 : ncov1 is incorrect
c
c
c $Header: dreps.f,v 2.100.1.1 86/10/07 12:55:24 lindstrom Exp $
c
      integer sw,oldsw,itemp,i,j,k,cont,dfrep
      double precision temp,diff,one,machpr,denom,wmin,wmax
c
      info = 0
      one = 1.0d0
      machpr = 1.0d0
   10 machpr = machpr/2.0d0
      if (one .lt. 1.0d0 + machpr) goto 10
      machpr = machpr*2.0d0
c
      do 20 i = 1,nobs
          order(i) = i
 	  xrep(i) = 0
   20 continue
      if (job .ne. 0) call dset(nobs,0.0d0,c1,1)
c			sort des and s 
      sw = nobs - 1
   30 if (sw .le. 0) goto 90 
          oldsw = sw
          sw = 0
          do 80 i = 1,oldsw 
	      cont = 1
	      k = 1
   40         if (cont .eq. 0) goto 70
                  if (k .le. dim) then
	              diff = des(i,k) - des(i+1,k)
                  else
	              diff = s(i,k-dim) - s(i+1,k-dim)
                  endif
		  if (diff .lt. 0.0d0) then
		      if (k .gt. dim) info = 1
		      cont = 0
		  else if (diff .gt. 0.0d0) then
		      if (k .gt. dim) info = 1
c		      switch the order of i and i+1
		      itemp = order(i)
		      order(i)=order(i+1)
		      order(i+1) = itemp
		      itemp = xrep(i)
		      xrep(i)= xrep(i+1)
		      xrep(i+1)= itemp
		      do 50 j = 1,dim
			  temp = des(i,j)
			  des(i,j) = des(i+1,j)
			  des(i+1,j) = temp
   50 		      continue
  		      do 60 j = 1,ncov1+ncov2
			  temp = s(i,j)
			  s(i,j) = s(i+1,j)
			  s(i+1,j) = temp
   60 		      continue
		      sw = i
		      cont = 0
                  else if (k .eq. dim + ncov1) then
		      xrep(i + 1) = 1
		      cont = 0
                  else 
	    	      k = k + 1
                  endif
              goto 40
   70         continue
   80     continue
      goto 30
   90 continue
c			compute range of design
      denom=0.0d0
      do 120 j=1,dim
          wmin = des(1,j)
	  wmax = des(1,j)
          do 110 i=1,nobs
	      if (des(i,j) .lt. wmin) wmin = des(i,j)
	      if (des(i,j) .gt. wmax) wmax = des(i,j)
  110     continue
	  denom = denom + (wmax-wmin)**2
  120 continue
	    
c			check for design points too close together
      do 140 i=1,nobs-1
	  if (xrep(i+1) .eq. 0) then
	     diff = 0.0d0
	     do 130 j=1,dim
	        diff = diff + (des(i,j)-des(i+1,j))**2
  130	     continue
	     if (abs(diff)/denom .lt. 100*machpr) xrep(i+1)=1
	  endif 
  140 continue
c			compute dfrep and c1
      dfrep = 0
      j = 0
       do 150 i = 1,nobs
	   j = j + 1 - xrep(i)
	   if (job .ne. 0) c1(j) = xrep(i)*c1(j) + 1.0d0
	   dfrep = dfrep + xrep(i)
  150 continue
      nuobs = nobs - dfrep
      if (job .eq. 0 ) return
      do 160 i = 1,nuobs
	  c1(i) = sqrt(c1(i))
  160 continue
      return 
      end
      subroutine drsap(fg,ldfg,nobs,nnull,fgaux,u,ldu,nmh,npsing,z,
     * ssqw2,addend,work)
      integer ldfg,nobs,nnull,ldu,nmh,npsing
      double precision fg(ldfg,nnull),fgaux(nnull),u(ldu,npsing),
     * z(nobs),ssqw2,addend,work(npsing)
c
c Purpose: apply Householder transformations to a response vector and
c	collect its inner product with u and the addend which are used 
c	to define the generalized cross validation function with a
c	semi-norm.
c
c On Entry:
c   fg(ldfg,nnull)	information on the Householder transformations 
c			that define f and g
c   ldfg		leading dimension of fg as declared in the 
c			calling program
c   nobs		number of rows in fg
c   nnull		number of columns in fg
c   fgaux(nnull)	auxiliary information on the fg	Householder 
c			transformations
c   u(ldu,npsing)		left singular vectors 
c   ldu			leading dimension of u as declared in the
c			calling program
c   nmh	     		number of rows in u. nmh = nobs - nnull
c   npsing		number of columns in u (maximum of npar - nnull)
c   z(nobs)		response vector
c
c On Exit:
c   z(nobs)		the first nnull positions contain w1 and the 
c			next npsing positions contain u'w2
c   ssqw2		the squared length of w2
c   addend		the squared length of z minus the squared length
c			of u'w2
c
c Work Arrays:
c   work(npsing) 	double precision work vector
c
c Subprograms Called Directly:
c      Linpack - dqrsl 
c      Blas    - ddot dcopy dgemv
c
c $Header: drsap.f,v 2.100.1.1 86/10/07 12:55:34 lindstrom Exp $
c
      integer locinf,hp1
      double precision dummy(1)
      double precision ddot
c			apply Householder transformations
c			which define f
      call dqrsl (fg,ldfg,nobs,nnull,fgaux,z,dummy,z,dummy,dummy,dummy,
     * 01000,locinf)
c			w1 in first nnull positions of z,w2 in
c			last nmh
      hp1 = nnull + 1
      addend = ddot(nmh,z(hp1),1,z(hp1),1)
      ssqw2=addend
      call dgemv('T',nmh,npsing,1.0d0,u,ldu,z(hp1),1,0.0d0,work,1)
c			u'w2 in positions nnull+1 to
c			nnull+npsing of z
      call dcopy (npsing,work,1,z(hp1),1)
      addend = addend - ddot(npsing,z(hp1),1,z(hp1),1)
      return
      end
      subroutine dsetup(des,lddes,su1,ldsu1,dim,m,ncov1,nuobs,c1,
     * tusu1,ldtu,ncts1,fgaux,ku,ldku,work,iwork,info)
      integer lddes,ldsu1,dim,m,ncov1,nuobs,ldtu,ncts1,ldku,
     * iwork(dim),info
      double precision des(lddes,dim),su1(ldsu1,*),c1(nuobs),
     * tusu1(ldtu,ncts1),fgaux(ncts1),ku(ldku,nuobs),work(ncts1)
c
c Purpose: set up [tu:su1] as f g and ku as f'c1 ku c1'f.
c
c On Entry:
c   des(lddes,dim)  	variables to be splined (unique rows)
c   lddes		leading dimension of des as declared in the 
c			calling	program 
c   su1(ldsu1,ncov1)	covariates (unique rows)
c   ldsu1		leading dimension of su1 as declared in the
c			calling	program 
c   dim 		dimension of the variables to be splined
c   m			order of the derivatives in the penalty
c   ncov1		number of covariates
c   nuobs		number of unique rows in des
c   c1(nuobs)		c1(i) contains the square root of the number of
c			replicates of the ith sorted design point
c   ldtu		leading dimension of tusu1 as declared in the 
c			calling	program 
c   ldku		leading dimension of ku as declared in the 
c			calling program 
c
c On Exit:
c   tusu1(ldtu,ncts1)	the qr decomposition of [tu:su1]
c   ncts1		number of columns in [tu:su1] = npoly + ncov1
c   fgaux(ncts1)	the auxiliary info on the qr decomposition of 
c			[tu:su1]
c   ku(p,p)  		f'ku f
c   info		error indicator
c			   0 : successful completion
c			   1 : error in dmaket
c
c Work Arrays:
c   work(ncts1)		double precision work vector
c   iwork(dim)		integer work vector
c
c Subroutines called directly
c	Gcvpack - dmaket dmakek
c	Linpack - dqrdc
c	Other   - dftkf mkpoly
c
c Subroutines called indirectly
c	Blas    - dcopy
c	Gcvpack - dqrsl
c	Other   - fact mkpoly
c
c $Header: dsetup.f,v 2.100.1.1 86/10/07 12:55:48 lindstrom Exp $
c
      integer npoly,i,j
      integer mkpoly,idummy(1)
	  double precision dummy(1)
c
      info = 0
      npoly=mkpoly(m,dim)
      ncts1=npoly+ncov1
c			make [tu:su1] and ku
      call dmaket(m,nuobs,dim,des,lddes,su1,ldsu1,ncov1,npoly,tusu1,
     * ldtu,iwork,info)
      if (info .ne. 0) then
	 return
      endif
      call dmakek(m,nuobs,dim,des,lddes,nuobs,des,lddes,ku,ldku)
      if (c1(1) .ne. 0) then
         do 30 i = 1,nuobs 
            do 10 j = 1,npoly+ncov1
               tusu1(i,j) = tusu1(i,j) * c1(i)
   10       continue
            do 20 j = 1,nuobs
               ku(i,j) = ku(i,j) * c1(i) * c1(j)
   20       continue
   30    continue
      endif
c			decompose [tu:su1] into fg
      call dqrdc(tusu1,ldtu,nuobs,ncts1,fgaux,idummy,dummy,0)
c      			calculate f'ku f
      call dftkf(tusu1,ldtu,nuobs,ncts1,fgaux,ku,ldku,work)
      return
      end
      subroutine dsgdc (sigma,ldsigm,npar,nnull,qraux,sgpvt,info)
      integer ldsigm,npar,nnull,sgpvt(npar),info
      double precision sigma(ldsigm,npar),qraux(npar)
c
c Purpose: decompose the semi-norm smoothing matrix into a QR 
c	decomposition of the transpose of the Cholesky factor.
c
c On Entry:
c   sigma(ldsigm,npar)	symmetric matrix which defines the semi-norm
c   ldsigm		leading dimension of sigma as declared in the
c			calling program
c   npar		number of parameters
c   nnull		dimension of the null space of sigma
c
c On Exit:
c   sigma(ldsigm,npar)	overwritten with the QR decomposition
c			of the Cholesky factor of sigma
c   nnull		if input nnull is too small it is replaced by
c			larger value such that sigma has rank npar-nnull
c   qraux(npar)		auxiliary information for the QR decomposition
c   sgpvt(npar) 	permuted indices from the Cholesky decomposition
c			with pivots of sigma
c   info		error indicator
c			   0 : successful completion
c			  -3 : nnull is too small (not fatal)
c			   1 : nnull is too large
c
c Subprograms Called Directly:
c	Linpack - dchdc dqrdc  
c	Blas    - dcopy
c	Other   - dset
c
c $Header: dsgdc.f,v 2.100.1.1 86/10/07 12:56:05 lindstrom Exp $
c
      integer locinf,i,j,nsm,idummy(1)
      double precision dummy(1),machpr,one
c
c
      info = 0
      one = 1.0d0
      machpr = 1.0d0
   10 machpr = machpr/2.0d0
      if (one .lt. 1.0d0 + machpr) goto 10
      machpr = machpr*2.0d0
c
      locinf = 0
c			Cholesky decomposition of sigma
      do 20 j = 1,npar
         sgpvt(j) = 0
   20 continue
      call dchdc (sigma,ldsigm,npar,qraux,sgpvt,1,locinf)
      do 30 i=1,locinf 
         if (((sigma(i,i)/sigma(1,1))**2) .gt. machpr) nsm = i
   30 continue
      if (nsm .gt. npar - nnull) then
         info = 1
         return
      endif
      if (nsm .lt. npar - nnull) then
         nnull = npar - nsm
         info = -3
      endif
c			copy transpose of Cholesky factor to sigma
      do 40 i = 1,nsm 
         call dcopy (npar,sigma(i,1),ldsigm,sigma(1,i),1)
         j = i - 1
         call dset (j,0.0d0,sigma(1,i),1)
   40 continue
c			QR decomposition of Cholesky transpose
      call dqrdc (sigma,ldsigm,npar,nsm,qraux,idummy,dummy,0)
      return
      end
      subroutine dsgdc1(f2kf2,ldfkf,p,svals,iwork,work,lwa,info)
      integer ldfkf,p,iwork(p),lwa,info
      double precision f2kf2(ldfkf,p),svals(*),work(lwa)
c 
c Purpose: form the singular value decomposition of the Cholesky factor 
c 	of f2'k f2.
c
c On Entry:
c   f2kf2(ldfkf,p)	f2'k f2
c   ldfkf		leading dimension of f2'k f2 as declared in the 
c			calling	program 
c   p			number of rows and columns in f2'k f2
c
c On Exit:
c   f2kf2(p,p)  	overwritten with singular value decomposition 
c 			of Cholesky factor of f2'k f2
c   svals(p) 		the singular values of the Cholesky factor of 
c			f2'k f2 if info = 0.
c			if info = 3 then svals is as it was returned 
c			from dsvdc.
c   info 	   	error indicator
c			  0 : successful completion
c			  1 : lwa too small
c			  2 : f2'k f2 is not of full rank 
c			  3 : error in dsvdc
c   p			if info = 3 p contains info as it was returned 
c			from dsvdc (otherwise unchanged)
c
c Work Arrays:
c   work(lwa)		double precision work vector
c   lwa			length of work as declared in the calling
c			program
c			must be at least 2*p
c   iwork(p)		integer work vector
c
c Subprograms Called:
c	Linpack - dchdc dsvdc
c	Blas    - dcopy
c	Other   - dset dprmut
c
c $Header: dsgdc1.f,v 2.100.1.1 86/10/07 12:56:12 lindstrom Exp $
c
      integer i,j,pp1,locinf,k
      double precision dummy(1,1),one,machpr
c
      info = 0
      if (lwa .lt. 2*p) then
	  info = 1
	  return
      endif
      pp1 = p + 1
      call dset(p,0.0d0,svals,1)
c
      one = 1.0d0
      machpr = 1.0d0
   10 machpr = machpr/2.0d0
      if (one .lt. 1.0d0 + machpr) goto 10
      machpr = machpr*2.0d0
c			Cholesky decomposition of f2'k f2
      do 20 j = 1,p
         iwork(j) = 0
   20 continue
      call dchdc (f2kf2,ldfkf,p,work,iwork,1,locinf)
      do 30 i=1,locinf 
         if ((f2kf2(i,i)/f2kf2(1,1))**2 .gt. machpr) k = i
   30 continue
      if (k .lt. p) then
	 info = 2
	 return
      endif
c    		copy f2kf2' into f2kf2 
c    		svd of f2 k f2' = udv' return only u
      do 40 j = 1,p 
         call dcopy(1+p-j,f2kf2(j,j),ldfkf,f2kf2(j,j),1)
         call dset(j-1,0.0d0,f2kf2(1,j),1)
   40 continue
      call dsvdc(f2kf2,ldfkf,p,p,svals,work,f2kf2,ldfkf,dummy,1,
     * work(pp1),20,info)
      if (info .ne. 0) then
	  p = info
	  info = 3
	  return
      endif
      do 50 j=1,p
         call dprmut(f2kf2(1,j),p,iwork,1)
   50 continue
      return
      end
      subroutine dsnsm (x,ldx,y,sigma,ldsigm,nobs,npar,nnull,adiag,
     * tau,lamlim,ntbl,dout,iout,coef,svals,tbl,ldtbl,auxtbl,
     * iwork,liwa,work,lwa,job,info)
      integer ldx,ldsigm,nobs,npar,nnull,ntbl,iout(3),ldtbl,liwa,
     * iwork(liwa),lwa,job,info
      double precision x(ldx,npar),y(nobs),sigma(ldsigm,npar),
     * adiag(nobs),tau,lamlim(2),dout(5),coef(npar),svals(*),
     * tbl(ldtbl,3),auxtbl(3,3),work(lwa)
c
c Purpose: determine the generalized cross validation estimate of the 
c	smoothing parameter and fit model parameters for a penalized 
c	least squares problem with a semi-norm smoothing matrix.
c
c On Entry:
c   x(ldx,npar)		design matrix
c   ldx			leading dimension of x as declared in the 
c			calling program, must be at least max(nobs,npar)
c   y(nobs)		response vector
c   sigma(ldsigm,npar)	symmetric matrix that defines the semi-norm
c   ldsigm		leading dimension of sigma as declared
c			in the calling program
c   nobs		number of observations
c   npar		number of parameters
c   nnull		dimension of the null space of sigma
c   adiag(nobs)		"true" y values on entry if computation of 
c			predictive mse is requested
c   lamlim(2)		limits on lambda hat search (in log10(nobs*
c			lambda) scale) if user input limits are 
c			requested if lamlim(1) = lamlim(2) then lamhat
c			is set to (10**lamlim(1))/nobs
c   tau			multiplier controlling the amount of truncation
c			if truncation is requested (try tau = 1 
c			to start then try 10 and 100)
c   ntbl		number of evenly spaced values for 
c			log10(nobs*lambda) to be used in the initial 
c			grid search for lambda hat
c			if ntbl = 0 only a golden ratio search will be 
c			done and tbl is not referenced, if ntbl > 0
c			there will be ntbl rows returned in tbl
c   ldtbl		leading dimension of tbl as declared in the 
c			calling program	
c   job			integer with decimal expansion abcd
c			if a is nonzero then truncation is used
c			if b is nonzero then predictive mse is computed
c			   using adiag as true y
c			if c is nonzero then user input limits on search
c			   for lambda hat are used
c			if d is nonzero then the diagonal of the hat 
c			   matrix is calculated
c
c On Exit:
c   x(ldx,npar)		overwritten with many intermediate results
c   y(nobs)		predicted values
c   sigma(ldsigm,npar)	overwritten with the QR decomposition of the 
c			Cholesky factor of sigma
c   adiag(nobs)		diagonal elements of the hat matrix if requested
c   lamlim(2)		limits on lambda hat search 
c			(in log10(nobs*lambda) scale)
c   dout(5)		contains:
c  			1 lamhat   generalized cross validation 
c				   estimate of the smoothing parameter
c			2 penlty   smoothing penalty
c			3 rss	   residual sum of squares
c			4 tr(I-A)  trace of I - A
c			5 truncation ratio = 1/(1+(normk/(nobs*lamhat)))
c				   where normk = norm(R - R sub k)**2
c   iout(3)		contains:
c			1  npsing   number of positive singular
c				    values
c				    if info indicates nonzero info in 
c				    dsvdc then iout(1) contains info as
c				    it was returned from dsvdc
c			2  npar	    number of parameters
c			3  nnull    size of the null space of sigma
c   coef(npar)		coefficient estimates
c   svals(npar-nnull)	first npsing entries contain singular values 
c			of the matrix j2 
c			if info indicates nonzero info in dsvdc then 
c			svals is as it was returned from dsvdc
c   tbl(ldtbl,3)	column	contains
c			  1 	grid of log10(nobs*lambda) 
c			  2  	V(lambda)
c			  3     R(lambda) if requested
c   auxtbl(3,3)		auxiliary table
c			1st row contains:
c			    log10(nobs*lamhat), V(lamhat) and  
c			    R(lamhat) if requested
c			    where lamhat is the gcv estimate of lambda
c			2nd row contains:
c			    0, V(0) and  R(0) if requested
c			3rd row contains:
c			    0, V(infinity) and R(infinity) if requested
c   info		error indicator
c			   0 : successful completion
c			  -3 : nnull is too small (not fatal)
c			  -2 : log10(nobs*lamhat) >= lamlim(2) 
c			       (not fatal)
c			  -1 : log10(nobs*lamhat) <= lamlim(1)
c			       (not fatal)
c			   1 : dimension error 	
c			   2 : lwa (length of work) is too small
c			   3 : liwa (length of iwork) is too small
c			   4 : error in ntbl or tau
c			  100< info <200 : 100 + nonzero info returned
c					   from ddcom
c			  200< info <300 : 200 + nonzero info returned
c					   from dgcv
c
c Work Arrays:
c   work(lwa)		double precision work vector
c   lwa			length of work as declared in the calling 
c			program 
c			must be at least 
c			(npar-nnull)*(npar-2*nnull+2+nobs)+npar+nobs
c   iwork(liwa)		integer work vector
c   liwa		length of iwork as declared in the calling 
c			program
c			must be at least 2*npar - nnull 
c
c Subprograms Called Directly:
c	Gcvpack - ddcom dgcv
c
c Subprograms Called Indirectly:
c	Gcvpack - dcrtz dsgdc dcfcr drsap dvlop dtsvdc
c		  dpmse dvmin dvl dzdc dpdcr ddiag
c	Linpack - dchdc dqrdc dqrsl dtrsl dsvdc dtrco
c	Blas    - dcopy ddot dgemv dswap
c	Other 	- dcpmut dprmut dset 
c
c $Header: dsnsm.f,v 2.100.1.2 86/11/19 09:24:39 lindstrom Exp $
c
      integer bigp,pmh,pp1,lwa2,wsize,ldcaux,job1,job2,npsing,sinfo
      sinfo = 0
      info = 0
      iout(2) = npar
c			check dimensions
      pmh = npar - nnull
      if ((nobs .le. 0) .or. (npar .le. 0) .or. (nnull .lt. 0) 
     * .or.(pmh .le. 0)) then
         info = 1
         return
      endif
      wsize = pmh**2 +2*npar - nnull + pmh*(nobs - nnull) + pmh + nobs
      if (lwa .lt. wsize) then
         info = 2
         return
      endif
      wsize = npar + pmh
      if (liwa .lt. wsize) then
         info = 3
         return
      endif
      if ((ntbl .lt. 0) .or. (ntbl .gt. ldtbl) .or. (tau .lt. 0)) then
         info = 4
         return
      endif
c 			first p positions of iwork contain sgpvt
c 			first bigp-1  positions of work contain dcaux 
c			iwork vector starts at pp1
      pp1 = npar + 1
c 	work vector starts at bigp
      bigp = pmh**2 + 2*npar - nnull + 1
      ldcaux = pmh**2 + 2*npar - nnull
      lwa2 = lwa - ldcaux
      job1=job/1000
      job2=mod(job,1000)
c			decompose sigma and design matrix
      call ddcom (x,ldx,sigma,ldsigm,nobs,npar,nnull,tau,npsing,
     * svals,iwork,work,ldcaux,dout(5),work(bigp),lwa2,iwork(pp1),pmh,
     * job1,info)
      iout(1) = npsing
      iout(3) = nnull
      if (info .gt. 0) then
         info = info + 100
         return
      endif
      if (info .lt. 0) sinfo = info
c			compute lambda hat and other parameters
      call dgcv(x,ldx,y,sigma,ldsigm,nobs,npar,nnull,iwork,work,ldcaux,
     * adiag, lamlim,ntbl,dout,coef,tbl,ldtbl,auxtbl,work(bigp),lwa2,
     * job2,info)
      if (info .gt. 0) then
         info = info + 200
         return
      endif
      dout(5) = 1.0d0/(1.0d0+dout(5)/nobs/dout(1))
      if (sinfo .lt. 0) info = sinfo
      return
      end
      subroutine dsuy(y,nobs,nuobs,ytrue,c1,order,xrep,ssqrep,
     * job)
      integer nobs,nuobs,order(nobs),xrep(nobs),job
      double precision y(nobs),ytrue(nobs),c1(nuobs),ssqrep
c
c Purpose: compute B1'y, B1'ytrue and ssq for replication.
c
c On Entry:
c   y(nobs)  		response vector
c   nobs		number of observations
c   nuobs		number of unique design points
c   ytrue(nobs)		"true" response, if job is nonzero 1
c			not referenced if job = 0
c   c1(nuobs)		c1(i) contains the square root of the number of
c			replicates of the ith sorted design point
c   order(nobs)		order of sorted des
c   xrep(nobs)		xrep(i) = 1 if the ith sorted design point is a
c			replicate, 0 if not
c   job 		job is nonzero if B1'ytrue should be calculated 
c			job = 0 otherwise
c
c On Exit:
c   y(nuobs)  		B1'y
c   ytrue(nuobs)	B1'ytrue if job is nonzero
c   ssqrep		sum of squares for replication
c
c $Header: dsuy.f,v 2.100.1.1 86/10/07 12:57:32 lindstrom Exp $
c
      integer first,i,j
      double precision accum
c
      accum = 0.0d0
      ssqrep = 0.0d0
      call dprmut (y,nobs,order,0)
      if (job .ne. 0) call dprmut (ytrue,nobs,order,0)
c			compute ssq for replication
      first = 0
      do 20 i = 1,nobs
	  if (xrep(i) .eq. 1) then
	      accum = accum + y(i)
	  else if (first .eq. i - 1) then
	      first = i
	      accum = y(i)
	  else
	      accum = accum/(i-first)
	      do 10 j = first,i-1
	          ssqrep = (y(j)-accum)**2 + ssqrep
   10         continue
	      first = i
	      accum = y(i)
          endif
   20 continue
      if (xrep(nobs) .eq. 1) then
          accum = accum/(nobs + 1 - first)
          do 30 j = first,nobs
              ssqrep = (y(j)-accum)**2 + ssqrep
   30     continue
      endif
c			compute B1'y and B1'ytrue
      j = 0
      do 40 i = 1,nobs
 	  if (xrep(i) .eq. 0) then
 	      if (j .ne. 0) then
 		  y(j) = y(j) / c1(j)
 		  if (job .ne. 0) ytrue(j) = ytrue(j) / c1(j)
              endif
 	      j = j + 1
	      y(j) = y(i)
	      if (job .ne. 0) ytrue(j) = ytrue(i)
	  else
	      y(j) = y(j) + y(i)
	      if (job .ne. 0) ytrue(j) = ytrue(j) + ytrue(i)
 	  endif
   40 continue
      y(j) = y(j) / c1(j)
      if (job .ne. 0) ytrue(j) = ytrue(j) / c1(j)
      return 
      end
      subroutine dtpss(des,lddes,nobs,dim,m,s,lds,ncov,y,ntbl,adiag,
     * lamlim,dout,iout,coef,svals,tbl,ldtbl,auxtbl,work,lwa,
     * iwork,liwa,job,info)
      integer lddes,nobs,dim,m,lds,ncov,ntbl,iout(4),ldtbl,lwa,
     * liwa,iwork(liwa),job,info
      double precision des(lddes,dim),s(lds,*),y(nobs),
     * adiag(nobs),lamlim(2),dout(5),coef(*),svals(*),
     * tbl(ldtbl,3),auxtbl(3,3),work(lwa)
c
c Purpose: determine the generalized cross validation estimate of the 
c 	smoothing parameter and fit model parameters for a thin plate
c 	smoothing spline.
c
c On Entry:
c   des(lddes,dim) 	design for the variables to be splined
c   lddes		leading dimension of des as declared in calling
c   			program 
c   nobs		number of observations
c   dim			number of columns in des
c   m			order of the derivatives in the penalty
c   s(lds,ncov) 	design for the covariates. The covariates
c			must duplicate the replication structure of des.
c			See dptpss to handle covariates which do not.
c   lds			leading dimension of s as declared in calling
c   			program 
c   ncov		number of covariates 
c   y(nobs)		response vector
c   ntbl		number of evenly spaced values for 
c			log10(nobs*lambda) to be used in the initial 
c			grid search for lambda hat 
c			if ntbl = 0 only a golden ratio search will be 
c			done and tbl is not referenced, if ntbl > 0
c			there will be ntbl rows returned in tbl
c   adiag(nobs)	 	"true" y values on entry if predictive mse is 
c			requested
c   lamlim(2)		limits on lambda hat search (in log10(nobs*
c			lambda)	scale) if user input limits are 
c			requested if lamlim(1) = lamlim(2) then lamhat
c			is set to (10**lamlim(1))/nobs
c   ldtbl		leading dimension of tbl as declared in the 
c			calling program	
c   job			integer with decimal expansion abdc
c			if a is nonzero then predictive mse is computed
c			   using adiag as true y
c			if b is nonzero then user input limits on search
c			   for lambda hat are used
c			if c is nonzero then adiag will be calculated
c			if d is nonzero then there are replicates in the
c			   design
c
c On Exit:
c   des(lddes,dim)	sorted unique rows of des if job indicates that
c			there are replicates otherwise not changed
c   s(lds,ncov)		unique rows of s sorted to correspond to des
c   y(nobs)		predicted values
c   adiag(nobs)		diagonal elements of the hat matrix if requested
c   lamlim(2)		limits on lambda hat search 
c			(in log10(nobs*lambda) scale)
c   dout(5)		contains:
c  			1  lamhat   generalized cross validation 
c				    estimate of the smoothing parameter
c			2  penlty   smoothing penalty
c			3  rss	    residual sum of squares
c			4  tr(I-A)  trace of I - A
c   			5  ssqrep   sum of squares for replication
c   iout(4)		contains:
c			1  npsing   number of positive singular
c				    values (npsing = nuobs - ncts).
c				    if info indicates nonzero info in 
c				    dsvdc then npsing contains info as 
c				    it was returned from dsvdc.
c			2  npar	    number of parameters
c				    (npar = nuobs + ncts)
c			3  ncts     dimension of the polynomial space 
c				    plus ncov	
c				    ((m+dim-1 choose dim) + ncov)
c			4  nuobs    number of unique rows in des
c   coef(npar)		coefficient estimates [beta':alpha':delta']'
c			coef must have a dimension of at least nuobs+
c			ncts
c   svals(npar-nnull)	singular values of the matrix j2 if info = 0
c			if info indicates nonzero info from dsvdc then 
c			svals is as it was returned from dsvdc.
c   tbl(ldtbl,3)	column	contains
c			  1 	grid of log10(nobs*lambda) 
c			  2  	V(lambda)
c			  3     R(lambda) if requested
c   auxtbl(3,3)		auxiliary table
c			1st row contains:
c			    log10(nobs*lamhat), V(lamhat) and  
c			    R(lamhat) if requested
c			    where lamhat is the gcv estimate of lambda
c			2nd row contains:
c			    0, V(0) and  R(0) if requested
c			3rd row contains:
c			    0, V(infinity) and R(infinity) if requested
c   info		error indicator
c			  0 : successful completion
c			 -1 : log10(nobs*lamhat) <= lamlim(1)
c			      (not fatal)
c			 -2 : log10(nobs*lamhat) >= lamlim(2)
c			      (not fatal)
c			  1 : dimension error 	
c			  2 : error in dreps, covariates do not 
c			      duplicate the replication structure of des
c			  3 : lwa (length of work) is too small
c			  4 : liwa (length of iwork) is too small
c			  10 < info < 20 : 10 + nonzero info returned 
c					   from dsetup
c			  100< info <200 : 100 + nonzero info returned
c					   from dsgdc1
c			  200< info <300 : 200 + nonzero info returned
c					   from dgcv1
c
c Work Arrays:
c   work(lwa)		double precision work vector
c   lwa			length of work as declared in the calling 
c			program 
c			Must be at least nuobs(2+ncts+nuobs)+nobs
c   iwork(liwa)		integer work vector
c   liwa		length of iwork as declared in the calling 
c			program
c			Must be at least 2*nobs + nuobs - ncts
c
c Subprograms Called Directly:
c       Gcvpack - dreps duni dsuy dsetup dsgdc1 dgcv1 
c       Other   - dprmut mkpoly
c
c Subprograms Called Indirectly:
c	Gcvpack - dcfcr1 drsap dvlop dsvtc dpdcr dpmse 
c		  dvmin dvl dmaket dmakek ddiag
c	Linpack - dchdc dqrdc dqrsl dtrsl dsvdc
c	Blas    - ddot dcopy dgemv
c	Other 	- dprmut dset dftkf fact mkpoly 
c
c $Header: dtpss.f,v 2.100.1.1 86/10/07 12:57:52 lindstrom Exp $
c
      integer ncts,jrep,jadiag,p1,p1p1,p2,p2p1,p3,p3p1,ip1,ip1p1,ip2,
     * ip2p1,nuobs,p4,p4p1,j,i,lwa2,wsize,q1,q1p1,npsing,
     * jobgcv
      double precision ssqrep
      integer mkpoly
c
c
      info = 0
      ncts = mkpoly(m,dim) + ncov
      iout(3) = ncts
      jrep = mod(job,10)
      jadiag = mod(job,100)/10
      jobgcv = job/10
c			check dimensions
      if((nobs .le. 0) .or. (m .le. 0) .or. (dim .le. 0) .or.
     * (ntbl .lt. 0) .or. (ntbl .gt. ldtbl) .or. 
     * (2*m-dim .le. 0)) then
         info = 1
         return
      endif
c			set up pointers for iwork vector
c  			first nobs positions of iwork contain order
      ip1 = nobs
c			next nobs positions of iwork contain xrep    
      ip1p1 = ip1 + 1
      ip2 = ip1 + nobs
c			rest of iwork is a integer work vector
      ip2p1 = ip2 + 1
      if (jrep .ne. 0) then
         call dreps(des,lddes,nobs,dim,s,lds,ncov,0,work,iwork,nuobs,
     *    iwork(ip1p1),1,info)
         if (info .ne. 0) then
            info = 2
            return
         endif
c			put unique row of des into des
         call duni(des,lddes,nobs,dim,iwork(ip1p1),des,lddes)
c			put unique row of s into s
         call duni(s,lds,nobs,ncov,iwork(ip1p1),s,lds)
      endif
      if (jrep .eq. 0) then
         work(1)= 0
         nuobs = nobs
         ssqrep = 0.0d0
      endif
      iout(2) = nuobs + ncts
      iout(4) = nuobs
c			check size of work vectors
      wsize = nuobs*(2+ncts+nuobs)+nobs
      if (lwa .lt. wsize) then
         info = 3
         return
      endif
      wsize = 2*nobs + nuobs - ncts
      if (liwa .lt. wsize) then
         info = 4
         return
      endif
c			set up pointers for c1,[tu:su1],ku,fgaux in work
c
c		c1  	 runs from 1    to p1,   p1 = nuobs
c		[tu:su1] runs from p1p1 to p2,   p2 = p1 + nuobs*ncts
c		ku	 runs from p2p1 to p3,   p3 = p2 + nuobs**2
c		fgaux    runs from p3p1 to p4,   p4 = p3 + ncts
c		the rest of work is a work vector

c  			after the call to dsetup 
c		f2'k f2  runs from q1p1 to p3, q1 = p2 + nuobs*ncts+ncts
c
      p1 = nuobs
      p1p1 = p1 + 1
      p2 = p1 + nuobs*ncts
      p2p1 = p2 + 1
      q1 = p2 + nuobs*ncts+ncts
      q1p1 = q1 + 1
      p3 = p2 + nuobs**2
      p3p1 = p3 + 1
      p4 = p3 + ncts
      p4p1 = p4 + 1
      lwa2 = lwa - (nuobs*(1+ncts+nuobs) + ncts)
c			set up structures needed for dsgdc1 and dgcv1
      call dsetup(des,lddes,s,lds,dim,m,ncov,nuobs,work,work(p1p1),
     * nuobs,ncts,work(p3p1),work(p2p1),nuobs,work(p4p1),iwork(ip2p1),
     * info)
      if (info .ne. 0) then
         info = info + 10
         return
      endif

c			decompose f2' k f2 
      npsing = nuobs - ncts
      call dsgdc1(work(q1p1),nuobs,npsing,svals,iwork(ip2p1),
     * work(p4p1),lwa2,info)
      iout(1) = npsing
      if (info .gt. 0) then
         info = info + 100
         return
      endif
c			setup y
      if (jrep .ne. 0) then
         call dsuy(y,nobs,nuobs,adiag,work,iwork,iwork(ip1p1),
     *    ssqrep,jadiag)
      endif
      dout(5) = ssqrep
c			compute lambda hat and other parameters
      call dgcv1(work(p2p1),nuobs,y,nuobs,nobs,work(p1p1),nuobs,ncts,
     * work(p3p1),svals,adiag,lamlim,ssqrep,ntbl,dout,coef,tbl,ldtbl,
     * auxtbl,work(p4p1),lwa2,jobgcv,info)
      if (info .gt. 0) then
         info = info + 200
         return
      endif
c			if there are replicates then rescale the coef. 
c			vector,	the predicted values, and diagonal of A
      if (nuobs .ne. nobs) then
         do 10 i = 1,nuobs
            coef(ncts+i) = coef(ncts+i) * work(i)
   10    continue
         j = nuobs
         do 20 i = nobs,1,-1 
            y(i)=y(j)/work(j)
            if (jadiag .ne. 0) then
               adiag(i) = adiag(j)/work(j)**2
            endif
            if (iwork(ip1 + i) .eq. 0) then
               j = j - 1
            endif
   20    continue
      endif
c			undo the permutation of the predicted values and
c			adiag if necessary
      if (jrep .ne. 0) then
         call dprmut (y,nobs,iwork,1)
         if (jadiag .ne. 0) call dprmut (adiag,nobs,iwork,1)
      endif
      return
      end
      subroutine dtsvdc(x,ldx,n,p,minrat,k,s,u,ldu,v,ldv,normk,
     * work,lwa,iwork,liwa,job,info)
      integer ldx,n,p,k,ldu,ldv,lwa,iwork(liwa),liwa,job,info
      double precision x(ldx,p),minrat,s(p),u(ldu,*),v(ldv,p),
     * normk,work(lwa)
c
c Purpose: form the singular value decomposition of a truncated matrix 
c	obtained by first taking a qr decomposition with pivoting.
c
c On Entry:
c   x(ldx,p)   		matrix to be decomposed, x is destroyed
c   ldx     		leading dimension of x in the calling program
c   n       		number of rows in x
c   p       		number of columns in x   (p <= n)
c   minrat  		minimum ratio to determine truncation
c   ldu     		leading dimension of u as declared in the 
c			calling program
c   ldv     		leading dimension of v as declared in the 
c			calling program
c   job     		controls the computation of the singular vectors
c               	it has the decimal expansion ab with the meaning
c                            a = 0 no left singular vectors
c                            a = 1 all n left singular vectors
c                            a > 1 first p left singular vectors
c                            b = 0 no right singular vectors
c                            b > 0 all p right singular vectors
c
c On Exit:
c   k       		number of positive singular values 
c   s(p)		an approximation to the first k singular values
c			in decreasing order
c   u(ldu,m)   		first k singular vectors if job > 1 or all n 
c			left singular vectors if joba = 1 otherwise it 
c			is not accessed it may be identified with x in 
c			the subroutine call if joba > 1
c   v(ldv,p)   		the first k right singular vectors if job b > 0
c               	otherwise it is not accessed
c   normk		norm of the k by k lower right sub matrix of r
c   info    		error indicator
c		   	   0 : successful completion
c		  	   info > 0 : info was returned from dsvdc
c
c Work Arrays:
c   iwork(liwa) 	integer work vector, holds the pivot vector from
c			dqrdc
c   liwa		must be at least p
c   work(lwa)		work vector
c   lwa			must be at least p**2 + p + n if joba > 0, 
c			otherwise it must be at least 2*p
c
c Subprograms Called Directly
c	Linpack - dsvdc dqrdc dqrsl
c   	Blas    - ddot dcopy
c	Other   - dset dprmut
c
c $Header: dtsvdc.f,v 2.100.1.1 86/10/07 12:58:56 lindstrom Exp $
c
      integer i,j,jj,pp1,ppnp1,nmk,sjob,locinf
      double precision trxpx,accum,dummy(1),mintr
      double precision ddot
c
      info = 0
      call dset(p,0.0d0,s,1)
      pp1 = p+1
      ppnp1 = p+n+1
c                       calculate trace of x' x
      trxpx = 0.0d0
      do 10 j = 1,p
         trxpx = trxpx+ddot(n,x(1,j),1,x(1,j),1)
   10 continue
      mintr = trxpx * minrat
c                       qr decomposition of x
      do 20 j = 1,p
         iwork(j) = 0
   20 continue
      call dqrdc(x,ldx,n,p,work,iwork,work(pp1),1)
c                       calculate ratios for the truncated matrix
      accum = 0.0d0
      k = p
      do 30 i = p,1,-1
         accum = accum+ddot(pp1-i,x(i,i),ldx,x(i,i),ldx)
         if (accum .lt. mintr) then
	     k = i
	     normk = accum
         endif
   30 continue
      if (job/10.le.0) then
c                       no left singular vectors
c                       copy rk' (transpose of the first k rows of r
c			from qr decomposition) into x
         do 40 j = 1,k 
            call dcopy(p,x(j,1),ldx,x(1,j),1)
            call dset(j-1,0.0d0,x(1,j),1)
   40    continue
c                       svd of rk'
         sjob = 0
         if (mod(job,10).ne.0) sjob = 20
         call dsvdc(x,ldx,p,k,s,work,v,ldv,dummy,0,work(pp1),sjob,info)
         if (info.ne.0) return
         if (mod(job,10).ne.0) then
            do 50 j = 1,k
               call dprmut(v(1,j),p,iwork,1)
   50       continue
         endif
      else
c                       u or u1 is to be created
         i = n+1
c                       copy rk' to work
         do 60 j = 1,k 
            i = i+p
            call dcopy(p,x(j,1),ldx,work(i),1)
            call dset(j-1,0.0d0,work(i),1)
   60    continue
c                       create u2 if requested
         if (job/10.eq.1) then
            nmk = n-k
            do 70 i = 1,nmk 
               j = n+1-i
               call dset(n,0.0d0,work(pp1),1)
               work(p+j) = 1.0d0
               call dqrsl(x,ldx,n,min(j,p),work,work(pp1),work(pp1),
     *          dummy,dummy,dummy,dummy,10000,locinf)
               call dcopy(n,work(pp1),1,u(1,j),1)
   70       continue
         endif
         do 80 j = 1,k 
            jj = k+1-j
            call dset(n,0.0d0,work(pp1),1)
            i = p+jj
            work(i) = 1.0d0
            call dqrsl(x,ldx,n,jj,work,work(pp1),work(pp1),dummy,dummy,
     *       dummy,dummy,10000,locinf)
            call dcopy(n,work(pp1),1,u(1,jj),1)
   80    continue
c			  svd of rk'
         sjob = 1
         if (mod(job,10).ne.0) then
            sjob = 21
         endif
         call dsvdc(work(ppnp1),p,p,k,s,work,v,ldv,work(ppnp1),p,
     *    work(pp1),sjob,info)
         if (info.ne.0) return
         do 100 i = 1,n 
            call dcopy(k,u(i,1),ldu,work,1)
            jj = n+1
            do 90 j = 1,k 
               jj = jj+p
               u(i,j) = ddot(k,work,1,work(jj),1)
   90       continue
  100    continue
         if (mod(job,10).ne.0) then
c		 undo pivots on right singular vectors
            do 110 j = 1,k
               call dprmut(v(1,j),p,iwork,1)
  110       continue
         endif
      endif
      end
      subroutine duni(x,ldx,nobs,ncx,xrep,xu,ldxu)
      integer ldx,nobs,ncx,ldxu,xrep(nobs)
      double precision x(ldx,*),xu(ldxu,*)
c
c Purpose: compute xu.
c
c On Entry:
c   x(ldx,ncx)		a matrix to be reduced to unique rows
c   ldx			leading dimension of x as declared in the
c			calling	program 
c   nobs		number of observations
c   ncx			number of columns in x
c   xrep(nobs)		xrep(i) contains 1 if ith row of x is a 
c			replicate row, 0 if not
c   ldxu		leading dimension of xu as declared in the 
c			calling	program 
c On Exit:
c   xu(ldxu,ncx) 	unique rows of x
c			may be identified with x in the calling sequence
c
c $Header: duni.f,v 2.100.1.1 86/10/07 12:59:15 lindstrom Exp $
c
      integer i,j,k
c
      j = 0
      do  20 i = 1,nobs
 	  if (xrep(i) .eq. 0) then
 	      j = j + 1
 	      do 10 k = 1,ncx
 	          xu(j,k) = x(i,k)
   10         continue
 	  endif
   20 continue
      return 
      end
      double precision function dvl(lgnlam,svals,z,npsing)
      integer npsing
      double precision lgnlam,svals(npsing),z(npsing)
c
c Purpose: evaluate the cross-validation function with a semi-norm.
c
c On Entry:
c   lgnlam		log10(nobs*lambda) where lambda is the value of
c			lambda for which V is evaluated
c   svals(npsing)	singular values 
c   z(npsing)		data vector in canonical coordinates
c   npsing		number of positive svals 
c
c On Exit:
c   dvl			V(lambda)
c
c $Header: dvl.f,v 2.100.1.1 86/10/07 12:59:22 lindstrom Exp $
c
      integer j
      double precision nlam,numrtr,denom,factor
c
      common / gcvcom / addend,rss,tria,n,h
      integer n,h
      double precision rss,tria,addend
c     			see dvlop for definition of common block 
c			variables
c
      nlam = 10**lgnlam
      numrtr = addend
      denom = dble(n - h - npsing)
      do 10 j = 1,npsing 
         factor = 1.0d0/(1.0d0 + (svals(j)**2)/nlam)
         numrtr = numrtr + (factor*z(j))**2
         denom = denom + factor
   10 continue
      rss=numrtr
      tria=denom
      dvl=dble(n)*numrtr/denom**2
      return
      end
      subroutine dvlop(z,svals,nobs,nnull,npsing,inadd,ssqw2,
     * lamlim,ntbl,nlamht,tbl,ldtbl,auxtbl,dout,job,info)
      integer nobs,nnull,npsing,ntbl,ldtbl,job,info
      double precision z(npsing),svals(npsing),inadd,ssqw2,lamlim(2),
     * nlamht,tbl(ldtbl,3),auxtbl(3,3),dout(2)
c
c Purpose: determine the optimal lambda for the generalized cross 
c	validation function given singular values and the data vector 
c	in canonical coordinates.
c
c On Entry:
c   z(npsing)		data vector in canonical coordinates
c   svals(npsing)	singular values 
c   nobs		number of observations
c   nnull		dimension of the null space of sigma
c   npsing		number of positive elements of svals 
c   inadd		constant term in expression for V
c   ssqw2		squared length of w2
c   lamlim(2)		limits on lambda hat search (in log10(nobs*
c			lambda) scale) if user input limits are 
c			requested. if lamlim(1) = lamlim(2) then nlamht
c			is set to 10**lamlim(1)
c   ntbl		number of evenly spaced values for 
c			log10(nobs*lambda) to be used in the initial 
c			grid search for lambda hat
c			if ntbl = 0 only a golden ratio search will be 
c			done and tbl is not referenced, if ntbl > 0
c			there will be ntbl rows returned in tbl
c   ldtbl		leading dimension of tbl as declared in the 
c			calling program	
c   job      		if job is nonzero then user input limits on 
c			lambda hat search are used
c
c On Exit:
c   lamlim(2)		limits on lambda hat search 
c			(in log10(nobs*lambda) scale)
c   nlamht		nobs*(lambda hat) where lambda hat is the gcv 
c			estimate of lambda
c   tbl(ldtbl,3)	column	contains
c			  1 	grid of log10(nobs*lambda) 
c			  2  	V(lambda)
c   auxtbl(3,3)		auxiliary table
c			1st row contains:
c			    log10(nobs*lambda hat), V(lambda hat)
c			2nd row contains:
c			    0, V(0) 
c			3rd row contains:
c			    0, V(infinity) 
c   dout(2)		contains:
c			1  rss
c			2  tr(I-A)
c   info		error indicator
c			  0 : successful completion
c			 -1 : log10(nlamht) <= lamlim(1) (not fatal)
c			 -2 : log10(nlamht) >= lamlim(2) (not fatal)
c			  1 : svals(1) = 0.0d0
c			  2 : npsing is incorrect
c			  3 : lamlim(1) > lamlim(2)
c
c Subprograms Called Directly:
c	Gcvpack - dvmin 
c
c Subprograms Called Indirectly:
c	Gcvpack - dvl
c
c $Header: dvlop.f,v 2.100.1.1 86/10/07 12:59:38 lindstrom Exp $
c
      integer i,k
      double precision vlamht,w
      double precision dvmin
c
      common / gcvcom / addend,rss,tria,n,h
      integer n,h
      double precision addend,rss,tria,machpr,one
c
      info = 0
      one = 1.0d0
      machpr = 1.0d0
   10 machpr = machpr/2.0d0
      if (one .lt. 1.0d0 + machpr) goto 10
      machpr = machpr*2.0d0
c
      n=nobs
      h=nnull
      addend = inadd
      if (svals(1) .eq. 0.0d0) then
         info = 1
         return
      endif
      k = 0
      do 20 i = 1,npsing
         if (svals(i) .gt. 0) then
            k = i
         endif
   20 continue
      if (k .ne. npsing) then
	 info = 2
    	 return
      endif
      if (job .ne. 0 .and. (lamlim(1) .gt. lamlim(2))) then
         info = 3
         return
      endif
      if (job .eq. 0) then
         lamlim(2) = 2.0d0*dlog10(svals(1))+2.0d0
         lamlim(1) = 2.0d0*dlog10(svals(npsing))-2.0d0
      endif
      nlamht = dvmin (lamlim(1),lamlim(2),svals,z,npsing,ntbl,tbl,
     * ldtbl,vlamht,info)
      dout(1) = rss
      dout(2) = tria
c			compute auxtbl
      auxtbl(1,1)=nlamht
      auxtbl(1,2)=vlamht
c			lambda = 0
      auxtbl(2,1)=0.0d0
      auxtbl(2,2)=0.0d0
      if ((nobs-nnull) .ne. npsing) then
         auxtbl(2,2)=inadd*(nobs)/(nobs-nnull-npsing)**2
      endif
      if ((nobs-nnull) .eq. npsing) then
         w=0.0d0
         do 30 i=npsing,1,-1
c           w=w+(z(i)*svals(npsing)**2/(svals(i)**2))**2
            w=w+(z(i)*(svals(npsing)/svals(i))**2)**2
   30    continue
         auxtbl(2,2)=nobs*w
         w=0.0d0
         do 40 i=npsing,1,-1
            w=w+(svals(npsing)/svals(i))**2
   40    continue
         auxtbl(2,2)=auxtbl(2,2)/(w**2)
      endif
c			lambda = infinity
      auxtbl(3,1)=0.0d0
      auxtbl(3,2)=ssqw2/(nobs - nnull)
      nlamht = 10**nlamht
      return
      end
      double precision function dvmin(lower,upper,svals,z,npsing,
     * ntbl,tbl,ldtbl,vlamht,info)
      integer npsing,ntbl,ldtbl,info
      double precision lower,upper,svals(npsing),z(npsing),
     * tbl(ldtbl,3),vlamht
c
c Purpose: evaluate V(lambda) for a grid of ln(nobs*lambda) values 
c	between	lower and upper, store these in the array tbl, and find
c	minimizer of v.
c
c On Entry:
c   lower		lower bound of interval (in nobs*ln(lambda)
c			scale) over which V(lambda) is to be minimized
c   upper		upper bound of interval (in nobs*ln(lambda)
c			scale) over which V(lambda) is to be minimized
c   svals(npsing)	singular values 
c   z(npsing)		data vector in canonical coordinates
c   npsing		number of positive elements of svals 
c   ntbl		number of evenly spaced values for 
c			ln(nobs*lambda)	to be used in the initial grid 
c			search for lambda hat
c			if ntbl = 0 only a golden ratio search will be 
c			done and tbl is not referenced, if ntbl > 0
c			there will be ntbl rows returned in tbl
c   ldtbl		leading dimension of tbl as declared in the 
c			calling program	
c   
c On Exit:
c   tbl(ldtbl,3)	column	contains
c			  1 	grid of ln(nobs*lambda) 
c			  2  	V(lambda)
c   vlamht		V(lambda hat)
c   dvmin		ln(nobs*lambda hat)
c   info		error indicator
c			  0 : successful completion
c			 -1 : dvmin <= lower (not fatal)
c			 -2 : dvmin >= upper (not fatal)
c
c Subprograms Called Directly:
c	Gcvpack - dvl
c
c $Header: dvmin.f,v 2.100.1.1 86/10/07 12:59:51 lindstrom Exp $
c

      double precision a,b,c,d,vc,vd,del,k1,k2,x,v
      integer j,jmin,k
      double precision dvl
c				null interval
      if (lower .eq. upper) then
	 dvmin = lower
	 info = -1
	 vlamht = dvl(lower,svals,z,npsing)
	 do 10 j = 1, ntbl
	    tbl(j,1) = lower
	    tbl(j,2) = vlamht
   10    continue
	 return
      end if
c				non-null interval
      info = 0
      a = lower
      b = upper
      if (ntbl .eq. 1) then
	 x = (a + b)/2
	 tbl(1,1) = x
	 tbl(1,2) = dvl(x,svals,z,npsing)
      else if (ntbl .ge. 2) then
c			do grid search
	 v=dvl(lower,svals,z,npsing)*2.0d0
	 del=(upper-lower)/(ntbl-1)
	 do 20 j = 1, ntbl
	    tbl(j,1) = lower + (j - 1) * del
	    tbl(j,2) = dvl(tbl(j,1),svals,z,npsing)
	    if (tbl(j,2) .le. v) then
	       jmin = j
	       v = tbl(j,2)
	    endif
   20    continue	      
	 a=tbl(jmin,1)-del
	 b=tbl(jmin,1)+del
      end if
c			do golden ratio search			
      k1=(3.0d0-dsqrt(5.0d0))/2.0d0
      k2=(dsqrt(5.0d0)-1)/2.0d0
      c = a + k1*(b - a)
      d = a + k2*(b - a)
      vc = dvl(c,svals,z,npsing)
      vd = dvl(d,svals,z,npsing)
      do 30 k=1,50
	 if (vd .lt. vc) then
	    a = c
	    c = d
	    d = a + k2*(b - a)
	    vc = vd
	    vd = dvl(d,svals,z,npsing)
	 else
	    b = d
	    d = c
	    c = a + k1*(b - a)
	    vd = vc
	    vc = dvl(c,svals,z,npsing)
	 end if
   30 continue
      x=(a+b)/2
      if (x .le. lower) info = -1
      if (x .ge. upper) info = -2
      vlamht=dvl(x,svals,z,npsing)
      dvmin = x
      return
      end
      subroutine dzdc(x,ldx,nobs,npar,pmh,tau,fgaux,svals,npsing,
     * v,ldv,normk,work,lwa,iwork,liwa,job,info)
      integer ldx,nobs,npar,pmh,npsing,ldv,lwa,iwork(liwa),liwa,info,job
      double precision x(ldx,npar),tau,fgaux(*),svals(*),v(ldv,pmh),
     * normk,work(lwa)
c
c Purpose: create and decompose j and decompose j2
c
c On Entry:
c   x(ldx,npar)		z = design matrix after sigma rotations 
c   ldx		     	leading dimension of x as declared in the 
c			calling	program
c   nobs   		number of observations
c   npar      		number of parameters
c   pmh    		npar minus the dimension of the null space of 
c			the semi-norm
c   tau			multiplier controlling the amount of truncation 
c			if truncation is requested
c   ldv    		leading dimension of v as declared in the 
c			calling	program
c   job			controls use of truncated singular value 
c			decomposition
c			   if job is nonzero then dtsvdc is used
c			   if job = 0 then dsvdc is called directly
c
c On Exit:
c   x(ldx,npar)		overwritten with many intermediate results
c   fgaux(npar-pmh)	auxiliary information on the qr decomposition 
c			of the matrix z2
c   svals(npsing)      	singular values of the matrix j2 if info = 0.
c			if info = 4,5,6 or 7 then svals is as 
c			it was returned from dsvdc.
c   npsing       	if info = 0 then npsing contains the number 
c			of singular values calculated 
c			if info = 4,5,6 or 7 then npsing contains info
c			as it was returned from dsvdc.
c   v(ldv,pmh) 		right singular vectors of the matrix j2
c   normk		frobenius norm of the k by k lower right sub 
c			matrix of j2
c   info    		error indicator
c			  0 : successful completion
c			  1 : lwa (length of work) is too small
c			  2 : tau < 0
c			  3 : transpose of j2 is necessary 
c				(npar .gt. nobs) and npar .gt. ldx
c			  4 : error in dtsvdc (using j2')
c			  5 : error in dsvdc (using j2')
c			  6 : error in dtsvdc
c			  7 : error in dsvdc
c
c Work Arrays:
c   work(lwa)		double precision work vector
c   lwa     		length of work as declared in the calling 
c			program	Must be at least 
c			pmh*(nobs - (npar-pmh) ) + pmh + nobs
c   iwork(liwa)		integer work vector
c   liwa		length of integer work vector, must be at least
c			pmh
c
c Subprograms Called Directly:
c	Gcvpack - dtsvdc 
c	Linpack - dqrdc dqrsl dsvdc
c	Blas    - dcopy dswap
c
c Subprograms Called Indirectly:
c	Linpack - dqrdc dsvdc dqrsl
c	Blas    - ddot dcopy
c	Other   - dset dprmut dswap
c
c $Header: dzdc.f,v 2.100.1.1 86/10/07 13:00:03 lindstrom Exp $
c
      double precision dummy(1),machpr,one,minrat
      integer i,j,idummy(1),nnull,hp1,nmh,pmhp1,wsize,locinf
c
      info = 0
c
      one = 1.0d0
      machpr = 1.0d0
   10 machpr = machpr/2.0d0
      if (one .lt. 1.0d0 + machpr) goto 10
      machpr = machpr*2.0d0
      minrat = machpr*tau
c			check dimensions
      wsize = pmh**2 + pmh + nobs
      if (lwa .lt. wsize) then
         info = 1
         return
      endif
      if (tau .lt. 0) then
	 info = 2
	 return
      endif
c			z1 is the first pmh columns of x
c			z2 is next h columns
      pmhp1 = pmh + 1
      nnull = npar - pmh
      hp1 = nnull + 1
      nmh = nobs - nnull
c			decompose z2 into f and g
      call dqrdc (x(1,pmhp1),ldx,nobs,nnull,fgaux,idummy,dummy,0)
c			create j as f' z1 
      do 30 j = 1,pmh
         call dqrsl (x(1,pmhp1),ldx,nobs,nnull,fgaux,x(1,j),dummy,
     *    x(1,j),dummy,dummy,dummy,01000,locinf)
   30 continue
c			decompose j2 (last npar - nnull rows of j)
c			as udv'
      if (pmh .gt. nmh) then
c			transpose the j2 matrix
         if (npar .gt. ldx) then
            info = 3
            return
         endif
         do 40 i = 0, nmh-1
            call dcopy(pmh,x(hp1+i,1),ldx,work(i*pmh),1)
   40    continue
         do 50 i = 0, nmh-1
            call dcopy (pmh, work(i*pmh),1,x(hp1,i+1),1)
   50    continue
         if (job .ne. 0) then
            call dtsvdc (x(hp1,1),ldx,pmh,nmh,minrat,npsing,svals,
     *       x(hp1,1),ldx,v,ldv,normk,work,lwa,iwork,liwa,21,info)
            if (info .ne. 0) then
	       npsing = info
               info = 4
	       return
            endif
         else
            call dsvdc (x(hp1,1),ldx,pmh,nmh,svals,work(pmhp1),x(hp1,1),
     *       ldx,v,ldv,work,21,info)
            if (info .ne. 0) then
	       npsing = info
               info = 5
	       return
            endif
            npsing=nmh
	    normk=0.0d0
         endif
         do 60 i = 1, npsing
            call dswap(pmh,x(hp1,i),1,v(1,i),1)
   60    continue
      else
         if (job .ne. 0) then
            call dtsvdc (x(hp1,1),ldx,nmh,pmh,minrat,npsing,svals,
     *       x(hp1,1),ldx,v,ldv,normk,work,lwa,iwork,liwa,21,info)
            if (info .ne. 0) then
	       npsing = info
               info = 6
	       return
            endif
         else
            call dsvdc (x(hp1,1),ldx,nmh,pmh,svals,work(nmh+1),x(hp1,1),
     *       ldx,v,ldv,work,21,info)
            if (info .ne. 0) then
	       npsing = info
               info = 7
	       return
            endif
            npsing=pmh
	    normk = 0.0d0
         endif
      endif
      return
      end
      integer function fact(i)
      integer i
c 
c Purpose: quick factorial function for the bspline routine
c	returns zero for negative i.
c
c On Entry:
c   i			a non-negative integer
c On Exit:
c   fact		i factorial
c
c $Header: fact.f,v 2.100.1.1 86/10/07 13:00:27 lindstrom Exp $
c
      integer j
      fact = 0
      if (i .ge. 0) fact = 1
      if (i .le. 1) return
      do 10 j = 2,i
	 fact = fact*j
   10	continue
      return
      end
      integer function mkpoly(m,dim)
      integer m,dim
c
c  Purpose: compute the binomial coefficient of m + dim - 1 choose dim.
c  	This is the dimension of the space of polynomials which are in
c  	the null space of the smoothing penalty. Uses Andy Jaworski's
c	binomial coefficient algorithm that only requires integer
c	arithmetic.
c
c  On Entry:
c   m			order of derivatives in the penalty
c   dim	 		dimension of the variables to be splined
c
c  On Exit:
c   mkploy		(m + dim - 1) choose dim
c
c $Header: mkpoly.f,v 2.100.1.1 86/10/07 13:00:35 lindstrom Exp $
c
      integer i,j,k,k1,kcoef,n
c 			compute binomial coefficient 
c			m + dim - 1 choose dim
      n = m + dim - 1
      k1 = dim
      if (k1 .gt. n .or. k1 .lt. 0) then
         mkpoly = 0
         return
      endif
      k = k1
      if ((n - k1) .lt. k) then
         k = n - k1
      endif
      kcoef = 1
      j = n - k
      do 10 i = 1, k
         j = j + 1
         kcoef = (kcoef * j) / i
   10 continue
      mkpoly = kcoef
      return
      end
cc      subroutine mktpar(bign,nobs,testno,bandwd,tpar,t,f0,
cc     *   k0,info)
cc      double precision bandwd,tpar(4,*),t(nobs),f0(nobs),k0(nobs)
cc      integer bign,nobs,testno,info

c
c  Purpose: To construct the "true" parameter vectors for an
c	integral equation test case.
c
c  On Entry:
c   bign		dimension of Hilbert space
c   nobs		number of observations
c   testno		test number (determines the spread of the peaks)
c   bandwd		band width for the filter function
c   t(nobs)		the location of the observations
c
c  On Exit:
c   tpar(4,bign/2+1)	values of the "true" parameters
c   f0(nobs)		the "true" values of f at locations t
c   k0(nobs)		the "true" values of k at locations t
c   info		error indicator. 0 indicates successful 
c			completion. Other values are:
c				1 - invalid testno
c				2 - bign <= 1 or nobs <= 1
c				3 - bign not even 
c  $Header: mktpar.f,v 2.100.1.2 86/10/13 14:09:36 lindstrom Exp $ 
cc      double precision rt2pi,twopi,sn,cs,fo,f00,ko,k00,fdiff,
cc     *	iovern
cc      integer j,i
c     
cc      twopi = 8.0d0 * atan(1.0d0)
cc      rt2pi = sqrt(twopi)
cc      info = 1
cc      if ((testno .le. 0) .or. (testno .gt. 4)) return
cc      info = 2
cc      if ((bign .le. 1) .or. (nobs .le. 1)) return
cc      info = 3
cc      if (int(bign/2) .ne. dble(bign)/2.0d0) return
cc      info = 0
cc      fdiff = f00(0,testno,rt2pi) - f00(1,testno,rt2pi)
cc      call dset(4*(bign/2+1),0.0d0,tpar,1)
cc      do 20 j = 0, bign/2 
cc	  do 10 i = 1,bign
cc	      iovern = dble(i)/dble(bign)
cc	      sn = sin(twopi*j*iovern)
cc	      cs = cos(twopi*j*iovern)
cc	      fo = f00(iovern,testno,rt2pi)+( iovern-0.5)*fdiff
cc	      tpar(1,j+1)= tpar(1,j+1)+cs*fo
cc	      tpar(2,j+1)= tpar(2,j+1)+sn*fo
cc	      ko = k00(iovern,bandwd,rt2pi)
cc	      tpar(3,j+1)= tpar(3,j+1)+cs*ko
cc	      tpar(4,j+1)= tpar(4,j+1)+sn*ko
cc   10     continue
cc   20 continue
cc      call dscal(4*(bign/2+1),1.0d0/dble(bign),tpar,1)
cc      do 30 j = 1,nobs
cc          f0(j) = f00(t(j),testno,rt2pi)+(t(j)-0.5)*fdiff
cc          k0(j) = k00(t(j),bandwd,rt2pi)
cc   30 continue
cc      return
cc      end
c
      double precision function f00(x,testno,rt2pi)
      double precision x,rt2pi
      integer testno
c
      double precision s1, s2, mu(4)
      data mu(1), mu(2), mu(3), mu(4) / 0.2, 0.15, 0.1, 0.05 /,
     *   s1, s2 / 0.015, 0.045 /
c
      f00 = (exp(-((x-0.3)/s1)**2/2.0d0)/s1+
     *    2.0d0*exp(-((x-0.3-mu(testno))/s2)**2/2.0d0)/s2)/(3.0d0*rt2pi)
      return
      end
c
      double precision function k00(x,bandwd,rt2pi)
      double precision x,bandwd,rt2pi
c
      k00 = (exp(-(x/bandwd)**2/2.0d0)+
     *    exp(-((1-x)/bandwd)**2/2.0d0))/(rt2pi*bandwd)
      return
      end
      subroutine mkxys(bign,nobs,p,m,t,tpar,x,ldx,truey,
     *  sigma,ldsig,ftrig,ktrig,info)
      double precision t(nobs),tpar(4,*),x(ldx,*),truey(nobs),
     *  sigma(ldsig,*),ftrig(nobs),ktrig(nobs)
      integer bign,nobs,p,m,ldx,ldsig,info
 
c  Purpose: to create the design matrix, the observation vector and the
c	sigma matrix for the integral equation example.
c
c  On Entry:
c   bign		dimension of Hilbert space
c   nobs		number of observations
c   p			the number of parameters
c   m			order of the derivatives in the penalty
c   t(nobs)		locations for observations
c   tpar		values of the "true" parameters
c   ldx			leading dimension of x
c   ldsig 		leading dimension of sigma
c
c  On Exit:
c   x(nobs,p)		design matrix for dsnsm
c   truey(nobs)		"true" values for y (no noise added)
c   sigma(ldsig,p)	sigma matrix
c   ftrig(nobs)         trigonometric interpolant to f0
c   ktrig(nobs)         trigonometric interpolant to k0
c   info		error indicator. 0 indicates successful 
c			completion. Other values are:
c				1 - invalid testno
c				2 - bign <= 1 or nobs <= 1
c				3 - bign not even 
c  $Header: mkxys.f,v 2.100.1.1 86/10/07 13:30:50 lindstrom Exp $ 
      double precision twopi,sn,cs,tmp
      integer i,j,kb,r,flag
      
      twopi = 8.0d0 * atan(1.0d0)
      info = 1
      if (bign .le. 1 .or. nobs .le. 1) return
      info = 2
      if (int(bign/2) .ne. dble(bign)/2.0d0) return
      info = 3
      flag = 0
      if (p .eq. bign) then
	  flag = 1
      else
          if (int(p/2) .eq. dble(p)/2.0d0) return
      endif
      info = 0
      call dset(nobs,tpar(3,1),x(1,1),1)
      if (flag .eq. 1) then
	 kb = bign/2 + 1
	 r = bign/2 - 1
      else
         kb = (p - 1)/2 + 1
	 r = (p - 1)/2 
      endif
      do 20 i = 1, nobs
	  cs = cos( twopi*bign*t(i)/2.0d0)
	  truey(i) = tpar(1,1)*tpar(3,1)+ tpar(1,bign/2+1)*
     *	    tpar(3,bign/2+1)*cs*2.0d0
	  ftrig(i)= tpar(1,1) + tpar(1,bign/2+1)*cs/2.0d0
	  ktrig(i)= tpar(3,1) + tpar(3,bign/2+1)*cs/2.0d0
	  do 10 j = 1,bign/2-1
	      sn = sin(twopi*j*t(i))
	      cs = cos(twopi*j*t(i))
	      truey(i) = truey(i) + 2*((tpar(1,j+1)*tpar(3,j+1)-
     *	       tpar(2,j+1)*tpar(4,j+1))*cs+
     *	       (tpar(1,j+1)*tpar(4,j+1)+tpar(2,j+1)*tpar(3,j+1))*sn)
	      ftrig(i) = ftrig(i) + 2.0d0*tpar(1,j+1)*cs+
     *		    2.0d0*tpar(2,j+1)*sn
	      ktrig(i) = ktrig(i) + 2.0d0*tpar(3,j+1)*cs+
     *		    2.0d0*tpar(4,j+1)*sn
	      if (j .le. r) then
	          x(i,1+j) = 2.0d0*(tpar(3,j+1)*cs+tpar(4,j+1)*sn)
	          x(i,kb+j) = 2.0d0*(tpar(3,j+1)*sn-tpar(4,j+1)*cs)
	      endif
   10     continue
      if (flag .eq. 1) then
	  x(i,kb) = tpar(3,bign/2+1)*cos(twopi*bign*t(i)/2.0d0)/2.0d0
      endif
   20 continue
c			calculate sigma
      do 30 i = 1,p
	  call dset(p,0.0d0,sigma(1,i),1)
   30 continue
      do 40 i = 1,r
	  tmp = (twopi*i)**(2*m)
	  sigma(1+i,1+i) = tmp
	  sigma(kb+i,kb+i) = tmp
   40 continue
      if (flag .eq. 1) then
	  sigma(kb,kb) = ((twopi*dble(bign/2.0d0))**(2*m))/2.0d0
      endif
      return
      end
      subroutine prmut (x,npar,jpvt,job)
      integer npar,x(npar),jpvt(npar),job
c
c Purpose: permute the elements of the array x according to the index 
c	vector jpvt (either forward or backward permutation).
c
c On Entry:
c   x(npar)		array to be permuted
c   npar		size of x (and jpvt)
c   jpvt		indices of the permutation
c   job			indicator of forward or backward permutation
c			if job = 0 forward permutation  
c				x(jpvt(i)) moved to x(i)
c			if job is nonzero backward permutation 
c				x(i) moved to x(jpvt(i))
c On Exit:
c   x(npar)		array with permuted entries
c
c   Written:	Yin Ling	U. of Maryland, August,1978
c
c $Header: prmut.f,v 2.100.1.1 86/10/07 13:00:47 lindstrom Exp $
c
      integer i,j,k,t
c
      if (npar .le. 1) then
         return
      endif
      do 10 j = 1,npar
         jpvt(j) = -jpvt(j)
   10 continue
      if (job .eq. 0) then
c		forward permutation
         do 30 i = 1,npar 
            if (jpvt(i) .gt. 0) then
               goto 30
            endif
            j = i
            jpvt(j) = -jpvt(j)
            k = jpvt(j)
c           while
   20       if (jpvt(k) .lt. 0) then
               t = x(j)
               x(j) = x(k)
               x(k) = t
               jpvt(k) = -jpvt(k)
               j = k
               k = jpvt(k)
               goto 20
c           endwhile
            endif
   30    continue
      endif
      if (job .ne. 0 ) then
c			backward permutation
         do 50 i = 1,npar 
            if (jpvt(i) .gt. 0) then
               goto 50
            endif
            jpvt(i) = -jpvt(i)
            j = jpvt(i)
c           while
   40       if (j .ne. i) then
               t = x(i)
               x(i) = x(j)
               x(j) = t
               jpvt(j) = -jpvt(j)
               j = jpvt(j)
               goto 40
c           endwhile
            endif
   50    continue
      endif
      return
      end
c
c     these are three routines not in /lib/libcmalib.a and
c     ../../rech/ARBR/For/rkpk/lib/lib.a, they are
c         dchdc.f, dsvdc.f, dtrco.f 
c     all from linpack of ftp site netlib.att.com
c
c
      subroutine dchdc(a,lda,p,work,jpvt,job,info)
      integer lda,p,jpvt(1),job,info
      double precision a(lda,1),work(1)
c
c     dchdc computes the cholesky decomposition of a positive definite
c     matrix.  a pivoting option allows the user to estimate the
c     condition of a positive definite matrix or determine the rank
c     of a positive semidefinite matrix.
c
c     on entry
c
c         a      double precision(lda,p).
c                a contains the matrix whose decomposition is to
c                be computed.  onlt the upper half of a need be stored.
c                the lower part of the array a is not referenced.
c
c         lda    integer.
c                lda is the leading dimension of the array a.
c
c         p      integer.
c                p is the order of the matrix.
c
c         work   double precision.
c                work is a work array.
c
c         jpvt   integer(p).
c                jpvt contains integers that control the selection
c                of the pivot elements, if pivoting has been requested.
c                each diagonal element a(k,k)
c                is placed in one of three classes according to the
c                value of jpvt(k).
c
c                   if jpvt(k) .gt. 0, then x(k) is an initial
c                                      element.
c
c                   if jpvt(k) .eq. 0, then x(k) is a free element.
c
c                   if jpvt(k) .lt. 0, then x(k) is a final element.
c
c                before the decomposition is computed, initial elements
c                are moved by symmetric row and column interchanges to
c                the beginning of the array a and final
c                elements to the end.  both initial and final elements
c                are frozen in place during the computation and only
c                free elements are moved.  at the k-th stage of the
c                reduction, if a(k,k) is occupied by a free element
c                it is interchanged with the largest free element
c                a(l,l) with l .ge. k.  jpvt is not referenced if
c                job .eq. 0.
c
c        job     integer.
c                job is an integer that initiates column pivoting.
c                if job .eq. 0, no pivoting is done.
c                if job .ne. 0, pivoting is done.
c
c     on return
c
c         a      a contains in its upper half the cholesky factor
c                of the matrix a as it has been permuted by pivoting.
c
c         jpvt   jpvt(j) contains the index of the diagonal element
c                of a that was moved into the j-th position,
c                provided pivoting was requested.
c
c         info   contains the index of the last positive diagonal
c                element of the cholesky factor.
c
c     for positive definite matrices info = p is the normal return.
c     for pivoting with positive semidefinite matrices info will
c     in general be less than p.  however, info may be greater than
c     the rank of a, since rounding error can cause an otherwise zero
c     element to be positive. indefinite systems will always cause
c     info to be less than p.
c
c     linpack. this version dated 08/14/78 .
c     j.j. dongarra and g.w. stewart, argonne national laboratory and
c     university of maryland.
c
c
c     blas daxpy,dswap
c     fortran dsqrt
c
c     internal variables
c
      integer pu,pl,plp1,i,j,jp,jt,k,kb,km1,kp1,l,maxl
      double precision temp
      double precision maxdia
      logical swapk,negk
c
      pl = 1
      pu = 0
      info = p
      if (job .eq. 0) go to 160
c
c        pivoting has been requested. rearrange the
c        the elements according to jpvt.
c
         do 70 k = 1, p
            swapk = jpvt(k) .gt. 0
            negk = jpvt(k) .lt. 0
            jpvt(k) = k
            if (negk) jpvt(k) = -jpvt(k)
            if (.not.swapk) go to 60
               if (k .eq. pl) go to 50
                  call dswap(pl-1,a(1,k),1,a(1,pl),1)
                  temp = a(k,k)
                  a(k,k) = a(pl,pl)
                  a(pl,pl) = temp
                  plp1 = pl + 1
                  if (p .lt. plp1) go to 40
                  do 30 j = plp1, p
                     if (j .ge. k) go to 10
                        temp = a(pl,j)
                        a(pl,j) = a(j,k)
                        a(j,k) = temp
                     go to 20
   10                continue
                     if (j .eq. k) go to 20
                        temp = a(k,j)
                        a(k,j) = a(pl,j)
                        a(pl,j) = temp
   20                continue
   30             continue
   40             continue
                  jpvt(k) = jpvt(pl)
                  jpvt(pl) = k
   50          continue
               pl = pl + 1
   60       continue
   70    continue
         pu = p
         if (p .lt. pl) go to 150
         do 140 kb = pl, p
            k = p - kb + pl
            if (jpvt(k) .ge. 0) go to 130
               jpvt(k) = -jpvt(k)
               if (pu .eq. k) go to 120
                  call dswap(k-1,a(1,k),1,a(1,pu),1)
                  temp = a(k,k)
                  a(k,k) = a(pu,pu)
                  a(pu,pu) = temp
                  kp1 = k + 1
                  if (p .lt. kp1) go to 110
                  do 100 j = kp1, p
                     if (j .ge. pu) go to 80
                        temp = a(k,j)
                        a(k,j) = a(j,pu)
                        a(j,pu) = temp
                     go to 90
   80                continue
                     if (j .eq. pu) go to 90
                        temp = a(k,j)
                        a(k,j) = a(pu,j)
                        a(pu,j) = temp
   90                continue
  100             continue
  110             continue
                  jt = jpvt(k)
                  jpvt(k) = jpvt(pu)
                  jpvt(pu) = jt
  120          continue
               pu = pu - 1
  130       continue
  140    continue
  150    continue
  160 continue
      do 270 k = 1, p
c
c        reduction loop.
c
         maxdia = a(k,k)
         kp1 = k + 1
         maxl = k
c
c        determine the pivot element.
c
         if (k .lt. pl .or. k .ge. pu) go to 190
            do 180 l = kp1, pu
               if (a(l,l) .le. maxdia) go to 170
                  maxdia = a(l,l)
                  maxl = l
  170          continue
  180       continue
  190    continue
c
c        quit if the pivot element is not positive.
c
         if (maxdia .gt. 0.0d0) go to 200
            info = k - 1
c     ......exit
            go to 280
  200    continue
         if (k .eq. maxl) go to 210
c
c           start the pivoting and update jpvt.
c
            km1 = k - 1
            call dswap(km1,a(1,k),1,a(1,maxl),1)
            a(maxl,maxl) = a(k,k)
            a(k,k) = maxdia
            jp = jpvt(maxl)
            jpvt(maxl) = jpvt(k)
            jpvt(k) = jp
  210    continue
c
c        reduction step. pivoting is contained across the rows.
c
         work(k) = dsqrt(a(k,k))
         a(k,k) = work(k)
         if (p .lt. kp1) go to 260
         do 250 j = kp1, p
            if (k .eq. maxl) go to 240
               if (j .ge. maxl) go to 220
                  temp = a(k,j)
                  a(k,j) = a(j,maxl)
                  a(j,maxl) = temp
               go to 230
  220          continue
               if (j .eq. maxl) go to 230
                  temp = a(k,j)
                  a(k,j) = a(maxl,j)
                  a(maxl,j) = temp
  230          continue
  240       continue
            a(k,j) = a(k,j)/work(k)
            work(j) = a(k,j)
            temp = -a(k,j)
            call daxpy(j-k,temp,work(kp1),1,a(kp1,j),1)
  250    continue
  260    continue
  270 continue
  280 continue
      return
      end
      subroutine dsvdc(x,ldx,n,p,s,e,u,ldu,v,ldv,work,job,info)
      integer ldx,n,p,ldu,ldv,job,info
      double precision x(ldx,1),s(1),e(1),u(ldu,1),v(ldv,1),work(1)
c
c
c     dsvdc is a subroutine to reduce a double precision nxp matrix x
c     by orthogonal transformations u and v to diagonal form.  the
c     diagonal elements s(i) are the singular values of x.  the
c     columns of u are the corresponding left singular vectors,
c     and the columns of v the right singular vectors.
c
c     on entry
c
c         x         double precision(ldx,p), where ldx.ge.n.
c                   x contains the matrix whose singular value
c                   decomposition is to be computed.  x is
c                   destroyed by dsvdc.
c
c         ldx       integer.
c                   ldx is the leading dimension of the array x.
c
c         n         integer.
c                   n is the number of rows of the matrix x.
c
c         p         integer.
c                   p is the number of columns of the matrix x.
c
c         ldu       integer.
c                   ldu is the leading dimension of the array u.
c                   (see below).
c
c         ldv       integer.
c                   ldv is the leading dimension of the array v.
c                   (see below).
c
c         work      double precision(n).
c                   work is a scratch array.
c
c         job       integer.
c                   job controls the computation of the singular
c                   vectors.  it has the decimal expansion ab
c                   with the following meaning
c
c                        a.eq.0    do not compute the left singular
c                                  vectors.
c                        a.eq.1    return the n left singular vectors
c                                  in u.
c                        a.ge.2    return the first min(n,p) singular
c                                  vectors in u.
c                        b.eq.0    do not compute the right singular
c                                  vectors.
c                        b.eq.1    return the right singular vectors
c                                  in v.
c
c     on return
c
c         s         double precision(mm), where mm=min(n+1,p).
c                   the first min(n,p) entries of s contain the
c                   singular values of x arranged in descending
c                   order of magnitude.
c
c         e         double precision(p), 
c                   e ordinarily contains zeros.  however see the
c                   discussion of info for exceptions.
c
c         u         double precision(ldu,k), where ldu.ge.n.  if
c                                   joba.eq.1 then k.eq.n, if joba.ge.2
c                                   then k.eq.min(n,p).
c                   u contains the matrix of left singular vectors.
c                   u is not referenced if joba.eq.0.  if n.le.p
c                   or if joba.eq.2, then u may be identified with x
c                   in the subroutine call.
c
c         v         double precision(ldv,p), where ldv.ge.p.
c                   v contains the matrix of right singular vectors.
c                   v is not referenced if job.eq.0.  if p.le.n,
c                   then v may be identified with x in the
c                   subroutine call.
c
c         info      integer.
c                   the singular values (and their corresponding
c                   singular vectors) s(info+1),s(info+2),...,s(m)
c                   are correct (here m=min(n,p)).  thus if
c                   info.eq.0, all the singular values and their
c                   vectors are correct.  in any event, the matrix
c                   b = trans(u)*x*v is the bidiagonal matrix
c                   with the elements of s on its diagonal and the
c                   elements of e on its super-diagonal (trans(u)
c                   is the transpose of u).  thus the singular
c                   values of x and b are the same.
c
c     linpack. this version dated 08/14/78 .
c              correction made to shift 2/84.
c     g.w. stewart, university of maryland, argonne national lab.
c
c     dsvdc uses the following functions and subprograms.
c
c     external drot
c     blas daxpy,ddot,dscal,dswap,dnrm2,drotg
c     fortran dabs,dmax1,max0,min0,mod,dsqrt
c
c     internal variables
c
      integer i,iter,j,jobu,k,kase,kk,l,ll,lls,lm1,lp1,ls,lu,m,maxit,
     *        mm,mm1,mp1,nct,nctp1,ncu,nrt,nrtp1
      double precision ddot,t,r
      double precision b,c,cs,el,emm1,f,g,dnrm2,scale,shift,sl,sm,sn,
     *                 smm1,t1,test,ztest
      logical wantu,wantv
c
c
c     set the maximum number of iterations.
c
      maxit = 30
c
c     determine what is to be computed.
c
      wantu = .false.
      wantv = .false.
      jobu = mod(job,100)/10
      ncu = n
      if (jobu .gt. 1) ncu = min0(n,p)
      if (jobu .ne. 0) wantu = .true.
      if (mod(job,10) .ne. 0) wantv = .true.
c
c     reduce x to bidiagonal form, storing the diagonal elements
c     in s and the super-diagonal elements in e.
c
      info = 0
      nct = min0(n-1,p)
      nrt = max0(0,min0(p-2,n))
      lu = max0(nct,nrt)
      if (lu .lt. 1) go to 170
      do 160 l = 1, lu
         lp1 = l + 1
         if (l .gt. nct) go to 20
c
c           compute the transformation for the l-th column and
c           place the l-th diagonal in s(l).
c
            s(l) = dnrm2(n-l+1,x(l,l),1)
            if (s(l) .eq. 0.0d0) go to 10
               if (x(l,l) .ne. 0.0d0) s(l) = dsign(s(l),x(l,l))
               call dscal(n-l+1,1.0d0/s(l),x(l,l),1)
               x(l,l) = 1.0d0 + x(l,l)
   10       continue
            s(l) = -s(l)
   20    continue
         if (p .lt. lp1) go to 50
         do 40 j = lp1, p
            if (l .gt. nct) go to 30
            if (s(l) .eq. 0.0d0) go to 30
c
c              apply the transformation.
c
               t = -ddot(n-l+1,x(l,l),1,x(l,j),1)/x(l,l)
               call daxpy(n-l+1,t,x(l,l),1,x(l,j),1)
   30       continue
c
c           place the l-th row of x into  e for the
c           subsequent calculation of the row transformation.
c
            e(j) = x(l,j)
   40    continue
   50    continue
         if (.not.wantu .or. l .gt. nct) go to 70
c
c           place the transformation in u for subsequent back
c           multiplication.
c
            do 60 i = l, n
               u(i,l) = x(i,l)
   60       continue
   70    continue
         if (l .gt. nrt) go to 150
c
c           compute the l-th row transformation and place the
c           l-th super-diagonal in e(l).
c
            e(l) = dnrm2(p-l,e(lp1),1)
            if (e(l) .eq. 0.0d0) go to 80
               if (e(lp1) .ne. 0.0d0) e(l) = dsign(e(l),e(lp1))
               call dscal(p-l,1.0d0/e(l),e(lp1),1)
               e(lp1) = 1.0d0 + e(lp1)
   80       continue
            e(l) = -e(l)
            if (lp1 .gt. n .or. e(l) .eq. 0.0d0) go to 120
c
c              apply the transformation.
c
               do 90 i = lp1, n
                  work(i) = 0.0d0
   90          continue
               do 100 j = lp1, p
                  call daxpy(n-l,e(j),x(lp1,j),1,work(lp1),1)
  100          continue
               do 110 j = lp1, p
                  call daxpy(n-l,-e(j)/e(lp1),work(lp1),1,x(lp1,j),1)
  110          continue
  120       continue
            if (.not.wantv) go to 140
c
c              place the transformation in v for subsequent
c              back multiplication.
c
               do 130 i = lp1, p
                  v(i,l) = e(i)
  130          continue
  140       continue
  150    continue
  160 continue
  170 continue
c
c     set up the final bidiagonal matrix or order m.
c
      m = min0(p,n+1)
      nctp1 = nct + 1
      nrtp1 = nrt + 1
      if (nct .lt. p) s(nctp1) = x(nctp1,nctp1)
      if (n .lt. m) s(m) = 0.0d0
      if (nrtp1 .lt. m) e(nrtp1) = x(nrtp1,m)
      e(m) = 0.0d0
c
c     if required, generate u.
c
      if (.not.wantu) go to 300
         if (ncu .lt. nctp1) go to 200
         do 190 j = nctp1, ncu
            do 180 i = 1, n
               u(i,j) = 0.0d0
  180       continue
            u(j,j) = 1.0d0
  190    continue
  200    continue
         if (nct .lt. 1) go to 290
         do 280 ll = 1, nct
            l = nct - ll + 1
            if (s(l) .eq. 0.0d0) go to 250
               lp1 = l + 1
               if (ncu .lt. lp1) go to 220
               do 210 j = lp1, ncu
                  t = -ddot(n-l+1,u(l,l),1,u(l,j),1)/u(l,l)
                  call daxpy(n-l+1,t,u(l,l),1,u(l,j),1)
  210          continue
  220          continue
               call dscal(n-l+1,-1.0d0,u(l,l),1)
               u(l,l) = 1.0d0 + u(l,l)
               lm1 = l - 1
               if (lm1 .lt. 1) go to 240
               do 230 i = 1, lm1
                  u(i,l) = 0.0d0
  230          continue
  240          continue
            go to 270
  250       continue
               do 260 i = 1, n
                  u(i,l) = 0.0d0
  260          continue
               u(l,l) = 1.0d0
  270       continue
  280    continue
  290    continue
  300 continue
c
c     if it is required, generate v.
c
      if (.not.wantv) go to 350
         do 340 ll = 1, p
            l = p - ll + 1
            lp1 = l + 1
            if (l .gt. nrt) go to 320
            if (e(l) .eq. 0.0d0) go to 320
               do 310 j = lp1, p
                  t = -ddot(p-l,v(lp1,l),1,v(lp1,j),1)/v(lp1,l)
                  call daxpy(p-l,t,v(lp1,l),1,v(lp1,j),1)
  310          continue
  320       continue
            do 330 i = 1, p
               v(i,l) = 0.0d0
  330       continue
            v(l,l) = 1.0d0
  340    continue
  350 continue
c
c     main iteration loop for the singular values.
c
      mm = m
      iter = 0
  360 continue
c
c        quit if all the singular values have been found.
c
c     ...exit
         if (m .eq. 0) go to 620
c
c        if too many iterations have been performed, set
c        flag and return.
c
         if (iter .lt. maxit) go to 370
            info = m
c     ......exit
            go to 620
  370    continue
c
c        this section of the program inspects for
c        negligible elements in the s and e arrays.  on
c        completion the variables kase and l are set as follows.
c
c           kase = 1     if s(m) and e(l-1) are negligible and l.lt.m
c           kase = 2     if s(l) is negligible and l.lt.m
c           kase = 3     if e(l-1) is negligible, l.lt.m, and
c                        s(l), ..., s(m) are not negligible (qr step).
c           kase = 4     if e(m-1) is negligible (convergence).
c
         do 390 ll = 1, m
            l = m - ll
c        ...exit
            if (l .eq. 0) go to 400
            test = dabs(s(l)) + dabs(s(l+1))
            ztest = test + dabs(e(l))
            if (ztest .ne. test) go to 380
               e(l) = 0.0d0
c        ......exit
               go to 400
  380       continue
  390    continue
  400    continue
         if (l .ne. m - 1) go to 410
            kase = 4
         go to 480
  410    continue
            lp1 = l + 1
            mp1 = m + 1
            do 430 lls = lp1, mp1
               ls = m - lls + lp1
c           ...exit
               if (ls .eq. l) go to 440
               test = 0.0d0
               if (ls .ne. m) test = test + dabs(e(ls))
               if (ls .ne. l + 1) test = test + dabs(e(ls-1))
               ztest = test + dabs(s(ls))
               if (ztest .ne. test) go to 420
                  s(ls) = 0.0d0
c           ......exit
                  go to 440
  420          continue
  430       continue
  440       continue
            if (ls .ne. l) go to 450
               kase = 3
            go to 470
  450       continue
            if (ls .ne. m) go to 460
               kase = 1
            go to 470
  460       continue
               kase = 2
               l = ls
  470       continue
  480    continue
         l = l + 1
c
c        perform the task indicated by kase.
c
         go to (490,520,540,570), kase
c
c        deflate negligible s(m).
c
  490    continue
            mm1 = m - 1
            f = e(m-1)
            e(m-1) = 0.0d0
            do 510 kk = l, mm1
               k = mm1 - kk + l
               t1 = s(k)
               call drotg(t1,f,cs,sn)
               s(k) = t1
               if (k .eq. l) go to 500
                  f = -sn*e(k-1)
                  e(k-1) = cs*e(k-1)
  500          continue
               if (wantv) call drot(p,v(1,k),1,v(1,m),1,cs,sn)
  510       continue
         go to 610
c
c        split at negligible s(l).
c
  520    continue
            f = e(l-1)
            e(l-1) = 0.0d0
            do 530 k = l, m
               t1 = s(k)
               call drotg(t1,f,cs,sn)
               s(k) = t1
               f = -sn*e(k)
               e(k) = cs*e(k)
               if (wantu) call drot(n,u(1,k),1,u(1,l-1),1,cs,sn)
  530       continue
         go to 610
c
c        perform one qr step.
c
  540    continue
c
c           calculate the shift.
c
            scale = dmax1(dabs(s(m)),dabs(s(m-1)),dabs(e(m-1)),
     *                    dabs(s(l)),dabs(e(l)))
            sm = s(m)/scale
            smm1 = s(m-1)/scale
            emm1 = e(m-1)/scale
            sl = s(l)/scale
            el = e(l)/scale
            b = ((smm1 + sm)*(smm1 - sm) + emm1**2)/2.0d0
            c = (sm*emm1)**2
            shift = 0.0d0
            if (b .eq. 0.0d0 .and. c .eq. 0.0d0) go to 550
               shift = dsqrt(b**2+c)
               if (b .lt. 0.0d0) shift = -shift
               shift = c/(b + shift)
  550       continue
            f = (sl + sm)*(sl - sm) + shift
            g = sl*el
c
c           chase zeros.
c
            mm1 = m - 1
            do 560 k = l, mm1
               call drotg(f,g,cs,sn)
               if (k .ne. l) e(k-1) = f
               f = cs*s(k) + sn*e(k)
               e(k) = cs*e(k) - sn*s(k)
               g = sn*s(k+1)
               s(k+1) = cs*s(k+1)
               if (wantv) call drot(p,v(1,k),1,v(1,k+1),1,cs,sn)
               call drotg(f,g,cs,sn)
               s(k) = f
               f = cs*e(k) + sn*s(k+1)
               s(k+1) = -sn*e(k) + cs*s(k+1)
               g = sn*e(k+1)
               e(k+1) = cs*e(k+1)
               if (wantu .and. k .lt. n)
     *            call drot(n,u(1,k),1,u(1,k+1),1,cs,sn)
  560       continue
            e(m-1) = f
            iter = iter + 1
         go to 610
c
c        convergence.
c
  570    continue
c
c           make the singular value  positive.
c
            if (s(l) .ge. 0.0d0) go to 580
               s(l) = -s(l)
               if (wantv) call dscal(p,-1.0d0,v(1,l),1)
  580       continue
c
c           order the singular value.
c
  590       if (l .eq. mm) go to 600
c           ...exit
               if (s(l) .ge. s(l+1)) go to 600
               t = s(l)
               s(l) = s(l+1)
               s(l+1) = t
               if (wantv .and. l .lt. p)
     *            call dswap(p,v(1,l),1,v(1,l+1),1)
               if (wantu .and. l .lt. n)
     *            call dswap(n,u(1,l),1,u(1,l+1),1)
               l = l + 1
            go to 590
  600       continue
            iter = 0
            m = m - 1
  610    continue
      go to 360
  620 continue
      return
      end
      subroutine dtrco(t,ldt,n,rcond,z,job)
      integer ldt,n,job
      double precision t(ldt,1),z(1)
      double precision rcond
c
c     dtrco estimates the condition of a double precision triangular
c     matrix.
c
c     on entry
c
c        t       double precision(ldt,n)
c                t contains the triangular matrix. the zero
c                elements of the matrix are not referenced, and
c                the corresponding elements of the array can be
c                used to store other information.
c
c        ldt     integer
c                ldt is the leading dimension of the array t.
c
c        n       integer
c                n is the order of the system.
c
c        job     integer
c                = 0         t  is lower triangular.
c                = nonzero   t  is upper triangular.
c
c     on return
c
c        rcond   double precision
c                an estimate of the reciprocal condition of  t .
c                for the system  t*x = b , relative perturbations
c                in  t  and  b  of size  epsilon  may cause
c                relative perturbations in  x  of size  epsilon/rcond .
c                if  rcond  is so small that the logical expression
c                           1.0 + rcond .eq. 1.0
c                is true, then  t  may be singular to working
c                precision.  in particular,  rcond  is zero  if
c                exact singularity is detected or the estimate
c                underflows.
c
c        z       double precision(n)
c                a work vector whose contents are usually unimportant.
c                if  t  is close to a singular matrix, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,dasum
c     fortran dabs,dmax1,dsign
c
c     internal variables
c
      double precision w,wk,wkm,ek
      double precision tnorm,ynorm,s,sm,dasum
      integer i1,j,j1,j2,k,kk,l
      logical lower
c
      lower = job .eq. 0
c
c     compute 1-norm of t
c
      tnorm = 0.0d0
      do 10 j = 1, n
         l = j
         if (lower) l = n + 1 - j
         i1 = 1
         if (lower) i1 = j
         tnorm = dmax1(tnorm,dasum(l,t(i1,j),1))
   10 continue
c
c     rcond = 1/(norm(t)*(estimate of norm(inverse(t)))) .
c     estimate = norm(z)/norm(y) where  t*z = y  and  trans(t)*y = e .
c     trans(t)  is the transpose of t .
c     the components of  e  are chosen to cause maximum local
c     growth in the elements of y .
c     the vectors are frequently rescaled to avoid overflow.
c
c     solve trans(t)*y = e
c
      ek = 1.0d0
      do 20 j = 1, n
         z(j) = 0.0d0
   20 continue
      do 100 kk = 1, n
         k = kk
         if (lower) k = n + 1 - kk
         if (z(k) .ne. 0.0d0) ek = dsign(ek,-z(k))
         if (dabs(ek-z(k)) .le. dabs(t(k,k))) go to 30
            s = dabs(t(k,k))/dabs(ek-z(k))
            call dscal(n,s,z,1)
            ek = s*ek
   30    continue
         wk = ek - z(k)
         wkm = -ek - z(k)
         s = dabs(wk)
         sm = dabs(wkm)
         if (t(k,k) .eq. 0.0d0) go to 40
            wk = wk/t(k,k)
            wkm = wkm/t(k,k)
         go to 50
   40    continue
            wk = 1.0d0
            wkm = 1.0d0
   50    continue
         if (kk .eq. n) go to 90
            j1 = k + 1
            if (lower) j1 = 1
            j2 = n
            if (lower) j2 = k - 1
            do 60 j = j1, j2
               sm = sm + dabs(z(j)+wkm*t(k,j))
               z(j) = z(j) + wk*t(k,j)
               s = s + dabs(z(j))
   60       continue
            if (s .ge. sm) go to 80
               w = wkm - wk
               wk = wkm
               do 70 j = j1, j2
                  z(j) = z(j) + w*t(k,j)
   70          continue
   80       continue
   90    continue
         z(k) = wk
  100 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
c
      ynorm = 1.0d0
c
c     solve t*z = y
c
      do 130 kk = 1, n
         k = n + 1 - kk
         if (lower) k = kk
         if (dabs(z(k)) .le. dabs(t(k,k))) go to 110
            s = dabs(t(k,k))/dabs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
  110    continue
         if (t(k,k) .ne. 0.0d0) z(k) = z(k)/t(k,k)
         if (t(k,k) .eq. 0.0d0) z(k) = 1.0d0
         i1 = 1
         if (lower) i1 = k + 1
         if (kk .ge. n) go to 120
            w = -z(k)
            call daxpy(n-kk,w,t(i1,k),1,z(i1),1)
  120    continue
  130 continue
c     make znorm = 1.0
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
c
      if (tnorm .ne. 0.0d0) rcond = ynorm/tnorm
      if (tnorm .eq. 0.0d0) rcond = 0.0d0
      return
      end
c
c     these are more supplemenatry routines,
c         drot.f, drotg.f
c     they are from blas
c

      subroutine  drot (n,dx,incx,dy,incy,c,s)
c
c     applies a plane rotation.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1),dtemp,c,s
      integer i,incx,incy,ix,iy,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = c*dx(ix) + s*dy(iy)
        dy(iy) = c*dy(iy) - s*dx(ix)
        dx(ix) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
   20 do 30 i = 1,n
        dtemp = c*dx(i) + s*dy(i)
        dy(i) = c*dy(i) - s*dx(i)
        dx(i) = dtemp
   30 continue
      return
      end
      subroutine drotg(da,db,c,s)
c
c     construct givens plane rotation.
c     jack dongarra, linpack, 3/11/78.
c                    modified 9/27/86.
c
      double precision da,db,c,s,roe,scale,r,z
c
      roe = db
      if( dabs(da) .gt. dabs(db) ) roe = da
      scale = dabs(da) + dabs(db)
      if( scale .ne. 0.0d0 ) go to 10
         c = 1.0d0
         s = 0.0d0
         r = 0.0d0
         go to 20
   10 r = scale*dsqrt((da/scale)**2 + (db/scale)**2)
      r = dsign(1.0d0,roe)*r
      c = da/r
      s = db/r
   20 z = s
      if( dabs(c) .gt. 0.0d0 .and. dabs(c) .le. s ) z = 1.0d0/c
      da = r
      db = z
      return
      end

c libcost2.f Updated Feb 23, 2009 by Junqing Wu, UCSB
c ----------------------------------------------------------------------
      SUBROUTINE LIBCOST2(NOBS, Y, BAS, NNULL, NBAS, MAXBAS,  
     & IDFHAT, nrep, GDFHAT, TAUHAT, ISB, RSS,IWK, YY, BASWK, WKA,
     & WKALL,WKAMAX,ICB,PTB,U,H,FITS,FITSWK,Z,ZK,YYPTB)

      INTEGER NOBS, NNULL, NBAS, MAXBAS, ISB(*),  
     & IWK(*), nrep, WKAMAX

      DOUBLE PRECISION Y(*), BAS(NOBS,*), IDFHAT(MAXBAS),  
     & RSS(*), YY(*), YYPTB(NOBS,*), BASWK(NOBS,*), WKA(NOBS,*), 
     & GDFHAT(MAXBAS), TAUHAT,tmp1,tmp2,
     & dbar(NOBS), YWK(NOBS), WKALL(NOBS,WKAMAX,nrep),
     & d(NOBS,nrep),ICB(MAXBAS,nrep),DISB(MAXBAS), FITSWK(NOBS,MAXBAS),
     & PTB(NOBS,nrep),U(*),H(*),FITS(NOBS,MAXBAS,nrep),Z(*),ZK(*)

      REAL RNOR

c########################################################################
c     NOBS:   number of observations
c     Y:      response.
c     BAS:    NOBSxNBAS, Bases. The first NNULL columns are bases of the 
c             null space
c     NNULL:  Number of basis functions already in the model
c     NBAS:   number of bases
c     MAXBAS: the maximum number of basis functions to be chosen
c     IDFHAT: The estimated IDF values
c     NREP:   Number of repetitions of fitting perturbed responses 
c     GDFHAT: Total Cost for all chosen basis functions at each step
c     TAUHAT: .5*(est of sigma)  
c     ISB:    Index of Selected Bases, (MAXBAS), the column index of
c             chosen bases in BAS. It always include the first NNULL
c             columns.
c     RSS:    residual sum of squares for each model
c     IWK:    working vector, (NBAS), each element indicates if a basis
c		has been chosen, 1-yes, 0-no
c     YY:     working vector for storing perturbed responses
c     BASWK:  working matrix (NOBSxNBAS) for basis functions
c     WKA:    working matrix (NOBSxNBAS) for "working vectors"
c     WKALL:  working matrix (NOBSxNBASxNREP) for all "working vectors"
c             for all perturbed responses.
c     WKAMAX: must be >MAXBAS, it ensures enough memory space for 
c		selecting up to MAXBAS basis functions.
c     ICB: Indices of Chosen Bases, similar to ISB, except for perturbed
c	    responses.
c     PTB: perturbations generated inside Fortran.
c     U: The norms of the original vectors from which the "working
c	  vectors" are derived with the opposite sign of the first
c	  element.
c     H: The first elements of the vectors corresponding to the
c        "working vectors" obtained from the transformed y's 
c     FITS:   array of fits, (NOBSx(MAXBAS-NNULL+1)xNREP)
c     FITSWK: working matrix to store fits temporarily
c     Z: The first k elements are obtained from y after k 
c        transformations, the rest are zeroes.
c     ZK: The first K elements are obtained from the (k+1)th chosen
c         basis after k transformations, the rest are zeroes.
c     YYPTB: perturbed responses after being transformed by the
c            "working vectors".
c Modified by Junqing Wu, UCSB on Sep 28, 2009
c########################################################################

      INTEGER I, J, K, MBAS, NB, IDXTMP(NBAS)
      DOUBLE PRECISION TMP
 
      MBAS = MAXBAS
	  DO 5 K = 1,NBAS
	    IDXTMP(K)=IWK(K)
5	  CONTINUE
c	   TMP = RSTART(123)
         DO 20 J = 1,nrep
           DO 10 i=1, NOBS
              d(i,j) = TAUHAT*dble(rnor())
              ywk(i) = dble(Y(i)) + dble(d(i,j))
 10        CONTINUE

        MAXBAS = MBAS
        CALL DCOPY(NOBS,YWK(1),1,YY(1),1)
        CALL DCOPY(NOBS*NBAS, BAS, 1, BASWK, 1)
        
        DO 12 I= 1,MAXBAS
        CALL DSET(NOBS,0.D0,FITSWK(1,I),1)
12      CONTINUE

          CALL PRESELE2(BASWK,NOBS,NBAS,YY,NNULL,MAXBAS,ISB,
     &             RSS,IWK,WKA,U,H,FITSWK,Z,ZK,BAS)
	    CALL DCOPY(NOBS,YY,1,YYPTB(1,J),1)
	    CALL DCOPY(NOBS*MAXBAS,WKA(1,1),1,WKALL(1,1,J),1)
          CALL DCOPY(NOBS*MAXBAS,FITSWK(1,1),1,FITS(1,1,J),1)

          DO 15 i=1,MAXBAS
            DISB(i)=dble(ISB(i))
 15   CONTINUE

          CALL DCOPY(MAXBAS,DISB,1,ICB(1,J),1)
	  DO 17 K = 1,NBAS
	    IWK(K)=IDXTMP(K)
 17   CONTINUE
 20   CONTINUE

          CALL DCOPY(NOBS*nrep,d,1,PTB,1)
          CALL DCOPY(NOBS,Y(1),1,YY(1),1)
          CALL DCOPY(NOBS*NBAS, BAS, 1, BASWK, 1)
          CALL PRESELE2(BASWK,NOBS,NBAS,YY,NNULL,MAXBAS,ISB,
     &             RSS,IWK,WKA,U,H,FITSWK,Z,ZK,BAS)



	DO 25 I=1, NOBS
	DBAR(I) = 0.D0
	  DO 23 J=1, nrep
	  DBAR(I) = DBAR(I) + D(I,J)
23	  CONTINUE
	DBAR(I) = DBAR(I)/DBLE(nrep)
25	CONTINUE

	DO 28 I=1, NNULL
	  GDFHAT(I) = DBLE(I)
	  IDFHAT(I) = 1.D0
28	CONTINUE

	DO 50 K=NNULL+1, MAXBAS
	  GDFHAT(K) = 0.D0
        IDFHAT(K) = 0.D0
	  DO 40 I=1, NOBS
	tmp1 = 0.d0
	tmp2 = 0.d0
	    DO 30 J=1, nrep
	    tmp1 = tmp1 + (D(I,J)-DBAR(I))*FITS(I,K,J)
	    tmp2 = tmp2 + (D(I,J)-DBAR(I))**2
30	    CONTINUE
	  GDFHAT(K) = GDFHAT(K) + tmp1/tmp2
40	  CONTINUE
	IDFHAT(K) = (GDFHAT(K)-NNULL)/(MAXBAS-NNULL)
	IF (IDFHAT(K).LE.0.D0) THEN
		IDFHAT(K) = 0.D0
	ENDIF
50	CONTINUE
      RETURN
      END



c suppl.f, (Aug 10 1995, 17:30)
c some supplementary programs used in HAS: (collected on July 9, 1995) 
c **********************************************************************
c     MISC      dset.f iset.f idamin.f                           
c     BLAS      daxpy.f dcopy.f ddot.f*  dnrm2.f* dscal.f dswap.f
c     LINPACK   dqrdc.f dqrsl.f dasum.f* dtrsl.f 
c
c Note: those with * are functions.
c **********************************************************************
c
c ----------------------------------------------------------------------
c misc.f
c
c     dset, iset, idamin
c
c **********************************************************************
      subroutine  dset(n,da,dx,incx)
      integer n,incx
      double precision da,dx(*)
c
c Purpose : set vector dx to constant da. Unrolled loops are used for 
c       increment equal to one.
c
c On Entry:
c   n                   length of dx
c   da                  any constant
c   incx                increment for dx
c
c On Exit:
c   dx(n)               vector with all n entries set to da
c
c $Header: dset.f,v 2.1 86/04/08 14:06:25 lindstrom Exp $
c
      integer i,m,mp1,nincx
c
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da
        dx(i + 1) = da
        dx(i + 2) = da
        dx(i + 3) = da
        dx(i + 4) = da
   50 continue
      return
      end
c
      subroutine  iset(n,ia,ix,incx)
      integer n,incx,ia,ix(*)
c (modified from dset of gcvpack)
c Purpose : set vector ix to constant ia. Unrolled loops are used for 
c        increment equal to one.
c
c On Entry:
c     n                        length of ix
c   ia                        any constant
c   incx                increment for ix
c
c On Exit:
c   ix(n)                vector with all n entries set to ia
c
c $Header: iset.f,v 2.1 86/04/08 14:06:25 lindstrom Exp $
c
      integer i,m,mp1,nincx
c
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        ix(i) = ia
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        ix(i) = ia
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        ix(i) = ia
        ix(i + 1) = ia
        ix(i + 2) = ia
        ix(i + 3) = ia
        ix(i + 4) = ia
   50 continue
      return
      end
C **********************************************************************
      INTEGER FUNCTION IDAMIN(N,DX,INCX)
C  modified further from IDAMAX in Zheng Lou's program by taking out DABS 
C     FINDS THE INDEX OF ELEMENT HAVING MIN. VALUE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DMIN
      INTEGER I,INCX,IX,N
C
      IDAMIN = 0
      IF( N .LT. 1 ) RETURN
      IDAMIN = 1
      IF(N.EQ.1)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      DMIN = DX(1)
      IX = IX + INCX
      DO 10 I = 2,N
         IF(DX(IX).GE.DMIN) GO TO 5
         IDAMIN = I
         DMIN = DX(IX)
    5    IX = IX + INCX
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 DMIN = DX(1)
      DO 30 I = 2,N
         IF(DX(I).GE.DMIN) GO TO 30
         IDAMIN = I
         DMIN = DX(I)
   30 CONTINUE
      RETURN
      END
c **********************************************************************
c
c ----------------------------------------------------------------------
c blas.f
c
c some blas routines used in HAS:
c  daxpy dcopy ddot dnrm2 dscal dswap
c **********************************************************************
      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c     overwrite double precision dy with double precision da*dx + dy.
c     for i = 0 to n-1, replace  dy(ly+i*incy) with da*dx(lx+i*incx) +
c       dy(ly+i*incy), where lx = 1 if incx .ge. 0, else lx = (-incx)*n,
c       and ly is defined in a similar way using incy.
c
      double precision dx(1),dy(1),da
      if(n.le.0.or.da.eq.0.d0) return
      if(incx.eq.incy) if(incx-1) 5,20,60
    5 continue
c
c        code for nonequal or nonpositive increments.
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop so remaining vector length is a multiple of 4.
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
c
c        code for equal, positive, nonunit increments.
c
   60 continue
      ns = n*incx
          do 70 i=1,ns,incx
          dy(i) = da*dx(i) + dy(i)
   70     continue
      return
      end
C **********************************************************************
      subroutine dcopy(n,dx,incx,dy,incy)
c
c     copy double precision dx to double precision dy.
c     for i = 0 to n-1, copy dx(lx+i*incx) to dy(ly+i*incy),
c     where lx = 1 if incx .ge. 0, else lx = (-incx)*n, and ly is
c     defined in a similar way using incy.
c
      double precision dx(1),dy(1)
      if(n.le.0)return
      if(incx.eq.incy) if(incx-1) 5,20,60
    5 continue
c
c        code for unequal or nonpositive increments.
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop so remaining vector length is a multiple of 7.
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 continue
      return
c
c        code for equal, positive, nonunit increments.
c
   60 continue
      ns=n*incx
          do 70 i=1,ns,incx
          dy(i) = dx(i)
   70     continue
      return
      end
C **********************************************************************
      double precision function ddot(n,dx,incx,dy,incy)
c
c     returns the dot product of double precision dx and dy.
c     ddot = sum for i = 0 to n-1 of  dx(lx+i*incx) * dy(ly+i*incy)
c     where lx = 1 if incx .ge. 0, else lx = (-incx)*n, and ly is
c     defined in a similar way using incy.
c
      double precision dx(1),dy(1)
      ddot = 0.d0
      if(n.le.0)return
      if(incx.eq.incy) if(incx-1) 5,20,60
    5 continue
c
c         code for unequal or nonpositive increments.
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
         ddot = ddot + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1.
c
c
c        clean-up loop so remaining vector length is a multiple of 5.
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
         ddot = ddot + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
         ddot = ddot + dx(i)*dy(i) + dx(i+1)*dy(i+1) +
     1   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
      return
c
c         code for positive equal increments .ne.1.
c
   60 continue
      ns = n*incx
          do 70 i=1,ns,incx
          ddot = ddot + dx(i)*dy(i)
   70     continue
      return
      end
C **********************************************************************
      double precision function dnrm2 ( n, dx, incx)
      integer          next
      double precision   dx(1), cutlo, cuthi, hitest, sum, xmax,zero,one
      data   zero, one /0.0d0, 1.0d0/
c
c     euclidean norm of the n-vector stored in dx() with storage
c     increment incx .
c     if    n .le. 0 return with result = 0.
c     if n .ge. 1 then incx must be .ge. 1
c
c           c.l.lawson, 1978 jan 08
c
c     four phase method     using two built-in constants that are
c     hopefully applicable to all machines.
c         cutlo = maximum of  dsqrt(u/eps)  over all known machines.
c         cuthi = minimum of  dsqrt(v)      over all known machines.
c     where
c         eps = smallest no. such that eps + 1. .gt. 1.
c         u   = smallest positive no.   (underflow limit)
c         v   = largest  no.            (overflow  limit)
c
c     brief outline of algorithm..
c
c     phase 1    scans zero components.
c     move to phase 2 when a component is nonzero and .le. cutlo
c     move to phase 3 when a component is .gt. cutlo
c     move to phase 4 when a component is .ge. cuthi/m
c     where m = n for x() real and m = 2*n for complex.
c
c     values for cutlo and cuthi..
c     from the environmental parameters listed in the imsl converter
c     document the limiting values are as follows..
c     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
c                   univac and dec at 2**(-103)
c                   thus cutlo = 2**(-51) = 4.44089e-16
c     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.
c                   thus cuthi = 2**(63.5) = 1.30438e19
c     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.
c                   thus cutlo = 2**(-33.5) = 8.23181d-11
c     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19
c     data cutlo, cuthi / 8.232d-11,  1.304d19 /
c     data cutlo, cuthi / 4.441e-16,  1.304e19 /
      data cutlo, cuthi / 8.232d-11,  1.304d19 /
c
      if(n .gt. 0) go to 10
         dnrm2  = zero
         go to 300
c
   10 next = 30
      sum = zero
      nn = n * incx
c                                                 begin main loop
      i = 1
   20    if (next.eq.30) go to 30
		 if (next.eq.50) go to 50
		 if (next.eq.70) go to 70
		 if (next.eq.110) go to 110
   30 if( dabs(dx(i)) .gt. cutlo) go to 85
      next = 50
      xmax = zero
c
c                        phase 1.  sum is zero
c
   50 if( dx(i) .eq. zero) go to 200
      if( dabs(dx(i)) .gt. cutlo) go to 85
c
c                                prepare for phase 2.
      next = 70
      go to 105
c
c                                prepare for phase 4.
c
  100 i = j
      next = 110
      sum = (sum / dx(i)) / dx(i)
  105 xmax = dabs(dx(i))
      go to 115
c
c                   phase 2.  sum is small.
c                             scale to avoid destructive underflow.
c
   70 if( dabs(dx(i)) .gt. cutlo ) go to 75
c
c                     common code for phases 2 and 4.
c                     in phase 4 sum is large.  scale to avoid overflow.
c
  110 if( dabs(dx(i)) .le. xmax ) go to 115
         sum = one + sum * (xmax / dx(i))**2
         xmax = dabs(dx(i))
         go to 200
c
  115 sum = sum + (dx(i)/xmax)**2
      go to 200
c
c
c                  prepare for phase 3.
c
   75 sum = (sum * xmax) * xmax
c
c
c     for real or d.p. set hitest = cuthi/n
c     for complex      set hitest = cuthi/(2*n)
c
   85 hitest = cuthi/float( n )
c
c                   phase 3.  sum is mid-range.  no scaling.
c
      do 95 j =i,nn,incx
      if(dabs(dx(j)) .ge. hitest) go to 100
   95    sum = sum + dx(j)**2
      dnrm2 = dsqrt( sum )
      go to 300
c
  200 continue
      i = i + incx
      if ( i .le. nn ) go to 20
c
c              end of main loop.
c
c              compute square root and adjust for scaling.
c
      dnrm2 = xmax * dsqrt(sum)
  300 continue
      return
      end
C **********************************************************************
      subroutine dscal(n,da,dx,incx)
c
c     replace double precision dx by double precision da*dx.
c     for i = 0 to n-1, replace dx(1+i*incx) with  da * dx(1+i*incx)
c
      double precision da,dx(1)
      if(n.le.0)return
      if(incx.eq.1)goto 20
c
c        code for increments not equal to 1.
c
      ns = n*incx
          do 10 i = 1,ns,incx
          dx(i) = da*dx(i)
   10     continue
      return
c
c        code for increments equal to 1.
c
c
c        clean-up loop so remaining vector length is a multiple of 5.
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
C **********************************************************************
      subroutine dswap(n,dx,incx,dy,incy)
c
c     interchange double precision dx and double precision dy.
c     for i = 0 to n-1, interchange  dx(lx+i*incx) and dy(ly+i*incy),
c     where lx = 1 if incx .ge. 0, else lx = (-incx)*n, and ly is
c     defined in a similar way using incy.
c
      double precision dx(1),dy(1),dtemp1,dtemp2,dtemp3
      if(n.le.0)return
      if(incx.eq.incy) if(incx-1) 5,20,60
    5 continue
c
c       code for unequal or nonpositive increments.
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp1 = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp1
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
c
c       clean-up loop so remaining vector length is a multiple of 3.
c
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp1 = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp1
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        dtemp1 = dx(i)
        dtemp2 = dx(i+1)
        dtemp3 = dx(i+2)
        dx(i) = dy(i)
        dx(i+1) = dy(i+1)
        dx(i+2) = dy(i+2)
        dy(i) = dtemp1
        dy(i+1) = dtemp2
        dy(i+2) = dtemp3
   50 continue
      return
   60 continue
c
c     code for equal, positive, nonunit increments.
c
      ns = n*incx
        do 70 i=1,ns,incx
        dtemp1 = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp1
   70   continue
      return
      end
c **********************************************************************
c
c ----------------------------------------------------------------------
c lpk.f 
c
c  some LINPACK subroutines used in HAS:
c  dqrdc dqrsl dasum dtrsl
c **********************************************************************
      subroutine dqrdc(x,ldx,n,p,qraux,jpvt,work,job)
      integer ldx,n,p,job
      integer jpvt(1)
      double precision x(ldx,1),qraux(1),work(1)
c
c     dqrdc uses householder transformations to compute the qr
c     factorization of an n by p matrix x.  column pivoting
c     based on the 2-norms of the reduced columns may be
c     performed at the users option.
c
c     on entry
c
c        x       double precision(ldx,p), where ldx .ge. n.
c                x contains the matrix whose decomposition is to be
c                computed.
c
c        ldx     integer.
c                ldx is the leading dimension of the array x.
c
c        n       integer.
c                n is the number of rows of the matrix x.
c
c        p       integer.
c                p is the number of columns of the matrix x.
c
c        jpvt    integer(p).
c                jpvt contains integers that control the selection
c                of the pivot columns.  the k-th column x(k) of x
c                is placed in one of three classes according to the
c                value of jpvt(k).
c
c                   if jpvt(k) .gt. 0, then x(k) is an initial
c                                      column.
c
c                   if jpvt(k) .eq. 0, then x(k) is a free column.
c
c                   if jpvt(k) .lt. 0, then x(k) is a final column.
c
c                before the decomposition is computed, initial columns
c                are moved to the beginning of the array x and final
c                columns to the end.  both initial and final columns
c                are frozen in place during the computation and only
c                free columns are moved.  at the k-th stage of the
c                reduction, if x(k) is occupied by a free column
c                it is interchanged with the free column of largest
c                reduced norm.  jpvt is not referenced if
c                job .eq. 0.
c
c        work    double precision(p).
c                work is a work array.  work is not referenced if
c                job .eq. 0.
c
c        job     integer.
c                job is an integer that initiates column pivoting.
c                if job .eq. 0, no pivoting is done.
c                if job .ne. 0, pivoting is done.
c
c     on return
c
c        x       x contains in its upper triangle the upper
c                triangular matrix r of the qr factorization.
c                below its diagonal x contains information from
c                which the orthogonal part of the decomposition
c                can be recovered.  note that if pivoting has
c                been requested, the decomposition is not that
c                of the original matrix x but that of x
c                with its columns permuted as described by jpvt.
c
c        qraux   double precision(p).
c                qraux contains further information required to recover
c                the orthogonal part of the decomposition.
c
c        jpvt    jpvt(k) contains the index of the column of the
c                original matrix that has been interchanged into
c                the k-th column, if pivoting was requested.
c
c     linpack. this version dated 08/14/78 .
c     g.w. stewart, university of maryland, argonne national lab.
c
c     dqrdc uses the following functions and subprograms.
c
c     blas daxpy,ddot,dscal,dswap,dnrm2
c     fortran dabs,dmax1,min0,dsqrt
c
c     internal variables
c
      integer j,jp,l,lp1,lup,maxj,pl,pu
      double precision maxnrm,dnrm2,tt
      double precision ddot,nrmxl,t
      logical negj,swapj
c
c
      pl = 1
      pu = 0
      if (job .eq. 0) go to 60
c
c        pivoting has been requested.  rearrange the columns
c        according to jpvt.
c
         do 20 j = 1, p
            swapj = jpvt(j) .gt. 0
            negj = jpvt(j) .lt. 0
            jpvt(j) = j
            if (negj) jpvt(j) = -j
            if (.not.swapj) go to 10
               if (j .ne. pl) call dswap(n,x(1,pl),1,x(1,j),1)
               jpvt(j) = jpvt(pl)
               jpvt(pl) = j
               pl = pl + 1
   10       continue
   20    continue
         pu = p
         do 50 jj = 1, p
            j = p - jj + 1
            if (jpvt(j) .ge. 0) go to 40
               jpvt(j) = -jpvt(j)
               if (j .eq. pu) go to 30
                  call dswap(n,x(1,pu),1,x(1,j),1)
                  jp = jpvt(pu)
                  jpvt(pu) = jpvt(j)
                  jpvt(j) = jp
   30          continue
               pu = pu - 1
   40       continue
   50    continue
   60 continue
c
c     compute the norms of the free columns.
c
      if (pu .lt. pl) go to 80
      do 70 j = pl, pu
         qraux(j) = dnrm2(n,x(1,j),1)
         work(j) = qraux(j)
   70 continue
   80 continue
c
c     perform the householder reduction of x.
c
      lup = min0(n,p)
      do 200 l = 1, lup
         if (l .lt. pl .or. l .ge. pu) go to 120
c
c           locate the column of largest norm and bring it
c           into the pivot position.
c
            maxnrm = 0.0d0
            maxj = l
            do 100 j = l, pu
               if (qraux(j) .le. maxnrm) go to 90
                  maxnrm = qraux(j)
                  maxj = j
   90          continue
  100       continue
            if (maxj .eq. l) go to 110
               call dswap(n,x(1,l),1,x(1,maxj),1)
               qraux(maxj) = qraux(l)
               work(maxj) = work(l)
               jp = jpvt(maxj)
               jpvt(maxj) = jpvt(l)
               jpvt(l) = jp
  110       continue
  120    continue
         qraux(l) = 0.0d0
         if (l .eq. n) go to 190
c
c           compute the householder transformation for column l.
c
            nrmxl = dnrm2(n-l+1,x(l,l),1)
            if (nrmxl .eq. 0.0d0) go to 180
               if (x(l,l) .ne. 0.0d0) nrmxl = dsign(nrmxl,x(l,l))
               call dscal(n-l+1,1.0d0/nrmxl,x(l,l),1)
               x(l,l) = 1.0d0 + x(l,l)
c
c              apply the transformation to the remaining columns,
c              updating the norms.
c
               lp1 = l + 1
               if (p .lt. lp1) go to 170
               do 160 j = lp1, p
                  t = -ddot(n-l+1,x(l,l),1,x(l,j),1)/x(l,l)
                  call daxpy(n-l+1,t,x(l,l),1,x(l,j),1)
                  if (j .lt. pl .or. j .gt. pu) go to 150
                  if (qraux(j) .eq. 0.0d0) go to 150
                     tt = 1.0d0 - (dabs(x(l,j))/qraux(j))**2
                     tt = dmax1(tt,0.0d0)
                     t = tt
                     tt = 1.0d0 + 0.05d0*tt*(qraux(j)/work(j))**2
                     if (tt .eq. 1.0d0) go to 130
                        qraux(j) = qraux(j)*dsqrt(t)
                     go to 140
  130                continue
                        qraux(j) = dnrm2(n-l,x(l+1,j),1)
                        work(j) = qraux(j)
  140                continue
  150             continue
  160          continue
  170          continue
c
c              save the transformation.
c
               qraux(l) = x(l,l)
               x(l,l) = -nrmxl
  180       continue
  190    continue
  200 continue
      return
      end
C **********************************************************************
      subroutine dqrsl(x,ldx,n,k,qraux,y,qy,qty,b,rsd,xb,job,info)
      integer ldx,n,k,job,info
      double precision x(ldx,1),qraux(1),y(1),qy(1),qty(1),b(1),rsd(1),
     *                 xb(1)
c
c     dqrsl applies the output of dqrdc to compute coordinate
c     transformations, projections, and least squares solutions.
c     for k .le. min(n,p), let xk be the matrix
c
c            xk = (x(jpvt(1)),x(jpvt(2)), ... ,x(jpvt(k)))
c
c     formed from columnns jpvt(1), ... ,jpvt(k) of the original
c     n x p matrix x that was input to dqrdc (if no pivoting was
c     done, xk consists of the first k columns of x in their
c     original order).  dqrdc produces a factored orthogonal matrix q
c     and an upper triangular matrix r such that
c
c              xk = q * (r)
c                       (0)
c
c     this information is contained in coded form in the arrays
c     x and qraux.
c
c     on entry
c
c        x      double precision(ldx,p).
c               x contains the output of dqrdc.
c
c        ldx    integer.
c               ldx is the leading dimension of the array x.
c
c        n      integer.
c               n is the number of rows of the matrix xk.  it must
c               have the same value as n in dqrdc.
c
c        k      integer.
c               k is the number of columns of the matrix xk.  k
c               must nnot be greater than min(n,p), where p is the
c               same as in the calling sequence to dqrdc.
c
c        qraux  double precision(p).
c               qraux contains the auxiliary output from dqrdc.
c
c        y      double precision(n)
c               y contains an n-vector that is to be manipulated
c               by dqrsl.
c
c        job    integer.
c               job specifies what is to be computed.  job has
c               the decimal expansion abcde, with the following
c               meaning.
c
c                    if a.ne.0, compute qy.
c                    if b,c,d, or e .ne. 0, compute qty.
c                    if c.ne.0, compute b.
c                    if d.ne.0, compute rsd.
c                    if e.ne.0, compute xb.
c
c               note that a request to compute b, rsd, or xb
c               automatically triggers the computation of qty, for
c               which an array must be provided in the calling
c               sequence.
c
c     on return
c
c        qy     double precision(n).
c               qy conntains q*y, if its computation has been
c               requested.
c
c        qty    double precision(n).
c               qty contains trans(q)*y, if its computation has
c               been requested.  here trans(q) is the
c               transpose of the matrix q.
c
c        b      double precision(k)
c               b contains the solution of the least squares problem
c
c                    minimize norm2(y - xk*b),
c
c               if its computation has been requested.  (note that
c               if pivoting was requested in dqrdc, the j-th
c               component of b will be associated with column jpvt(j)
c               of the original matrix x that was input into dqrdc.)
c
c        rsd    double precision(n).
c               rsd contains the least squares residual y - xk*b,
c               if its computation has been requested.  rsd is
c               also the orthogonal projection of y onto the
c               orthogonal complement of the column space of xk.
c
c        xb     double precision(n).
c               xb contains the least squares approximation xk*b,
c               if its computation has been requested.  xb is also
c               the orthogonal projection of y onto the column space
c               of x.
c
c        info   integer.
c               info is zero unless the computation of b has
c               been requested and r is exactly singular.  in
c               this case, info is the index of the first zero
c               diagonal element of r and b is left unaltered.
c
c     the parameters qy, qty, b, rsd, and xb are not referenced
c     if their computation is not requested and in this case
c     can be replaced by dummy variables in the calling program.
c     to save storage, the user may in some cases use the same
c     array for different parameters in the calling sequence.  a
c     frequently occuring example is when one wishes to compute
c     any of b, rsd, or xb and does not need y or qty.  in this
c     case one may identify y, qty, and one of b, rsd, or xb, while
c     providing separate arrays for anything else that is to be
c     computed.  thus the calling sequence
c
c          call dqrsl(x,ldx,n,k,qraux,y,dum,y,b,y,dum,110,info)
c
c     will result in the computation of b and rsd, with rsd
c     overwriting y.  more generally, each item in the following
c     list contains groups of permissible identifications for
c     a single callinng sequence.
c
c          1. (y,qty,b) (rsd) (xb) (qy)
c
c          2. (y,qty,rsd) (b) (xb) (qy)
c
c          3. (y,qty,xb) (b) (rsd) (qy)
c
c          4. (y,qy) (qty,b) (rsd) (xb)
c
c          5. (y,qy) (qty,rsd) (b) (xb)
c
c          6. (y,qy) (qty,xb) (b) (rsd)
c
c     in any group the value returned in the array allocated to
c     the group corresponds to the last member of the group.
c
c     linpack. this version dated 08/14/78 .
c     g.w. stewart, university of maryland, argonne national lab.
c
c     dqrsl uses the following functions and subprograms.
c
c     blas daxpy,dcopy,ddot
c     fortran dabs,min0,mod
c
c     internal variables
c
      integer i,j,jj,ju,kp1
      double precision ddot,t,temp
      logical cb,cqy,cqty,cr,cxb
c
c
c     set info flag.
c
      info = 0
c
c     determine what is to be computed.
c
      cqy = job/10000 .ne. 0
      cqty = mod(job,10000) .ne. 0
      cb = mod(job,1000)/100 .ne. 0
      cr = mod(job,100)/10 .ne. 0
      cxb = mod(job,10) .ne. 0
      ju = min0(k,n-1)
c
c     special action when n=1.
c
      if (ju .ne. 0) go to 40
         if (cqy) qy(1) = y(1)
         if (cqty) qty(1) = y(1)
         if (cxb) xb(1) = y(1)
         if (.not.cb) go to 30
            if (x(1,1) .ne. 0.0d0) go to 10
               info = 1
            go to 20
   10       continue
               b(1) = y(1)/x(1,1)
   20       continue
   30    continue
         if (cr) rsd(1) = 0.0d0
      go to 250
   40 continue
c
c        set up to compute qy or qty.
c
         if (cqy) call dcopy(n,y,1,qy,1)
         if (cqty) call dcopy(n,y,1,qty,1)
         if (.not.cqy) go to 70
c
c           compute qy.
c
            do 60 jj = 1, ju
               j = ju - jj + 1
               if (qraux(j) .eq. 0.0d0) go to 50
                  temp = x(j,j)
                  x(j,j) = qraux(j)
                  t = -ddot(n-j+1,x(j,j),1,qy(j),1)/x(j,j)
                  call daxpy(n-j+1,t,x(j,j),1,qy(j),1)
                  x(j,j) = temp
   50          continue
   60       continue
   70    continue
         if (.not.cqty) go to 100
c
c           compute trans(q)*y.
c
            do 90 j = 1, ju
               if (qraux(j) .eq. 0.0d0) go to 80
                  temp = x(j,j)
                  x(j,j) = qraux(j)
                  t = -ddot(n-j+1,x(j,j),1,qty(j),1)/x(j,j)
                  call daxpy(n-j+1,t,x(j,j),1,qty(j),1)
                  x(j,j) = temp
   80          continue
   90       continue
  100    continue
c
c        set up to compute b, rsd, or xb.
c
         if (cb) call dcopy(k,qty,1,b,1)
         kp1 = k + 1
         if (cxb) call dcopy(k,qty,1,xb,1)
         if (cr .and. k .lt. n) call dcopy(n-k,qty(kp1),1,rsd(kp1),1)
         if (.not.cxb .or. kp1 .gt. n) go to 120
            do 110 i = kp1, n
               xb(i) = 0.0d0
  110       continue
  120    continue
         if (.not.cr) go to 140
            do 130 i = 1, k
               rsd(i) = 0.0d0
  130       continue
  140    continue
         if (.not.cb) go to 190
c
c           compute b.
c
            do 170 jj = 1, k
               j = k - jj + 1
               if (x(j,j) .ne. 0.0d0) go to 150
                  info = j
c           ......exit
                  go to 180
  150          continue
               b(j) = b(j)/x(j,j)
               if (j .eq. 1) go to 160
                  t = -b(j)
                  call daxpy(j-1,t,x(1,j),1,b,1)
  160          continue
  170       continue
  180       continue
  190    continue
         if (.not.cr .and. .not.cxb) go to 240
c
c           compute rsd or xb as required.
c
            do 230 jj = 1, ju
               j = ju - jj + 1
               if (qraux(j) .eq. 0.0d0) go to 220
                  temp = x(j,j)
                  x(j,j) = qraux(j)
                  if (.not.cr) go to 200
                     t = -ddot(n-j+1,x(j,j),1,rsd(j),1)/x(j,j)
                     call daxpy(n-j+1,t,x(j,j),1,rsd(j),1)
  200             continue
                  if (.not.cxb) go to 210
                     t = -ddot(n-j+1,x(j,j),1,xb(j),1)/x(j,j)
                     call daxpy(n-j+1,t,x(j,j),1,xb(j),1)
  210             continue
                  x(j,j) = temp
  220          continue
  230       continue
  240    continue
  250 continue
      return
      end
C***********************************************************************
      double precision function dasum(n,dx,incx)
c
c     takes the sum of the absolute values.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      double precision dx(1),dtemp
      integer i,incx,m,mp1,n,nincx
c
      dasum = 0.0d0
      dtemp = 0.0d0
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dtemp = dtemp + dabs(dx(i))
   10 continue
      dasum = dtemp
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,6)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dabs(dx(i))
   30 continue
      if( n .lt. 6 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,6
        dtemp = dtemp + dabs(dx(i)) + dabs(dx(i + 1)) + dabs(dx(i + 2))
     *  + dabs(dx(i + 3)) + dabs(dx(i + 4)) + dabs(dx(i + 5))
   50 continue
   60 dasum = dtemp
      return
      end
C***********************************************************************
      subroutine dtrsl(t,ldt,n,b,job,info)
      integer ldt,n,job,info
      double precision t(ldt,1),b(1)
c
c
c     dtrsl solves systems of the form
c
c                   t * x = b
c     or
c                   trans(t) * x = b
c
c     where t is a triangular matrix of order n. here trans(t)
c     denotes the transpose of the matrix t.
c
c     on entry
c
c         t         double precision(ldt,n)
c                   t contains the matrix of the system. the zero
c                   elements of the matrix are not referenced, and
c                   the corresponding elements of the array can be
c                   used to store other information.
c
c         ldt       integer
c                   ldt is the leading dimension of the array t.
c
c         n         integer
c                   n is the order of the system.
c
c         b         double precision(n).
c                   b contains the right hand side of the system.
c
c         job       integer
c                   job specifies what kind of system is to be solved.
c                   if job is
c
c                        00   solve t*x=b, t lower triangular,
c                        01   solve t*x=b, t upper triangular,
c                        10   solve trans(t)*x=b, t lower triangular,
c                        11   solve trans(t)*x=b, t upper triangular.
c
c     on return
c
c         b         b contains the solution, if info .eq. 0.
c                   otherwise b is unaltered.
c
c         info      integer
c                   info contains zero if the system is nonsingular.
c                   otherwise info contains the index of
c                   the first zero diagonal element of t.
c
c     linpack. this version dated 08/14/78 .
c     g. w. stewart, university of maryland, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c     fortran mod
c
c     internal variables
c
      double precision ddot,temp
      integer case,j,jj
c
c     begin block permitting ...exits to 150
c
c        check for zero diagonal elements.
c
         do 10 info = 1, n
c     ......exit
            if (t(info,info) .eq. 0.0d0) go to 150
   10    continue
         info = 0
c
c        determine the task and go to it.
c
         case = 1
         if (mod(job,10) .ne. 0) case = 2
         if (mod(job,100)/10 .ne. 0) case = case + 2
         go to (20,50,80,110), case
c
c        solve t*x=b for t lower triangular
c
   20    continue
            b(1) = b(1)/t(1,1)
            if (n .lt. 2) go to 40
            do 30 j = 2, n
               temp = -b(j-1)
               call daxpy(n-j+1,temp,t(j,j-1),1,b(j),1)
               b(j) = b(j)/t(j,j)
   30       continue
   40       continue
         go to 140
c
c        solve t*x=b for t upper triangular.
c
   50    continue
            b(n) = b(n)/t(n,n)
            if (n .lt. 2) go to 70
            do 60 jj = 2, n
               j = n - jj + 1
               temp = -b(j+1)
               call daxpy(j,temp,t(1,j+1),1,b(1),1)
               b(j) = b(j)/t(j,j)
   60       continue
   70       continue
         go to 140
c
c        solve trans(t)*x=b for t lower triangular.
c
   80    continue
            b(n) = b(n)/t(n,n)
            if (n .lt. 2) go to 100
            do 90 jj = 2, n
               j = n - jj + 1
               b(j) = b(j) - ddot(jj-1,t(j+1,j),1,b(j+1),1)
               b(j) = b(j)/t(j,j)
   90       continue
  100       continue
         go to 140
c
c        solve trans(t)*x=b for t upper triangular.
c
  110    continue
            b(1) = b(1)/t(1,1)
            if (n .lt. 2) go to 130
            do 120 j = 2, n
               b(j) = b(j) - ddot(j-1,t(1,j),1,b(1),1)
               b(j) = b(j)/t(j,j)
  120       continue
  130       continue
  140    continue
  150 continue
      return
      end
C***********************************************************************
C
C ----------------------------------------------------------------------

c nbsele.f, Yuedong Wang, Fri Nov  8 13:54:11 PST 2002
c ----------------------------------------------------------------------
      SUBROUTINE NBSELE (NOBS,NNULL,RSS,MAXBAS,CRIT,IDF,NB,SCORE)
c########################################################################c
c Use GCV(modified), AIC, BIC and Cp to decide how many basis to be used.
c 
c Input:
c     NOBS,NNULL,MAXBAS,IDF: same as in has1.f
c     RSS(MAXBAS): output from presele.
c     CRIT: criteria for deciding NB. CRIT=1,2,3 and 4 corresponds
c             to GCV (with IDF), AIC, BIC and Cp respectively
c Output:
c     NB: the number of basis functions chosen by modified GCV criterion
c     SCORE(MAXBAS): scores correponding to the selected criteria
c Functions used: idamin
c########################################################################c
      INTEGER NOBS,NNULL,MAXBAS,NB,MBAS, CRIT, NTEMP
      DOUBLE PRECISION RSS(*),SCORE(*),IDF
      INTEGER IDAMIN
      MBAS = MAXBAS
      NTEMP = NNULL
      IF(NNULL.EQ.0) NNULL = 1
      if(CRIT.eq.1) MBAS = MIN(MAXBAS,INT(DBLE(NOBS-NNULL)/IDF)+NNULL)

      if(CRIT.eq.1) go to 10
      if(CRIT.eq.2) go to 20
      if(CRIT.eq.3) go to 30
      if(CRIT.eq.4) go to 40

 10   DO 15 I=NNULL,MBAS
          SCORE(I) = (RSS(I)/DBLE(NOBS))/
     &      (1.D0-(IDF*DBLE(I-NNULL)+DBLE(NNULL))/DBLE(NOBS))**2
 15   CONTINUE
      NB = IDAMIN(MBAS-(NNULL-1),SCORE(NNULL),1)+(NNULL-1)
      NNULL = NTEMP
      RETURN

 20   DO 25 I=NNULL,MBAS
          SCORE(I) = DBLE(NOBS)*dlog(RSS(I)/DBLE(NOBS))
     & +2.D0*DBLE(I)
 25   CONTINUE
      NB = IDAMIN(MBAS-(NNULL-1),SCORE(NNULL),1)+(NNULL-1)
      NNULL = NTEMP
      RETURN

 30   DO 35 I=NNULL,MBAS
          SCORE(I) = DBLE(NOBS)*dlog(RSS(I)/DBLE(NOBS))+
     &               dlog(DBLE(NOBS))*DBLE(I)
 35   CONTINUE
      NB = IDAMIN(MBAS-(NNULL-1),SCORE(NNULL),1)+(NNULL-1)
      NNULL = NTEMP
      RETURN

 40   DO 45 I=NNULL,MBAS
          SCORE(I) = RSS(I)+2.D0*DBLE(I)*RSS(MBAS)/
     &               (DBLE(NOBS)-DBLE(MBAS))
 45   CONTINUE
      NB = IDAMIN(MBAS-(NNULL-1),SCORE(NNULL),1)+(NNULL-1)
      NNULL = NTEMP
      RETURN
      END

c presele2.f updated on Sep 28, 2009 by Junqing Wu, UCSB
c ---------------------------------------------------------------------
      SUBROUTINE PRESELE2(BAS,N,NBAS,Y,NNULL,MAXBAS,ISB,RSS,IWK,
     & WKA,U,H,FITS,Z,ZK,BASO)


c#####################################################################
c A subroutine for PREliminary SELEcting basis functions.
c     Select, from NBAS candidates, MAXBAS basis functions
c     which have the largest RSS reduction by LS fit in turn.
c     The fits with 0~MAXBAS (depending on NNULL and MAXBAS) 
c     basis functions in the model will also be calculated (new for
c     PRESELE2).
c	
c Use QR decomposition done by a sequence of Householder reflections
c     to compute each RSS in LS fit (without pivoting).
c     (ref. SEBER (1977) pp 338-340 )
c     QR decomposition and the application of it are updated;
c     BAS and Y are destroyed after calling.
c
c Subroutines used: 
c     those in suppl.f: DCOPY, DSET, ISET, (indirect: DAXPY),
c     those at the end of this file: HOUSEHOLDER, APPLYHOUSEHOLDER 
c Functions used: DDOT, DNRM2 (they are in suppl.f too)
c
c     BAS: Matrix of basis functions, NOBSxNBAS, the first NNULL 
c          columns are bases of the null space. Its
c          content will be modified due to Householder transformations
c          being applied throughout the basis selection process
c     N:      the number of observations
c     NBAS:   the number of candidate basis functions
c     Y:      (N), the response
c     NNULL:  the number of basis functions which will be forced
c             into the selected subset
c     MAXBAS: the maximum number of basis functions to be chosen
c     ISB:    (MAXBAS), the column index of chosen basis functions
c             in BAS
c     RSS:    (MAXBAS), RSS(I) is the residual sum of squares of the LS 
c             fit when regressing Y on the first I chosen basis
c             functions.
c     IWK:    (NBAS), work vector,
c             IWK(I) is the index of whether the I-th basis function
c             in BAS has been chosen or not;  0 -- no, 1 -- yes
c     WKA:    (N,MAXBAS), work array, stores the "work vectors" for
c		contructing Householder matrices
c     U: The norms of the original vectors from which the "working
c	  vectors" are derived with the opposite sign of the first
c	  element.
c     H: The first elements of the vectors corresponding to the
c        "working vectors" obtained from the transformed y's 
c     Z: The first k elements are obtained from y after k 
c        transformations, the rest are zeroes.
c     ZK: The first K elements are obtained from the (k+1)th chosen
c         basis after k transformations, the rest are zeroes.
c     BASO: Same as the initial content of BAS, but its content 
c           will not be changed.
c
c Sep 28, 2009 modified by Junqing Wu, UCSB
c######################################################################

      INTEGER N, NBAS, NNULL, MAXBAS,  IWK(NBAS), ISB(*), FLAG
      DOUBLE PRECISION BAS(N,*), Y(*), RSS(*),WKA(N,*),PREV,CURR
      DOUBLE PRECISION U(*),H(*),FITS(N,*),Z(*),ZK(*),FF(N)
      INTEGER I, J, IB, JB
      DOUBLE PRECISION NU, DE, A, BNEW, BASO(N,*)
      DOUBLE PRECISION DDOT, ZER, DIF


ccccccccccc    added for the case of NNULL=0
      CALL DSET(N,0.d0,WKA(1,1),1)
      IF(NNULL.eq.0) go to 11
ccccccccccc
 
c     first NNULL basis functions in BAS are forced into selected set.
      ISB(1) = 1
      IWK(1) = 1
      CALL HOUSEHOLDER(N,1,BAS(1,1),WKA(1,1),U)
      NU = DDOT(N,BAS(1,1),1,Y(1),1)
      DE = DDOT(N,BAS(1,1),1,BAS(1,1),1)
      A = NU/DE
      DO 5 I = 1, N
      FITS(I,1) = A * BAS(I,1)
5     CONTINUE

cccccccccccccccccccccccccccccccccccccccc
      IF(NNULL.EQ.1) GOTO 11

      DO 10 IB = 2,NNULL
         CALL APPLYHOUSEHOLDER(N,WKA(1,IB-1),Y)
         H(IB-1) = Y(IB-1)

      IF (IB.LT.3) GOTO 8
        CALL DSET(N-IB+2,0.d0,Z(IB-1),1)
        CALL DSET(N-IB+2,0.d0,ZK(IB-1),1)
      DO 7 I = 1, (IB-2)
        CALL APPLYHOUSEHOLDER(N,WKA(1,IB-I-1),Z(1))
        CALL APPLYHOUSEHOLDER(N,WKA(1,IB-I-1),ZK(1))
 7    CONTINUE
      BNEW = H(IB-1)/U(IB-1)*(-1.d0)
        CALL DAXPY(N,-1.D0,BASO(1,IB-1),1,ZK(1),1)
        CALL DAXPY(N,BNEW,ZK(1),1,Z(1),1)
        CALL DCOPY(N,Z(1),1,FITS(1,IB-1),1)

 8       CALL DCOPY(IB-1,Y(1),1,Z(1),1)
 
       DO 9 J = IB,NBAS
            IF(IWK(J).EQ.1) GOTO 9
            CALL APPLYHOUSEHOLDER(N,WKA(1,IB-1),BAS(1,J))
 9       CONTINUE
         CALL DCOPY(IB-1,BAS(1,IB),1,ZK(1),1)
         ISB(IB) = IB
         IWK(IB) = 1
         CALL HOUSEHOLDER(N,IB,BAS(1,IB),WKA(1,IB),U)
 10   CONTINUE
cccccccccccccccccccccccccccccccccccccccc        
11      DO 20 IB = NNULL+1, MAXBAS+1
	  PREV = 9.D30
	  IF (IB.EQ.1) GOTO 12
         CALL APPLYHOUSEHOLDER(N,WKA(1,IB-1),Y)
         H(IB-1) = Y(IB-1)

         RSS(IB-1) = DDOT(N-IB+1,Y(IB),1,Y(IB),1)
         PREV = 0.D0
12       ISB(IB) = 0 
         CURR = 0.D0

      IF (IB.LT.3) GOTO 15
        CALL DSET(N-IB+2,0.d0,Z(IB-1),1)
        CALL DSET(N-IB+2,0.d0,ZK(IB-1),1)
      DO 13 I = 1, (IB-2)
        CALL APPLYHOUSEHOLDER(N,WKA(1,IB-I-1),Z(1))
        CALL APPLYHOUSEHOLDER(N,WKA(1,IB-I-1),ZK(1))
 13   CONTINUE
      BNEW = H(IB-1)/U(IB-1)*(-1.d0)
        CALL DAXPY(N,-1.D0,BASO(1,ISB(IB-1)),1,ZK(1),1)
        CALL DAXPY(N,BNEW,ZK(1),1,Z(1),1)
        CALL DCOPY(N,Z(1),1,FITS(1,IB-1),1)

 15     IF (IB.GT.MAXBAS) GOTO 20
        CALL DCOPY(IB-1,Y(1),1,Z(1),1)
 
        DO 18 JB = NNULL+1,NBAS
          IF(IWK(JB).EQ.1) GOTO 18
		IF (IB.EQ.1) GOTO 16
            CALL APPLYHOUSEHOLDER(N,WKA(1,IB-1),BAS(1,JB))
            CURR = (DDOT(N-IB+1,BAS(IB,JB),1,Y(IB),1))**2/
     &           DDOT(N-IB+1,BAS(IB,JB),1,BAS(IB,JB),1)
            IF(CURR.LE.PREV) GOTO 18
            PREV = CURR
            ISB(IB) = JB

16	   IF (IB.GT.1) GOTO 18
	   NU = DDOT(N,BAS(1,JB),1,Y(1),1)
         DE = DDOT(N,BAS(1,JB),1,BAS(1,JB),1)
         A = NU/DE
         DO 17 I = 1, N
         FF(I) = A * BAS(I,JB)
17       CONTINUE
         CALL DAXPY(N,-1.D0,Y(1),1,FF(1),1)
	   CURR = DDOT(N,FF(1),1,FF(1),1)
	   IF (CURR.GT.PREV) GOTO 18
	   PREV = CURR
	   ISB(IB) = JB

 18      CONTINUE
         CALL DCOPY(IB-1,BAS(1,ISB(IB)),1,ZK(1),1)

         IWK(ISB(IB)) = 1
         CALL HOUSEHOLDER(N,IB,BAS(1,ISB(IB)),WKA(1,IB),U)
         IF (IB.GT.1) GOTO 20
      NU = DDOT(N,BAS(1,ISB(IB)),1,Y(1),1)
      DE = DDOT(N,BAS(1,ISB(IB)),1,BAS(1,ISB(IB)),1)
      A = NU/DE
      DO 19 I = 1, N
      FITS(I,1) = A * BAS(I,ISB(IB))
 19   CONTINUE
 20   CONTINUE
      RETURN
      END 


      SUBROUTINE HOUSEHOLDER(N,J,A,W,U)
c     Compute the Householder matrix for eliminating the last
c     (N-J) elements of N-vector A. Store the transformation matrix
c     information in N-vector W (the first (J-1) elements of W
c     are zeros). The Householder matrix H = I - 2*W*W'.
c     U is defined in subroutine PRESELE2.
c     New version modified by Junqing Wu, UCSB. Sep 28, 2009

      INTEGER N,J, II
      DOUBLE PRECISION A(*), W(*), S, DNRM2, U(*)
      S = DNRM2(N-J+1,A(J),1)
      U(J) = DSIGN(S,A(J))*(-1.d0)
      CALL DSET(J-1,0.D0,W(1),1)
      W(J) = A(J) + DSIGN(S,A(J))
      S = DSQRT(2.D0*S*(S+DABS(A(J))))
      W(J) = W(J)/S
      DO 1 II = J+1,N
         W(II) = A(II)/S
 1    CONTINUE
      RETURN
      END


      SUBROUTINE APPLYHOUSEHOLDER(N,W,V)
c     Apply a Householder reflection, i.e. multiply a Householder
c     matrix H = I - 2*W*W', to the N-vector V, where W is a N-vector
c     and the first (J-1) elements of it are zeros.
      INTEGER N
      DOUBLE PRECISION W(*), V(*), T, DDOT
      T = - DDOT(N,W(1),1,V(1),1) * 2.D0
      CALL DAXPY(N,T,W(1),1,V(1),1)
      RETURN
      END

      REAL FUNCTION RNOR()
C***BEGIN PROLOGUE  RNOR
C***DATE WRITTEN   810915 (YYMMDD)
C***REVISION DATE  870419 (YYMMDD)
C***CATEGORY NO.  L6A14
C***KEYWORDS  RANDOM NUMBERS, NORMAL DEVIATES
C***AUTHOR    KAHANER, DAVID, SCIENTIFIC COMPUTING DIVISION, NBS
C             MARSAGLIA, GEORGE, SUPERCOMPUTER RES. INST., FLORIDA ST. U.
C
C***PURPOSE  GENERATES NORMAL RANDOM NUMBERS, WITH MEAN ZERO AND
C             UNIT STANDARD DEVIATION, OFTEN DENOTED N(0,1).
C***DESCRIPTION
C
C       RNOR generates normal random numbers with zero mean and
C       unit standard deviation, often denoted N(0,1).
C           From the book, "Numerical Methods and Software" by
C                D. Kahaner, C. Moler, S. Nash
C                Prentice Hall, 1988
C   Use 
C       First time....
C                   Z = RSTART(ISEED)
C                     Here ISEED is any  n o n - z e r o  integer.
C                     This causes initialization of the program.
C                     RSTART returns a real (single precision) echo of ISEED.
C
C       Subsequent times...
C                   Z = RNOR()
C                     Causes the next real (single precision) random number
C                           to be returned as Z.
C
C.....................................................................
C                 Typical usage
C
C                    REAL RSTART,RNOR,Z
C                    INTEGER ISEED,I
C                    ISEED = 305
C                    Z = RSTART(ISEED)
C                    DO 1 I = 1,10
C                       Z = RNOR()
C                       WRITE(*,*) Z
C                 1  CONTINUE
C                    END
C
C
C***REFERENCES  MARSAGLIA & TSANG, "A FAST, EASILY IMPLEMENTED
C                 METHOD FOR SAMPLING FROM DECREASING OR
C                 SYMMETRIC UNIMODAL DENSITY FUNCTIONS", TO BE
C                 PUBLISHED IN SIAM J SISC 1983.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  RNOR
      REAL AA,B,C,C1,C2,PC,X,Y,XN,V(65),RSTART,U(17),S,T,UN 
      INTEGER J,IA,IB,IC,II,JJ,ID,III,JJJ
      SAVE U,II,JJ 
C
      DATA AA,B,C/12.37586,.4878992,12.67706/
      DATA C1,C2,PC,XN/.9689279,1.301198,.1958303E-1,2.776994/
      DATA V/ .3409450, .4573146, .5397793, .6062427, .6631691
     +, .7136975, .7596125, .8020356, .8417227, .8792102, .9148948
     +, .9490791, .9820005, 1.0138492, 1.0447810, 1.0749254, 1.1043917
     +,1.1332738, 1.1616530, 1.1896010, 1.2171815, 1.2444516, 1.2714635
     +,1.2982650, 1.3249008, 1.3514125, 1.3778399, 1.4042211, 1.4305929
     +,1.4569915, 1.4834526, 1.5100121, 1.5367061, 1.5635712, 1.5906454
     +,1.6179680, 1.6455802, 1.6735255, 1.7018503, 1.7306045, 1.7598422
     +,1.7896223, 1.8200099, 1.8510770, 1.8829044, 1.9155830, 1.9492166
     +,1.9839239, 2.0198430, 2.0571356, 2.0959930, 2.1366450, 2.1793713
     +,2.2245175, 2.2725185, 2.3239338, 2.3795007, 2.4402218, 2.5075117
     +,2.5834658, 2.6713916, 2.7769943, 2.7769943, 2.7769943, 2.7769943/

C      Load data array in case user forgets to initialize.
C      This array is the result of calling UNI 100000 times
C         with seed 305.
      DATA U/
     *0.8668672834288,  0.3697986366357,  0.8008968294805,
     *0.4173889774680,  0.8254561579836,  0.9640965269077,
     *0.4508667414265,  0.6451309529668,  0.1645456024730,
     *0.2787901807898,  0.06761531340295, 0.9663226330820,
     *0.01963343943798, 0.02947398211399, 0.1636231515294,
     *0.3976343250467,  0.2631008574685/
C
      DATA II,JJ / 17, 5 /
C
C***FIRST EXECUTABLE STATEMENT  RNOR
C
C Fast part...
C
C 
C   Basic generator is Fibonacci
C 
      UN = U(II)-U(JJ)
      IF(UN.LT.0.0) UN = UN+1.
      U(II) = UN
C           U(II) and UN are uniform on [0,1)
C           VNI is uniform on [-1,1)
      VNI = UN + UN -1.
      II = II-1
      IF(II.EQ.0)II = 17
      JJ = JJ-1
      IF(JJ.EQ.0)JJ = 17
C        INT(UN(II)*128) in range [0,127],  J is in range [1,64]
      J = MOD(INT(U(II)*128),64)+1
C        Pick sign as VNI is positive or negative
      RNOR = VNI*V(J+1)
      IF(ABS(RNOR).LE.V(J))RETURN      
C
C Slow part; AA is a*f(0)
      X = (ABS(RNOR)-V(J))/(V(J+1)-V(J))
C          Y is uniform on [0,1)
      Y = U(II)-U(JJ)
      IF(Y.LT.0.0) Y = Y+1.
      U(II) = Y
      II = II-1
      IF(II.EQ.0)II = 17
      JJ = JJ-1
      IF(JJ.EQ.0)JJ = 17
C
      S = X+Y
      IF(S.GT.C2)GO TO 11
      IF(S.LE.C1)RETURN
      IF(Y.GT.C-AA*EXP(-.5*(B-B*X)**2))GO TO 11
      IF(EXP(-.5*V(J+1)**2)+Y*PC/V(J+1).LE.EXP(-.5*RNOR**2))RETURN
C
C Tail part; .3601016 is 1./XN
C       Y is uniform on [0,1)
   22 Y = U(II)-U(JJ)
      IF(Y.LE.0.0) Y = Y+1.
      U(II) = Y
      II = II-1
      IF(II.EQ.0)II = 17
      JJ = JJ-1
      IF(JJ.EQ.0)JJ = 17
C 
      X = 0.3601016*LOG(Y)
C       Y is uniform on [0,1)
      Y = U(II)-U(JJ)
      IF(Y.LE.0.0) Y = Y+1.
      U(II) = Y
      II = II-1
      IF(II.EQ.0)II = 17
      JJ = JJ-1
      IF(JJ.EQ.0)JJ = 17
      IF( -2.*LOG(Y).LE.X**2 )GO TO 22
      RNOR = SIGN(XN-X,RNOR)
      RETURN
   11 RNOR = SIGN(B-B*X,RNOR)
      RETURN
C
C
C  Fill
      ENTRY RSTART(ISEED)
      IF(ISEED.NE.0) THEN
C 
C          Set up ...
C              Generate random bit pattern in array based on given seed
C 
        II = 17
        JJ = 5
        IA = MOD(ABS(ISEED),32707)
        IB = 1111
        IC = 1947
        DO 2 III = 1,17
          S = 0.0
          T = .50
C             Do for each of the bits of mantissa of word 
C             Loop  over 64 bits, enough for all known machines
C                   in single precision
          DO 3 JJJ = 1,64
                  ID = IC-IA
                  IF(ID.GE.0)GOTO 4
                  ID = ID+32707
                  S = S+T
    4             IA = IB
                  IB = IC
                  IC = ID
    3     T = .5*T
    2   U(III) = S
      ENDIF
C       Return floating echo of ISEED
      RSTART=ISEED
      RETURN
      END 
c Updated (corrected) 3-08-06
c has1 (Modified by Yuedong Wang, Thu Nov  7 08:06:23 PST 2002)
c modified by jeff sklar 1-19-2003 to allow for smoothing parameter 
c selection options.  Now dsnsm1 is called (altered from dsnsm)
c ----------------------------------------------------------------------
      SUBROUTINE has1(NOBS, Y, BAS, NNULL, NBAS, MAXBAS, CRIT, 
     & IDF, DELTA, COND,
     & NB, ISB, COEF, HASFIT, HASDEV, HASDF, SCORE, RSS, INFO,
     & LAMLIM, MAXTBL, NTBL, JOB, TAU,
     & IWK, IWORK, YY, BASWK, WKA, DES, SIGMA, ADIAG,
     & SVALS, TBL,WORK,VMU,U,H,FITSWK,Z,ZK)

c     & SVALS, TBL, WORK,PREDES,BASWORK,VMU)

c########################################################################
c A driver program for estimating function f by HAS procedure,
c     using generic input of basis functions, from data:
c     y_i = f(x_i) + noise, i = 1,2,...,n                                
c     where x_i is in [0,1].
c                                                                        
c Subroutines called:                                                    
c     presele, nbsele, suppl.f (a group of subroutines from BLAS and 
c     LINPACK), and gcvpack.f (a group of subroutines from GCVPACK).
c                                                                        
c########################################################################
      INTEGER NOBS, NNULL, NBAS, MAXBAS, CRIT, NB, ISB(*), INFO, 
     & MAXTBL, NTBL, JOB, IWK(*), IWORK(*), COND, VMU
      DOUBLE PRECISION Y(*), BAS(NOBS,*), IDF, DELTA, COEF(*), 
     & HASFIT(*), HASDEV, HASDF, SCORE(*), RSS(*), TAU, LAMLIM(2), 
     & YY(*), BASWK(NOBS,*), WKA(NOBS,*), DES(NOBS,*), 
     & SIGMA(MAXBAS,*), ADIAG(*), SVALS(*), TBL(MAXTBL,3), WORK(*),
     & U(*),H(*),FITSWK(NOBS,MAXBAS),Z(*),ZK(*)
c########################################################################
c Inputs:
c     NOBS:   number of observations
c     Y:      response.
c     BAS:    NOBSxNBAS, Bases. The first NNULL columns are bases of the null space
c     NNULL:  dimension of null space.
c     NBAS:   number of bases
c     MAXBAS: the maximum number of basis functions to be chosen
c     CRIT: criteria for deciding NB. CRIT=1,2,3 and 4 corresponds
c             to GCV (with IDF), AIC, BIC and Cp respectively
c     IDF:    the inflated degrees of freedom for each adaptively 
c             selected basis function. (1.2 is recommended for general 
c             purpose, it should be, however, between 1 and 2).
c     DELTA:  a constant to be added to the diagnol of the SIGMA matrix
c 
c Outputs:
c     NB:     number of basis functions chosen 
c     ISB:    (MAXBAS), the column index of chosen bases in BAS. It
c             always include the first NNULL columns
c     COEF:   coefficients of the penalized regression on the 
c             selected basis functions. (first two are for the two
c             basis functions in null space, i.e., 1 and (x-.5).
c     HASFIT: HAS estimates of the noiseless response at original
c             data points.
c     HASDEV: estimate of standard deviation of noise from final HAS fit
c             (only for reference).
c     HASDF:  degrees of freedom from final HAS fit (only for reference).
c     RSS:    (MAXBAS), RSS(I) is the residual sum of squares of the LS 
c             fit when regressing Y on the first I chosen basis functions.
c     SCORE:  (MAXBAS), scores correponding to RSS
c########################################################################
c Variables used only in dsnsm1 (altered dsnsm from GCVPACK):
c     MAXTBL,LWA,LIWA,NTBL,IOUT,IWORK,DES,SIGMA,ADIAG,TAU,LAMLIM,
c     DOUT,SVALS,TBL,AUXTBL,WORK.
c     The meanings of them can be found in dsnsm.f in GCVPACK.
c     VMU:  integer indicating smoothing parameter selection method
c           1 = GCV, 0 = GML, 2 = UBR
c########################################################################
      INTEGER LWA, LIWA, IOUT(3), I, J
      DOUBLE PRECISION DUM, DOUT(5),AUXTBL(3,3)
c      NBAS = NNULL+NOBS  # blocked to allow more felxible bases
      LIWA = 2*MAXBAS-NNULL
      LWA = (MAXBAS-NNULL)*(MAXBAS-2*NNULL+2+NOBS)+MAXBAS+NOBS 
      
c call presele to select basis functions one by one (i.e. order them).
      CALL DCOPY(NOBS,Y(1),1,YY(1),1)
      call dcopy(NOBS*NBAS, BAS, 1, BASWK, 1)

ccc eliminate 3-8-06 causing a bug in ridge step
c      CALL PRESELE(BASWK,NOBS,NBAS,YY,NNULL,MAXBAS,ISB,RSS,IWK,
c     & WKA,FITS,PREDES,BASWORK)

      CALL PRESELE2(BASWK,NOBS,NBAS,YY,NNULL,MAXBAS,ISB,RSS,IWK,
     & WKA,U,H,FITSWK,Z,ZK,BAS)

c      if(cond.eq.1) return
c      if(nb.eq.nnull) return      
 
c call nbsele to decide the number of basis functions to be used
      CALL NBSELE(NOBS,NNULL,RSS,MAXBAS,CRIT,IDF,NB,SCORE)

      if(cond.eq.1) return
      if(nb.eq.nnull) return      
 
 

c   compute the first NB selected basis
      DO 43 I = 1,NOBS
         DO 42 J = 1,NB
            DES(I,J) = BAS(I,ISB(J))
 42      CONTINUE
 43   CONTINUE

c   compute SIGMA.
c   Sometimes a matrix DELTA*IDENTITY is added to SIGMA to ensure 
c   the solution of the linear system from dsnsm is numerically stable,
c   otherwise dsnsm may give the error message 221.
      DO 45 I=1,NB
         CALL DSET(NB,0.D0,SIGMA(1,I),1)
 45   CONTINUE
      DO 47 I=NNULL+1,NB
         DO 46 J=NNULL+1,NB
            SIGMA(I,J) = BAS(ISB(I)-NNULL,ISB(J))
 46      CONTINUE
         SIGMA(I,I) = SIGMA(I,I) + DELTA 
 47   CONTINUE
c   
c call dsnsm1 from GCVPACK to compute the penalized regression problem
      CALL DCOPY(NOBS,Y,1,HASFIT,1)
      CALL DSNSM1(DES,NOBS,HASFIT,SIGMA,MAXBAS,NOBS,NB,NNULL,ADIAG,
     &     TAU,LAMLIM,NTBL,DOUT,IOUT,COEF,SVALS,TBL,MAXTBL,AUXTBL,
     &     IWORK,LIWA,WORK,LWA,JOB,INFO,VMU)
      IF(INFO .gt. 0) RETURN 

      HASDEV = SQRT(DOUT(3)/(DBLE(NOBS)-(DBLE(NOBS)-DOUT(4)-
     & DBLE(NNULL))*IDF-DBLE(NNULL)))
      HASDF = DBLE(NOBS)-DOUT(4)

      RETURN
      END

c  pickbas.f updated Feb 23, 2009 by Junqing Wu, UCSB
c---------------------------------------------------------------
      SUBROUTINE PICKBAS(NOBS, Y, BAS, NNULL, NBAS,  
     & IDFHAT, nrep, GDFHAT, TAUHAT, ISB, ICB, RSS, IWK, YY, BASWK, WKA,
     & WKALL,WKAMAX,PTB,U,H,FITS,FITSWK,Z,ZK,YYPTB)

      INTEGER NOBS, NNULL, NBAS, ISB, nrep, WKAMAX, IWK(*), ICB(nrep)

      DOUBLE PRECISION Y(*), BAS(NOBS,*), IDFHAT,  
     & RSS, YY(*), YYPTB(NOBS,*), BASWK(NOBS,*), WKA(NOBS,*), 
     & GDFHAT, TAUHAT,tmp1,tmp2, WKALL(NOBS,WKAMAX,nrep),
     & dbar(NOBS), YWK(NOBS), WKAA(NOBS,WKAMAX),
     & d(NOBS,nrep),FITSWK(NOBS),
     & PTB(NOBS,nrep),U(*),H(*),FITS(NOBS,nrep),Z(*),ZK(*)

c########################################################################
c     NOBS:   number of observations
c     Y:      response.
c     BAS:    NOBSxNBAS, Bases. The first NNULL columns are bases of the 
c             null space
c     NNULL:  Number of basis functions already in the model
c     NBAS:   number of bases
c     IDFHAT: The estimated IDF values
c     NREP:   Number of repetitions of fitting perturbed responses 
c     GDFHAT: Total Cost for all chosen basis functions at each step
c     TAUHAT: .5*(est of sigma)  
c     ISB:    Index of Selected Bases, (MAXBAS), the column index of
c             chosen bases in BAS. It always include the first NNULL
c             columns.
c     RSS:    residual sum of squares for each model
c		has been chosen, 1-yes, 0-no
c     YY:     working vector for storing perturbed responses
c     BASWK:  working matrix (NOBSxNBAS) for basis functions
c     WKA:    working matrix (NOBSxNBAS) for "working vectors"
c     WKALL:  working matrix (NOBSxNBASxNREP) for all "working vectors"
c             for all perturbed responses.
c     WKAMAX: must be >MAXBAS, it ensures enough memory space for 
c		selecting up to MAXBAS basis functions.
c     PTB: perturbations generated inside Fortran.
c     U: The norms of the original vectors from which the "working
c	  vectors" are derived with the opposite sign of the first
c	  element.
c     H: The first elements of the vectors corresponding to the
c        "working vectors" obtained from the transformed y's 
c     FITS:   array of fits, (NOBSx(MAXBAS-NNULL+1)xNREP)
c     FITSWK: working matrix to store fits temporarily
c     Z: The first k elements are obtained from y after k 
c        transformations, the rest are zeroes.
c     ZK: The first K elements are obtained from the (k+1)th chosen
c         basis after k transformations, the rest are zeroes.
c     YYPTB: perturbed responses, updated upon the selection of a basis
c Modified by Junqing Wu, UCSB on Sep 28, 2009
c########################################################################


      INTEGER I, J, K, NB, IDXTMP(NBAS)
      DOUBLE PRECISION TMP
	  DO 5 K = 1,NBAS
	    IDXTMP(K)=IWK(K)
5	  CONTINUE 

         DO 20 J = 1,nrep
           DO 10 i=1, NOBS
              d(i,j) = PTB(i,j)
              ywk(i) = YYPTB(I,J)
 10        CONTINUE

        CALL DCOPY(NOBS,YWK(1),1,YY(1),1)
        CALL DCOPY(NOBS*NBAS, BAS, 1, BASWK, 1)
        
        CALL DSET(NOBS,0.D0,FITSWK,1)
	    CALL DCOPY(NOBS*WKAMAX,WKALL(1,1,J),1,WKAA(1,1),1)
          CALL PRESELE3(BASWK,NOBS,NBAS,YY,NNULL,ISB,
     &             RSS,IWK,WKAA,U,H,FITSWK,Z,ZK,BAS)
	    CALL DCOPY(NOBS*WKAMAX,WKAA(1,1),1,WKALL(1,1,J),1)
	    CALL DCOPY(NOBS,YY,1,YYPTB(1,J),1)
          CALL DCOPY(NOBS,FITSWK,1,FITS(1,J),1)

          ICB(J) = ISB

	  DO 17 K = 1,NBAS
	    IWK(K)=IDXTMP(K)
 17   CONTINUE

 20   CONTINUE

          CALL DCOPY(NOBS,Y(1),1,YY(1),1)
          CALL DCOPY(NOBS*NBAS, BAS, 1, BASWK, 1)
          CALL DSET(NOBS,0.D0,FITSWK,1)
          CALL PRESELE3(BASWK,NOBS,NBAS,YY,NNULL,ISB,
     &             RSS,IWK,WKA,U,H,FITSWK,Z,ZK,BAS)



	DO 25 I=1, NOBS
	DBAR(I) = 0.D0
	  DO 23 J=1, nrep
	  DBAR(I) = DBAR(I) + D(I,J)
23	  CONTINUE
	DBAR(I) = DBAR(I)/DBLE(nrep)
25	CONTINUE


	  GDFHAT = 0.D0
        IDFHAT = 0.D0
	  DO 40 I=1, NOBS
	tmp1 = 0.d0
	tmp2 = 0.d0
	    DO 30 J=1, nrep
	    tmp1 = tmp1 + (D(I,J)-DBAR(I))*FITS(I,J)
	    tmp2 = tmp2 + (D(I,J)-DBAR(I))**2
30	    CONTINUE
	  GDFHAT = GDFHAT + tmp1/tmp2
40	  CONTINUE
	IDFHAT = (GDFHAT-NNULL)
	IF (IDFHAT.LE.0.D0) THEN
		IDFHAT = 0.D0
	ENDIF
      RETURN
      END

c  pickbas2.f updated Jul 4, 2010 by Junqing Wu, UCSB
c---------------------------------------------------------------
      SUBROUTINE PICKBAS2(NOBS, Y, BAS, NNULL, NBAS,  
     & IDFHAT, nrep, GDFHAT, TAUHAT, ISB, ICB, RSS, IWK, YY, BASWK, WKA,
     & WKALL,WKAMAX,PTB,U,H,FITS,FITSWK,Z,ZK,YYPTB)

      INTEGER NOBS, NNULL, NBAS, ISB, nrep, WKAMAX, IWK(*), ICB(nrep)

      DOUBLE PRECISION Y(*), BAS(NOBS,*), IDFHAT,  
     & RSS, YY(*), YYPTB(NOBS,*), BASWK(NOBS,*), WKA(NOBS,*), 
     & GDFHAT, TAUHAT,tmp1,tmp2, WKALL(NOBS,WKAMAX,nrep),
     & dbar(NOBS), YWK(NOBS), WKAA(NOBS,WKAMAX), YBAR(NOBS),
     & d(NOBS,nrep),FITSWK(NOBS),
     & PTB(NOBS,nrep),U(*),H(*),FITS(NOBS,nrep),Z(*),ZK(*)

c########################################################################
c     NOBS:   number of observations
c     Y:      response.
c     BAS:    NOBSxNBAS, Bases. The first NNULL columns are bases of the 
c             null space
c     NNULL:  Number of basis functions already in the model
c     NBAS:   number of bases
c     IDFHAT: The estimated IDF values
c     NREP:   Number of repetitions of fitting perturbed responses 
c     GDFHAT: Total Cost for all chosen basis functions at each step
c     TAUHAT: .5*(est of sigma)  
c     ISB:    Index of Selected Bases, (MAXBAS), the column index of
c             chosen bases in BAS. It always include the first NNULL
c             columns.
c     RSS:    residual sum of squares for each model
c		has been chosen, 1-yes, 0-no
c     YY:     working vector for storing perturbed responses
c     BASWK:  working matrix (NOBSxNBAS) for basis functions
c     WKA:    working matrix (NOBSxNBAS) for "working vectors"
c     WKALL:  working matrix (NOBSxNBASxNREP) for all "working vectors"
c             for all perturbed responses.
c     WKAMAX: must be >MAXBAS, it ensures enough memory space for 
c		selecting up to MAXBAS basis functions.
c     PTB: perturbations generated inside Fortran.
c     U: The norms of the original vectors from which the "working
c	  vectors" are derived with the opposite sign of the first
c	  element.
c     H: The first elements of the vectors corresponding to the
c        "working vectors" obtained from the transformed y's 
c     FITS:   array of fits, (NOBSx(MAXBAS-NNULL+1)xNREP)
c     FITSWK: working matrix to store fits temporarily
c     Z: The first k elements are obtained from y after k 
c        transformations, the rest are zeroes.
c     ZK: The first K elements are obtained from the (k+1)th chosen
c         basis after k transformations, the rest are zeroes.
c     YYPTB: perturbed responses, updated upon the selection of a basis
c Modified by Junqing Wu, UCSB on Sep 28, 2009
c########################################################################


      INTEGER I, J, K, NB, IDXTMP(NBAS)
      DOUBLE PRECISION TMP
	  DO 5 K = 1,NBAS
	    IDXTMP(K)=IWK(K)
5	  CONTINUE 

         DO 20 J = 1,nrep
           DO 10 i=1, NOBS
              d(i,j) = PTB(i,j)
              ywk(i) = YYPTB(I,J)
 10        CONTINUE

        CALL DCOPY(NOBS,YWK(1),1,YY(1),1)
        CALL DCOPY(NOBS*NBAS, BAS, 1, BASWK, 1)
        
        CALL DSET(NOBS,0.D0,FITSWK,1)
	    CALL DCOPY(NOBS*WKAMAX,WKALL(1,1,J),1,WKAA(1,1),1)
          CALL PRESELE3(BASWK,NOBS,NBAS,YY,NNULL,ISB,
     &             RSS,IWK,WKAA,U,H,FITSWK,Z,ZK,BAS)
	    CALL DCOPY(NOBS*WKAMAX,WKAA(1,1),1,WKALL(1,1,J),1)
	    CALL DCOPY(NOBS,YY,1,YYPTB(1,J),1)
          CALL DCOPY(NOBS,FITSWK,1,FITS(1,J),1)

		  ICB(J)=ISB
		  
	  DO 17 K = 1,NBAS
	    IWK(K)=IDXTMP(K)
 17   CONTINUE

 20   CONTINUE

          CALL DCOPY(NOBS,Y(1),1,YY(1),1)
          CALL DCOPY(NOBS*NBAS, BAS, 1, BASWK, 1)
          CALL DSET(NOBS,0.D0,FITSWK,1)
          CALL PRESELE3(BASWK,NOBS,NBAS,YY,NNULL,ISB,
     &             RSS,IWK,WKA,U,H,FITSWK,Z,ZK,BAS)


	DO 25 J= 1, NOBS
	  DBAR(J) = 0.D0
	  YBAR(J) = 0.D0
	  DO 23 I = 1, nrep
	    DBAR(J) = DBAR(J) + D(J,I)
	    YBAR(J) = YBAR(J) + FITS(J,I)
23	  CONTINUE
	  DBAR(J)=DBAR(J)/nrep
	  YBAR(J)=YBAR(J)/nrep
25	CONTINUE

	GDFHAT = 0.D0
c    IDFHAT = 0.D0
	DO 40 J = 1, NOBS
	  DO 30 I = 1,nrep
	    GDFHAT = GDFHAT+(D(J,I)-DBAR(J))*
     & FITS(J,I)/(nrep-1)
30	  CONTINUE
40	CONTINUE

C	IDFHAT = (GDFHAT-NNULL)
C	IF (IDFHAT.LE.0.D0) THEN
C		IDFHAT = 0.D0
C	ENDIF
      RETURN
      END
	  
c presele3.f updated Sep 28, 2009 by Junqing Wu, UCSB
c ---------------------------------------------------------------------
      SUBROUTINE PRESELE3(BAS,N,NBAS,Y,NNULL,ISB,RSS,IWK,
     & WKA,U,H,FITS,Z,ZK,BASO)

      INTEGER N, NBAS, NNULL, ISB, IWK(NBAS)
      DOUBLE PRECISION BAS(N,*), Y(*), RSS,WKA(N,*),PREV,CURR
      DOUBLE PRECISION U(*),H(*),FITS(*),Z(*),ZK(*)
      INTEGER I, J, JB
      DOUBLE PRECISION BNEW, BASO(N,*)
      DOUBLE PRECISION DDOT

c#####################################################################
c     BAS: Matrix of basis functions, NOBSxNBAS, the first NNULL 
c          columns are bases of the null space. Its
c          content will be modified due to Householder transformations
c          being applied throughout the basis selection process
c     N:      the number of observations
c     NBAS:   the number of candidate basis functions
c     Y:      (N), the response
c     NNULL:  the number of basis functions which will be forced
c             into the selected subset
c     ISB:    (MAXBAS), the column index of chosen basis functions
c             in BAS
c     RSS:    (MAXBAS), RSS(I) is the residual sum of squares of the LS 
c             fit when regressing Y on the first I chosen basis
c             functions.
c             IWK(I) is the index of whether the I-th basis function
c             in BAS has been chosen or not;  0 -- no, 1 -- yes
c     WKA:    (N,MAXBAS), work array, stores the "work vectors" for
c		contructing Householder matrices
c     U: The norms of the original vectors from which the "working
c	  vectors" are derived with the opposite sign of the first
c	  element.
c     H: The first elements of the vectors corresponding to the
c        "working vectors" obtained from the transformed y's 
c     FITS: (NOBS), vector of fits
c     Z: The first k elements are obtained from y after k 
c        transformations, the rest are zeroes.
c     ZK: The first K elements are obtained from the (k+1)th chosen
c         basis after k transformations, the rest are zeroes.
c     BASO: Same as the initial content of BAS, but its content 
c           will not be changed.
c
c Sep 28, 2009 modified by Junqing Wu, UCSB
c######################################################################


	DO 3 I = NNULL+1, NBAS
      IF(IWK(I).EQ.1) GOTO 3
	  DO 2 J = 1, NNULL
	  CALL APPLYHOUSEHOLDER(N,WKA(1,J),BAS(1,I))
2	  CONTINUE
3	CONTINUE

         ISB = 0 
         PREV = 0.D0
         CURR = 0.D0
	   RSS = 0.D0

      CALL DCOPY(NNULL,Y(1),1,Z(1),1)

         DO 18 JB = NNULL+1,NBAS          
          IF(IWK(JB).EQ.1) GOTO 18
            CURR = (DDOT(N-NNULL,BAS(NNULL+1,JB),1,Y(NNULL+1),1))**2/
     &           DDOT(N-NNULL,BAS(NNULL+1,JB),1,BAS(NNULL+1,JB),1)
            IF(CURR.LE.PREV) GOTO 18
            PREV = CURR
            ISB = JB
 18      CONTINUE
         CALL DCOPY(NNULL,BAS(1,ISB),1,ZK(1),1)
         IWK(ISB) = 1
         CALL HOUSEHOLDER(N,NNULL+1,BAS(1,ISB),WKA(1,NNULL+1),U)
         CALL APPLYHOUSEHOLDER(N,WKA(1,NNULL+1),Y)
         H(NNULL+1) = Y(NNULL+1)
         RSS = DDOT(N-NNULL-1,Y(NNULL+2),1,Y(NNULL+2),1)

        CALL DSET(N-NNULL,0.d0,Z(NNULL+1),1)
        CALL DSET(N-NNULL,0.d0,ZK(NNULL+1),1)
      DO 20 I = 1, NNULL
        CALL APPLYHOUSEHOLDER(N,WKA(1,NNULL+1-I),Z(1))
        CALL APPLYHOUSEHOLDER(N,WKA(1,NNULL+1-I),ZK(1))
 20   CONTINUE
      BNEW = - H(NNULL+1)/U(NNULL+1)
        CALL DAXPY(N,-1.D0,BASO(1,ISB),1,ZK(1),1)
        CALL DAXPY(N,BNEW,ZK(1),1,Z(1),1)
        CALL DCOPY(N,Z(1),1,FITS,1)

      RETURN
      END 

c  libcost3.f updated Mar 9, 2009 by Junqing Wu, UCSB
c----------------------------------------------------------------------
      SUBROUTINE LIBCOST3(NOBS, Y, BAS, NNULL, NBAS, MAXBAS,  
     & IDFHAT, nrep, GDFHAT, TAUHAT, ISB, RSS,IWK, YY, BASWK, WKA,
     & WKALL,WKAMAX,ICB,PTB,U,H,FITS,FITSWK,Z,ZK,YYPTB)

      INTEGER NOBS, NNULL, NBAS, MAXBAS, ISB(*), 
     & IWK(*), nrep, WKAMAX

      DOUBLE PRECISION Y(*), BAS(NOBS,*), IDFHAT(MAXBAS),  
     & RSS(*), YY(*), YYPTB(NOBS,*), BASWK(NOBS,*), WKA(NOBS,*), 
     & GDFHAT(MAXBAS), TAUHAT,tmp1,tmp2,
     & dbar(NOBS), YWK(NOBS), WKALL(NOBS,WKAMAX,nrep), YBAR(NOBS),
     & d(NOBS,nrep),ICB(MAXBAS,nrep),DISB(MAXBAS), FITSWK(NOBS,MAXBAS),
     & PTB(NOBS,nrep),U(*),H(*),FITS(NOBS,MAXBAS,nrep),Z(*),ZK(*)

      REAL RNOR, UNI

c########################################################################
c     NOBS:   number of observations
c     Y:      response.
c     BAS:    NOBSxNBAS, Bases. The first NNULL columns are bases of the 
c             null space
c     NNULL:  Number of basis functions already in the model
c     NBAS:   number of bases
c     MAXBAS: the maximum number of basis functions to be chosen
c     IDFHAT: The estimated IDF values
c     NREP:   Number of repetitions of fitting perturbed responses 
c     GDFHAT: Total Cost for all chosen basis functions at each step
c     TAUHAT: .5*(est of sigma)  
c     ISB:    Index of Selected Bases, (MAXBAS), the column index of
c             chosen bases in BAS. It always include the first NNULL
c             columns.
c     RSS:    residual sum of squares for each model
c     IWK:    working vector, (NBAS), each element indicates if a basis
c		has been chosen, 1-yes, 0-no
c     YY:     working vector for storing perturbed responses
c     BASWK:  working matrix (NOBSxNBAS) for basis functions
c     WKA:    working matrix (NOBSxNBAS) for "working vectors"
c     WKALL:  working matrix (NOBSxNBASxNREP) for all "working vectors"
c             for all perturbed responses.
c     WKAMAX: must be >MAXBAS, it ensures enough memory space for 
c		selecting up to MAXBAS basis functions.
c     ICB: Indices of Chosen Bases, similar to ISB, except for perturbed
c	    responses.
c     PTB: perturbations generated inside Fortran.
c     U: The norms of the original vectors from which the "working
c	  vectors" are derived with the opposite sign of the first
c	  element.
c     H: The first elements of the vectors corresponding to the
c        "working vectors" obtained from the transformed y's 
c     FITS:   array of fits, (NOBSx(MAXBAS-NNULL+1)xNREP)
c     FITSWK: working matrix to store fits temporarily
c     Z: The first k elements are obtained from y after k 
c        transformations, the rest are zeroes.
c     ZK: The first K elements are obtained from the (k+1)th chosen
c         basis after k transformations, the rest are zeroes.
c     YYPTB: perturbed responses after being transformed by the
c	      "working vectors".
c Modified by Junqing Wu, UCSB on Sep 28, 2009
c########################################################################

      INTEGER I, J, K, KK, START, DRAW, INDX(NOBS), IDXTMP(NBAS)
      DOUBLE PRECISION TMP, FITT(NOBS,MAXBAS), DRW(NOBS)
	DOUBLE PRECISION RESID(NOBS,MAXBAS)

	DO 1 K = 1,NBAS
	    IDXTMP(K)=IWK(K)
1	  CONTINUE
 
	START = NNULL
	IF (NNULL.EQ.0) THEN
	  START=1
	ENDIF

      CALL DCOPY(NOBS,Y(1),1,YY(1),1)
      CALL DCOPY(NOBS*NBAS, BAS, 1, BASWK, 1)
      CALL PRESELE2(BASWK,NOBS,NBAS,YY,NNULL,MAXBAS,ISB,
     & RSS,IWK,WKA,U,H,FITSWK,Z,ZK,BAS)
	CALL DCOPY(NOBS*MAXBAS, FITSWK, 1, FITT, 1)
	  DO 2 K = 1,NBAS
	    IWK(K)=IDXTMP(K)
2     CONTINUE

	DO 5 J = 1, MAXBAS
	  DO 3 I = 1, NOBS
	    RESID(I,J) = Y(I) - FITT(I,J)
3	CONTINUE
5	CONTINUE

	DO 21 K = START, MAXBAS
        DO 16 J = 1,nrep
          DO 10 I = 1, NOBS
c		DRAW = INDX(I)
c            d(I,J) = RESID(DRAW,MAXBAS)
		d(I,J) = TAUHAT*dble(rnor())
c		 YWK(I) = FITT(DRAW,K)
            YWK(I) = dble(FITT(I,K)) + dble(d(I,J))
 	 	d(I,J) = YWK(I)
 10        CONTINUE

        CALL DCOPY(NOBS,YWK(1),1,YY(1),1)
        CALL DCOPY(NOBS*NBAS, BAS, 1, BASWK, 1)
        
        DO 12 I= 1,MAXBAS
        CALL DSET(NOBS,0.D0,FITSWK(1,I),1)
12      CONTINUE

        CALL PRESELE2(BASWK,NOBS,NBAS,YY,NNULL,K,ISB,
     & RSS,IWK,WKA,U,H,FITSWK,Z,ZK,BAS)
	  CALL DCOPY(NOBS,YY,1,YYPTB(1,J),1)
	  CALL DCOPY(NOBS*MAXBAS,WKA(1,1),1,WKALL(1,1,J),1)
        CALL DCOPY(NOBS,FITSWK(1,K),1,FITS(1,K,J),1)

        DO 15 i=1,MAXBAS
          DISB(i)=dble(ISB(i))
15      CONTINUE
        CALL DCOPY(MAXBAS,DISB,1,ICB(1,J),1)

	  DO 77 KK = 1,NBAS
	    IWK(KK)=IDXTMP(KK)
77    CONTINUE
		
16	CONTINUE

	DO 18 J= 1, NOBS
	  DBAR(J) = 0.D0
	  YBAR(J) = 0.D0
	  DO 17 I = 1, nrep
	    DBAR(J) = DBAR(J) + D(J,I)
	    YBAR(J) = YBAR(J) + FITS(J,K,I)
17	  CONTINUE
	  DBAR(J)=DBAR(J)/nrep
	  YBAR(J)=YBAR(J)/nrep
18	CONTINUE

	GDFHAT(K) = 0.D0
	DO 20 J = 1, NOBS
	  DO 19 I = 1,nrep
	    GDFHAT(K) = GDFHAT(K)+(D(J,I)-DBAR(J))*
     & FITS(J,K,I)/(nrep-1)
19	  CONTINUE
20	CONTINUE
21    CONTINUE

      CALL DCOPY(NOBS*nrep,d,1,PTB,1)

      CALL DCOPY(NOBS,Y(1),1,YY(1),1)
      CALL DCOPY(NOBS*NBAS, BAS, 1, BASWK, 1)
      CALL PRESELE2(BASWK,NOBS,NBAS,YY,NNULL,MAXBAS,ISB,
     & RSS,IWK,WKA,U,H,FITSWK,Z,ZK,BAS)

      RETURN
      END

	        SUBROUTINE BASCHECK(BASEL,SELEX,N,NUMCOL)
      INTEGER N, NUMCOL, I, J, SELEX(NUMCOL)
      DOUBLE PRECISION BASEL(N,*),S(N),DD,DNRM2
      DO 30 I = 1,NUMCOL-1
        IF (SELEX(I).EQ.1) GOTO 30
        DO 20 J = (I+1),NUMCOL
          IF (SELEX(J).EQ.1) GOTO 20
          CALL DCOPY(N,BASEL(1,J),1,S(1),1)
          CALL DAXPY(N,-1.D0,BASEL(1,I),1,S(1),1)
          DD = DNRM2(N,S(1),1)
          IF (DABS(DD).LT.(1.D-60)) SELEX(J)=1
          CALL DCOPY(N,BASEL(1,J),1,S(1),1)
          CALL DAXPY(N,1.D0,BASEL(1,I),1,S(1),1)
          DD = DNRM2(N,S(1),1)
          IF (DABS(DD).LT.(1.D-60)) SELEX(J)=1
20      CONTINUE
30    CONTINUE
      RETURN
      END

      SUBROUTINE STDZ(X,N,NUMCOL)
      INTEGER N, NUMCOL, I, J
      DOUBLE PRECISION X(N,*), XTEMP(N), S,DNRM2
      DO 30 I = 1,NUMCOL
        CALL DCOPY(N,X(1,I),1,XTEMP(1),1)
	    S = DNRM2(N,XTEMP(1),1)
		IF (DABS(S).LT.(1.D-60)) GOTO 30
        DO 20 J = 1,N
          XTEMP(J)=XTEMP(J)/S
20      CONTINUE
        CALL DCOPY(N,XTEMP(1),1,X(1,I),1)
30    CONTINUE
      RETURN
      END
