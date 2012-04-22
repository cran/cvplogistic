cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     This is a complete package for following 
!     1. SCAD: MM only since adaptive rescaling is not applicable
!     2. MCP: MM(kappa<1/4) and adaptive rescaling
!     
!     Modified from mcpv13.f & scadv11.f
!     Sep 28, 2011
!
!     Version 2
!     SCAD: only the MMCD, along kappa, along lambda and hybrid
!     MCP:  MMCD and Adaptive rescaling, three solution surface
!     March 15, 2012
!
!     Version 3
!     LLA1: using adaptive lasso approach to compute concave solutions
!     LLA2: strictly implement Zou,H., Li,R., 2008 Annals Stat paper
!     For both LLA, only solution path along kappa is implemented.
!     May 23, 2012
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     LEVEL I FUNCTIONS
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C     perform standardization
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine standard(as,sz,z,n,p)
      integer n,p
      double precision as(p),sz(n,p),z(n,p)
      integer j
      do 10000 j=1,p
         as(j)=sqrt( dble(n)/dot_product(z(:,j),z(:,j)) )
         sz(:,j)=as(j)*z(:,j)
10000 continue
      end

C     soft threshold
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine soft(out,m,lambda)
      double precision out,m,lambda
      if (abs(m) .le. lambda) then
         out=0.d0
      else if (m .gt. lambda) then
         out=m-lambda
      else if (m .lt. -lambda) then
         out=m+lambda
      endif
      end
	  
C     convergence check, using L2 norm
c***********************************************************************
      subroutine converge(tag,newz,oldz,lgz,epsilon)
      integer tag,lgz
      double precision newz(lgz),oldz(lgz),epsilon
c     local var
      integer i
      double precision nz(lgz),oz(lgz),sn,so
      tag=1
      do 00006 i=1,lgz
         nz(i)=newz(i)*newz(i)
         oz(i)=oldz(i)*oldz(i)
00006 continue
      sn=sqrt( sum(nz(:)) )
      so=sqrt( sum(oz(:)) )
      if ( abs(sn-so)/(so+0.01) .gt. epsilon) tag=0
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine converge1(tag,newx,oldx,lng,epsilon)
      integer tag,lng
      double precision newx(lng),oldx(lng), epsilon
      integer i
      double precision nnorm(lng),onorm(lng),nsum,osum
      tag=1
      do 00003 i=1,lng
         nnorm(i)=newx(i)*newx(i)
         onorm(i)=oldx(i)*oldx(i)
00003 continue
      nsum=sqrt(sum(nnorm(:)))
      osum=sqrt(sum(onorm(:)))
      if ( abs(nsum-osum)/(osum+0.01) .gt. epsilon ) tag=0
      end

C     convergence check, using L2 norm of both alpha and beta
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine converge2(tag,newx,oldx,lgx,newz,oldz,lgz,epsilon)
      integer tag,lgx,lgz
      double precision newx(lgx),oldx(lgx),newz(lgz),oldz(lgz),epsilon
      integer i
      double precision nx(lgx),ox(lgx),nz(lgz),oz(lgz),sn,so
      tag=1
      do 00005 i=1,lgx
         nx(i)=newx(i)*newx(i)
         ox(i)=oldx(i)*oldx(i)
00005 continue
      do 00006 i=1,lgz
         nz(i)=newz(i)*newz(i)
         oz(i)=oldz(i)*oldz(i)
00006 continue
      sn=sqrt( sum(nx(:)) + sum(nz(:)) )
      so=sqrt( sum(ox(:)) + sum(oz(:)) )
      if ( abs(sn-so)/(so+0.01) .gt. epsilon) tag=0
      end


c     calculate the objective function given solution, kapa,lmda
C     MCP only
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fobj(obj,ab,y,x,z,n,q,p,ka,lmda)
      integer n,q,p
      double precision obj,ab(q+p),y(n),x(n,q),z(n,p),ka,lmda
      integer j,i
      double precision alpha(q),beta(p),eta(n),logl,pen(p),sumpen,t
      alpha(:)=ab(1:q)
      beta(:)=ab((q+1):(q+p))
      call dgemv("N",n,q,1.d0,x,n,alpha,1,0.d0,eta,1)
      call dgemv("N",n,p,1.d0,z,n, beta,1,1.d0,eta,1)
      logl=0.d0
      do 11001 j=1,n
         logl=logl+ y(j)*eta(j) - log( 1.d0+exp(eta(j)) )
11001 continue
c     penalty
      if (ka .eq. 0.d0) then
         pen(:)=abs(beta(:))*lmda
      else
         do 1000 i=1,p
            t=abs(beta(i))
            if (t .ge. lmda/ka) then
               pen(i)=lmda*lmda/(2.d0*ka)
            else
               pen(i)=lmda*t-ka*t*t/2.d0
            endif
 1000    continue
      endif
      sumpen=sum(pen(:))
      obj= -logl/dble(n)+sumpen
      end

c     Given solution, kappa,lmda,
c     determine model size, eigenvalue, aic, bic, objective value
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fmeabo(ms,evgx,aic,bic,obj,ab,y,x,z,n,q,p,ka,lmda)
      integer ms,evgx,n,q,p
      double precision aic,bic,obj,ab(q+p),y(n),x(n,q),z(n,p),ka,lmda
!     .. local arguments
      integer msz,j,nzidx(p),lwork,i,info
      double precision alpha(q),beta(p),eta(n),logl,pi(n),work1(p),
     +     sumpen,t
      alpha(:)=ab(1:q)
      beta(:)=ab((q+1):(q+p))
c     model size
      msz=0
      do 1000 j=1,p
         if (beta(j) .ne. 0.d0) then
            msz=msz+1
            nzidx(msz)=j
         endif
 1000 continue
      ms=q+msz
      end

c     Given solution, kappa,lmda,
c     determine model size, eigenvalue, aic, bic, objective value
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fmeabo2(ms,evgx,aic,bic,obj,ab,y,x,z,n,q,p,ka,lmda)
      integer ms,evgx,n,q,p
      double precision aic,bic,obj,ab(q+p),y(n),x(n,q),z(n,p),ka,lmda
!     .. local arguments
      integer msz,j,nzidx(p),lwork,i,info
      double precision alpha(q),beta(p),eta(n),logl,pi(n),work1(p),
     +     sumpen,t
      double precision, dimension(:,:), allocatable:: xnz
      double precision, dimension(:), allocatable:: anb
      double precision, dimension(:,:), allocatable::qxnz
      double precision, dimension(:),allocatable::w
      double precision, dimension(:),allocatable::work2 
      double precision, dimension(:),allocatable::pen
      alpha(:)=ab(1:q)
      beta(:)=ab((q+1):(q+p))
c     model size
      msz=0
      do 1000 j=1,p
         if (beta(j) .ne. 0.d0) then
            msz=msz+1
            nzidx(msz)=j
         endif
 1000 continue
      ms=q+msz
c     assign non-zero matrix and coefficients
      allocate(xnz(n,ms))
      xnz(:,1:q)=x(:,:)
      xnz(:,(q+1):ms)=z(:,nzidx(1:msz))
      allocate(anb(ms))
      anb(1:q)=alpha(:)
      anb((q+1):ms)=beta(nzidx(1:msz))
!     convexty diagnosis
      if (ms .gt. n) then 
         evgx=0
      else
!     .... X'X 
         allocate(qxnz(ms,ms))
         call dsyrk("U","T",ms,n,1.d0,xnz,n,0.d0,qxnz,ms)
         allocate(w(ms))
         lwork=-1
         call dsyev("N","U",ms,qxnz,ms,w,work1,lwork,info)
         lwork=int(work1(1))
         allocate(work2(lwork))
         call dsyev("N","U",ms,qxnz,ms,w,work2,lwork,info)
         if (w(1) .ge. dble(n)*ka) then 
            evgx=1
         else
            evgx=0
         endif
         deallocate(qxnz)
         deallocate(w)
         deallocate(work2)
      endif
!     .. loglikelihood, pi
      call dgemv("N",n,ms,1.d0,xnz,n,anb,1,0.d0,eta,1)
      logl=0.d0
      do 11001 i=1,n
         pi(i)=1.d0/(1.d0+exp(-eta(i)))
         logl=logl+ y(i)*eta(i) - log( 1.d0+exp(eta(i)) )
11001 continue
c     aic/bic
      aic= -2.d0*logl+ ms*2.d0
      bic= -2.d0*logl+ ms*log(dble(n))
c     objective function
      allocate(pen(msz))
      if (ka .eq. 0.d0) then
         pen(:)=abs( anb((q+1):ms) )*lmda
      else
         do 1002 i=1,msz
            t=abs(anb(i+q))
            if (t .ge. lmda/ka) then
               pen(i)=lmda*lmda/(2.d0*ka)
            else
               pen(i)=lmda*t-ka*t*t/2.d0
            endif
 1002       continue
      endif
      sumpen=sum(pen(:))
      obj= -logl/dble(n)+sumpen
c     .. deallocate 
      deallocate(xnz)
      deallocate(anb)
      deallocate(pen) 
      end

c     Given solution, kappa,lmda,
c     determine model size, eigenvalue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fme(full,evgx,alpha,beta,y,x,z,n,q,p)
      integer full,evgx,n,q,p
      double precision alpha(q),beta(p),y(n),x(n,q),z(n,p)
!     .. local arguments
      integer ms,msz,j
c     model size
      msz=0
      do 1000 j=1,p
         if (beta(j) .ne. 0.d0) then
            msz=msz+1
         endif
 1000 continue
      ms=q+msz
      if (ms .le. n) then
         full=0
      else
         full=1
      endif
      end

!     Calculate the AUC of ROC given cvpy,cvpi
!***********************************************************************
      subroutine aucroc(auc,y,pi,n)
      integer n
      double precision auc,y(n),pi(n)
!     local vars
      integer i,j,n0,grv(n),ltv(n),eqv(n),tmp,l
      double precision rankv(n),u0
      n0=0
!     get the group size, n0, n1
      do 1000 i=1,n
         if (y(i) .eq. 0.d0) n0=n0+1
 1000 continue
!     get the ranks of data
      do 1001 i=1,n
         grv(i)=0
         ltv(i)=0
         eqv(i)=0
         do 1002 j=1,n
            if ( pi(i) .lt. pi(j) ) then
               grv(i)=grv(i)+1
            else if ( pi(i) .eq. pi(j) ) then
               eqv(i)=eqv(i)+1
            else
               ltv(i)=ltv(i)+1
            end if
 1002    continue
         grv(i)=n-grv(i)
         ltv(i)=1+ltv(i)
         if ( eqv(i) .eq. 1) then
            rankv(i)=dble( grv(i) )
         else
            tmp=0
            do 1003 l=ltv(i),grv(i)
               tmp=tmp+l
 1003       continue
            rankv(i)=dble(tmp)/dble(eqv(i))
         end if
 1001 continue
!     calculate u statistics
      u0= - dble( n0*(n0+1)/2 )
      do 1004 i=1,n
         if ( y(i) .eq. 0.d0) u0=u0+rankv(i)
1004  continue
      auc=dble(u0)/dble( n0*(n-n0) )
      if (auc .lt. 0.5) auc=1.d0-auc
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     LEVEL II FUNCTIONS
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c     calculate initial value of alpha,beta, byproduct:lambda_max
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine maxbi(lmdamax,alpha,beta,eta,y,x,z,n,q,p,
     +     epsilon,maxit)
      integer n,q,p,maxit
      double precision lmdamax,alpha(q),beta(p), eta(n),
     +    y(n),x(n,q),z(n,p),epsilon
      integer count,i,tag
      double precision estp,alphaold(q),etax(n),pi(n),r(n),tmp(p)
c     ... assign initial value of alpha, beta
      estp=sum(y(:))/dble(n)
      alpha(1)=log(estp)-log(1.d0-estp)
      if (q .lt. 1) alpha(2:q)=0.d0
      beta(:)=0.d0
c     ... start of iteration
      count=0
 1204 continue
      alphaold(:)=alpha(:)
      call dgemv("N",n,q,1.d0,x,n,alpha,1,0.d0,etax,1)
      do 1205 i=1,n
         pi(i)=1.d0/(1.d0+exp(-etax(i)))
c     if (pi(i) .lt. 0.0001) pi(i)=0.d0
c     if (pi(i) .gt. 0.9999) pi(i)=1.d0
         r(i)=y(i)-pi(i)
 1205 continue
      call dgemv("T",n,q,4.d0/dble(n),x,n,r,1,1.d0,alpha,1)
      count=count+1
      call converge1(tag,alpha,alphaold,q,epsilon)
      if (tag .eq. 1) go to 1206
      if (count .ge. maxit) call rexit("diverge of alpha initials!\n")
      if (count .lt. maxit) go to 1204
 1206 continue
c     ... end of iteration
c     ... compute lambda_max
      call dgemv("n",n,q,1.d0,x,n,alpha,1,0.d0,eta,1)
      do 1207 i=1,n
         pi(i)=1.d0/(1.d0+exp(-eta(i)))
c     if (pi(i) .lt. 0.0001) pi(i)=0.d0
c     if (pi(i) .gt. 0.9999) pi(i)=1.d0
         r(i)=y(i)-pi(i)
 1207 continue
      do 1208 i=1,p
         tmp(i)=abs(dot_product(z(:,i),r))/dble(n)
 1208 continue
      lmdamax=maxval(tmp)
      end

c     Lasso solution
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine lasso(alpha,beta,eta,lmda,y,x,z,n,q,p,
     +     epsilon,maxit)
      integer n,q,p,maxit
      double precision alpha(q),beta(p),eta(n),lmda,y(n),x(n,q),
     +     z(n,p),epsilon
      integer count,j,i,tag
      double precision alphaold(q),betaold(p),pi(n),r(n),m,numor
c     ... start of iteration
      count=0
 1210 continue
      alphaold(:)=alpha(:)
      betaold(:)=beta(:)
!     ... update beta part
      do 1211 j=1,p
         do 1212 i=1,n
            pi(i)=1.d0/(1.d0+exp(-eta(i)))
            if (pi(i) .lt. 0.0001) pi(i)=0.d0
            if (pi(i) .gt. 0.9999) pi(i)=1.d0
            r(i)=y(i)-pi(i)
 1212    continue
!     ...... m=beta_{j}/4+Z'_{j}r/n
         m=beta(j)/4.d0 + dot_product(z(:,j),r)/dble(n)
         call soft(numor,m,lmda)
         beta(j)=4.d0*numor
!     ...... eta^{s+1}=eta^{s}+z_{j}*(b^{s+1}-b^{s})
         eta(:)=eta(:) + z(:,j)*( beta(j) - betaold(j) )
 1211 continue
!     ... update alpha part
      do 1215 i=1,n
         pi(i)=1.d0/(1.d0+exp(-eta(i)))
         if (pi(i) .lt. 0.0001) pi(i)=0.d0
         if (pi(i) .gt. 0.9999) pi(i)=1.d0
         r(i)=y(i)-pi(i)
 1215 continue
      call dgemv("T",n,q,4.d0/dble(n),x,n,r,1,1.d0,alpha,1)
!     ...... eta^{s+1}=eta^{s}+z_{j}*(b^{s+1}-b^{s})
      do 1217 j=1,q
         eta(:)=eta(:) + x(:,j)*( alpha(j) - alphaold(j) )
 1217 continue
      count=count+1
      call converge2(tag,alpha,alphaold,q,beta,betaold,p,epsilon)
      if (tag .eq. 1) go to 1216
      if (count .ge. maxit) call rexit("Lasso path diverges! \n")
      if (count .lt. maxit) go to 1210
 1216 continue
c     ... end of iteration
      end


!***********************************************************************
!     SCAD part: MMCD
!***********************************************************************

c     SCAD solution
c***********************************************************************
      subroutine scad(alpha,beta,eta,lmda,ka,y,x,z,n,q,p,
     +     epsilon,maxit)
      integer n,q,p,maxit
      double precision alpha(q),beta(p),eta(n),lmda,ka,y(n),x(n,q),
     +     z(n,p),epsilon
      integer count,j,i,tag
      double precision alphaold(q),betaold(p),pi(n),r(n),m,numor,gamma,
     +     newlmda
c     ... begin of iteration
      gamma=1.d0/ka
      count=0
 1000 continue
      alphaold(:)=alpha(:)
      betaold(:)=beta(:)
!     ... update beta part
      do 1001 j=1,p
         do 1003 i=1,n
            pi(i)=1.d0/(1.d0+exp(-eta(i)))
            if (pi(i) .lt. 0.0001) pi(i)=0.d0
            if (pi(i) .gt. 0.9999) pi(i)=1.d0
            r(i)=y(i)-pi(i)
 1003    continue
!     ...... m=beta_{j}/4+Z'_{j}r/n
         m=beta(j)/4.d0 + dot_product(z(:,j),r)/dble(n)
         if (abs(m) .lt. 1.25*lmda) then
            call soft(numor,m,lmda)
            beta(j)=4.d0*numor
         else if ( abs(m) .ge. 0.25*lmda*gamma ) then
            beta(j)=4.d0*m
         else 
            newlmda=lmda*gamma/(gamma-1.d0)
            call soft(numor,m,newlmda)
            beta(j)=4.d0*(gamma-1.d0)*numor/(gamma-5.d0)
         end if
!     ...... eta^{s+1}=eta^{s}+z_{j}*(b^{s+1}-b^{s})
         eta(:)=eta(:) + z(:,j)*( beta(j) - betaold(j) )
 1001 continue
!     ... update alpha part
      do 1215 i=1,n
         pi(i)=1.d0/(1.d0+exp(-eta(i)))
         if (pi(i) .lt. 0.0001) pi(i)=0.d0
         if (pi(i) .gt. 0.9999) pi(i)=1.d0
         r(i)=y(i)-pi(i)
 1215 continue
      call dgemv("T",n,q,4.d0/dble(n),x,n,r,1,1.d0,alpha,1)
!     ...... eta^{s+1}=eta^{s}+z_{jl}*(b^{s+1}-b^{s})
      do 1217 j=1,q
         eta(:)=eta(:) + x(:,j)*( alpha(j) - alphaold(j) )
 1217 continue
      count=count+1
      call converge2(tag,alpha,alphaold,q,beta,betaold,p,epsilon)
      if (tag .eq. 1) go to 1006
      if (count .ge. maxit) call rexit("Diverge at kappa & lambda! \n")
      if (count .lt. maxit) go to 1000
 1006 continue
c     ... end of iteration
      end



c     SCAD solution, update non-zero only
c***********************************************************************
      subroutine scad2(alpha,beta,eta,upidx,lmda,ka,y,x,z,n,q,p,
     +     epsilon,maxit)
      integer n,q,p,maxit,upidx(p)
      double precision alpha(q),beta(p),eta(n),lmda,ka,y(n),x(n,q),
     +     z(n,p),epsilon
      integer count,j,i,tag
      double precision alphaold(q),betaold(p),pi(n),r(n),m,numor,gamma,
     +     newlmda
c     ... begin of iteration
      gamma=1.d0/ka
      count=0
 1000 continue
      alphaold(:)=alpha(:)
      betaold(:)=beta(:)
!     ... update beta part, non-zero lasso only
      do 1001 j=1,p
         if (upidx(j) .eq. 1) then 
            do 1003 i=1,n
               pi(i)=1.d0/(1.d0+exp(-eta(i)))
               if (pi(i) .lt. 0.0001) pi(i)=0.d0
               if (pi(i) .gt. 0.9999) pi(i)=1.d0
               r(i)=y(i)-pi(i)
 1003       continue
!     ...... m=beta_{j}/4+Z'_{j}r/n
            m=beta(j)/4.d0 + dot_product(z(:,j),r)/dble(n)
            if (abs(m) .lt. 1.25*lmda) then
               call soft(numor,m,lmda)
               beta(j)=4.d0*numor
            else if ( abs(m) .ge. 0.25*lmda*gamma ) then
               beta(j)=4.d0*m
            else
               newlmda=lmda*gamma/(gamma-1)
               call soft(numor,m,newlmda)
               beta(j)=4.d0*(gamma-1)*numor/(gamma-5)
            end if
!     ...... eta^{s+1}=eta^{s}+z_{j}*(b^{s+1}-b^{s})
            eta(:)=eta(:) + z(:,j)*( beta(j) - betaold(j) )
         endif
 1001 continue
!     ... update alpha part
      do 1215 i=1,n
         pi(i)=1.d0/(1.d0+exp(-eta(i)))
         if (pi(i) .lt. 0.0001) pi(i)=0.d0
         if (pi(i) .gt. 0.9999) pi(i)=1.d0
         r(i)=y(i)-pi(i)
 1215 continue
      call dgemv("T",n,q,4.d0/dble(n),x,n,r,1,1.d0,alpha,1)
!     ...... eta^{s+1}=eta^{s}+z_{jl}*(b^{s+1}-b^{s})
      do 1217 j=1,q
         eta(:)=eta(:) + x(:,j)*( alpha(j) - alphaold(j) )
 1217 continue
      count=count+1
      call converge2(tag,alpha,alphaold,q,beta,betaold,p,epsilon)
      if (tag .eq. 1) go to 1006
      if (count .ge. maxit) call rexit("Diverge at kappa & lambda! \n")
      if (count .lt. maxit) go to 1000
 1006 continue
c     ... end of iteration
      end

     

!***********************************************************************
!     MCP part: MMCD
!***********************************************************************

c     MCP solution
c***********************************************************************
      subroutine mcp(alpha,beta,eta,lmda,ka,y,x,z,n,q,p,
     +     epsilon,maxit)
      integer n,q,p,maxit
      double precision alpha(q),beta(p),eta(n),lmda,ka,y(n),x(n,q),
     +     z(n,p),epsilon
      integer count,j,i,tag
      double precision alphaold(q),betaold(p),pi(n),r(n),m,numor
c     ... begin of iteration
      count=0
 1000 continue
      alphaold(:)=alpha(:)
      betaold(:)=beta(:)
!     ... update beta part
      do 1001 j=1,p
         do 1003 i=1,n
            pi(i)=1.d0/(1.d0+exp(-eta(i)))
            if (pi(i) .lt. 0.0001) pi(i)=0.d0
            if (pi(i) .gt. 0.9999) pi(i)=1.d0
            r(i)=y(i)-pi(i)
 1003    continue
!     ...... m=beta_{j}/4+Z'_{j}r/n
         m=beta(j)/4.d0 + dot_product(z(:,j),r)/dble(n)
         if (abs(m) .lt. 0.25*lmda/ka) then
            call soft(numor,m,lmda)
            beta(j)=numor/(0.25-ka)
         else
            beta(j)=4.d0*m
         end if
!     ...... eta^{s+1}=eta^{s}+z_{j}*(b^{s+1}-b^{s})
         eta(:)=eta(:) + z(:,j)*( beta(j) - betaold(j) )
 1001 continue
!     ... update alpha part
      do 1215 i=1,n
         pi(i)=1.d0/(1.d0+exp(-eta(i)))
         if (pi(i) .lt. 0.0001) pi(i)=0.d0
         if (pi(i) .gt. 0.9999) pi(i)=1.d0
         r(i)=y(i)-pi(i)
 1215 continue
      call dgemv("T",n,q,4.d0/dble(n),x,n,r,1,1.d0,alpha,1)
!     ...... eta^{s+1}=eta^{s}+z_{jl}*(b^{s+1}-b^{s})
      do 1217 j=1,q
         eta(:)=eta(:) + x(:,j)*( alpha(j) - alphaold(j) )
 1217 continue
      count=count+1
      call converge2(tag,alpha,alphaold,q,beta,betaold,p,epsilon)
      if (tag .eq. 1) go to 1006
      if (count .ge. maxit) call rexit("Diverge at kappa & lambda! \n")
      if (count .lt. maxit) go to 1000
 1006 continue
c     ... end of iteration
      end


c     MCP solution, update non-zero only
c***********************************************************************
      subroutine mcp2(alpha,beta,eta,upidx,lmda,ka,y,x,z,n,q,p,
     +     epsilon,maxit)
      integer n,q,p,maxit,upidx(p)
      double precision alpha(q),beta(p),eta(n),lmda,ka,y(n),x(n,q),
     +     z(n,p),epsilon
      integer count,j,i,tag
      double precision alphaold(q),betaold(p),pi(n),r(n),m,numor
c     ... begin of iteration
      count=0
 1000 continue
      alphaold(:)=alpha(:)
      betaold(:)=beta(:)
!     ... update beta part, non-zero lasso only
      do 1001 j=1,p
         if (upidx(j) .eq. 1) then 
            do 1003 i=1,n
               pi(i)=1.d0/(1.d0+exp(-eta(i)))
               if (pi(i) .lt. 0.0001) pi(i)=0.d0
               if (pi(i) .gt. 0.9999) pi(i)=1.d0
               r(i)=y(i)-pi(i)
 1003       continue
!     ...... m=beta_{j}/4+Z'_{j}r/n
            m=beta(j)/4.d0 + dot_product(z(:,j),r)/dble(n)
            if (abs(m) .lt. 0.25*lmda/ka) then
               call soft(numor,m,lmda)
               beta(j)=numor/(0.25-ka)
            else
               beta(j)=4.d0*m
            end if
!     ...... eta^{s+1}=eta^{s}+z_{j}*(b^{s+1}-b^{s})
            eta(:)=eta(:) + z(:,j)*( beta(j) - betaold(j) )
         endif
 1001 continue
!     ... update alpha part
      do 1215 i=1,n
         pi(i)=1.d0/(1.d0+exp(-eta(i)))
         if (pi(i) .lt. 0.0001) pi(i)=0.d0
         if (pi(i) .gt. 0.9999) pi(i)=1.d0
         r(i)=y(i)-pi(i)
 1215 continue
      call dgemv("T",n,q,4.d0/dble(n),x,n,r,1,1.d0,alpha,1)
!     ...... eta^{s+1}=eta^{s}+z_{jl}*(b^{s+1}-b^{s})
      do 1217 j=1,q
         eta(:)=eta(:) + x(:,j)*( alpha(j) - alphaold(j) )
 1217 continue
      count=count+1
      call converge2(tag,alpha,alphaold,q,beta,betaold,p,epsilon)
      if (tag .eq. 1) go to 1006
      if (count .ge. maxit) call rexit("Diverge at kappa & lambda! \n")
      if (count .lt. maxit) go to 1000
 1006 continue
c     ... end of iteration
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     LEVEL III FUNCTIONS
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!***********************************************************************
!     SCAD
!***********************************************************************

c     solution path along kappa
c***********************************************************************
      subroutine scadkapa(olmdas,okas,ocoef,oaic,obic,oobj,odf,oevidx,
     +     y,x,z,n,q,p,nka,maxka,nlmda,minlmda,epsilon,maxit)
      integer n,q,p,nka,nlmda,maxit,odf(nka*nlmda),oevidx(nka*nlmda)
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),oaic(nka*nlmda),obic(nka*nlmda),
     +     oobj(nka*nlmda),y(n),x(n,q),z(n,p),
     +     maxka,minlmda,epsilon
      integer i,j
      double precision as(p),sz(n,p),kas(nka),unitka,unitlmda,
     +     lmdas(nlmda),alpha(q),beta(p),eta(n),
     +     inia(q,nlmda),inib(p,nlmda),etamat(n,nlmda)
c     standardization of Z
      call standard(as,sz,z,n,p)
c     calculate kappas
      if (nka .eq. 1) then
         kas(1)=0.d0
      else
         unitka=maxka/dble(nka-1)
         do 10000 i=1,nka
            kas(i)=dble(i-1)*unitka
10000    continue
      endif
c     calculate lambdas
      unitlmda=log(minlmda)/dble(nlmda-1)
      call maxbi(lmdas(1),alpha,beta,eta,y,x,sz,n,q,p,
     +     epsilon,maxit)
      do 10001 i=2,nlmda
         lmdas(i)=lmdas(1)*exp(unitlmda*dble(i-1))
10001 continue
c     Lasso solution path as initials
      inia(:,1)=alpha(:)
      inib(:,1)=beta(:)
      etamat(:,1)=eta(:)
      do 10002 i=2,nlmda
         call lasso(alpha,beta,eta,lmdas(i),y,x,sz,n,q,p,
     +        epsilon,maxit)
         inia(:,i)=alpha(:)
         inib(:,i)=beta(:)
         etamat(:,i)=eta(:)
10002 continue
c     SCAD solution path along kappa
      if (nka .eq. 1) then
         olmdas(:)=lmdas(:)
         okas(1:nlmda)=kas(1)
         ocoef(1:q,:)=inia(:,:)
         ocoef((q+1):(q+p),:)=inib(:,:)
      else
         do 1003 i=1,nlmda
            olmdas(((i-1)*nka+1):(i*nka))=lmdas(i)
            okas(((i-1)*nka+1):(i*nka))=kas(:)
 1003    continue
!     ... i=1 case
         i=1
         alpha(:)=inia(:,i)
         beta(:) =inib(:,i)
         eta(:)  =etamat(:,i)
         do 1005 j=1,nka
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1005    continue
         do 1004 i=2,nlmda
            alpha(:)=inia(:,i)
            beta(:) =inib(:,i)
            eta(:)  =etamat(:,i)
c     ...... j=1 case
            j=1
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
            do 1006 j=2,nka
               call scad(alpha,beta,eta,lmdas(i),kas(j),y,x,sz,
     +              n,q,p,epsilon,maxit)
               ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
               ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1006       continue
 1004    continue
      endif
c     compute model size,eigenvalue,aic,bic, objective function
      do 10 i=1,(nka*nlmda)
         call fmeabo(odf(i),oevidx(i),oaic(i),obic(i),oobj(i),
     +        ocoef(:,i),y,x,sz,n,q,p,okas(i),olmdas(i))
 10   continue
!     change the penalized ocoefficients back
      do 10005 i=1,p
         ocoef(q+i,:)=as(i)*ocoef(q+i,:)
10005 continue
      end



c     hybrid path, update the non-zero lasso only
c***********************************************************************
      subroutine scadkapa2(olmdas,okas,ocoef,oaic,obic,oobj,odf,oevidx,
     +     y,x,z,n,q,p,nka,maxka,nlmda,minlmda,epsilon,maxit)
      integer n,q,p,nka,nlmda,maxit,odf(nka*nlmda),oevidx(nka*nlmda)
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),oaic(nka*nlmda),obic(nka*nlmda),
     +     oobj(nka*nlmda),y(n),x(n,q),z(n,p),
     +     maxka,minlmda,epsilon
      integer i,j,k,upidx(p)
      double precision as(p),sz(n,p),kas(nka),unitka,unitlmda,
     +     lmdas(nlmda),alpha(q),beta(p),eta(n),
     +     inia(q,nlmda),inib(p,nlmda),etamat(n,nlmda)
c     standardization of Z
      call standard(as,sz,z,n,p)
c     calculate kappas
      if (nka .eq. 1) then
         kas(1)=0.d0
      else
         unitka=maxka/dble(nka-1)
         do 10000 i=1,nka
            kas(i)=dble(i-1)*unitka
10000    continue
      endif
c     calculate lambdas
      unitlmda=log(minlmda)/dble(nlmda-1)
      call maxbi(lmdas(1),alpha,beta,eta,y,x,sz,n,q,p,
     +     epsilon,maxit)
      do 10001 i=2,nlmda
         lmdas(i)=lmdas(1)*exp(unitlmda*dble(i-1))
10001 continue
c     Lasso solution path as initials
      inia(:,1)=alpha(:)
      inib(:,1)=beta(:)
      etamat(:,1)=eta(:)
      do 10002 i=2,nlmda
         call lasso(alpha,beta,eta,lmdas(i),y,x,sz,n,q,p,
     +        epsilon,maxit)
         inia(:,i)=alpha(:)
         inib(:,i)=beta(:)
         etamat(:,i)=eta(:)
10002 continue
c     SCAD solution path along kappa
      if (nka .eq. 1) then
         olmdas(:)=lmdas(:)
         okas(1:nlmda)=kas(1)
         ocoef(1:q,:)=inia(:,:)
         ocoef((q+1):(q+p),:)=inib(:,:)
      else
         do 1003 i=1,nlmda
            olmdas(((i-1)*nka+1):(i*nka))=lmdas(i)
            okas(((i-1)*nka+1):(i*nka))=kas(:)
 1003    continue
!     ... i=1 case
         i=1
         alpha(:)=inia(:,i)
         beta(:) =inib(:,i)
         eta(:)  =etamat(:,i)
         do 1005 j=1,nka
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1005    continue
         do 1004 i=2,nlmda
            alpha(:)=inia(:,i)
            beta(:) =inib(:,i)
            eta(:)  =etamat(:,i)
c     ...... j=1 case
            j=1
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
            upidx(:)=1
            do 9000 k=1,p
               if (beta(k) .eq. 0.d0) upidx(k)=0
 9000       continue
            do 1006 j=2,nka
               call scad2(alpha,beta,eta,upidx,lmdas(i),kas(j),y,x,sz,
     +              n,q,p,epsilon,maxit)
               ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
               ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1006       continue
 1004    continue
      endif
c     compute model size,eigenvalue,aic,bic, objective function
      do 10 i=1,(nka*nlmda)
         call fmeabo(odf(i),oevidx(i),oaic(i),obic(i),oobj(i),
     +        ocoef(:,i),y,x,sz,n,q,p,okas(i),olmdas(i))
 10   continue
!     change the penalized ocoefficients back
      do 10005 i=1,p
         ocoef(q+i,:)=as(i)*ocoef(q+i,:)
10005 continue
      end


c     solution path along lambda
c***********************************************************************
      subroutine scadlmda(olmdas,okas,ocoef,oaic,obic,oobj,odf,oevidx,
     +     y,x,z,n,q,p,nka,maxka,nlmda,minlmda,epsilon,maxit)
      integer n,q,p,nka,nlmda,maxit,odf(nka*nlmda),oevidx(nka*nlmda)
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),oaic(nka*nlmda),obic(nka*nlmda),
     +     oobj(nka*nlmda),y(n),x(n,q),z(n,p),
     +     maxka,minlmda,epsilon
      integer i,j
      double precision as(p),sz(n,p),kas(nka),unitka,unitlmda,
     +     lmdas(nlmda),inia(q),inib(p),inieta(n),
     +     alpha(q),beta(p),eta(n)
c     standardization of Z
      call standard(as,sz,z,n,p)
c     calculate kappas
      if (nka .eq. 1) then
         kas(1)=0.d0
      else
         unitka=maxka/dble(nka-1)
         do 10000 i=1,nka
            kas(i)=dble(i-1)*unitka
10000    continue
      endif
c     calculate lambdas
      unitlmda=log(minlmda)/dble(nlmda-1)
      call maxbi(lmdas(1),inia,inib,inieta,y,x,sz,n,q,p,
     +     epsilon,maxit)
      do 10001 i=2,nlmda
         lmdas(i)=lmdas(1)*exp(unitlmda*dble(i-1))
10001 continue
c     output lambdas and kappas
      do 10002 i=1,nlmda
         olmdas( ((i-1)*nka+1):(i*nka) )=lmdas(i)
         okas( ((i-1)*nka+1):(i*nka) )=kas(:)
10002 continue
c     Lasso solution path
      if (nka .eq. 1) then
         alpha(:)=inia(:)
         beta(:)=inib(:)
         eta(:)=inieta(:)
         ocoef(1:q,1)=alpha(:)
         ocoef((q+1):(q+p),1)=beta(:)
         do 10004 i=2,nlmda
            call lasso(alpha,beta,eta,lmdas(i),y,x,sz,n,q,p,
     +           epsilon,maxit)
            ocoef(1:q,i)=alpha(:)
            ocoef((q+1):(q+p),i)=beta(:)
10004    continue
c     SCAD solution path along lambda given kappa
      else
         alpha(:)=inia(:)
         beta(:)=inib(:)
         eta(:)=inieta(:)
!     i=1 case
!     .. j=1 case
         j=1
         ocoef(1:q,        1+(j-1)*nka)=alpha(:)
         ocoef((q+1):(q+p),1+(j-1)*nka)=beta(:)
         do 1005 j=2,nlmda
            call lasso(alpha,beta,eta,lmdas(j),y,x,sz,n,q,p,
     +           epsilon,maxit)
            ocoef(1:q,        1+(j-1)*nka)=alpha(:)
            ocoef((q+1):(q+p),1+(j-1)*nka)=beta(:)
 1005    continue
         do 10005 i=2,nka
            alpha(:)=inia(:)
            beta(:)=inib(:)
            eta(:)=inieta(:)
            ocoef(1:q,        i)=alpha(:)
            ocoef((q+1):(q+p),i)=beta(:)
            do 10006 j=2,nlmda
               call scad(alpha,beta,eta,lmdas(j),kas(i),y,x,sz,
     +              n,q,p,epsilon,maxit)
               ocoef(1:q,        i+(j-1)*nka)=alpha(:)
               ocoef((q+1):(q+p),i+(j-1)*nka)=beta(:)
10006       continue
10005    continue
      endif
c     compute model size,eigenvalue,aic,bic, objective function
      do 10 i=1,(nka*nlmda)
         call fmeabo(odf(i),oevidx(i),oaic(i),obic(i),oobj(i),
     +        ocoef(:,i),y,x,sz,n,q,p,okas(i),olmdas(i))
 10   continue
!     change the penalized ocoefficients back
      do 10008 i=1,p
         ocoef(q+i,:)=as(i)*ocoef(q+i,:)
10008 continue
      end

!***********************************************************************
!     MCP
!***********************************************************************

c     solution path along kappa
c***********************************************************************
      subroutine mcpkapa(olmdas,okas,ocoef,oaic,obic,oobj,odf,oevidx,
     +     y,x,z,n,q,p,nka,maxka,nlmda,minlmda,epsilon,maxit)
      integer n,q,p,nka,nlmda,maxit,odf(nka*nlmda),oevidx(nka*nlmda)
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),oaic(nka*nlmda),obic(nka*nlmda),
     +     oobj(nka*nlmda),y(n),x(n,q),z(n,p),
     +     maxka,minlmda,epsilon
      integer i,j
      double precision as(p),sz(n,p),kas(nka),unitka,unitlmda,
     +     lmdas(nlmda),alpha(q),beta(p),eta(n),
     +     inia(q,nlmda),inib(p,nlmda),etamat(n,nlmda)
c     standardization of Z
      call standard(as,sz,z,n,p)
c     calculate kappas
      if (nka .eq. 1) then
         kas(1)=0.d0
      else
         unitka=maxka/dble(nka-1)
         do 10000 i=1,nka
            kas(i)=dble(i-1)*unitka
10000    continue
      endif
c     calculate lambdas
      unitlmda=log(minlmda)/dble(nlmda-1)
      call maxbi(lmdas(1),alpha,beta,eta,y,x,sz,n,q,p,
     +     epsilon,maxit)
      do 10001 i=2,nlmda
         lmdas(i)=lmdas(1)*exp(unitlmda*dble(i-1))
10001 continue
c     Lasso solution path as initials
      inia(:,1)=alpha(:)
      inib(:,1)=beta(:)
      etamat(:,1)=eta(:)
      do 10002 i=2,nlmda
         call lasso(alpha,beta,eta,lmdas(i),y,x,sz,n,q,p,
     +        epsilon,maxit)
         inia(:,i)=alpha(:)
         inib(:,i)=beta(:)
         etamat(:,i)=eta(:)
10002 continue
c     MCP solution path along kappa
      if (nka .eq. 1) then
         olmdas(:)=lmdas(:)
         okas(1:nlmda)=kas(1)
         ocoef(1:q,:)=inia(:,:)
         ocoef((q+1):(q+p),:)=inib(:,:)
      else
         do 1003 i=1,nlmda
            olmdas(((i-1)*nka+1):(i*nka))=lmdas(i)
            okas(((i-1)*nka+1):(i*nka))=kas(:)
 1003    continue
!     ... i=1 case
         i=1
         alpha(:)=inia(:,i)
         beta(:) =inib(:,i)
         eta(:)  =etamat(:,i)
         do 1005 j=1,nka
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1005    continue
         do 1004 i=2,nlmda
            alpha(:)=inia(:,i)
            beta(:) =inib(:,i)
            eta(:)  =etamat(:,i)
c     ...... j=1 case
            j=1
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
            do 1006 j=2,nka
               call mcp(alpha,beta,eta,lmdas(i),kas(j),y,x,sz,
     +              n,q,p,epsilon,maxit)
               ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
               ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1006       continue
 1004    continue
      endif
c     compute model size,eigenvalue,aic,bic, objective function
      do 10 i=1,(nka*nlmda)
         call fmeabo(odf(i),oevidx(i),oaic(i),obic(i),oobj(i),
     +        ocoef(:,i),y,x,sz,n,q,p,okas(i),olmdas(i))
 10   continue
!     change the penalized ocoefficients back
      do 10005 i=1,p
         ocoef(q+i,:)=as(i)*ocoef(q+i,:)
10005 continue
      end


c     Hybrid path, update the non-zero lasso only
c***********************************************************************
      subroutine mcpkapa2(olmdas,okas,ocoef,oaic,obic,oobj,odf,oevidx,
     +     y,x,z,n,q,p,nka,maxka,nlmda,minlmda,epsilon,maxit)
      integer n,q,p,nka,nlmda,maxit,odf(nka*nlmda),oevidx(nka*nlmda)
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),oaic(nka*nlmda),obic(nka*nlmda),
     +     oobj(nka*nlmda),y(n),x(n,q),z(n,p),
     +     maxka,minlmda,epsilon
      integer i,j,k,upidx(p)
      double precision as(p),sz(n,p),kas(nka),unitka,unitlmda,
     +     lmdas(nlmda),alpha(q),beta(p),eta(n),
     +     inia(q,nlmda),inib(p,nlmda),etamat(n,nlmda)
c     standardization of Z
      call standard(as,sz,z,n,p)
c     calculate kappas
      if (nka .eq. 1) then
         kas(1)=0.d0
      else
         unitka=maxka/dble(nka-1)
         do 10000 i=1,nka
            kas(i)=dble(i-1)*unitka
10000    continue
      endif
c     calculate lambdas
      unitlmda=log(minlmda)/dble(nlmda-1)
      call maxbi(lmdas(1),alpha,beta,eta,y,x,sz,n,q,p,
     +     epsilon,maxit)
      do 10001 i=2,nlmda
         lmdas(i)=lmdas(1)*exp(unitlmda*dble(i-1))
10001 continue
c     Lasso solution path as initials
      inia(:,1)=alpha(:)
      inib(:,1)=beta(:)
      etamat(:,1)=eta(:)
      do 10002 i=2,nlmda
         call lasso(alpha,beta,eta,lmdas(i),y,x,sz,n,q,p,
     +        epsilon,maxit)
         inia(:,i)=alpha(:)
         inib(:,i)=beta(:)
         etamat(:,i)=eta(:)
10002 continue
c     MCP solution path along kappa
      if (nka .eq. 1) then
         olmdas(:)=lmdas(:)
         okas(1:nlmda)=kas(1)
         ocoef(1:q,:)=inia(:,:)
         ocoef((q+1):(q+p),:)=inib(:,:)
      else
         do 1003 i=1,nlmda
            olmdas(((i-1)*nka+1):(i*nka))=lmdas(i)
            okas(((i-1)*nka+1):(i*nka))=kas(:)
 1003    continue
!     ... i=1 case
         i=1
         alpha(:)=inia(:,i)
         beta(:) =inib(:,i)
         eta(:)  =etamat(:,i)
         do 1005 j=1,nka
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1005    continue
         do 1004 i=2,nlmda
            alpha(:)=inia(:,i)
            beta(:) =inib(:,i)
            eta(:)  =etamat(:,i)
c     ...... j=1 case
            j=1
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
            upidx(:)=1
            do 9000 k=1,p
               if (beta(k) .eq. 0.d0) upidx(k)=0
 9000       continue
            do 1006 j=2,nka
               call mcp2(alpha,beta,eta,upidx,lmdas(i),kas(j),y,x,sz,
     +              n,q,p,epsilon,maxit)
               ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
               ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1006       continue
 1004    continue
      endif
c     compute model size,eigenvalue,aic,bic, objective function
      do 10 i=1,(nka*nlmda)
         call fmeabo(odf(i),oevidx(i),oaic(i),obic(i),oobj(i),
     +        ocoef(:,i),y,x,sz,n,q,p,okas(i),olmdas(i))
 10   continue
!     change the penalized ocoefficients back
      do 10005 i=1,p
         ocoef(q+i,:)=as(i)*ocoef(q+i,:)
10005 continue
      end


c     solution path along lambda
c***********************************************************************
      subroutine mcplmda(olmdas,okas,ocoef,oaic,obic,oobj,odf,oevidx,
     +     y,x,z,n,q,p,nka,maxka,nlmda,minlmda,epsilon,maxit)
      integer n,q,p,nka,nlmda,maxit,odf(nka*nlmda),oevidx(nka*nlmda)
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),oaic(nka*nlmda),obic(nka*nlmda),
     +     oobj(nka*nlmda),y(n),x(n,q),z(n,p),
     +     maxka,minlmda,epsilon
      integer i,j
      double precision as(p),sz(n,p),kas(nka),unitka,unitlmda,
     +     lmdas(nlmda),inia(q),inib(p),inieta(n),
     +     alpha(q),beta(p),eta(n)
c     standardization of Z
      call standard(as,sz,z,n,p)
c     calculate kappas
      if (nka .eq. 1) then
         kas(1)=0.d0
      else
         unitka=maxka/dble(nka-1)
         do 10000 i=1,nka
            kas(i)=dble(i-1)*unitka
10000    continue
      endif
c     calculate lambdas
      unitlmda=log(minlmda)/dble(nlmda-1)
      call maxbi(lmdas(1),inia,inib,inieta,y,x,sz,n,q,p,
     +     epsilon,maxit)
      do 10001 i=2,nlmda
         lmdas(i)=lmdas(1)*exp(unitlmda*dble(i-1))
10001 continue
c     output lambdas and kappas
      do 10002 i=1,nlmda
         olmdas( ((i-1)*nka+1):(i*nka) )=lmdas(i)
         okas( ((i-1)*nka+1):(i*nka) )=kas(:)
10002 continue
c     Lasso solution path
      if (nka .eq. 1) then
         alpha(:)=inia(:)
         beta(:)=inib(:)
         eta(:)=inieta(:)
         ocoef(1:q,1)=alpha(:)
         ocoef((q+1):(q+p),1)=beta(:)
         do 10004 i=2,nlmda
            call lasso(alpha,beta,eta,lmdas(i),y,x,sz,n,q,p,
     +           epsilon,maxit)
            ocoef(1:q,i)=alpha(:)
            ocoef((q+1):(q+p),i)=beta(:)
10004    continue
c     MCP solution path along lambda given kappa
      else
         alpha(:)=inia(:)
         beta(:)=inib(:)
         eta(:)=inieta(:)
!     i=1 case
!     .. j=1 case
         j=1
         ocoef(1:q,        1+(j-1)*nka)=alpha(:)
         ocoef((q+1):(q+p),1+(j-1)*nka)=beta(:)
         do 1005 j=2,nlmda
            call lasso(alpha,beta,eta,lmdas(j),y,x,sz,n,q,p,
     +           epsilon,maxit)
            ocoef(1:q,        1+(j-1)*nka)=alpha(:)
            ocoef((q+1):(q+p),1+(j-1)*nka)=beta(:)
 1005    continue
         do 10005 i=2,nka
            alpha(:)=inia(:)
            beta(:)=inib(:)
            eta(:)=inieta(:)
            ocoef(1:q,        i)=alpha(:)
            ocoef((q+1):(q+p),i)=beta(:)
            do 10006 j=2,nlmda
               call mcp(alpha,beta,eta,lmdas(j),kas(i),y,x,sz,
     +              n,q,p,epsilon,maxit)
               ocoef(1:q,        i+(j-1)*nka)=alpha(:)
               ocoef((q+1):(q+p),i+(j-1)*nka)=beta(:)
10006       continue
10005    continue
      endif
c     compute model size,eigenvalue,aic,bic, objective function
      do 10 i=1,(nka*nlmda)
         call fmeabo(odf(i),oevidx(i),oaic(i),obic(i),oobj(i),
     +        ocoef(:,i),y,x,sz,n,q,p,okas(i),olmdas(i))
 10   continue
!     change the penalized ocoefficients back
      do 10008 i=1,p
         ocoef(q+i,:)=as(i)*ocoef(q+i,:)
10008 continue
      end




cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Level IV function
c     CV component for cross valiation 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!***********************************************************************
!     SCAD
!***********************************************************************

c     path along kapa 
c***********************************************************************
      subroutine cvsdka(pauc,full,cvx,lmdas,kas,nlmda,nka,cvtn,cvpn,q,p,
     +     cvty,cvtx,cvtz,epsilon,maxit,cvpy,cvpx,cvpz)
      integer nlmda,nka,cvtn,cvpn,q,p,
     +     full(nka*nlmda),cvx(nka*nlmda)
      double precision pauc(nlmda*nka),lmdas(nlmda),kas(nka),
     +     cvty(cvtn),cvtx(cvtn,q),cvtz(cvtn,p),epsilon,cvpy(cvpn),
     +     cvpx(cvpn,q),cvpz(cvpn,p)
c     local vars
      integer j,i
      double precision cvtas(p),cvtsz(cvtn,p),
     +     nulllmda,cvalpha(q),cvbeta(p),cvteta(cvtn),
     +     cvtetam(cvtn,nlmda),
     +     inicva(q,nlmda),inicvb(p,nlmda),
     +     cvas(q,nlmda*nka),cvbs(p,nlmda*nka),
     +     cvpeta(cvpn),cvpi(cvpn),tmp
!     column-wise standardization
      call standard(cvtas,cvtsz,cvtz,cvtn,p)
c     solution path along lambdas, initial values for along kas
      call maxbi(nulllmda,cvalpha,cvbeta,cvteta,cvty,cvtx,cvtsz,
     +     cvtn,q,p,epsilon,maxit)
      do 1004 i=1,nlmda
         if (lmdas(i) .ge. nulllmda) then
            inicva(:,i)=cvalpha(:)
            inicvb(:,i)=cvbeta(:)
            cvtetam(:,i)=cvteta(:)
         else
            call lasso(cvalpha,cvbeta,cvteta,lmdas(i),cvty,
     +           cvtx,cvtsz,cvtn,q,p,epsilon,maxit)
            inicva(:,i)=cvalpha(:)
            inicvb(:,i)=cvbeta(:)
            cvtetam(:,i)=cvteta(:)
         endif
 1004 continue
c     solution path along kappa
      if (nka .eq. 1) then
         cvas(:,:)=inicva(:,:)
         cvbs(:,:)=inicvb(:,:)
      else
         do 1005 i=1,nlmda
            if (lmdas(i) .ge. nulllmda) then
               do 1006 j=1,nka
                  cvas(:,(i-1)*nka+j)=inicva(:,i)
                  cvbs(:,(i-1)*nka+j)=inicvb(:,i)
 1006          continue
            else
               cvalpha(:)=inicva(:,i)
               cvbeta(:) =inicvb(:,i)
               cvteta(:)=cvtetam(:,i)
c     ....     j=1 case 
               j=1
               cvas(:,(i-1)*nka+j)=cvalpha(:)
               cvbs(:,(i-1)*nka+j)=cvbeta(:)
               do 1007 j=2,nka
                  call scad(cvalpha,cvbeta,cvteta,lmdas(i),kas(j),
     +                 cvty,cvtx,cvtsz,cvtn,q,p,
     +                 epsilon,maxit)
                  cvas(:,(i-1)*nka+j)=cvalpha(:)
                  cvbs(:,(i-1)*nka+j)=cvbeta(:)
 1007          continue
            endif
 1005    continue
      endif
!     determin whether satured or not 0: Not, 1: satured
      do 10006 i=1,nka*nlmda
         call fme(full(i),cvx(i),cvas(:,i),cvbs(:,i),
     +        cvty,cvtx,cvtsz,cvtn,q,p)
10006 continue
!     change the penalized coefficients back
      do 10005 i=1,p
         cvbs(i,:)=cvtas(i)*cvbs(i,:)
10005 continue
c     get the paucs for the prediction set
      do 1011 i=1,(nka*nlmda)
         do 1012 j=1,cvpn
            cvpeta(j)=dot_product(cvpx(j,:),cvas(:,i))+dot_product(
     +           cvpz(j,:),cvbs(:,i))
 1012    continue
         do 1013 j=1,cvpn
            cvpi(j)=1.d0/(1.d0+exp(-cvpeta(j)))
1013     continue
         call aucroc(tmp,cvpy,cvpi,cvpn)
         pauc(i)=tmp
1011  continue
      end


c     hybrid path Lasso-SCAD hybrid algorithm along kapa 
c***********************************************************************
      subroutine cvsdka2(pauc,full,cvx,lmdas,kas,nlmda,nka,cvtn,cvpn,
     +     q,p,cvty,cvtx,cvtz,epsilon,maxit,cvpy,cvpx,cvpz)
      integer nlmda,nka,cvtn,cvpn,q,p,
     +     full(nka*nlmda),cvx(nka*nlmda)
      double precision pauc(nlmda*nka),lmdas(nlmda),kas(nka),
     +     cvty(cvtn),cvtx(cvtn,q),cvtz(cvtn,p),epsilon,cvpy(cvpn),
     +     cvpx(cvpn,q),cvpz(cvpn,p)
c     local vars
      integer j,i,k,cvupidx(p)
      double precision cvtas(p),cvtsz(cvtn,p),
     +     nulllmda,cvalpha(q),cvbeta(p),cvteta(cvtn),
     +     cvtetam(cvtn,nlmda),
     +     inicva(q,nlmda),inicvb(p,nlmda),
     +     cvas(q,nlmda*nka),cvbs(p,nlmda*nka),
     +     cvpeta(cvpn),cvpi(cvpn),tmp
!     column-wise standardization
      call standard(cvtas,cvtsz,cvtz,cvtn,p)
c     solution path along lambdas, initial values for along kas
      call maxbi(nulllmda,cvalpha,cvbeta,cvteta,cvty,cvtx,cvtsz,
     +     cvtn,q,p,epsilon,maxit)
      do 1004 i=1,nlmda
         if (lmdas(i) .ge. nulllmda) then
            inicva(:,i)=cvalpha(:)
            inicvb(:,i)=cvbeta(:)
            cvtetam(:,i)=cvteta(:)
         else
            call lasso(cvalpha,cvbeta,cvteta,lmdas(i),cvty,
     +           cvtx,cvtsz,cvtn,q,p,epsilon,maxit)
            inicva(:,i)=cvalpha(:)
            inicvb(:,i)=cvbeta(:)
            cvtetam(:,i)=cvteta(:)
         endif
 1004 continue
c     solution path along kappa
      if (nka .eq. 1) then
         cvas(:,:)=inicva(:,:)
         cvbs(:,:)=inicvb(:,:)
      else
         do 1005 i=1,nlmda
            if (lmdas(i) .ge. nulllmda) then
               do 1006 j=1,nka
                  cvas(:,(i-1)*nka+j)=inicva(:,i)
                  cvbs(:,(i-1)*nka+j)=inicvb(:,i)
 1006          continue
            else
               cvalpha(:)=inicva(:,i)
               cvbeta(:) =inicvb(:,i)
               cvteta(:)=cvtetam(:,i)
c     ....     j=1 case 
               j=1
               cvas(:,(i-1)*nka+j)=cvalpha(:)
               cvbs(:,(i-1)*nka+j)=cvbeta(:)
               cvupidx(:)=1
               do 9000 k=1,p
                  if (cvbeta(k) .eq. 0.d0) cvupidx(k)=0
 9000          continue
               do 1007 j=2,nka
                  call scad2(cvalpha,cvbeta,cvteta,cvupidx,
     +                 lmdas(i),kas(j),cvty,cvtx,cvtsz,cvtn,q,p,
     +                 epsilon,maxit)
                  cvas(:,(i-1)*nka+j)=cvalpha(:)
                  cvbs(:,(i-1)*nka+j)=cvbeta(:)
 1007          continue
            endif
 1005    continue
      endif
!     determin whether satured or not 0: Not, 1: satured
      do 10006 i=1,nka*nlmda
         call fme(full(i),cvx(i),cvas(:,i),cvbs(:,i),
     +        cvty,cvtx,cvtsz,cvtn,q,p)
10006 continue
!     change the penalized coefficients back
      do 10005 i=1,p
         cvbs(i,:)=cvtas(i)*cvbs(i,:)
10005 continue
c     get the paucs for the prediction set
      do 1011 i=1,(nka*nlmda)
         do 1012 j=1,cvpn
            cvpeta(j)=dot_product(cvpx(j,:),cvas(:,i))+dot_product(
     +           cvpz(j,:),cvbs(:,i))
 1012    continue
         do 1013 j=1,cvpn
            cvpi(j)=1.d0/(1.d0+exp(-cvpeta(j)))
1013     continue
         call aucroc(tmp,cvpy,cvpi,cvpn)
         pauc(i)=tmp
1011  continue
      end


c     Along the Lambda for SCAD
c***********************************************************************
      subroutine cvsdlm(pauc,full,cvx,lmdas,kas,nlmda,nka,cvtn,cvpn,q,p,
     +     cvty,cvtx,cvtz,epsilon,maxit,cvpy,cvpx,cvpz)
      integer nlmda,nka,cvtn,cvpn,q,p,
     +     full(nka*nlmda),cvx(nka*nlmda)
      double precision pauc(nlmda*nka),lmdas(nlmda),kas(nka),
     +     cvty(cvtn),cvtx(cvtn,q),cvtz(cvtn,p),epsilon,cvpy(cvpn),
     +     cvpx(cvpn,q),cvpz(cvpn,p)
c     local vars
      integer j,i,k,cvupidx(p)
      double precision cvtas(p),cvtsz(cvtn,p),
     +     nulllmda,cvinia(q),cvinib(p),cvinieta(cvtn),
     +     cvalpha(q),cvbeta(p),cvteta(cvtn),
     +     cvas(q,nlmda*nka),cvbs(p,nlmda*nka),
     +     cvpeta(cvpn),cvpi(cvpn),tmp
!     column-wise standardization
      call standard(cvtas,cvtsz,cvtz,cvtn,p)
c     solution path along lambdas, initial values for along kas
      call maxbi(nulllmda,cvinia,cvinib,cvinieta,cvty,cvtx,cvtsz,
     +     cvtn,q,p,epsilon,maxit)
c     Lasso solution path
      if (nka .eq. 1) then 
         cvalpha(:)=cvinia(:)
         cvbeta(:)=cvinib(:)
         cvteta(:)=cvinieta(:)
         do 1004 i=1,nlmda
            if (lmdas(i) .ge. nulllmda) then
               cvas(:,i)=cvalpha(:)
               cvbs(:,i)=cvbeta(:)
            else
               call lasso(cvalpha,cvbeta,cvteta,lmdas(i),cvty,
     +              cvtx,cvtsz,cvtn,q,p,epsilon,maxit)
               cvas(:,i)=cvalpha(:)
               cvbs(:,i)=cvbeta(:)
            endif
 1004    continue
c     SCAD solution along lambda given kappa
      else
!     i=1
         cvalpha(:)=cvinia(:)
         cvbeta(:)=cvinib(:)
         cvteta(:)=cvinieta(:)
         do 10004 j=1,nlmda
            if (lmdas(j) .ge. nulllmda) then
               cvas(:,1+(j-1)*nka)=cvalpha(:)
               cvbs(:,1+(j-1)*nka)=cvbeta(:)
            else
               call lasso(cvalpha,cvbeta,cvteta,lmdas(j),cvty,
     +              cvtx,cvtsz,cvtn,q,p,epsilon,maxit)
               cvas(:,1+(j-1)*nka)=cvalpha(:)
               cvbs(:,1+(j-1)*nka)=cvbeta(:)
            endif
10004    continue
!     i=2...
         do 11005 i=2,nka
            cvalpha(:)=cvinia(:)
            cvbeta(:)=cvinib(:)
            cvteta(:)=cvinieta(:)
            do 11006 j=1,nlmda
               if (lmdas(j) .ge. nulllmda) then
                  cvas(:,i+(j-1)*nka)=cvalpha(:)
                  cvbs(:,i+(j-1)*nka)=cvbeta(:)
               else
                  call scad(cvalpha,cvbeta,cvteta,lmdas(j),kas(i),cvty,
     +                 cvtx,cvtsz,cvtn,q,p,epsilon,maxit)
                  cvas(:,i+(j-1)*nka)=cvalpha(:)
                  cvbs(:,i+(j-1)*nka)=cvbeta(:)
               endif
11006       continue   
11005    continue
      endif
!     determin whether satured or not 0: Not, 1: satured
      do 10006 i=1,nka*nlmda
         call fme(full(i),cvx(i),cvas(:,i),cvbs(:,i),
     +        cvty,cvtx,cvtsz,cvtn,q,p)
10006 continue
!     change the penalized coefficients back
      do 10005 i=1,p
         cvbs(i,:)=cvtas(i)*cvbs(i,:)
10005 continue
c     get the paucs for the prediction set
      do 1011 i=1,(nka*nlmda)
         do 1012 j=1,cvpn
            cvpeta(j)=dot_product(cvpx(j,:),cvas(:,i))+dot_product(
     +           cvpz(j,:),cvbs(:,i))
 1012    continue
         do 1013 j=1,cvpn
            cvpi(j)=1.d0/(1.d0+exp(-cvpeta(j)))
1013     continue
         call aucroc(tmp,cvpy,cvpi,cvpn)
         pauc(i)=tmp
1011  continue
      end



!***********************************************************************
!     MCP
!***********************************************************************

c     along kapa 
c***********************************************************************
      subroutine cvkapa(pauc,full,cvx,lmdas,kas,nlmda,nka,cvtn,cvpn,q,p,
     +     cvty,cvtx,cvtz,epsilon,maxit,cvpy,cvpx,cvpz)
      integer nlmda,nka,cvtn,cvpn,q,p,
     +     full(nka*nlmda),cvx(nka*nlmda)
      double precision pauc(nlmda*nka),lmdas(nlmda),kas(nka),
     +     cvty(cvtn),cvtx(cvtn,q),cvtz(cvtn,p),epsilon,cvpy(cvpn),
     +     cvpx(cvpn,q),cvpz(cvpn,p)
c     local vars
      integer j,i
      double precision cvtas(p),cvtsz(cvtn,p),
     +     nulllmda,cvalpha(q),cvbeta(p),cvteta(cvtn),
     +     cvtetam(cvtn,nlmda),
     +     inicva(q,nlmda),inicvb(p,nlmda),
     +     cvas(q,nlmda*nka),cvbs(p,nlmda*nka),
     +     cvpeta(cvpn),cvpi(cvpn),tmp
!     column-wise standardization
      call standard(cvtas,cvtsz,cvtz,cvtn,p)
c     solution path along lambdas, initial values for along kas
      call maxbi(nulllmda,cvalpha,cvbeta,cvteta,cvty,cvtx,cvtsz,
     +     cvtn,q,p,epsilon,maxit)
      do 1004 i=1,nlmda
         if (lmdas(i) .ge. nulllmda) then
            inicva(:,i)=cvalpha(:)
            inicvb(:,i)=cvbeta(:)
            cvtetam(:,i)=cvteta(:)
         else
            call lasso(cvalpha,cvbeta,cvteta,lmdas(i),cvty,
     +           cvtx,cvtsz,cvtn,q,p,epsilon,maxit)
            inicva(:,i)=cvalpha(:)
            inicvb(:,i)=cvbeta(:)
            cvtetam(:,i)=cvteta(:)
         endif
 1004 continue
c     solution path along kappa
      if (nka .eq. 1) then
         cvas(:,:)=inicva(:,:)
         cvbs(:,:)=inicvb(:,:)
      else
         do 1005 i=1,nlmda
            if (lmdas(i) .ge. nulllmda) then
               do 1006 j=1,nka
                  cvas(:,(i-1)*nka+j)=inicva(:,i)
                  cvbs(:,(i-1)*nka+j)=inicvb(:,i)
 1006          continue
            else
               cvalpha(:)=inicva(:,i)
               cvbeta(:) =inicvb(:,i)
               cvteta(:)=cvtetam(:,i)
c     ....     j=1 case 
               j=1
               cvas(:,(i-1)*nka+j)=cvalpha(:)
               cvbs(:,(i-1)*nka+j)=cvbeta(:)
               do 1007 j=2,nka
                  call mcp(cvalpha,cvbeta,cvteta,lmdas(i),kas(j),
     +                 cvty,cvtx,cvtsz,cvtn,q,p,
     +                 epsilon,maxit)
                  cvas(:,(i-1)*nka+j)=cvalpha(:)
                  cvbs(:,(i-1)*nka+j)=cvbeta(:)
 1007          continue
            endif
 1005    continue
      endif
!     determin whether satured or not 0: Not, 1: satured
      do 10006 i=1,nka*nlmda
         call fme(full(i),cvx(i),cvas(:,i),cvbs(:,i),
     +        cvty,cvtx,cvtsz,cvtn,q,p)
10006 continue
!     change the penalized coefficients back
      do 10005 i=1,p
         cvbs(i,:)=cvtas(i)*cvbs(i,:)
10005 continue
c     get the paucs for the prediction set
      do 1011 i=1,(nka*nlmda)
         do 1012 j=1,cvpn
            cvpeta(j)=dot_product(cvpx(j,:),cvas(:,i))+dot_product(
     +           cvpz(j,:),cvbs(:,i))
 1012    continue
         do 1013 j=1,cvpn
            cvpi(j)=1.d0/(1.d0+exp(-cvpeta(j)))
1013     continue
         call aucroc(tmp,cvpy,cvpi,cvpn)
         pauc(i)=tmp
1011  continue
      end



c     hybrid path: only update the non-zero lasso solution
c***********************************************************************
      subroutine cvkapa2(pauc,full,cvx,lmdas,kas,nlmda,nka,cvtn,cvpn,
     +     q,p,cvty,cvtx,cvtz,epsilon,maxit,cvpy,cvpx,cvpz)
      integer nlmda,nka,cvtn,cvpn,q,p,
     +     full(nka*nlmda),cvx(nka*nlmda)
      double precision pauc(nlmda*nka),lmdas(nlmda),kas(nka),
     +     cvty(cvtn),cvtx(cvtn,q),cvtz(cvtn,p),epsilon,cvpy(cvpn),
     +     cvpx(cvpn,q),cvpz(cvpn,p)
c     local vars
      integer j,i,k,cvupidx(p)
      double precision cvtas(p),cvtsz(cvtn,p),
     +     nulllmda,cvalpha(q),cvbeta(p),cvteta(cvtn),
     +     cvtetam(cvtn,nlmda),
     +     inicva(q,nlmda),inicvb(p,nlmda),
     +     cvas(q,nlmda*nka),cvbs(p,nlmda*nka),
     +     cvpeta(cvpn),cvpi(cvpn),tmp
!     column-wise standardization
      call standard(cvtas,cvtsz,cvtz,cvtn,p)
c     solution path along lambdas, initial values for along kas
      call maxbi(nulllmda,cvalpha,cvbeta,cvteta,cvty,cvtx,cvtsz,
     +     cvtn,q,p,epsilon,maxit)
      do 1004 i=1,nlmda
         if (lmdas(i) .ge. nulllmda) then
            inicva(:,i)=cvalpha(:)
            inicvb(:,i)=cvbeta(:)
            cvtetam(:,i)=cvteta(:)
         else
            call lasso(cvalpha,cvbeta,cvteta,lmdas(i),cvty,
     +           cvtx,cvtsz,cvtn,q,p,epsilon,maxit)
            inicva(:,i)=cvalpha(:)
            inicvb(:,i)=cvbeta(:)
            cvtetam(:,i)=cvteta(:)
         endif
 1004 continue
c     solution path along kappa
      if (nka .eq. 1) then
         cvas(:,:)=inicva(:,:)
         cvbs(:,:)=inicvb(:,:)
      else
         do 1005 i=1,nlmda
            if (lmdas(i) .ge. nulllmda) then
               do 1006 j=1,nka
                  cvas(:,(i-1)*nka+j)=inicva(:,i)
                  cvbs(:,(i-1)*nka+j)=inicvb(:,i)
 1006          continue
            else
               cvalpha(:)=inicva(:,i)
               cvbeta(:) =inicvb(:,i)
               cvteta(:)=cvtetam(:,i)
c     ....     j=1 case 
               j=1
               cvas(:,(i-1)*nka+j)=cvalpha(:)
               cvbs(:,(i-1)*nka+j)=cvbeta(:)
               cvupidx(:)=1
               do 9000 k=1,p
                  if (cvbeta(k) .eq. 0.d0) cvupidx(k)=0
 9000          continue
               do 1007 j=2,nka
                  call mcp2(cvalpha,cvbeta,cvteta,cvupidx,
     +                 lmdas(i),kas(j),cvty,cvtx,cvtsz,cvtn,q,p,
     +                 epsilon,maxit)
                  cvas(:,(i-1)*nka+j)=cvalpha(:)
                  cvbs(:,(i-1)*nka+j)=cvbeta(:)
 1007          continue
            endif
 1005    continue
      endif
!     determin whether satured or not 0: Not, 1: satured
      do 10006 i=1,nka*nlmda
         call fme(full(i),cvx(i),cvas(:,i),cvbs(:,i),
     +        cvty,cvtx,cvtsz,cvtn,q,p)
10006 continue
!     change the penalized coefficients back
      do 10005 i=1,p
         cvbs(i,:)=cvtas(i)*cvbs(i,:)
10005 continue
c     get the paucs for the prediction set
      do 1011 i=1,(nka*nlmda)
         do 1012 j=1,cvpn
            cvpeta(j)=dot_product(cvpx(j,:),cvas(:,i))+dot_product(
     +           cvpz(j,:),cvbs(:,i))
 1012    continue
         do 1013 j=1,cvpn
            cvpi(j)=1.d0/(1.d0+exp(-cvpeta(j)))
1013     continue
         call aucroc(tmp,cvpy,cvpi,cvpn)
         pauc(i)=tmp
1011  continue
      end


c     Along the Lambda for MCP
c***********************************************************************
      subroutine cvlmda(pauc,full,cvx,lmdas,kas,nlmda,nka,cvtn,cvpn,q,p,
     +     cvty,cvtx,cvtz,epsilon,maxit,cvpy,cvpx,cvpz)
      integer nlmda,nka,cvtn,cvpn,q,p,
     +     full(nka*nlmda),cvx(nka*nlmda)
      double precision pauc(nlmda*nka),lmdas(nlmda),kas(nka),
     +     cvty(cvtn),cvtx(cvtn,q),cvtz(cvtn,p),epsilon,cvpy(cvpn),
     +     cvpx(cvpn,q),cvpz(cvpn,p)
c     local vars
      integer j,i,k,cvupidx(p)
      double precision cvtas(p),cvtsz(cvtn,p),
     +     nulllmda,cvinia(q),cvinib(p),cvinieta(cvtn),
     +     cvalpha(q),cvbeta(p),cvteta(cvtn),
     +     cvas(q,nlmda*nka),cvbs(p,nlmda*nka),
     +     cvpeta(cvpn),cvpi(cvpn),tmp
!     column-wise standardization
      call standard(cvtas,cvtsz,cvtz,cvtn,p)
c     solution path along lambdas, initial values for along kas
      call maxbi(nulllmda,cvinia,cvinib,cvinieta,cvty,cvtx,cvtsz,
     +     cvtn,q,p,epsilon,maxit)
c     Lasso solution path
      if (nka .eq. 1) then 
         cvalpha(:)=cvinia(:)
         cvbeta(:)=cvinib(:)
         cvteta(:)=cvinieta(:)
         do 1004 i=1,nlmda
            if (lmdas(i) .ge. nulllmda) then
               cvas(:,i)=cvalpha(:)
               cvbs(:,i)=cvbeta(:)
            else
               call lasso(cvalpha,cvbeta,cvteta,lmdas(i),cvty,
     +              cvtx,cvtsz,cvtn,q,p,epsilon,maxit)
               cvas(:,i)=cvalpha(:)
               cvbs(:,i)=cvbeta(:)
            endif
 1004    continue
c     MCP solution along lambda given kappa
      else
!     i=1
         cvalpha(:)=cvinia(:)
         cvbeta(:)=cvinib(:)
         cvteta(:)=cvinieta(:)
         do 10004 j=1,nlmda
            if (lmdas(j) .ge. nulllmda) then
               cvas(:,1+(j-1)*nka)=cvalpha(:)
               cvbs(:,1+(j-1)*nka)=cvbeta(:)
            else
               call lasso(cvalpha,cvbeta,cvteta,lmdas(j),cvty,
     +              cvtx,cvtsz,cvtn,q,p,epsilon,maxit)
               cvas(:,1+(j-1)*nka)=cvalpha(:)
               cvbs(:,1+(j-1)*nka)=cvbeta(:)
            endif
10004    continue
!     i=2...
         do 11005 i=2,nka
            cvalpha(:)=cvinia(:)
            cvbeta(:)=cvinib(:)
            cvteta(:)=cvinieta(:)
            do 11006 j=1,nlmda
               if (lmdas(j) .ge. nulllmda) then
                  cvas(:,i+(j-1)*nka)=cvalpha(:)
                  cvbs(:,i+(j-1)*nka)=cvbeta(:)
               else
                  call mcp(cvalpha,cvbeta,cvteta,lmdas(j),kas(i),cvty,
     +                 cvtx,cvtsz,cvtn,q,p,epsilon,maxit)
                  cvas(:,i+(j-1)*nka)=cvalpha(:)
                  cvbs(:,i+(j-1)*nka)=cvbeta(:)
               endif
11006       continue   
11005    continue
      endif
!     determin whether satured or not 0: Not, 1: satured
      do 10006 i=1,nka*nlmda
         call fme(full(i),cvx(i),cvas(:,i),cvbs(:,i),
     +        cvty,cvtx,cvtsz,cvtn,q,p)
10006 continue
!     change the penalized coefficients back
      do 10005 i=1,p
         cvbs(i,:)=cvtas(i)*cvbs(i,:)
10005 continue
c     get the paucs for the prediction set
      do 1011 i=1,(nka*nlmda)
         do 1012 j=1,cvpn
            cvpeta(j)=dot_product(cvpx(j,:),cvas(:,i))+dot_product(
     +           cvpz(j,:),cvbs(:,i))
 1012    continue
         do 1013 j=1,cvpn
            cvpi(j)=1.d0/(1.d0+exp(-cvpeta(j)))
1013     continue
         call aucroc(tmp,cvpy,cvpi,cvpn)
         pauc(i)=tmp
1011  continue
      end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Level V function
c     Tuning parameter selection based on CV-AUC
c     only select on df, do not consider convexity
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      

!***********************************************************************
!     SCAD
!***********************************************************************


c     Aong kappa
c***********************************************************************
      subroutine cvaucsdka(out,pauc,
     +     olmdas,okas,ocoef,oaic,obic,oobj,odf,ocvx,ofull,
     +     cvcvx,cvfull,
     +     nindex,cvk,y,x,z,n,q,p,nka,maxka,nlmda,minlmda,epsilon,maxit)
      integer cvk,n,q,p,nka,nlmda,nindex(n),maxit,
     +     odf(nka*nlmda),ocvx(nka*nlmda),ofull(nka*nlmda),
     +     cvcvx(nka*nlmda),cvfull(nka*nlmda)
      double precision out(3+q+p),pauc(nka*nlmda),olmdas(nka*nlmda),
     +     okas(nka*nlmda),ocoef(q+p,nka*nlmda),oaic(nka*nlmda),
     +     obic(nka*nlmda),oobj(nka*nlmda),
     +     y(n),x(n,q),z(n,p),maxka,minlmda,epsilon
c     local vars
      integer i,cv,cvtn,cvpn,k,
     +     cvfullm(nka*nlmda,cvk),cvcvxm(nka*nlmda,cvk),
     +     use,loc
      double precision lmdas(nlmda),kas(nka),
     +     cvpy(n),cvpx(n,q),cvpz(n,p),cvty(n),cvtx(n,q),cvtz(n,p),
     +     cvpauc(nka*nlmda,cvk),
     +     omlmdas(nka*nlmda),omkas(nka*nlmda),omcoef(q+p,nka*nlmda),
     +     ompauc(nka*nlmda)
c     solution path of original dataset
      call scadkapa(olmdas,okas,ocoef,oaic,obic,oobj,odf,ocvx,
     +     y,x,z,n,q,p,nka,maxka,nlmda,minlmda,epsilon,maxit)
c     Determine sature model, 0:not, 1: satured
      do 10000 i=1,nka*nlmda
         if ( odf(i) .le. n) then
            ofull(i)=0
         else
            ofull(i)=1
         endif
10000 continue
c     get lambdas,kappas
      do 10005 i=1,nlmda
         lmdas(i)= olmdas((i-1)*nka+1)
10005 continue
      kas(:)=okas(1:nka)
c     start cross validation
      do 10006 cv=1,cvk
         cvtn=0
         cvpn=0
         do 10007 k=1,n
            if (nindex(k) .eq. cv) then
               cvpn=1+cvpn
               cvpy(cvpn)=y(k)
               cvpx(cvpn,:)=x(k,:)
               cvpz(cvpn,:)=z(k,:)
            else
               cvtn=1+cvtn
               cvty(cvtn)=y(k)
               cvtx(cvtn,:)=x(k,:)
               cvtz(cvtn,:)=z(k,:)
            end if
10007    continue
c     .... solve alpha, beta, and pauc for cv datasets
         call cvsdka(cvpauc(:,cv),cvfullm(:,cv),cvcvxm(:,cv),lmdas,kas,
     +        nlmda,nka,cvtn,cvpn,q,p,
     +        cvty(1:cvtn),cvtx(1:cvtn,:),cvtz(1:cvtn,:),
     +        epsilon,maxit,
     +        cvpy(1:cvpn),cvpx(1:cvpn,:),cvpz(1:cvpn,:))
10006 continue
!     get average pauc, satured, convex info  in cross-validation
      do 10008 i=1,(nlmda*nka)
         pauc(i)=sum(cvpauc(i,:))/dble(cvk)
         if ( sum(cvfullm(i,:)) .eq. 0) then
            cvfull(i)=0
         else
            cvfull(i)=1
         endif
!     .....if ( sum(cvcvx(i,:)) .eq. 5) cvcvx(i)=1
10008 continue
!     assign useful output
      use=0
      do 10009 i=1,nka*nlmda
         if ( (ofull(i) .eq. 0) .and. (cvfull(i) .eq. 0) ) then
            use=use+1
            ompauc(use) =pauc(i)
            omlmdas(use)=olmdas(i)
            omkas(use)  =okas(i)
            omcoef(:,use)=ocoef(:,i)
         endif
10009 continue
!     get maximum
      loc=maxloc( ompauc(1:use), use )
      out(1)=ompauc(loc)
      out(2)=omlmdas(loc)
      out(3)=omkas(loc)
      out(4:(3+q+p))=omcoef(:,loc)
      end


c     hybrid path: update the non-zero variables
c***********************************************************************
      subroutine cvaucsdka2(out,pauc,
     +     olmdas,okas,ocoef,oaic,obic,oobj,odf,ocvx,ofull,
     +     cvcvx,cvfull,
     +     nindex,cvk,y,x,z,n,q,p,nka,maxka,nlmda,minlmda,epsilon,maxit)
      integer cvk,n,q,p,nka,nlmda,nindex(n),maxit,
     +     odf(nka*nlmda),ocvx(nka*nlmda),ofull(nka*nlmda),
     +     cvcvx(nka*nlmda),cvfull(nka*nlmda)
      double precision out(3+q+p),pauc(nka*nlmda),olmdas(nka*nlmda),
     +     okas(nka*nlmda),ocoef(q+p,nka*nlmda),oaic(nka*nlmda),
     +     obic(nka*nlmda),oobj(nka*nlmda),
     +     y(n),x(n,q),z(n,p),maxka,minlmda,epsilon
c     local vars
      integer i,cv,cvtn,cvpn,k,
     +     cvfullm(nka*nlmda,cvk),cvcvxm(nka*nlmda,cvk),
     +     use,loc
      double precision lmdas(nlmda),kas(nka),
     +     cvpy(n),cvpx(n,q),cvpz(n,p),cvty(n),cvtx(n,q),cvtz(n,p),
     +     cvpauc(nka*nlmda,cvk),
     +     omlmdas(nka*nlmda),omkas(nka*nlmda),omcoef(q+p,nka*nlmda),
     +     ompauc(nka*nlmda)
c     solution path of original dataset
      call scadkapa2(olmdas,okas,ocoef,oaic,obic,oobj,odf,ocvx,
     +     y,x,z,n,q,p,nka,maxka,nlmda,minlmda,epsilon,maxit)
c     Determine sature model, 0:not, 1: satured
      do 10000 i=1,nka*nlmda
         if ( odf(i) .le. n) then
            ofull(i)=0
         else
            ofull(i)=1
         endif
10000 continue
c     get lambdas,kappas
      do 10005 i=1,nlmda
         lmdas(i)= olmdas((i-1)*nka+1)
10005 continue
      kas(:)=okas(1:nka)
c     start cross validation
      do 10006 cv=1,cvk
         cvtn=0
         cvpn=0
         do 10007 k=1,n
            if (nindex(k) .eq. cv) then
               cvpn=1+cvpn
               cvpy(cvpn)=y(k)
               cvpx(cvpn,:)=x(k,:)
               cvpz(cvpn,:)=z(k,:)
            else
               cvtn=1+cvtn
               cvty(cvtn)=y(k)
               cvtx(cvtn,:)=x(k,:)
               cvtz(cvtn,:)=z(k,:)
            end if
10007    continue
c     .... solve alpha, beta, and pauc for cv datasets
         call cvsdka2(cvpauc(:,cv),cvfullm(:,cv),cvcvxm(:,cv),lmdas,kas,
     +        nlmda,nka,cvtn,cvpn,q,p,
     +        cvty(1:cvtn),cvtx(1:cvtn,:),cvtz(1:cvtn,:),
     +        epsilon,maxit,
     +        cvpy(1:cvpn),cvpx(1:cvpn,:),cvpz(1:cvpn,:))
10006 continue
!     get average pauc, satured, convex info  in cross-validation
      do 10008 i=1,(nlmda*nka)
         pauc(i)=sum(cvpauc(i,:))/dble(cvk)
         if ( sum(cvfullm(i,:)) .eq. 0) then
            cvfull(i)=0
         else
            cvfull(i)=1
         endif
!     .....if ( sum(cvcvx(i,:)) .eq. 5) cvcvx(i)=1
10008 continue
!     assign useful output
      use=0
      do 10009 i=1,nka*nlmda
         if ( (ofull(i) .eq. 0) .and. (cvfull(i) .eq. 0) ) then
            use=use+1
            ompauc(use) =pauc(i)
            omlmdas(use)=olmdas(i)
            omkas(use)  =okas(i)
            omcoef(:,use)=ocoef(:,i)
         endif
10009 continue
!     get maximum
      loc=maxloc( ompauc(1:use), use )
      out(1)=ompauc(loc)
      out(2)=omlmdas(loc)
      out(3)=omkas(loc)
      out(4:(3+q+p))=omcoef(:,loc)
      end


c     along lambda
c***********************************************************************
      subroutine cvaucsdlm(out,pauc,
     +     olmdas,okas,ocoef,oaic,obic,oobj,odf,ocvx,ofull,
     +     cvcvx,cvfull,
     +     nindex,cvk,y,x,z,n,q,p,nka,maxka,nlmda,minlmda,epsilon,maxit)
      integer cvk,n,q,p,nka,nlmda,nindex(n),maxit,
     +     odf(nka*nlmda),ocvx(nka*nlmda),ofull(nka*nlmda),
     +     cvcvx(nka*nlmda),cvfull(nka*nlmda)
      double precision out(3+q+p),pauc(nka*nlmda),olmdas(nka*nlmda),
     +     okas(nka*nlmda),ocoef(q+p,nka*nlmda),oaic(nka*nlmda),
     +     obic(nka*nlmda),oobj(nka*nlmda),
     +     y(n),x(n,q),z(n,p),maxka,minlmda,epsilon
c     local vars
      integer i,cv,cvtn,cvpn,k,
     +     cvfullm(nka*nlmda,cvk),cvcvxm(nka*nlmda,cvk),
     +     use,loc
      double precision lmdas(nlmda),kas(nka),
     +     cvpy(n),cvpx(n,q),cvpz(n,p),cvty(n),cvtx(n,q),cvtz(n,p),
     +     cvpauc(nka*nlmda,cvk),
     +     omlmdas(nka*nlmda),omkas(nka*nlmda),omcoef(q+p,nka*nlmda),
     +     ompauc(nka*nlmda)
c     solution path of original dataset
      call scadlmda(olmdas,okas,ocoef,oaic,obic,oobj,odf,ocvx,
     +     y,x,z,n,q,p,nka,maxka,nlmda,minlmda,epsilon,maxit)
c     Determine sature model, 0:not, 1: satured
      do 10000 i=1,nka*nlmda
         if ( odf(i) .le. n) then
            ofull(i)=0
         else
            ofull(i)=1
         endif
10000 continue
c     get lambdas,kappas
      do 10005 i=1,nlmda
         lmdas(i)= olmdas((i-1)*nka+1)
10005 continue
      kas(:)=okas(1:nka)
c     start cross validation
      do 10006 cv=1,cvk
         cvtn=0
         cvpn=0
         do 10007 k=1,n
            if (nindex(k) .eq. cv) then
               cvpn=1+cvpn
               cvpy(cvpn)=y(k)
               cvpx(cvpn,:)=x(k,:)
               cvpz(cvpn,:)=z(k,:)
            else
               cvtn=1+cvtn
               cvty(cvtn)=y(k)
               cvtx(cvtn,:)=x(k,:)
               cvtz(cvtn,:)=z(k,:)
            end if
10007    continue
c     .... solve alpha, beta, and pauc for cv datasets
         call cvsdlm(cvpauc(:,cv),cvfullm(:,cv),cvcvxm(:,cv),lmdas,kas,
     +        nlmda,nka,cvtn,cvpn,q,p,
     +        cvty(1:cvtn),cvtx(1:cvtn,:),cvtz(1:cvtn,:),
     +        epsilon,maxit,
     +        cvpy(1:cvpn),cvpx(1:cvpn,:),cvpz(1:cvpn,:))
10006 continue
!     get average pauc, satured, convex info  in cross-validation
      do 10008 i=1,(nlmda*nka)
         pauc(i)=sum(cvpauc(i,:))/dble(cvk)
         if ( sum(cvfullm(i,:)) .eq. 0) then
            cvfull(i)=0
         else
            cvfull(i)=1
         endif
!     .....if ( sum(cvcvx(i,:)) .eq. 5) cvcvx(i)=1
10008 continue
!     assign useful output
      use=0
      do 10009 i=1,nka*nlmda
         if ( (ofull(i) .eq. 0) .and. (cvfull(i) .eq. 0) ) then
            use=use+1
            ompauc(use) =pauc(i)
            omlmdas(use)=olmdas(i)
            omkas(use)  =okas(i)
            omcoef(:,use)=ocoef(:,i)
         endif
10009 continue
!     get maximum
      loc=maxloc( ompauc(1:use), use )
      out(1)=ompauc(loc)
      out(2)=omlmdas(loc)
      out(3)=omkas(loc)
      out(4:(3+q+p))=omcoef(:,loc)
      end

!***********************************************************************
!     MCP 
!***********************************************************************

c     along kappa
c***********************************************************************
      subroutine cvauckapa(out,pauc,
     +     olmdas,okas,ocoef,oaic,obic,oobj,odf,ocvx,ofull,
     +     cvcvx,cvfull,
     +     nindex,cvk,y,x,z,n,q,p,nka,maxka,nlmda,minlmda,epsilon,maxit)
      integer cvk,n,q,p,nka,nlmda,nindex(n),maxit,
     +     odf(nka*nlmda),ocvx(nka*nlmda),ofull(nka*nlmda),
     +     cvcvx(nka*nlmda),cvfull(nka*nlmda)
      double precision out(3+q+p),pauc(nka*nlmda),olmdas(nka*nlmda),
     +     okas(nka*nlmda),ocoef(q+p,nka*nlmda),oaic(nka*nlmda),
     +     obic(nka*nlmda),oobj(nka*nlmda),
     +     y(n),x(n,q),z(n,p),maxka,minlmda,epsilon
c     local vars
      integer i,cv,cvtn,cvpn,k,
     +     cvfullm(nka*nlmda,cvk),cvcvxm(nka*nlmda,cvk),
     +     use,loc
      double precision lmdas(nlmda),kas(nka),
     +     cvpy(n),cvpx(n,q),cvpz(n,p),cvty(n),cvtx(n,q),cvtz(n,p),
     +     cvpauc(nka*nlmda,cvk),
     +     omlmdas(nka*nlmda),omkas(nka*nlmda),omcoef(q+p,nka*nlmda),
     +     ompauc(nka*nlmda)
c     solution path of original dataset
      call mcpkapa(olmdas,okas,ocoef,oaic,obic,oobj,odf,ocvx,
     +     y,x,z,n,q,p,nka,maxka,nlmda,minlmda,epsilon,maxit)
c     Determine sature model, 0:not, 1: satured
      do 10000 i=1,nka*nlmda
         if ( odf(i) .le. n) then
            ofull(i)=0
         else
            ofull(i)=1
         endif
10000 continue
c     get lambdas,kappas
      do 10005 i=1,nlmda
         lmdas(i)= olmdas((i-1)*nka+1)
10005 continue
      kas(:)=okas(1:nka)
c     start cross validation
      do 10006 cv=1,cvk
         cvtn=0
         cvpn=0
         do 10007 k=1,n
            if (nindex(k) .eq. cv) then
               cvpn=1+cvpn
               cvpy(cvpn)=y(k)
               cvpx(cvpn,:)=x(k,:)
               cvpz(cvpn,:)=z(k,:)
            else
               cvtn=1+cvtn
               cvty(cvtn)=y(k)
               cvtx(cvtn,:)=x(k,:)
               cvtz(cvtn,:)=z(k,:)
            end if
10007    continue
c     .... solve alpha, beta, and pauc for cv datasets
         call cvkapa(cvpauc(:,cv),cvfullm(:,cv),cvcvxm(:,cv),lmdas,kas,
     +        nlmda,nka,cvtn,cvpn,q,p,
     +        cvty(1:cvtn),cvtx(1:cvtn,:),cvtz(1:cvtn,:),
     +        epsilon,maxit,
     +        cvpy(1:cvpn),cvpx(1:cvpn,:),cvpz(1:cvpn,:))
10006 continue
!     get average pauc, satured, convex info  in cross-validation
      do 10008 i=1,(nlmda*nka)
         pauc(i)=sum(cvpauc(i,:))/dble(cvk)
         if ( sum(cvfullm(i,:)) .eq. 0) then
            cvfull(i)=0
         else
            cvfull(i)=1
         endif
!     .....if ( sum(cvcvx(i,:)) .eq. 5) cvcvx(i)=1
10008 continue
!     assign useful output
      use=0
      do 10009 i=1,nka*nlmda
         if ( (ofull(i) .eq. 0) .and. (cvfull(i) .eq. 0) ) then
            use=use+1
            ompauc(use) =pauc(i)
            omlmdas(use)=olmdas(i)
            omkas(use)  =okas(i)
            omcoef(:,use)=ocoef(:,i)
         endif
10009 continue
!     get maximum
      loc=maxloc( ompauc(1:use), use )
      out(1)=ompauc(loc)
      out(2)=omlmdas(loc)
      out(3)=omkas(loc)
      out(4:(3+q+p))=omcoef(:,loc)
      end


c     hybrid path: update the non-zero only
c***********************************************************************
      subroutine cvauckapa2(out,pauc,
     +     olmdas,okas,ocoef,oaic,obic,oobj,odf,ocvx,ofull,
     +     cvcvx,cvfull,
     +     nindex,cvk,y,x,z,n,q,p,nka,maxka,nlmda,minlmda,epsilon,maxit)
      integer cvk,n,q,p,nka,nlmda,nindex(n),maxit,
     +     odf(nka*nlmda),ocvx(nka*nlmda),ofull(nka*nlmda),
     +     cvcvx(nka*nlmda),cvfull(nka*nlmda)
      double precision out(3+q+p),pauc(nka*nlmda),olmdas(nka*nlmda),
     +     okas(nka*nlmda),ocoef(q+p,nka*nlmda),oaic(nka*nlmda),
     +     obic(nka*nlmda),oobj(nka*nlmda),
     +     y(n),x(n,q),z(n,p),maxka,minlmda,epsilon
c     local vars
      integer i,cv,cvtn,cvpn,k,
     +     cvfullm(nka*nlmda,cvk),cvcvxm(nka*nlmda,cvk),
     +     use,loc
      double precision lmdas(nlmda),kas(nka),
     +     cvpy(n),cvpx(n,q),cvpz(n,p),cvty(n),cvtx(n,q),cvtz(n,p),
     +     cvpauc(nka*nlmda,cvk),
     +     omlmdas(nka*nlmda),omkas(nka*nlmda),omcoef(q+p,nka*nlmda),
     +     ompauc(nka*nlmda)
c     solution path of original dataset
      call mcpkapa2(olmdas,okas,ocoef,oaic,obic,oobj,odf,ocvx,
     +     y,x,z,n,q,p,nka,maxka,nlmda,minlmda,epsilon,maxit)
c     Determine sature model, 0:not, 1: satured
      do 10000 i=1,nka*nlmda
         if ( odf(i) .le. n) then
            ofull(i)=0
         else
            ofull(i)=1
         endif
10000 continue
c     get lambdas,kappas
      do 10005 i=1,nlmda
         lmdas(i)= olmdas((i-1)*nka+1)
10005 continue
      kas(:)=okas(1:nka)
c     start cross validation
      do 10006 cv=1,cvk
         cvtn=0
         cvpn=0
         do 10007 k=1,n
            if (nindex(k) .eq. cv) then
               cvpn=1+cvpn
               cvpy(cvpn)=y(k)
               cvpx(cvpn,:)=x(k,:)
               cvpz(cvpn,:)=z(k,:)
            else
               cvtn=1+cvtn
               cvty(cvtn)=y(k)
               cvtx(cvtn,:)=x(k,:)
               cvtz(cvtn,:)=z(k,:)
            end if
10007    continue
c     .... solve alpha, beta, and pauc for cv datasets
         call cvkapa2(cvpauc(:,cv),cvfullm(:,cv),cvcvxm(:,cv),lmdas,kas,
     +        nlmda,nka,cvtn,cvpn,q,p,
     +        cvty(1:cvtn),cvtx(1:cvtn,:),cvtz(1:cvtn,:),
     +        epsilon,maxit,
     +        cvpy(1:cvpn),cvpx(1:cvpn,:),cvpz(1:cvpn,:))
10006 continue
!     get average pauc, satured, convex info  in cross-validation
      do 10008 i=1,(nlmda*nka)
         pauc(i)=sum(cvpauc(i,:))/dble(cvk)
         if ( sum(cvfullm(i,:)) .eq. 0) then
            cvfull(i)=0
         else
            cvfull(i)=1
         endif
!     .....if ( sum(cvcvx(i,:)) .eq. 5) cvcvx(i)=1
10008 continue
!     assign useful output
      use=0
      do 10009 i=1,nka*nlmda
         if ( (ofull(i) .eq. 0) .and. (cvfull(i) .eq. 0) ) then
            use=use+1
            ompauc(use) =pauc(i)
            omlmdas(use)=olmdas(i)
            omkas(use)  =okas(i)
            omcoef(:,use)=ocoef(:,i)
         endif
10009 continue
!     get maximum
      loc=maxloc( ompauc(1:use), use )
      out(1)=ompauc(loc)
      out(2)=omlmdas(loc)
      out(3)=omkas(loc)
      out(4:(3+q+p))=omcoef(:,loc)
      end


c     along lambda
c***********************************************************************
      subroutine cvauclmda(out,pauc,
     +     olmdas,okas,ocoef,oaic,obic,oobj,odf,ocvx,ofull,
     +     cvcvx,cvfull,
     +     nindex,cvk,y,x,z,n,q,p,nka,maxka,nlmda,minlmda,epsilon,maxit)
      integer cvk,n,q,p,nka,nlmda,nindex(n),maxit,
     +     odf(nka*nlmda),ocvx(nka*nlmda),ofull(nka*nlmda),
     +     cvcvx(nka*nlmda),cvfull(nka*nlmda)
      double precision out(3+q+p),pauc(nka*nlmda),olmdas(nka*nlmda),
     +     okas(nka*nlmda),ocoef(q+p,nka*nlmda),oaic(nka*nlmda),
     +     obic(nka*nlmda),oobj(nka*nlmda),
     +     y(n),x(n,q),z(n,p),maxka,minlmda,epsilon
c     local vars
      integer i,cv,cvtn,cvpn,k,
     +     cvfullm(nka*nlmda,cvk),cvcvxm(nka*nlmda,cvk),
     +     use,loc
      double precision lmdas(nlmda),kas(nka),
     +     cvpy(n),cvpx(n,q),cvpz(n,p),cvty(n),cvtx(n,q),cvtz(n,p),
     +     cvpauc(nka*nlmda,cvk),
     +     omlmdas(nka*nlmda),omkas(nka*nlmda),omcoef(q+p,nka*nlmda),
     +     ompauc(nka*nlmda)
c     solution path of original dataset
      call mcplmda(olmdas,okas,ocoef,oaic,obic,oobj,odf,ocvx,
     +     y,x,z,n,q,p,nka,maxka,nlmda,minlmda,epsilon,maxit)
c     Determine sature model, 0:not, 1: satured
      do 10000 i=1,nka*nlmda
         if ( odf(i) .le. n) then
            ofull(i)=0
         else
            ofull(i)=1
         endif
10000 continue
c     get lambdas,kappas
      do 10005 i=1,nlmda
         lmdas(i)= olmdas((i-1)*nka+1)
10005 continue
      kas(:)=okas(1:nka)
c     start cross validation
      do 10006 cv=1,cvk
         cvtn=0
         cvpn=0
         do 10007 k=1,n
            if (nindex(k) .eq. cv) then
               cvpn=1+cvpn
               cvpy(cvpn)=y(k)
               cvpx(cvpn,:)=x(k,:)
               cvpz(cvpn,:)=z(k,:)
            else
               cvtn=1+cvtn
               cvty(cvtn)=y(k)
               cvtx(cvtn,:)=x(k,:)
               cvtz(cvtn,:)=z(k,:)
            end if
10007    continue
c     .... solve alpha, beta, and pauc for cv datasets
         call cvlmda(cvpauc(:,cv),cvfullm(:,cv),cvcvxm(:,cv),lmdas,kas,
     +        nlmda,nka,cvtn,cvpn,q,p,
     +        cvty(1:cvtn),cvtx(1:cvtn,:),cvtz(1:cvtn,:),
     +        epsilon,maxit,
     +        cvpy(1:cvpn),cvpx(1:cvpn,:),cvpz(1:cvpn,:))
10006 continue
!     get average pauc, satured, convex info  in cross-validation
      do 10008 i=1,(nlmda*nka)
         pauc(i)=sum(cvpauc(i,:))/dble(cvk)
         if ( sum(cvfullm(i,:)) .eq. 0) then
            cvfull(i)=0
         else
            cvfull(i)=1
         endif
!     .....if ( sum(cvcvx(i,:)) .eq. 5) cvcvx(i)=1
10008 continue
!     assign useful output
      use=0
      do 10009 i=1,nka*nlmda
         if ( (ofull(i) .eq. 0) .and. (cvfull(i) .eq. 0) ) then
            use=use+1
            ompauc(use) =pauc(i)
            omlmdas(use)=olmdas(i)
            omkas(use)  =okas(i)
            omcoef(:,use)=ocoef(:,i)
         endif
10009 continue
!     get maximum
      loc=maxloc( ompauc(1:use), use )
      out(1)=ompauc(loc)
      out(2)=omlmdas(loc)
      out(3)=omkas(loc)
      out(4:(3+q+p))=omcoef(:,loc)
      end



c***********************************************************************
C     Adaptive rescaling for MCP only
c***********************************************************************


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     parallel to level II function
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     MCP: adaptive resaling 
c***********************************************************************
      subroutine adpmcp(alpha,beta,eta,lmda,ka,y,x,z,n,q,p,
     +     epsilon,maxit)
      integer n,q,p,maxit
      double precision alpha(q),beta(p),eta(n),lmda,ka,y(n),x(n,q),
     +     z(n,p),epsilon
      integer count,j,i,tag
      double precision alphaold(q),betaold(p),pi(n),w(n),r(n),v,m,numor
c     ... begin of iteration
      count=0
 1000 continue
      alphaold(:)=alpha(:)
      betaold(:)=beta(:)
!     ... update beta part
      do 1001 j=1,p
         v=0.d0
         do 1003 i=1,n
            pi(i)=1.d0/(1.d0+exp(-eta(i)))
            if (pi(i) .lt. 0.0001) pi(i)=0.d0
            if (pi(i) .gt. 0.9999) pi(i)=1.d0
            w(i)=pi(i)*(1.d0-pi(i))
            r(i)=y(i)-pi(i)
            v=v + z(i,j)*z(i,j)*w(i)
 1003    continue
         v=v/dble(n)
!     call dblepr("v",-1,v,1)
!     ...... m=beta_{j}v_{j}+Z'_{j}r/n
         m=beta(j)*v + dot_product(z(:,j),r)/dble(n)
         if (abs(m) .lt. lmda/ka) then
            call soft(numor,m,lmda)
            beta(j)=numor/(v*(1.d0-ka))
         else
            beta(j)=m/v
         end if
!     ...... eta^{s+1}=eta^{s}+z_{j}*(b^{s+1}-b^{s})
         eta(:)=eta(:) + z(:,j)*( beta(j) - betaold(j) )
 1001 continue
!     ... update alpha part
      do 11001 j=1,q
         v=0.d0
         do 11003 i=1,n
            pi(i)=1.d0/(1.d0+exp(-eta(i)))
            if (pi(i) .lt. 0.0001) pi(i)=0.d0
            if (pi(i) .gt. 0.9999) pi(i)=1.d0
            w(i)=pi(i)*(1.d0-pi(i))
            r(i)=y(i)-pi(i)
            v=v + x(i,j)*x(i,j)*w(i)
11003    continue
         v=v/dble(n)
         m=alpha(j)*v + dot_product(x(:,j),r)/dble(n)
         alpha(j)=m/v
!     ...... eta^{s+1}=eta^{s}+z_{j}*(b^{s+1}-b^{s})
         eta(:)=eta(:) + x(:,j)*( alpha(j) - alphaold(j) )      
11001 continue
!      call dblepr("beta",-1,beta,p)
!      call dblepr("alpha",-1,alpha,q)
!     do 1215 i=1,n
!         pi(i)=1.d0/(1.d0+exp(-eta(i)))
!         if (pi(i) .lt. 0.0001) pi(i)=0.d0
!         if (pi(i) .gt. 0.9999) pi(i)=1.d0
!         r(i)=y(i)-pi(i)
! 1215 continue
!      call dgemv("T",n,q,4.d0/dble(n),x,n,r,1,1.d0,alpha,1)
!     ...... eta^{s+1}=eta^{s}+z_{jl}*(b^{s+1}-b^{s})
!      do 1217 j=1,q
!     eta(:)=eta(:) + x(:,j)*( alpha(j) - alphaold(j) )
! 1217 continue
      count=count+1
      call converge2(tag,alpha,alphaold,q,beta,betaold,p,epsilon)
      if (tag .eq. 1) go to 1006
      if (count .ge. maxit) call rexit("Diverge at kappa & lambda! \n")
      if (count .lt. maxit) go to 1000
 1006 continue
c     ... end of iteration
      end



c     MCP: adaptive resaling, update the non-zero lasso solution only
c***********************************************************************
      subroutine adpmcp2(alpha,beta,eta,upidx,lmda,ka,y,x,z,n,q,p,
     +     epsilon,maxit)
      integer n,q,p,maxit,upidx(p)
      double precision alpha(q),beta(p),eta(n),lmda,ka,y(n),x(n,q),
     +     z(n,p),epsilon
      integer count,j,i,tag
      double precision alphaold(q),betaold(p),pi(n),w(n),r(n),v,m,numor
c     ... begin of iteration
      count=0
 1000 continue
      alphaold(:)=alpha(:)
      betaold(:)=beta(:)
!     ... update beta part
      do 1001 j=1,p
         if (upidx(j) .eq. 1) then 
            v=0.d0
            do 1003 i=1,n
               pi(i)=1.d0/(1.d0+exp(-eta(i)))
               if (pi(i) .lt. 0.0001) pi(i)=0.d0
               if (pi(i) .gt. 0.9999) pi(i)=1.d0
               w(i)=pi(i)*(1.d0-pi(i))
               r(i)=y(i)-pi(i)
               v=v + z(i,j)*z(i,j)*w(i)
 1003       continue
            v=v/dble(n)
!     call dblepr("v",-1,v,1)
!     ...... m=beta_{j}v_{j}+Z'_{j}r/n
            m=beta(j)*v + dot_product(z(:,j),r)/dble(n)
            if (abs(m) .lt. lmda/ka) then
               call soft(numor,m,lmda)
               beta(j)=numor/(v*(1.d0-ka))
            else
               beta(j)=m/v
            end if
!     ...... eta^{s+1}=eta^{s}+z_{j}*(b^{s+1}-b^{s})
            eta(:)=eta(:) + z(:,j)*( beta(j) - betaold(j) )
         endif
 1001 continue
!     ... update alpha part
      do 11001 j=1,q
         v=0.d0
         do 11003 i=1,n
            pi(i)=1.d0/(1.d0+exp(-eta(i)))
            if (pi(i) .lt. 0.0001) pi(i)=0.d0
            if (pi(i) .gt. 0.9999) pi(i)=1.d0
            w(i)=pi(i)*(1.d0-pi(i))
            r(i)=y(i)-pi(i)
            v=v + x(i,j)*x(i,j)*w(i)
11003    continue
         v=v/dble(n)
         m=alpha(j)*v + dot_product(x(:,j),r)/dble(n)
         alpha(j)=m/v
!     ...... eta^{s+1}=eta^{s}+z_{j}*(b^{s+1}-b^{s})
         eta(:)=eta(:) + x(:,j)*( alpha(j) - alphaold(j) )      
11001 continue
!     call dblepr("beta",-1,beta,p)
!     call dblepr("alpha",-1,alpha,q)
!     do 1215 i=1,n
!     pi(i)=1.d0/(1.d0+exp(-eta(i)))
!     if (pi(i) .lt. 0.0001) pi(i)=0.d0
!     if (pi(i) .gt. 0.9999) pi(i)=1.d0
!     r(i)=y(i)-pi(i)
!     1215 continue
!      call dgemv("T",n,q,4.d0/dble(n),x,n,r,1,1.d0,alpha,1)
!     ...... eta^{s+1}=eta^{s}+z_{jl}*(b^{s+1}-b^{s})
!     do 1217 j=1,q
!     eta(:)=eta(:) + x(:,j)*( alpha(j) - alphaold(j) )
!     1217 continue
      count=count+1
      call converge2(tag,alpha,alphaold,q,beta,betaold,p,epsilon)
      if (tag .eq. 1) go to 1006
      if (count .ge. maxit) call rexit("Diverge at kappa & lambda! \n")
      if (count .lt. maxit) go to 1000
 1006 continue
c     ... end of iteration
      end
      
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     parallel to LEVEL III FUNCTIONS
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     solution path along kappa
c***********************************************************************
      subroutine adpmcpkp(olmdas,okas,ocoef,oaic,obic,oobj,odf,oevidx,
     +     y,x,z,n,q,p,nka,maxka,nlmda,minlmda,epsilon,maxit)
      integer n,q,p,nka,nlmda,maxit,odf(nka*nlmda),oevidx(nka*nlmda)
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),oaic(nka*nlmda),obic(nka*nlmda),
     +     oobj(nka*nlmda),y(n),x(n,q),z(n,p),
     +     maxka,minlmda,epsilon
      integer i,j
      double precision as(p),sz(n,p),kas(nka),unitka,unitlmda,
     +     lmdas(nlmda),alpha(q),beta(p),eta(n),
     +     inia(q,nlmda),inib(p,nlmda),etamat(n,nlmda)
c     standardization of Z
      call standard(as,sz,z,n,p)
c     calculate kappas
      if (nka .eq. 1) then
         kas(1)=0.d0
      else
         unitka=maxka/dble(nka-1)
         do 10000 i=1,nka
            kas(i)=dble(i-1)*unitka
10000    continue
      endif
c     calculate lambdas
      unitlmda=log(minlmda)/dble(nlmda-1)
      call maxbi(lmdas(1),alpha,beta,eta,y,x,sz,n,q,p,
     +     epsilon,maxit)
      do 10001 i=2,nlmda
         lmdas(i)=lmdas(1)*exp(unitlmda*dble(i-1))
10001 continue
c     Lasso solution path as initials
      inia(:,1)=alpha(:)
      inib(:,1)=beta(:)
      etamat(:,1)=eta(:)
      do 10002 i=2,nlmda
         call lasso(alpha,beta,eta,lmdas(i),y,x,sz,n,q,p,
     +        epsilon,maxit)
         inia(:,i)=alpha(:)
         inib(:,i)=beta(:)
         etamat(:,i)=eta(:)
10002 continue
c     MCP solution path along kappa
      if (nka .eq. 1) then
         olmdas(:)=lmdas(:)
         okas(1:nlmda)=kas(1)
         ocoef(1:q,:)=inia(:,:)
         ocoef((q+1):(q+p),:)=inib(:,:)
      else
         do 1003 i=1,nlmda
            olmdas(((i-1)*nka+1):(i*nka))=lmdas(i)
            okas(((i-1)*nka+1):(i*nka))=kas(:)
 1003    continue
!     ... i=1 case
         i=1
         alpha(:)=inia(:,i)
         beta(:) =inib(:,i)
         eta(:)  =etamat(:,i)
         do 1005 j=1,nka
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1005    continue
         do 1004 i=2,nlmda
            alpha(:)=inia(:,i)
            beta(:) =inib(:,i)
            eta(:)  =etamat(:,i)
c     ...... j=1 case
            j=1
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
            do 1006 j=2,nka
!               call intpr("i",-1,i,1)
!               call intpr("j",-1,j,1)
               call adpmcp(alpha,beta,eta,lmdas(i),kas(j),y,x,sz,
     +              n,q,p,epsilon,maxit)
!     call dblepr("beta",-1,beta,p)
               ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
               ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1006       continue
 1004    continue
      endif
c     compute model size,eigenvalue,aic,bic, objective function
      do 10 i=1,(nka*nlmda)
         call fmeabo(odf(i),oevidx(i),oaic(i),obic(i),oobj(i),
     +        ocoef(:,i),y,x,sz,n,q,p,okas(i),olmdas(i))
 10   continue
!     change the penalized ocoefficients back
      do 10005 i=1,p
         ocoef(q+i,:)=as(i)*ocoef(q+i,:)
10005 continue
      end




c     Hybrid path, update the non-zero lasso only
c***********************************************************************
      subroutine adpmcpkp2(olmdas,okas,ocoef,oaic,obic,oobj,odf,oevidx,
     +     y,x,z,n,q,p,nka,maxka,nlmda,minlmda,epsilon,maxit)
      integer n,q,p,nka,nlmda,maxit,odf(nka*nlmda),oevidx(nka*nlmda)
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),oaic(nka*nlmda),obic(nka*nlmda),
     +     oobj(nka*nlmda),y(n),x(n,q),z(n,p),
     +     maxka,minlmda,epsilon
      integer i,j,k,upidx(p)
      double precision as(p),sz(n,p),kas(nka),unitka,unitlmda,
     +     lmdas(nlmda),alpha(q),beta(p),eta(n),
     +     inia(q,nlmda),inib(p,nlmda),etamat(n,nlmda)
c     standardization of Z
      call standard(as,sz,z,n,p)
c     calculate kappas
      if (nka .eq. 1) then
         kas(1)=0.d0
      else
         unitka=maxka/dble(nka-1)
         do 10000 i=1,nka
            kas(i)=dble(i-1)*unitka
10000    continue
      endif
c     calculate lambdas
      unitlmda=log(minlmda)/dble(nlmda-1)
      call maxbi(lmdas(1),alpha,beta,eta,y,x,sz,n,q,p,
     +     epsilon,maxit)
      do 10001 i=2,nlmda
         lmdas(i)=lmdas(1)*exp(unitlmda*dble(i-1))
10001 continue
c     Lasso solution path as initials
      inia(:,1)=alpha(:)
      inib(:,1)=beta(:)
      etamat(:,1)=eta(:)
      do 10002 i=2,nlmda
         call lasso(alpha,beta,eta,lmdas(i),y,x,sz,n,q,p,
     +        epsilon,maxit)
         inia(:,i)=alpha(:)
         inib(:,i)=beta(:)
         etamat(:,i)=eta(:)
10002 continue
c     MCP solution path along kappa
      if (nka .eq. 1) then
         olmdas(:)=lmdas(:)
         okas(1:nlmda)=kas(1)
         ocoef(1:q,:)=inia(:,:)
         ocoef((q+1):(q+p),:)=inib(:,:)
      else
         do 1003 i=1,nlmda
            olmdas(((i-1)*nka+1):(i*nka))=lmdas(i)
            okas(((i-1)*nka+1):(i*nka))=kas(:)
 1003    continue
!     ... i=1 case
         i=1
         alpha(:)=inia(:,i)
         beta(:) =inib(:,i)
         eta(:)  =etamat(:,i)
         do 1005 j=1,nka
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1005    continue
         do 1004 i=2,nlmda
            alpha(:)=inia(:,i)
            beta(:) =inib(:,i)
            eta(:)  =etamat(:,i)
c     ...... j=1 case
            j=1
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
            upidx(:)=1
            do 9000 k=1,p
               if (beta(k) .eq. 0.d0) upidx(k)=0
 9000       continue
            do 1006 j=2,nka
               call adpmcp2(alpha,beta,eta,upidx,lmdas(i),kas(j),y,x,sz,
     +              n,q,p,epsilon,maxit)
               ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
               ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1006       continue
 1004    continue
      endif
c     compute model size,eigenvalue,aic,bic, objective function
      do 10 i=1,(nka*nlmda)
         call fmeabo(odf(i),oevidx(i),oaic(i),obic(i),oobj(i),
     +        ocoef(:,i),y,x,sz,n,q,p,okas(i),olmdas(i))
 10   continue
!     change the penalized ocoefficients back
      do 10005 i=1,p
         ocoef(q+i,:)=as(i)*ocoef(q+i,:)
10005 continue
      end



c     solution path along lambda
c***********************************************************************
      subroutine adpmcplm(olmdas,okas,ocoef,oaic,obic,oobj,odf,oevidx,
     +     y,x,z,n,q,p,nka,maxka,nlmda,minlmda,epsilon,maxit)
      integer n,q,p,nka,nlmda,maxit,odf(nka*nlmda),oevidx(nka*nlmda)
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),oaic(nka*nlmda),obic(nka*nlmda),
     +     oobj(nka*nlmda),y(n),x(n,q),z(n,p),
     +     maxka,minlmda,epsilon
      integer i,j
      double precision as(p),sz(n,p),kas(nka),unitka,unitlmda,
     +     lmdas(nlmda),inia(q),inib(p),inieta(n),
     +     alpha(q),beta(p),eta(n)
c     standardization of Z
      call standard(as,sz,z,n,p)
c     calculate kappas
      if (nka .eq. 1) then
         kas(1)=0.d0
      else
         unitka=maxka/dble(nka-1)
         do 10000 i=1,nka
            kas(i)=dble(i-1)*unitka
10000    continue
      endif
c     calculate lambdas
      unitlmda=log(minlmda)/dble(nlmda-1)
      call maxbi(lmdas(1),inia,inib,inieta,y,x,sz,n,q,p,
     +     epsilon,maxit)
      do 10001 i=2,nlmda
         lmdas(i)=lmdas(1)*exp(unitlmda*dble(i-1))
10001 continue
c     output lambdas and kappas
      do 10002 i=1,nlmda
         olmdas( ((i-1)*nka+1):(i*nka) )=lmdas(i)
         okas( ((i-1)*nka+1):(i*nka) )=kas(:)
10002 continue
c     Lasso solution path
      if (nka .eq. 1) then
         alpha(:)=inia(:)
         beta(:)=inib(:)
         eta(:)=inieta(:)
         ocoef(1:q,1)=alpha(:)
         ocoef((q+1):(q+p),1)=beta(:)
         do 10004 i=2,nlmda
            call lasso(alpha,beta,eta,lmdas(i),y,x,sz,n,q,p,
     +           epsilon,maxit)
            ocoef(1:q,i)=alpha(:)
            ocoef((q+1):(q+p),i)=beta(:)
10004    continue
c     MCP solution path along lambda given kappa
      else
         alpha(:)=inia(:)
         beta(:)=inib(:)
         eta(:)=inieta(:)
!     i=1 case
!     .. j=1 case
         j=1
         ocoef(1:q,        1+(j-1)*nka)=alpha(:)
         ocoef((q+1):(q+p),1+(j-1)*nka)=beta(:)
         do 1005 j=2,nlmda
            call lasso(alpha,beta,eta,lmdas(j),y,x,sz,n,q,p,
     +           epsilon,maxit)
            ocoef(1:q,        1+(j-1)*nka)=alpha(:)
            ocoef((q+1):(q+p),1+(j-1)*nka)=beta(:)
 1005    continue
         do 10005 i=2,nka
            alpha(:)=inia(:)
            beta(:)=inib(:)
            eta(:)=inieta(:)
            ocoef(1:q,        i)=alpha(:)
            ocoef((q+1):(q+p),i)=beta(:)
            do 10006 j=2,nlmda
               call adpmcp(alpha,beta,eta,lmdas(j),kas(i),y,x,sz,
     +              n,q,p,epsilon,maxit)
               ocoef(1:q,        i+(j-1)*nka)=alpha(:)
               ocoef((q+1):(q+p),i+(j-1)*nka)=beta(:)
10006       continue
10005    continue
      endif
c     compute model size,eigenvalue,aic,bic, objective function
      do 10 i=1,(nka*nlmda)
         call fmeabo(odf(i),oevidx(i),oaic(i),obic(i),oobj(i),
     +        ocoef(:,i),y,x,sz,n,q,p,okas(i),olmdas(i))
 10   continue
!     change the penalized ocoefficients back
      do 10008 i=1,p
         ocoef(q+i,:)=as(i)*ocoef(q+i,:)
10008 continue
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     parallel to LEVEL IV FUNCTIONS
c     CV component
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     along kapa 
c***********************************************************************
      subroutine adpcvkp(pauc,full,cvx,lmdas,kas,nlmda,nka,cvtn,cvpn,
     +     q,p,cvty,cvtx,cvtz,epsilon,maxit,cvpy,cvpx,cvpz)
      integer nlmda,nka,cvtn,cvpn,q,p,
     +     full(nka*nlmda),cvx(nka*nlmda)
      double precision pauc(nlmda*nka),lmdas(nlmda),kas(nka),
     +     cvty(cvtn),cvtx(cvtn,q),cvtz(cvtn,p),epsilon,cvpy(cvpn),
     +     cvpx(cvpn,q),cvpz(cvpn,p)
c     local vars
      integer j,i
      double precision cvtas(p),cvtsz(cvtn,p),
     +     nulllmda,cvalpha(q),cvbeta(p),cvteta(cvtn),
     +     cvtetam(cvtn,nlmda),
     +     inicva(q,nlmda),inicvb(p,nlmda),
     +     cvas(q,nlmda*nka),cvbs(p,nlmda*nka),
     +     cvpeta(cvpn),cvpi(cvpn),tmp
!     column-wise standardization
      call standard(cvtas,cvtsz,cvtz,cvtn,p)
c     solution path along lambdas, initial values for along kas
      call maxbi(nulllmda,cvalpha,cvbeta,cvteta,cvty,cvtx,cvtsz,
     +     cvtn,q,p,epsilon,maxit)
      do 1004 i=1,nlmda
         if (lmdas(i) .ge. nulllmda) then
            inicva(:,i)=cvalpha(:)
            inicvb(:,i)=cvbeta(:)
            cvtetam(:,i)=cvteta(:)
         else
            call lasso(cvalpha,cvbeta,cvteta,lmdas(i),cvty,
     +           cvtx,cvtsz,cvtn,q,p,epsilon,maxit)
            inicva(:,i)=cvalpha(:)
            inicvb(:,i)=cvbeta(:)
            cvtetam(:,i)=cvteta(:)
         endif
 1004 continue
c     solution path along kappa
      if (nka .eq. 1) then
         cvas(:,:)=inicva(:,:)
         cvbs(:,:)=inicvb(:,:)
      else
         do 1005 i=1,nlmda
            if (lmdas(i) .ge. nulllmda) then
               do 1006 j=1,nka
                  cvas(:,(i-1)*nka+j)=inicva(:,i)
                  cvbs(:,(i-1)*nka+j)=inicvb(:,i)
 1006          continue
            else
               cvalpha(:)=inicva(:,i)
               cvbeta(:) =inicvb(:,i)
               cvteta(:)=cvtetam(:,i)
c     ....     j=1 case 
               j=1
               cvas(:,(i-1)*nka+j)=cvalpha(:)
               cvbs(:,(i-1)*nka+j)=cvbeta(:)
               do 1007 j=2,nka
                  call adpmcp(cvalpha,cvbeta,cvteta,lmdas(i),kas(j),
     +                 cvty,cvtx,cvtsz,cvtn,q,p,
     +                 epsilon,maxit)
                  cvas(:,(i-1)*nka+j)=cvalpha(:)
                  cvbs(:,(i-1)*nka+j)=cvbeta(:)
 1007          continue
            endif
 1005    continue
      endif
!     determin whether satured or not 0: Not, 1: satured
      do 10006 i=1,nka*nlmda
         call fme(full(i),cvx(i),cvas(:,i),cvbs(:,i),
     +        cvty,cvtx,cvtsz,cvtn,q,p)
10006 continue
!     change the penalized coefficients back
      do 10005 i=1,p
         cvbs(i,:)=cvtas(i)*cvbs(i,:)
10005 continue
c     get the paucs for the prediction set
      do 1011 i=1,(nka*nlmda)
         do 1012 j=1,cvpn
            cvpeta(j)=dot_product(cvpx(j,:),cvas(:,i))+dot_product(
     +           cvpz(j,:),cvbs(:,i))
 1012    continue
         do 1013 j=1,cvpn
            cvpi(j)=1.d0/(1.d0+exp(-cvpeta(j)))
1013     continue
         call aucroc(tmp,cvpy,cvpi,cvpn)
         pauc(i)=tmp
1011  continue
      end



c     hybrid path: only update the non-zero lasso solution
c***********************************************************************
      subroutine adpcvkp2(pauc,full,cvx,lmdas,kas,nlmda,nka,cvtn,cvpn,
     +     q,p,cvty,cvtx,cvtz,epsilon,maxit,cvpy,cvpx,cvpz)
      integer nlmda,nka,cvtn,cvpn,q,p,
     +     full(nka*nlmda),cvx(nka*nlmda)
      double precision pauc(nlmda*nka),lmdas(nlmda),kas(nka),
     +     cvty(cvtn),cvtx(cvtn,q),cvtz(cvtn,p),epsilon,cvpy(cvpn),
     +     cvpx(cvpn,q),cvpz(cvpn,p)
c     local vars
      integer j,i,k,cvupidx(p)
      double precision cvtas(p),cvtsz(cvtn,p),
     +     nulllmda,cvalpha(q),cvbeta(p),cvteta(cvtn),
     +     cvtetam(cvtn,nlmda),
     +     inicva(q,nlmda),inicvb(p,nlmda),
     +     cvas(q,nlmda*nka),cvbs(p,nlmda*nka),
     +     cvpeta(cvpn),cvpi(cvpn),tmp
!     column-wise standardization
      call standard(cvtas,cvtsz,cvtz,cvtn,p)
c     solution path along lambdas, initial values for along kas
      call maxbi(nulllmda,cvalpha,cvbeta,cvteta,cvty,cvtx,cvtsz,
     +     cvtn,q,p,epsilon,maxit)
      do 1004 i=1,nlmda
         if (lmdas(i) .ge. nulllmda) then
            inicva(:,i)=cvalpha(:)
            inicvb(:,i)=cvbeta(:)
            cvtetam(:,i)=cvteta(:)
         else
            call lasso(cvalpha,cvbeta,cvteta,lmdas(i),cvty,
     +           cvtx,cvtsz,cvtn,q,p,epsilon,maxit)
            inicva(:,i)=cvalpha(:)
            inicvb(:,i)=cvbeta(:)
            cvtetam(:,i)=cvteta(:)
         endif
 1004 continue
c     solution path along kappa
      if (nka .eq. 1) then
         cvas(:,:)=inicva(:,:)
         cvbs(:,:)=inicvb(:,:)
      else
         do 1005 i=1,nlmda
            if (lmdas(i) .ge. nulllmda) then
               do 1006 j=1,nka
                  cvas(:,(i-1)*nka+j)=inicva(:,i)
                  cvbs(:,(i-1)*nka+j)=inicvb(:,i)
 1006          continue
            else
               cvalpha(:)=inicva(:,i)
               cvbeta(:) =inicvb(:,i)
               cvteta(:)=cvtetam(:,i)
c     ....     j=1 case 
               j=1
               cvas(:,(i-1)*nka+j)=cvalpha(:)
               cvbs(:,(i-1)*nka+j)=cvbeta(:)
               cvupidx(:)=1
               do 9000 k=1,p
                  if (cvbeta(k) .eq. 0.d0) cvupidx(k)=0
 9000          continue
               do 1007 j=2,nka
                  call adpmcp2(cvalpha,cvbeta,cvteta,cvupidx,
     +                 lmdas(i),kas(j),cvty,cvtx,cvtsz,cvtn,q,p,
     +                 epsilon,maxit)
                  cvas(:,(i-1)*nka+j)=cvalpha(:)
                  cvbs(:,(i-1)*nka+j)=cvbeta(:)
 1007          continue
            endif
 1005    continue
      endif
!     determin whether satured or not 0: Not, 1: satured
      do 10006 i=1,nka*nlmda
         call fme(full(i),cvx(i),cvas(:,i),cvbs(:,i),
     +        cvty,cvtx,cvtsz,cvtn,q,p)
10006 continue
!     change the penalized coefficients back
      do 10005 i=1,p
         cvbs(i,:)=cvtas(i)*cvbs(i,:)
10005 continue
c     get the paucs for the prediction set
      do 1011 i=1,(nka*nlmda)
         do 1012 j=1,cvpn
            cvpeta(j)=dot_product(cvpx(j,:),cvas(:,i))+dot_product(
     +           cvpz(j,:),cvbs(:,i))
 1012    continue
         do 1013 j=1,cvpn
            cvpi(j)=1.d0/(1.d0+exp(-cvpeta(j)))
1013     continue
         call aucroc(tmp,cvpy,cvpi,cvpn)
         pauc(i)=tmp
1011  continue
      end


c     Along the Lambda for MCP
c***********************************************************************
      subroutine adpcvlm(pauc,full,cvx,lmdas,kas,nlmda,nka,cvtn,cvpn,
     +     q,p,cvty,cvtx,cvtz,epsilon,maxit,cvpy,cvpx,cvpz)
      integer nlmda,nka,cvtn,cvpn,q,p,
     +     full(nka*nlmda),cvx(nka*nlmda)
      double precision pauc(nlmda*nka),lmdas(nlmda),kas(nka),
     +     cvty(cvtn),cvtx(cvtn,q),cvtz(cvtn,p),epsilon,cvpy(cvpn),
     +     cvpx(cvpn,q),cvpz(cvpn,p)
c     local vars
      integer j,i,k,cvupidx(p)
      double precision cvtas(p),cvtsz(cvtn,p),
     +     nulllmda,cvinia(q),cvinib(p),cvinieta(cvtn),
     +     cvalpha(q),cvbeta(p),cvteta(cvtn),
     +     cvas(q,nlmda*nka),cvbs(p,nlmda*nka),
     +     cvpeta(cvpn),cvpi(cvpn),tmp
!     column-wise standardization
      call standard(cvtas,cvtsz,cvtz,cvtn,p)
c     solution path along lambdas, initial values for along kas
      call maxbi(nulllmda,cvinia,cvinib,cvinieta,cvty,cvtx,cvtsz,
     +     cvtn,q,p,epsilon,maxit)
c     Lasso solution path
      if (nka .eq. 1) then 
         cvalpha(:)=cvinia(:)
         cvbeta(:)=cvinib(:)
         cvteta(:)=cvinieta(:)
         do 1004 i=1,nlmda
            if (lmdas(i) .ge. nulllmda) then
               cvas(:,i)=cvalpha(:)
               cvbs(:,i)=cvbeta(:)
            else
               call lasso(cvalpha,cvbeta,cvteta,lmdas(i),cvty,
     +              cvtx,cvtsz,cvtn,q,p,epsilon,maxit)
               cvas(:,i)=cvalpha(:)
               cvbs(:,i)=cvbeta(:)
            endif
 1004    continue
c     MCP solution along lambda given kappa
      else
!     i=1
         cvalpha(:)=cvinia(:)
         cvbeta(:)=cvinib(:)
         cvteta(:)=cvinieta(:)
         do 10004 j=1,nlmda
            if (lmdas(j) .ge. nulllmda) then
               cvas(:,1+(j-1)*nka)=cvalpha(:)
               cvbs(:,1+(j-1)*nka)=cvbeta(:)
            else
               call lasso(cvalpha,cvbeta,cvteta,lmdas(j),cvty,
     +              cvtx,cvtsz,cvtn,q,p,epsilon,maxit)
               cvas(:,1+(j-1)*nka)=cvalpha(:)
               cvbs(:,1+(j-1)*nka)=cvbeta(:)
            endif
10004    continue
!     i=2...
         do 11005 i=2,nka
            cvalpha(:)=cvinia(:)
            cvbeta(:)=cvinib(:)
            cvteta(:)=cvinieta(:)
            do 11006 j=1,nlmda
               if (lmdas(j) .ge. nulllmda) then
                  cvas(:,i+(j-1)*nka)=cvalpha(:)
                  cvbs(:,i+(j-1)*nka)=cvbeta(:)
               else
                  call adpmcp(cvalpha,cvbeta,cvteta,lmdas(j),kas(i),
     +                 cvty,cvtx,cvtsz,cvtn,q,p,epsilon,maxit)
                  cvas(:,i+(j-1)*nka)=cvalpha(:)
                  cvbs(:,i+(j-1)*nka)=cvbeta(:)
               endif
11006       continue   
11005    continue
      endif
!     determin whether satured or not 0: Not, 1: satured
      do 10006 i=1,nka*nlmda
         call fme(full(i),cvx(i),cvas(:,i),cvbs(:,i),
     +        cvty,cvtx,cvtsz,cvtn,q,p)
10006 continue
!     change the penalized coefficients back
      do 10005 i=1,p
         cvbs(i,:)=cvtas(i)*cvbs(i,:)
10005 continue
c     get the paucs for the prediction set
      do 1011 i=1,(nka*nlmda)
         do 1012 j=1,cvpn
            cvpeta(j)=dot_product(cvpx(j,:),cvas(:,i))+dot_product(
     +           cvpz(j,:),cvbs(:,i))
 1012    continue
         do 1013 j=1,cvpn
            cvpi(j)=1.d0/(1.d0+exp(-cvpeta(j)))
1013     continue
         call aucroc(tmp,cvpy,cvpi,cvpn)
         pauc(i)=tmp
1011  continue
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Parallel to Level V function
c     Tuning parameter selection based on CV-AUC
c     only select on df, do not consider convexity
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     


c     along kappa
c***********************************************************************
      subroutine adpcvauckp(out,pauc,
     +     olmdas,okas,ocoef,oaic,obic,oobj,odf,ocvx,ofull,
     +     cvcvx,cvfull,
     +     nindex,cvk,y,x,z,n,q,p,nka,maxka,nlmda,minlmda,epsilon,maxit)
      integer cvk,n,q,p,nka,nlmda,nindex(n),maxit,
     +     odf(nka*nlmda),ocvx(nka*nlmda),ofull(nka*nlmda),
     +     cvcvx(nka*nlmda),cvfull(nka*nlmda)
      double precision out(3+q+p),pauc(nka*nlmda),olmdas(nka*nlmda),
     +     okas(nka*nlmda),ocoef(q+p,nka*nlmda),oaic(nka*nlmda),
     +     obic(nka*nlmda),oobj(nka*nlmda),
     +     y(n),x(n,q),z(n,p),maxka,minlmda,epsilon
c     local vars
      integer i,cv,cvtn,cvpn,k,
     +     cvfullm(nka*nlmda,cvk),cvcvxm(nka*nlmda,cvk),
     +     use,loc
      double precision lmdas(nlmda),kas(nka),
     +     cvpy(n),cvpx(n,q),cvpz(n,p),cvty(n),cvtx(n,q),cvtz(n,p),
     +     cvpauc(nka*nlmda,cvk),
     +     omlmdas(nka*nlmda),omkas(nka*nlmda),omcoef(q+p,nka*nlmda),
     +     ompauc(nka*nlmda)
c     solution path of original dataset
      call adpmcpkp(olmdas,okas,ocoef,oaic,obic,oobj,odf,ocvx,
     +     y,x,z,n,q,p,nka,maxka,nlmda,minlmda,epsilon,maxit)
c     Determine sature model, 0:not, 1: satured
      do 10000 i=1,nka*nlmda
         if ( odf(i) .le. n) then
            ofull(i)=0
         else
            ofull(i)=1
         endif
10000 continue
c     get lambdas,kappas
      do 10005 i=1,nlmda
         lmdas(i)= olmdas((i-1)*nka+1)
10005 continue
      kas(:)=okas(1:nka)
c     start cross validation
      do 10006 cv=1,cvk
         cvtn=0
         cvpn=0
         do 10007 k=1,n
            if (nindex(k) .eq. cv) then
               cvpn=1+cvpn
               cvpy(cvpn)=y(k)
               cvpx(cvpn,:)=x(k,:)
               cvpz(cvpn,:)=z(k,:)
            else
               cvtn=1+cvtn
               cvty(cvtn)=y(k)
               cvtx(cvtn,:)=x(k,:)
               cvtz(cvtn,:)=z(k,:)
            end if
10007    continue
c     .... solve alpha, beta, and pauc for cv datasets
         call adpcvkp(cvpauc(:,cv),cvfullm(:,cv),cvcvxm(:,cv),lmdas,kas,
     +        nlmda,nka,cvtn,cvpn,q,p,
     +        cvty(1:cvtn),cvtx(1:cvtn,:),cvtz(1:cvtn,:),
     +        epsilon,maxit,
     +        cvpy(1:cvpn),cvpx(1:cvpn,:),cvpz(1:cvpn,:))
10006 continue
!     get average pauc, satured, convex info  in cross-validation
      do 10008 i=1,(nlmda*nka)
         pauc(i)=sum(cvpauc(i,:))/dble(cvk)
         if ( sum(cvfullm(i,:)) .eq. 0) then
            cvfull(i)=0
         else
            cvfull(i)=1
         endif
!     .....if ( sum(cvcvx(i,:)) .eq. 5) cvcvx(i)=1
10008 continue
!     assign useful output
      use=0
      do 10009 i=1,nka*nlmda
         if ( (ofull(i) .eq. 0) .and. (cvfull(i) .eq. 0) ) then
            use=use+1
            ompauc(use) =pauc(i)
            omlmdas(use)=olmdas(i)
            omkas(use)  =okas(i)
            omcoef(:,use)=ocoef(:,i)
         endif
10009 continue
!     get maximum
      loc=maxloc( ompauc(1:use), use )
      out(1)=ompauc(loc)
      out(2)=omlmdas(loc)
      out(3)=omkas(loc)
      out(4:(3+q+p))=omcoef(:,loc)
      end


c     hybrid path: update the non-zero only
c***********************************************************************
      subroutine adpcvauckp2(out,pauc,
     +     olmdas,okas,ocoef,oaic,obic,oobj,odf,ocvx,ofull,
     +     cvcvx,cvfull,
     +     nindex,cvk,y,x,z,n,q,p,nka,maxka,nlmda,minlmda,epsilon,maxit)
      integer cvk,n,q,p,nka,nlmda,nindex(n),maxit,
     +     odf(nka*nlmda),ocvx(nka*nlmda),ofull(nka*nlmda),
     +     cvcvx(nka*nlmda),cvfull(nka*nlmda)
      double precision out(3+q+p),pauc(nka*nlmda),olmdas(nka*nlmda),
     +     okas(nka*nlmda),ocoef(q+p,nka*nlmda),oaic(nka*nlmda),
     +     obic(nka*nlmda),oobj(nka*nlmda),
     +     y(n),x(n,q),z(n,p),maxka,minlmda,epsilon
c     local vars
      integer i,cv,cvtn,cvpn,k,
     +     cvfullm(nka*nlmda,cvk),cvcvxm(nka*nlmda,cvk),
     +     use,loc
      double precision lmdas(nlmda),kas(nka),
     +     cvpy(n),cvpx(n,q),cvpz(n,p),cvty(n),cvtx(n,q),cvtz(n,p),
     +     cvpauc(nka*nlmda,cvk),
     +     omlmdas(nka*nlmda),omkas(nka*nlmda),omcoef(q+p,nka*nlmda),
     +     ompauc(nka*nlmda)
c     solution path of original dataset
      call adpmcpkp2(olmdas,okas,ocoef,oaic,obic,oobj,odf,ocvx,
     +     y,x,z,n,q,p,nka,maxka,nlmda,minlmda,epsilon,maxit)
c     Determine sature model, 0:not, 1: satured
      do 10000 i=1,nka*nlmda
         if ( odf(i) .le. n) then
            ofull(i)=0
         else
            ofull(i)=1
         endif
10000 continue
c     get lambdas,kappas
      do 10005 i=1,nlmda
         lmdas(i)= olmdas((i-1)*nka+1)
10005 continue
      kas(:)=okas(1:nka)
c     start cross validation
      do 10006 cv=1,cvk
         cvtn=0
         cvpn=0
         do 10007 k=1,n
            if (nindex(k) .eq. cv) then
               cvpn=1+cvpn
               cvpy(cvpn)=y(k)
               cvpx(cvpn,:)=x(k,:)
               cvpz(cvpn,:)=z(k,:)
            else
               cvtn=1+cvtn
               cvty(cvtn)=y(k)
               cvtx(cvtn,:)=x(k,:)
               cvtz(cvtn,:)=z(k,:)
            end if
10007    continue
c     .... solve alpha, beta, and pauc for cv datasets
         call adpcvkp2(cvpauc(:,cv),cvfullm(:,cv),cvcvxm(:,cv),lmdas,
     +        kas,nlmda,nka,cvtn,cvpn,q,p,
     +        cvty(1:cvtn),cvtx(1:cvtn,:),cvtz(1:cvtn,:),
     +        epsilon,maxit,
     +        cvpy(1:cvpn),cvpx(1:cvpn,:),cvpz(1:cvpn,:))
10006 continue
!     get average pauc, satured, convex info  in cross-validation
      do 10008 i=1,(nlmda*nka)
         pauc(i)=sum(cvpauc(i,:))/dble(cvk)
         if ( sum(cvfullm(i,:)) .eq. 0) then
            cvfull(i)=0
         else
            cvfull(i)=1
         endif
!     .....if ( sum(cvcvx(i,:)) .eq. 5) cvcvx(i)=1
10008 continue
!     assign useful output
      use=0
      do 10009 i=1,nka*nlmda
         if ( (ofull(i) .eq. 0) .and. (cvfull(i) .eq. 0) ) then
            use=use+1
            ompauc(use) =pauc(i)
            omlmdas(use)=olmdas(i)
            omkas(use)  =okas(i)
            omcoef(:,use)=ocoef(:,i)
         endif
10009 continue
!     get maximum
      loc=maxloc( ompauc(1:use), use )
      out(1)=ompauc(loc)
      out(2)=omlmdas(loc)
      out(3)=omkas(loc)
      out(4:(3+q+p))=omcoef(:,loc)
      end


c     along lambda
c***********************************************************************
      subroutine adpcvauclm(out,pauc,
     +     olmdas,okas,ocoef,oaic,obic,oobj,odf,ocvx,ofull,
     +     cvcvx,cvfull,
     +     nindex,cvk,y,x,z,n,q,p,nka,maxka,nlmda,minlmda,epsilon,maxit)
      integer cvk,n,q,p,nka,nlmda,nindex(n),maxit,
     +     odf(nka*nlmda),ocvx(nka*nlmda),ofull(nka*nlmda),
     +     cvcvx(nka*nlmda),cvfull(nka*nlmda)
      double precision out(3+q+p),pauc(nka*nlmda),olmdas(nka*nlmda),
     +     okas(nka*nlmda),ocoef(q+p,nka*nlmda),oaic(nka*nlmda),
     +     obic(nka*nlmda),oobj(nka*nlmda),
     +     y(n),x(n,q),z(n,p),maxka,minlmda,epsilon
c     local vars
      integer i,cv,cvtn,cvpn,k,
     +     cvfullm(nka*nlmda,cvk),cvcvxm(nka*nlmda,cvk),
     +     use,loc
      double precision lmdas(nlmda),kas(nka),
     +     cvpy(n),cvpx(n,q),cvpz(n,p),cvty(n),cvtx(n,q),cvtz(n,p),
     +     cvpauc(nka*nlmda,cvk),
     +     omlmdas(nka*nlmda),omkas(nka*nlmda),omcoef(q+p,nka*nlmda),
     +     ompauc(nka*nlmda)
c     solution path of original dataset
      call adpmcplm(olmdas,okas,ocoef,oaic,obic,oobj,odf,ocvx,
     +     y,x,z,n,q,p,nka,maxka,nlmda,minlmda,epsilon,maxit)
c     Determine sature model, 0:not, 1: satured
      do 10000 i=1,nka*nlmda
         if ( odf(i) .le. n) then
            ofull(i)=0
         else
            ofull(i)=1
         endif
10000 continue
c     get lambdas,kappas
      do 10005 i=1,nlmda
         lmdas(i)= olmdas((i-1)*nka+1)
10005 continue
      kas(:)=okas(1:nka)
c     start cross validation
      do 10006 cv=1,cvk
         cvtn=0
         cvpn=0
         do 10007 k=1,n
            if (nindex(k) .eq. cv) then
               cvpn=1+cvpn
               cvpy(cvpn)=y(k)
               cvpx(cvpn,:)=x(k,:)
               cvpz(cvpn,:)=z(k,:)
            else
               cvtn=1+cvtn
               cvty(cvtn)=y(k)
               cvtx(cvtn,:)=x(k,:)
               cvtz(cvtn,:)=z(k,:)
            end if
10007    continue
c     .... solve alpha, beta, and pauc for cv datasets
         call adpcvlm(cvpauc(:,cv),cvfullm(:,cv),cvcvxm(:,cv),lmdas,kas,
     +        nlmda,nka,cvtn,cvpn,q,p,
     +        cvty(1:cvtn),cvtx(1:cvtn,:),cvtz(1:cvtn,:),
     +        epsilon,maxit,
     +        cvpy(1:cvpn),cvpx(1:cvpn,:),cvpz(1:cvpn,:))
10006 continue
!     get average pauc, satured, convex info  in cross-validation
      do 10008 i=1,(nlmda*nka)
         pauc(i)=sum(cvpauc(i,:))/dble(cvk)
         if ( sum(cvfullm(i,:)) .eq. 0) then
            cvfull(i)=0
         else
            cvfull(i)=1
         endif
!     .....if ( sum(cvcvx(i,:)) .eq. 5) cvcvx(i)=1
10008 continue
!     assign useful output
      use=0
      do 10009 i=1,nka*nlmda
         if ( (ofull(i) .eq. 0) .and. (cvfull(i) .eq. 0) ) then
            use=use+1
            ompauc(use) =pauc(i)
            omlmdas(use)=olmdas(i)
            omkas(use)  =okas(i)
            omcoef(:,use)=ocoef(:,i)
         endif
10009 continue
!     get maximum
      loc=maxloc( ompauc(1:use), use )
      out(1)=ompauc(loc)
      out(2)=omlmdas(loc)
      out(3)=omkas(loc)
      out(4:(3+q+p))=omcoef(:,loc)
      end










cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     LLA2 approach,
c     strictly follow the approach by Zou,H., Li,R., 2008 Ann Stat
c     only applies to n>p data
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Calculate lambda_max, byproduct:alpha,beta,xinvxtx,r
c***********************************************************************
      subroutine maxga(lmdamax,beta,r,y,z,n,p)
      integer n,p
      double precision lmdamax,beta(p),r(n),y(n),z(n,p)
c     local vars
      integer i
      double precision tmp(p)
      beta(:)=0.d0
      r(:)=y(:)
c     tmp=z'_{ij}r/n
      do 1202 i=1,p
         tmp(i)=abs( (dot_product(z(:,i),r))/dble(n) )
 1202 continue
      lmdamax=maxval(tmp)
      end

c     Lasso solution
c***********************************************************************
      subroutine iniga(beta,r,lmda,y,z,n,p,epsilon,maxit)
      integer n,p,maxit
      double precision beta(p),r(n),lmda,y(n),z(n,p),epsilon
c     local vars
      integer count,j,tag
      double precision betaold(p),m,numor
c     begin of iteration
      count=0
 1204 continue
      betaold(:)=beta(:)
c     .. update of beta
      do 1205 j=1,p
         m=beta(j) + dot_product(z(:,j),r)/dble(n)
         call soft(numor,m,lmda)
         beta(j)=numor
c     ...update residuals, r(s+1)=r(s)+z(b(s)-b(s+1))
         r(:)=r(:)+z(:,j)*( betaold(j)-beta(j) )
 1205 continue
c     check convergence
      count=count+1
      call converge(tag,beta,betaold,p,epsilon)
      if (tag .eq. 1) go to 1210
      if (count .ge. maxit) call rexit("Lasso solution diverges! \n")
      if (count .lt. maxit) go to 1204
 1210 continue
c     end of iteration loop
      end

c***********************************************************************
      subroutine lassoga(ocoef,y,z,n,p,lmda,epsilon,maxit)
      integer n,p,maxit
      double precision ocoef(p),y(n),z(n,p),lmda,epsilon
c     local vars
      integer i,j
      double precision as(p),sz(n,p),lmdamax,beta(p),r(n)
c     Column-wise standardization of Z
      call standard(as,sz,z,n,p)
c     calculate lambda_max
      call maxga(lmdamax,beta,r,y,sz,n,p)
C     compute the solution for a given lmda
      if (lmda .ge. lmdamax) then
         beta(:)=0.d0
      else
         call iniga(beta,r,lmda,y,sz,n,p,epsilon,maxit)
      endif
!     change the penalized coefficients back
      do 10005 i=1,p
         ocoef(i)=as(i)*beta(i)
10005 continue
      end

c     lla approach for compute mcp for a given kappa and lambda
c***********************************************************************
      subroutine lla(alpha,beta,eta,lmda,ka,y,x,z,n,q,p,
     +     epsilon,maxit)
      integer n,q,p,maxit
      double precision alpha(q),beta(p),eta(n),lmda,ka,y(n),x(n,q),
     +     z(n,p),epsilon
      integer q1,qp,count,j,i,uvidx(q+p),u,v,s,t,info,tag
      double precision xz(n,q+p),ab(q+p),abold(q+p),pi(n),sqrtd(n),
     +     y1(n),xz1(n,q+p),hu(n,n),iu(n,n),y2(n)
      integer, dimension(:), allocatable:: umk
      integer, dimension(:), allocatable:: vmk
      double precision, dimension(:), allocatable:: abu
      double precision, dimension(:), allocatable:: abv
      double precision, dimension(:,:), allocatable:: xu1
      double precision, dimension(:,:), allocatable:: xv1
      double precision, dimension(:,:), allocatable:: rm
      double precision, dimension(:,:), allocatable:: sm
      double precision, dimension(:,:), allocatable:: pu
      double precision, dimension(:,:), allocatable:: xv2
      double precision, dimension(:), allocatable:: consv
c     new forms
      q1=q+1
      qp=q+p
      xz(:,1:q)=x(:,:)
      xz(:,q1:qp)=z(:,:)
      ab(1:q)=alpha(:)
      ab(q1:qp)=beta(:)
c     ... begin of iteration
      count=0
 1000 continue
      abold(:)=ab(:)
c     determine U and V parts
      uvidx(1:q)=1
      do 1004 j=q1,qp
         if ( abs(ab(j)) .ge. lmda/ka ) then
            uvidx(j)=1
         else
            uvidx(j)=0
         endif
 1004 continue
      u=sum(uvidx(:))
      v=qp-u
      if (v .eq. 0) then
         go to 1210
      else
c     compute the eta=xz%*%ab
         call dgemv("N",n,qp,1.d0,xz,n,ab,1,0.d0,eta,1)
c     compute new data, 
         do 1003 i=1,n
            pi(i)=1.d0/( 1.d0 + exp(-eta(i)) )
c     if (pi(i) .lt. 0.0001) pi(i)=0.d0
c     if (pi(i) .gt. 0.9999) pi(i)=1.d0
            sqrtd(i)=sqrt( pi(i)*(1.d0-pi(i)) )
            y1(i)=sqrtd(i)*eta(i)
            xz1(i,:)=sqrtd(i)*xz(i,:)
 1003    continue
c     divide into U and V parts   
         allocate(umk(u))
         allocate(vmk(v))
         s=1
         t=1
         do 1005 j=1,qp
            if (uvidx(j) .eq. 1) then 
               umk(s)=j
               s=s+1
            else
               vmk(t)=j
               t=t+1
            endif
 1005    continue
c     divide the ab, xz1 into U and V part
         allocate(abu(u))
         allocate(abv(v))
         allocate(xu1(n,u))
         allocate(xv1(n,v))
         abu(:)=ab(umk)
         abv(:)=ab(vmk)
         xu1(:,:)=xz1(:,umk)
         xv1(:,:)=xz1(:,vmk)
c     compute X**T * X=R**T * R
         allocate(rm(u,u))
         call dsyrk("U","T",u,n,1.d0,xu1,n,0.d0,rm,u)
         call dpotrf("U",u,rm,u,info)
c     S=R**(-T) * X**T
         allocate(sm(u,n))
         do 1006 i=1,u
            sm(i,:)=xu1(:,i)
 1006    continue
         call dtrtrs("U","T","N",u,n,rm,u,sm,u,info)
c     hu= X * (X**T * X)**(-1) * X**T=S**(T) * S
         call dsyrk("U","T",n,u,1.d0,sm,u,0.d0,hu,n)
c     (X**T * X)**(-1)=R**(-1) * R**(-T)
         call dtrtri("U","N",u,rm,u,info)
c     pu=(X**T * X)**(-1) * X**T=R**(-1) * R**(-T) * X**T=R**(-1) * S
         allocate(pu(u,n))
c     call dgemm("N","N",u,n,u,1.d0,rm,u,sm,u,0.d0,pu,u)
         pu(:,:)=sm(:,:)
         call dtrmm("L","U","N","N",u,n,1.d0,rm,u,pu,u)
c     iu = I - hu
         do 10007 i=1,n
            j=i
            iu(i,j)=1.d0-hu(i,j)
            do 10008 j=(i+1),n
               iu(i,j)= -hu(i,j)
10008       continue
10007    continue
c     compute the new data form, y2=(I-hu)*y=iu*y
         call dsymv("U",n,1.d0,iu,n,y1,1,0.d0,y2,1)
c     xv2=(I-hu)*xv1=iu*xv1
         allocate(xv2(n,v))
         call dsymm("L","U",n,v,1.d0,iu,n,xv1,n,0.d0,xv2,n)
c     xv2=xv2 %*% diag (consv)
         allocate(consv(v))
         do 1008 j=1,v
            consv(j)=1.d0-ka*abs(abv(j))/lmda
            xv2(:,j)=xv2(:,j)/consv(j)
1008     continue
c     use lasso for the solution
         call lassoga(abv,y2,xv2,n,v,lmda,epsilon,maxit)
         do 10009 j=1,v
            abv(j)=abv(j)/consv(j)
10009    continue
c     abu=pu*(y1-xv1*abv)
         call dgemv("N",n,v,-1.d0,xv1 ,n,abv,1,1.d0,y1 ,1)
         call dgemv("N",u,n, 1.d0,pu  ,u,y1 ,1,0.d0,abu,1)
c     reput into original order
         do 10010 j=1,u
            ab(umk(j))=abu(j)
10010    continue
         do 10011 j=1,v
            ab(vmk(j))=abv(j)
10011    continue
         deallocate(umk)
         deallocate(vmk)
         deallocate(abu)
         deallocate(abv)
         deallocate(xu1)
         deallocate(xv1)
         deallocate(rm)
         deallocate(sm)
         deallocate(pu)
         deallocate(xv2)
         deallocate(consv)
c     check the convergence
         count=count+1
         call converge(tag,ab,abold,qp,epsilon)
         if (tag .eq. 1) go to 1210
         if (count .ge. maxit) call rexit("LLA diverges! \n")
         if (count .lt. maxit) go to 1000
      endif
 1210 continue
      alpha(:)=ab(1:q)
      beta(:)=ab(q1:qp)
      end
      
      
c***********************************************************************
c     Main function using LLA1 approach                                    
c***********************************************************************
      subroutine flla(olmdas,okas,ocoef,oaic,obic,oobj,odf,oevidx,
     +     y,x,z,n,q,p,nka,maxka,nlmda,minlmda,epsilon,maxit)
      integer n,q,p,nka,nlmda,maxit,odf(nka*nlmda),oevidx(nka*nlmda)
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),oaic(nka*nlmda),obic(nka*nlmda),
     +     oobj(nka*nlmda),y(n),x(n,q),z(n,p),
     +     maxka,minlmda,epsilon
      integer i,j
      double precision as(p),sz(n,p),kas(nka),unitka,unitlmda,
     +     lmdas(nlmda),alpha(q),beta(p),eta(n),
     +     inia(q,nlmda),inib(p,nlmda),etamat(n,nlmda)
c     standardization of Z
      call standard(as,sz,z,n,p)
c     calculate kappas
      if (nka .eq. 1) then
         kas(1)=0.d0
      else
         unitka=maxka/dble(nka-1)
         do 10000 i=1,nka
            kas(i)=dble(i-1)*unitka
10000    continue
      endif
c     calculate lambdas
      unitlmda=log(minlmda)/dble(nlmda-1)
      call maxbi(lmdas(1),alpha,beta,eta,y,x,sz,n,q,p,
     +     epsilon,maxit)
      do 10001 i=2,nlmda
         lmdas(i)=lmdas(1)*exp(unitlmda*dble(i-1))
10001 continue
c     Lasso solution path as initials
      inia(:,1)=alpha(:)
      inib(:,1)=beta(:)
      etamat(:,1)=eta(:)
      do 10002 i=2,nlmda
         call lasso(alpha,beta,eta,lmdas(i),y,x,sz,n,q,p,
     +        epsilon,maxit)
         inia(:,i)=alpha(:)
         inib(:,i)=beta(:)
         etamat(:,i)=eta(:)
10002 continue
c     MCP solution path along kappa
      if (nka .eq. 1) then
         olmdas(:)=lmdas(:)
         okas(1:nlmda)=kas(1)
         ocoef(1:q,:)=inia(:,:)
         ocoef((q+1):(q+p),:)=inib(:,:)
      else
         do 1003 i=1,nlmda
            olmdas(((i-1)*nka+1):(i*nka))=lmdas(i)
            okas(((i-1)*nka+1):(i*nka))=kas(:)
 1003    continue
!     ... i=1 case
         i=1
         alpha(:)=inia(:,i)
         beta(:) =inib(:,i)
         eta(:)  =etamat(:,i)
         do 1005 j=1,nka
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1005    continue
         do 1004 i=2,nlmda
            alpha(:)=inia(:,i)
            beta(:) =inib(:,i)
            eta(:)  =etamat(:,i)
c     ...... j=1 case
            j=1
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
            do 1006 j=2,nka
               call lla(alpha,beta,eta,lmdas(i),kas(j),y,x,sz,
     +              n,q,p,epsilon,maxit)
               ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
               ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1006       continue
 1004    continue
      endif
c     compute model size,eigenvalue,aic,bic, objective function
      do 10 i=1,(nka*nlmda)
         call fmeabo(odf(i),oevidx(i),oaic(i),obic(i),oobj(i),
     +        ocoef(:,i),y,x,sz,n,q,p,okas(i),olmdas(i))
 10   continue
!     change the penalized ocoefficients back
      do 10005 i=1,p
         ocoef(q+i,:)=as(i)*ocoef(q+i,:)
10005 continue
      end





cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     LLA1 approach
c     Quadratic approximation to loss and LLA approximation to  penalty
c     using the adaptive Lasso to compute concave solution
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



c     along kapa 
c***********************************************************************
      subroutine cvllabi(pauc,full,cvx,lmdas,kas,nlmda,nka,cvtn,cvpn,
     +     q,p,cvty,cvtx,cvtz,epsilon,maxit,cvpy,cvpx,cvpz)
      integer nlmda,nka,cvtn,cvpn,q,p,
     +     full(nka*nlmda),cvx(nka*nlmda)
      double precision pauc(nlmda*nka),lmdas(nlmda),kas(nka),
     +     cvty(cvtn),cvtx(cvtn,q),cvtz(cvtn,p),epsilon,cvpy(cvpn),
     +     cvpx(cvpn,q),cvpz(cvpn,p)
c     local vars
      integer j,i
      double precision cvtas(p),cvtsz(cvtn,p),
     +     nulllmda,cvalpha(q),cvbeta(p),cvteta(cvtn),
     +     cvtetam(cvtn,nlmda),
     +     inicva(q,nlmda),inicvb(p,nlmda),
     +     cvas(q,nlmda*nka),cvbs(p,nlmda*nka),
     +     cvpeta(cvpn),cvpi(cvpn),tmp
!     column-wise standardization
      call standard(cvtas,cvtsz,cvtz,cvtn,p)
c     solution path along lambdas, initial values for along kas
      call maxbi(nulllmda,cvalpha,cvbeta,cvteta,cvty,cvtx,cvtsz,
     +     cvtn,q,p,epsilon,maxit)
      do 1004 i=1,nlmda
         if (lmdas(i) .ge. nulllmda) then
            inicva(:,i)=cvalpha(:)
            inicvb(:,i)=cvbeta(:)
            cvtetam(:,i)=cvteta(:)
         else
            call lasso(cvalpha,cvbeta,cvteta,lmdas(i),cvty,
     +           cvtx,cvtsz,cvtn,q,p,epsilon,maxit)
            inicva(:,i)=cvalpha(:)
            inicvb(:,i)=cvbeta(:)
            cvtetam(:,i)=cvteta(:)
         endif
 1004 continue
c     solution path along kappa
      if (nka .eq. 1) then
         cvas(:,:)=inicva(:,:)
         cvbs(:,:)=inicvb(:,:)
      else
         do 1005 i=1,nlmda
            if (lmdas(i) .ge. nulllmda) then
               do 1006 j=1,nka
                  cvas(:,(i-1)*nka+j)=inicva(:,i)
                  cvbs(:,(i-1)*nka+j)=inicvb(:,i)
 1006          continue
            else
               cvalpha(:)=inicva(:,i)
               cvbeta(:) =inicvb(:,i)
               cvteta(:)=cvtetam(:,i)
c     ....     j=1 case 
               j=1
               cvas(:,(i-1)*nka+j)=cvalpha(:)
               cvbs(:,(i-1)*nka+j)=cvbeta(:)
               do 1007 j=2,nka
                  call llabi(cvalpha,cvbeta,cvteta,lmdas(i),kas(j),
     +                 cvty,cvtx,cvtsz,cvtn,q,p,
     +                 epsilon,maxit)
                  cvas(:,(i-1)*nka+j)=cvalpha(:)
                  cvbs(:,(i-1)*nka+j)=cvbeta(:)
 1007          continue
            endif
 1005    continue
      endif
!     determin whether satured or not 0: Not, 1: satured
      do 10006 i=1,nka*nlmda
         call fme(full(i),cvx(i),cvas(:,i),cvbs(:,i),
     +        cvty,cvtx,cvtsz,cvtn,q,p)
10006 continue
!     change the penalized coefficients back
      do 10005 i=1,p
         cvbs(i,:)=cvtas(i)*cvbs(i,:)
10005 continue
c     get the paucs for the prediction set
      do 1011 i=1,(nka*nlmda)
         do 1012 j=1,cvpn
            cvpeta(j)=dot_product(cvpx(j,:),cvas(:,i))+dot_product(
     +           cvpz(j,:),cvbs(:,i))
 1012    continue
         do 1013 j=1,cvpn
            cvpi(j)=1.d0/(1.d0+exp(-cvpeta(j)))
1013     continue
         call aucroc(tmp,cvpy,cvpi,cvpn)
         pauc(i)=tmp
1011  continue
      end





c     along kappa
c***********************************************************************
      subroutine fcvllabi(out,pauc,
     +     olmdas,okas,ocoef,oaic,obic,oobj,odf,ocvx,ofull,
     +     cvcvx,cvfull,
     +     nindex,cvk,y,x,z,n,q,p,nka,maxka,nlmda,minlmda,epsilon,maxit)
      integer cvk,n,q,p,nka,nlmda,nindex(n),maxit,
     +     odf(nka*nlmda),ocvx(nka*nlmda),ofull(nka*nlmda),
     +     cvcvx(nka*nlmda),cvfull(nka*nlmda)
      double precision out(3+q+p),pauc(nka*nlmda),olmdas(nka*nlmda),
     +     okas(nka*nlmda),ocoef(q+p,nka*nlmda),oaic(nka*nlmda),
     +     obic(nka*nlmda),oobj(nka*nlmda),
     +     y(n),x(n,q),z(n,p),maxka,minlmda,epsilon
c     local vars
      integer i,cv,cvtn,cvpn,k,
     +     cvfullm(nka*nlmda,cvk),cvcvxm(nka*nlmda,cvk),
     +     use,loc
      double precision lmdas(nlmda),kas(nka),
     +     cvpy(n),cvpx(n,q),cvpz(n,p),cvty(n),cvtx(n,q),cvtz(n,p),
     +     cvpauc(nka*nlmda,cvk),
     +     omlmdas(nka*nlmda),omkas(nka*nlmda),omcoef(q+p,nka*nlmda),
     +     ompauc(nka*nlmda)
c     solution path of original dataset
      call fllabi(olmdas,okas,ocoef,oaic,obic,oobj,odf,ocvx,
     +     y,x,z,n,q,p,nka,maxka,nlmda,minlmda,epsilon,maxit)
c     Determine sature model, 0:not, 1: satured
      do 10000 i=1,nka*nlmda
         if ( odf(i) .le. n) then
            ofull(i)=0
         else
            ofull(i)=1
         endif
10000 continue
c     get lambdas,kappas
      do 10005 i=1,nlmda
         lmdas(i)= olmdas((i-1)*nka+1)
10005 continue
      kas(:)=okas(1:nka)
c     start cross validation
      do 10006 cv=1,cvk
         cvtn=0
         cvpn=0
         do 10007 k=1,n
            if (nindex(k) .eq. cv) then
               cvpn=1+cvpn
               cvpy(cvpn)=y(k)
               cvpx(cvpn,:)=x(k,:)
               cvpz(cvpn,:)=z(k,:)
            else
               cvtn=1+cvtn
               cvty(cvtn)=y(k)
               cvtx(cvtn,:)=x(k,:)
               cvtz(cvtn,:)=z(k,:)
            end if
10007    continue
c     .... solve alpha, beta, and pauc for cv datasets
         call cvllabi(cvpauc(:,cv),cvfullm(:,cv),cvcvxm(:,cv),lmdas,kas,
     +        nlmda,nka,cvtn,cvpn,q,p,
     +        cvty(1:cvtn),cvtx(1:cvtn,:),cvtz(1:cvtn,:),
     +        epsilon,maxit,
     +        cvpy(1:cvpn),cvpx(1:cvpn,:),cvpz(1:cvpn,:))
10006 continue
!     get average pauc, satured, convex info  in cross-validation
      do 10008 i=1,(nlmda*nka)
         pauc(i)=sum(cvpauc(i,:))/dble(cvk)
         if ( sum(cvfullm(i,:)) .eq. 0) then
            cvfull(i)=0
         else
            cvfull(i)=1
         endif
!     .....if ( sum(cvcvx(i,:)) .eq. 5) cvcvx(i)=1
10008 continue
!     assign useful output
      use=0
      do 10009 i=1,nka*nlmda
         if ( (ofull(i) .eq. 0) .and. (cvfull(i) .eq. 0) ) then
            use=use+1
            ompauc(use) =pauc(i)
            omlmdas(use)=olmdas(i)
            omkas(use)  =okas(i)
            omcoef(:,use)=ocoef(:,i)
         endif
10009 continue
!     get maximum
      loc=maxloc( ompauc(1:use), use )
      out(1)=ompauc(loc)
      out(2)=omlmdas(loc)
      out(3)=omkas(loc)
      out(4:(3+q+p))=omcoef(:,loc)
      end

















c***********************************************************************
c     Main function using LLA1 approach                                    
c***********************************************************************
      subroutine fllabi(olmdas,okas,ocoef,oaic,obic,oobj,odf,oevidx,
     +     y,x,z,n,q,p,nka,maxka,nlmda,minlmda,epsilon,maxit)
      integer n,q,p,nka,nlmda,maxit,odf(nka*nlmda),oevidx(nka*nlmda)
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),oaic(nka*nlmda),obic(nka*nlmda),
     +     oobj(nka*nlmda),y(n),x(n,q),z(n,p),
     +     maxka,minlmda,epsilon
      integer i,j
      double precision as(p),sz(n,p),kas(nka),unitka,unitlmda,
     +     lmdas(nlmda),alpha(q),beta(p),eta(n),
     +     inia(q,nlmda),inib(p,nlmda),etamat(n,nlmda)
c     standardization of Z
      call standard(as,sz,z,n,p)
c     calculate kappas
      if (nka .eq. 1) then
         kas(1)=0.d0
      else
         unitka=maxka/dble(nka-1)
         do 10000 i=1,nka
            kas(i)=dble(i-1)*unitka
10000    continue
      endif
c     calculate lambdas
      unitlmda=log(minlmda)/dble(nlmda-1)
      call maxbi(lmdas(1),alpha,beta,eta,y,x,sz,n,q,p,
     +     epsilon,maxit)
      do 10001 i=2,nlmda
         lmdas(i)=lmdas(1)*exp(unitlmda*dble(i-1))
10001 continue
c     Lasso solution path as initials
      inia(:,1)=alpha(:)
      inib(:,1)=beta(:)
      etamat(:,1)=eta(:)
      do 10002 i=2,nlmda
         call lasso(alpha,beta,eta,lmdas(i),y,x,sz,n,q,p,
     +        epsilon,maxit)
         inia(:,i)=alpha(:)
         inib(:,i)=beta(:)
         etamat(:,i)=eta(:)
10002 continue
c     MCP solution path along kappa
      if (nka .eq. 1) then
         olmdas(:)=lmdas(:)
         okas(1:nlmda)=kas(1)
         ocoef(1:q,:)=inia(:,:)
         ocoef((q+1):(q+p),:)=inib(:,:)
      else
         do 1003 i=1,nlmda
            olmdas(((i-1)*nka+1):(i*nka))=lmdas(i)
            okas(((i-1)*nka+1):(i*nka))=kas(:)
 1003    continue
!     ... i=1 case
         i=1
         alpha(:)=inia(:,i)
         beta(:) =inib(:,i)
         eta(:)  =etamat(:,i)
         do 1005 j=1,nka
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1005    continue
         do 1004 i=2,nlmda
            alpha(:)=inia(:,i)
            beta(:) =inib(:,i)
            eta(:)  =etamat(:,i)
c     ...... j=1 case
            j=1
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
            do 1006 j=2,nka
               call llabi(alpha,beta,eta,lmdas(i),kas(j),y,x,sz,
     +              n,q,p,epsilon,maxit)
               ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
               ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1006       continue
 1004    continue
      endif
c     compute model size,eigenvalue,aic,bic, objective function
      do 10 i=1,(nka*nlmda)
         call fmeabo(odf(i),oevidx(i),oaic(i),obic(i),oobj(i),
     +        ocoef(:,i),y,x,sz,n,q,p,okas(i),olmdas(i))
 10   continue
!     change the penalized ocoefficients back
      do 10005 i=1,p
         ocoef(q+i,:)=as(i)*ocoef(q+i,:)
10005 continue
      end



c     lla using adaptive lasso method
c***********************************************************************
      subroutine llabi(alpha, beta, eta, lmda, ka, y, x, z, n, q, p,
     +     epsilon, maxit)
      integer n, q, p, maxit
      double precision alpha(q), beta(p), eta(n), lmda, ka, 
     +     y(n), x(n,q), z(n,p), epsilon
      integer q1, qp, count, j, tag
      double precision xz(n,q+p), ab(q+p), abold(q+p), wt(q+p)
c     new forms of data
      q1 = q + 1
      qp = q + p
      xz(:,1:q) = x(:,:)
      xz(:,q1:qp) = z(:,:)
      ab(1:q) = alpha(:)
      ab(q1:qp) = beta(:)
c     ... begin of iteration
      count = 0
 9000 continue
      abold(:) = ab(:)
c     compute the weight for adaptive lasso
      do 1004 j = 1, q
         wt(j) = 0.d0
 1004 continue
      do 1005 j = q1, qp
         if (abs(ab(j)) .ge. (lmda/ka)) then 
            wt(j) = 0.d0
         else
            wt(j) = 1.d0 - (abs(ab(j))*ka/lmda)
         endif
 1005 continue
c     apply the adaptive lasso 
      call adplasbi(ab, y, xz, n, q, p, q1, qp, lmda, wt, 
     +     epsilon, maxit)
c     check the convergence
      count = count + 1
      call converge(tag, ab, abold, qp, epsilon)
      if (tag .eq. 1) go to 10000
      if (count .ge. maxit) call rexit("LLA diverges! \n")
      if (count .lt. maxit) go to 9000
10000 continue
      alpha(:)=ab(1:q)
      beta(:)=ab(q1:qp)
      end


c     adaptive lasso
c***********************************************************************
      subroutine adplasbi(ab, y, xz, n, q, p, q1, qp, lmda, wt, epsilon,
     +     maxit)
      integer n, q, p, q1, qp, maxit
      double precision ab(qp), y(n), xz(n,qp), lmda, wt(qp), epsilon
      integer count, j, i,  tag
      double precision eta(n), abold(qp), pi(n), r(n), var(n), cons(qp),
     +     m, newlmda, numor      
      call dgemv("N", n, qp, 1.d0, xz, n, ab, 1, 0.d0, eta, 1)
      count = 0
c     ... begin of iteration
 9000 continue
      abold(:) = ab(:)
c     update the alpha and beta
      do 1000 j = 1, qp
         cons(j) = 0.d0
         do 1001 i = 1, n
            pi(i) = 1.d0 /( 1.d0 + exp(-eta(i)) )
            if (pi(i) .lt. 0.0001) pi(i)=0.d0
            if (pi(i) .gt. 0.9999) pi(i)=1.d0
            r(i) = y(i) - pi(i)
            var(i) =  pi(i)*(1.d0 - pi(i))
            cons(j) = cons(j) + var(i)*xz(i,j)*xz(i,j)
 1001    continue
         cons(j) = cons(j)/dble(n)
         m = cons(j)*ab(j) + dot_product(xz(:,j), r)/dble(n)
         if (wt(j) .eq. 0.d0) then
            ab(j) = m/cons(j)
         else
            newlmda = lmda*wt(j)
            call soft(numor, m, newlmda)
            ab(j) = numor/cons(j)
         endif
         eta(:) = eta(:) + xz(:,j)* (ab(j) - abold(j))
 1000 continue   
c     check convergence
      count = count + 1
      call converge(tag, ab, abold, qp, epsilon)
      if (tag .eq. 1) go to 10000
      if (count .ge. maxit) call rexit("Adaptive Lasso diverges! \n")
      if (count .lt. maxit) go to 9000
10000 continue
c     end of iteration loop
      end      
















cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     this is a test of 
c     LLA using one-step estimator in loss function and LLA
c     using the adaptive Lasso to compute concave solution
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine adplasso(ab, y, xz, n, q, p, q1, qp, lmda, w, epsilon,
     +     maxit)
      integer n, q, p, q1, qp, maxit
      double precision ab(qp), y(n), xz(n,qp), lmda, w(qp), epsilon
      integer j, count, tag
      double precision cons(qp), r(n), abold(qp), m, newlmda, numor
      do 1000 j = 1, qp
         cons(j) = dot_product(xz(:,j), xz(:,j))/dble(n)
 1000 continue 
c     r = y - xb
      r(:) = y(:)
      call dgemv("N", n, qp, -1.d0, xz, n, ab, 1, 1.d0, r, 1)
      count = 0
c     begin of iteration
 9000 continue
      abold(:) = ab(:)
c     .. update alpha
      do 1001 j = 1, q
         m = cons(j)*ab(j) + dot_product(xz(:,j), r) / dble(n)
         ab(j) = m/cons(j)
         r(:) = r(:) + xz(:,j)*(abold(j) - ab(j))
 1001 continue
c     .. update beta
      do 1002 j = q1, qp
         m = cons(j)*ab(j) + dot_product(xz(:,j), r) / dble(n)
         if (w(j) .gt. 0.d0) then
            newlmda = lmda*w(j)
            call soft(numor, m, newlmda)
            ab(j) =  numor/cons(j)
         else
            ab(j) = m/cons(j)
         endif
         r(:) = r(:) + xz(:,j)*(abold(j) - ab(j))
1002  continue
c     check convergence
      count = count + 1
      call converge(tag, ab, abold, qp, epsilon)
      if (tag .eq. 1) go to 10000
      if (count .ge. maxit) call rexit("Adaptive Lasso diverges! \n")
      if (count .lt. maxit) go to 9000
10000 continue
c     end of iteration loop
      end      



c     lla using adaptive lasso method
c***********************************************************************
      subroutine llaadp(alpha, beta, eta, lmda, ka, y, x, z, n, q, p,
     +     epsilon, maxit)
      integer n, q, p, maxit
      double precision alpha(q), beta(p), eta(n), lmda, ka, 
     +     y(n), x(n,q), z(n,p), epsilon
      integer q1, qp, count, j, i, tag
      double precision xz(n,q+p), ab(q+p), abold(q+p), pi(n), sqrtd(n),
     +     y1(n), xz1(n,q+p), w(q+p)
c     new forms of data
      q1 = q + 1
      qp = q + p
      xz(:,1:q) = x(:,:)
      xz(:,q1:qp) = z(:,:)
      ab(1:q) = alpha(:)
      ab(q1:qp) = beta(:)
c     ... begin of iteration
      count = 0
 9000 continue
      abold(:) = ab(:)
c     compute the eta=xz%*%ab
      call dgemv("N", n, qp, 1.d0, xz, n, ab, 1, 0.d0, eta, 1)
c     compute new data, 
      do 1003 i = 1 ,n
         pi(i) = 1.d0 /( 1.d0 + exp(-eta(i)) )
c     if (pi(i) .lt. 0.0001) pi(i)=0.d0
c     if (pi(i) .gt. 0.9999) pi(i)=1.d0
         sqrtd(i) = sqrt( pi(i)*(1.d0-pi(i)) )
         y1(i) = sqrtd(i)*eta(i)
         xz1(i,:) = sqrtd(i)*xz(i,:)
 1003 continue
c     start adaptive Lasso 
      do 1004 j = 1, q
         w(j) = 0.d0
 1004 continue
      do 1005 j = q1, qp
         if (abs(ab(j)) .ge. (lmda/ka)) then 
            w(j) = 0.d0
         else
            w(j) = 1.d0 - (abs(ab(j))*ka/lmda)
         endif
 1005 continue
      call adplasso(ab, y1, xz1, n, q, p, q1, qp, lmda, w,
     +     epsilon, maxit)
c     check the convergence
      count = count + 1
      call converge(tag, ab, abold, qp, epsilon)
      if (tag .eq. 1) go to 10000
      if (count .ge. maxit) call rexit("LLA diverges! \n")
      if (count .lt. maxit) go to 9000
10000 continue
      alpha(:)=ab(1:q)
      beta(:)=ab(q1:qp)
      end





c***********************************************************************
c     Main function using LLA1 approach                                    
c***********************************************************************
      subroutine fllaadp(olmdas,okas,ocoef,oaic,obic,oobj,odf,oevidx,
     +     y,x,z,n,q,p,nka,maxka,nlmda,minlmda,epsilon,maxit)
      integer n,q,p,nka,nlmda,maxit,odf(nka*nlmda),oevidx(nka*nlmda)
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),oaic(nka*nlmda),obic(nka*nlmda),
     +     oobj(nka*nlmda),y(n),x(n,q),z(n,p),
     +     maxka,minlmda,epsilon
      integer i,j
      double precision as(p),sz(n,p),kas(nka),unitka,unitlmda,
     +     lmdas(nlmda),alpha(q),beta(p),eta(n),
     +     inia(q,nlmda),inib(p,nlmda),etamat(n,nlmda)
c     standardization of Z
      call standard(as,sz,z,n,p)
c     calculate kappas
      if (nka .eq. 1) then
         kas(1)=0.d0
      else
         unitka=maxka/dble(nka-1)
         do 10000 i=1,nka
            kas(i)=dble(i-1)*unitka
10000    continue
      endif
c     calculate lambdas
      unitlmda=log(minlmda)/dble(nlmda-1)
      call maxbi(lmdas(1),alpha,beta,eta,y,x,sz,n,q,p,
     +     epsilon,maxit)
      do 10001 i=2,nlmda
         lmdas(i)=lmdas(1)*exp(unitlmda*dble(i-1))
10001 continue
c     Lasso solution path as initials
      inia(:,1)=alpha(:)
      inib(:,1)=beta(:)
      etamat(:,1)=eta(:)
      do 10002 i=2,nlmda
         call lasso(alpha,beta,eta,lmdas(i),y,x,sz,n,q,p,
     +        epsilon,maxit)
         inia(:,i)=alpha(:)
         inib(:,i)=beta(:)
         etamat(:,i)=eta(:)
10002 continue
c     MCP solution path along kappa
      if (nka .eq. 1) then
         olmdas(:)=lmdas(:)
         okas(1:nlmda)=kas(1)
         ocoef(1:q,:)=inia(:,:)
         ocoef((q+1):(q+p),:)=inib(:,:)
      else
         do 1003 i=1,nlmda
            olmdas(((i-1)*nka+1):(i*nka))=lmdas(i)
            okas(((i-1)*nka+1):(i*nka))=kas(:)
 1003    continue
!     ... i=1 case
         i=1
         alpha(:)=inia(:,i)
         beta(:) =inib(:,i)
         eta(:)  =etamat(:,i)
         do 1005 j=1,nka
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1005    continue
         do 1004 i=2,nlmda
            alpha(:)=inia(:,i)
            beta(:) =inib(:,i)
            eta(:)  =etamat(:,i)
c     ...... j=1 case
            j=1
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
            do 1006 j=2,nka
               call llaadp(alpha,beta,eta,lmdas(i),kas(j),y,x,sz,
     +              n,q,p,epsilon,maxit)
               ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
               ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1006       continue
 1004    continue
      endif
c     compute model size,eigenvalue,aic,bic, objective function
      do 10 i=1,(nka*nlmda)
         call fmeabo(odf(i),oevidx(i),oaic(i),obic(i),oobj(i),
     +        ocoef(:,i),y,x,sz,n,q,p,okas(i),olmdas(i))
 10   continue
!     change the penalized ocoefficients back
      do 10005 i=1,p
         ocoef(q+i,:)=as(i)*ocoef(q+i,:)
10005 continue
      end
