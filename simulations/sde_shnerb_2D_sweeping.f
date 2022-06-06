      Program sde_Shnerb_2D
      Implicit none
      Double precision a,b,D,sigma,dt,dx,beta,beta_aux,coeff1,coeff2,
     &coeff3,smallz,tmax,ro_sum,ro,ro_av,points,t,ro_c,xsumj,E,Elap,
     &E_l,E_r,E_c,E_0,exp_aux,lambda,xtmp,x,tlog,tlogaux,dran_g,gamma,
     &xpoisson,ro_l,ro_r,ro_0,x2,sigma2,ro_aux,ro_aux2,E_aux,E_aux2,ro_u
     &,ro_d,E_u,E_d,dran_u,w,w2,D2,kappa,R,coeff4,coeff5,coeff6,coeff7,
     &E_av,E_sum
      Integer N,L,itermax,Nmax,k,Kmax,i,j,iv,iter,m,nactivs,Nbin,iv2,up,
     &down,left,right,j1,j2
      PARAMETER (N=2**08,L=N*N,itermax=10**7,Nbin=10000,Kmax=2000000)
      DIMENSION ro(L),ro_av(Nbin),points(Nbin),E(L),E_av(Nbin),iv(0:L,4)
     &,iv2(0:N+1),ro_aux(0:L+1),ro_aux2(0:L+1),E_aux(0:L+1),
     &E_aux2(0:L+1)


      CALL dran_ini(95713)
      call sgrnd(3111)

      a=0.2d0
      b=0d0
      D=0.25d0
      w=1.d0
      D2=0.25d0
      w2=1.2d0
      R=0.40825d0
      kappa=1d0
      sigma=dsqrt(2.d0)

      dt=0.1d0
      dx=1d0

      beta=-a-4d0*D/(dx*dx)

      coeff1=2d0*D/(sigma*dx)**2
      coeff2=2d0/(sigma*sigma)
      coeff3=b*dt
      coeff4=D2*dt/(dx*dx)
      coeff5=kappa*dt
      coeff6=w2*dt
      coeff7=R*dt
      smallz=2d0/(sigma*sigma*dt)

      tmax=itermax*dt

      Nmax=int(dlog10(dble(tmax))*2d2)+1

      DO i=1,N

         iv2(i)=i

      END DO

      DO i=1,L

         ro_aux(i)=0d0
         ro_aux2(i)=0d0

         E_aux(i)=0d0
         E_aux2(i)=0d0

      END DO

      ro_aux(0)=0d0
      ro_aux2(0)=0d0

      E_aux(0)=0d0
      E_aux2(0)=0d0

      iv2(0)=N
      iv2(N+1)=1

      DO i=1,L

c
c     Ahora existe un sitio 0, vecino de todos los sitios de los bordes
c     (dos veces unido con las esquinas), cuyo potencial pondremos a
c     cero siempre.
c
         up=i-N
         down=i+N
         left=iv2(mod(i,N))-1
         right=iv2(mod(i,N))+1

c         IF (up.lt.0) up=0
         IF (up.lt.1) up=i+L-N
c         IF (down.gt.L) down=0
         IF (down.gt.L) down=i-L+N 
         IF (left.le.0) THEN
c            left=0
            left=i+N-1
         ELSE
            left=i-1
         END IF
         
         IF (right.ge.N+1) THEN
c            right=0
            right=i-N+1
         ELSE
            right=i+1
         END IF

         iv(i,1)=up
         iv(i,2)=right
         iv(i,3)=down
         iv(i,4)=left

      END DO

c      iv(0,1)=0
c      iv(0,2)=0
c      iv(0,3)=0
c      iv(0,4)=0


      Do m=1,Nbin

         ro_av(m)=0d0
         E_av(m)=0d0
         points(m)=0d0

      End Do

      DO k=1,Kmax
         
         write(*,*) 'Run',k

         ro_sum=0d0
         E_sum=0d0

         do i=1,L

            ro(i)=1d0
            E(i)=1d0
            ro_sum=ro_sum+ro(i)
            E_sum=E_sum+E(i)

         end do

         nactivs=L

         t=1d0

         DO iter=1,itermax

c            IF (mod(iter,100).eq.0) write(*,*) iter,ro_sum/L,E_sum/L

c            IF (ro_sum.eq.0) GO TO 10

            m=int(dlog10(t)*2d2)+1

c            IF (m.le.Nmax) THEN

            ro_av(m)=ro_av(m)+dble(ro_sum)/dble(L)
            E_av(m)=E_av(m)+dble(E_sum)/dble(L)

            points(m)=points(m)+1d0

c            END IF

            nactivs=0
            ro_sum=0d0
            E_sum=0d0
  
c
c     The upper neighbours for the first row is just the last row
c     sites. The lower neighbour for the last row is just the first row.
c
            DO i=1,N
               ro_aux(i)=ro(N*(N-1)+i)
               ro_aux2(i)=ro(i)

               E_aux(i)=E(N*(N-1)+i)
               E_aux2(i)=E(i)
            END DO
            
c
c     j1 is the number for the row and j2 for the column.
c
            DO j1=1,N-1

            ro_0=ro((j1-1)*N+1)
            ro_c=ro((j1-1)*N+N)
            ro_r=ro((j1-1)*N+1)

            E_0=E((j1-1)*N+1)
            E_c=E((j1-1)*N+N)
            E_r=E((j1-1)*N+1)
               
            DO j2=1,N-1

               j=(j1-1)*N+j2

               ro_l=ro_c
               ro_c=ro_r
               ro_u=ro_aux(j2)
               ro_d=ro(iv(j,3))
               ro_r=ro(iv(j,2))
               
               ro_aux(j2)=ro(j)

               E_l=E_c
               E_c=E_r
               E_u=E_aux(j2)
               E_d=E(iv(j,3))
               E_r=E(iv(j,2))
               
               E_aux(j2)=E(j)

               xsumj=ro_l+ro_r+ro_u+ro_d
               Elap=E_l+E_r+E_u+E_d-4d0*E_c

               IF ((xsumj+ro_c).ne.0) THEN

               beta_aux=beta+w*E_c
               exp_aux=dexp(-beta_aux*dt)


               if (exp_aux.lt.1.d-5) then
                  
                  lambda=smallz
                  
               else
                  
                  lambda=coeff2*exp_aux*beta_aux/(1.d0-exp_aux)
c                  lambda=coeff2*beta_aux/(1.d0-exp_aux)
                  
               endif
               
               xtmp=gamma(coeff1*xsumj+
     &              xpoisson(lambda*ro_c/exp_aux))/lambda
               
c               ro(j)=xtmp/(1.d0+coeff3*xtmp)
               ro(j)=xtmp              
               IF (ro(j).gt.0d0) nactivs=nactivs+1

               ro_sum=ro_sum+ro(j)

               END IF

               E(j)=E_c+coeff4*Elap+coeff7-coeff5*E_c-coeff6*E_c*ro_c
               IF (E(j).lt.0d0) E(j)=0d0
               
               E_sum=E_sum+E(j)

            END DO  


            j=(j1-1)*N+N

            ro_l=ro_c
            ro_c=ro_r
            ro_u=ro_aux(j2)
            ro_d=ro(iv(j,3))
            ro_r=ro_0
            
            ro_aux(j2)=ro(j)
            
            E_l=E_c
            E_c=E_r
            E_u=E_aux(j2)
            E_d=E(iv(j,3))
            E_r=E_0
            
            E_aux(j2)=E(j)

            xsumj=ro_l+ro_r+ro_u+ro_d
            Elap=E_l+E_r+E_u+E_d-4d0*E_c

            IF ((xsumj+ro_c).ne.0d0) THEN
            
            beta_aux=beta+w*E_c
            exp_aux=dexp(-beta_aux*dt)
            
            if (exp_aux.lt.1.d-5) then
               
               lambda=smallz
               
            else
               
               lambda=coeff2*exp_aux*beta_aux/(1.d0-exp_aux)
c               lambda=coeff2*beta_aux/(1.d0-exp_aux)
               
            endif
            
            xtmp=gamma(coeff1*xsumj+
     &           xpoisson(lambda*ro_c/exp_aux))/lambda
            
c            ro(j)=xtmp/(1.d0+coeff3*xtmp)
            ro(j)=xtmp
            IF (ro(j).gt.0d0) nactivs=nactivs+1

            ro_sum=ro_sum+ro(j)

            END IF

            E(j)=E_c+coeff4*Elap+coeff7-coeff5*E_c-coeff6*E_c*ro_c
            IF (E(j).lt.0d0) E(j)=0d0
            
            E_sum=E_sum+E(j)

            END DO

            ro_0=ro((N-1)*N+1)
            ro_c=ro((N-1)*N+N)
            ro_r=ro((N-1)*N+1)

            E_0=E((N-1)*N+1)
            E_c=E((N-1)*N+N)
            E_r=E((N-1)*N+1)
               
            DO j2=1,N-1

               j=(N-1)*N+j2
                  
               ro_l=ro_c
               ro_c=ro_r
               ro_u=ro_aux(j2)
               ro_d=ro_aux2(j2)
               ro_r=ro(iv(j,2))

               ro_aux(j2)=ro(j)

               E_l=E_c
               E_c=E_r
               E_u=E_aux(j2)
               E_d=E_aux2(j2)
               E_r=E(iv(j,2))

               E_aux(j2)=E(j)

               xsumj=ro_l+ro_r+ro_u+ro_d
               Elap=E_l+E_r+E_u+E_d-4d0*E_c

               IF ((xsumj+ro_c).ne.0) THEN

               beta_aux=beta+w*E_c
               exp_aux=dexp(-beta_aux*dt)


               if (exp_aux.lt.1.d-5) then
                  
                  lambda=smallz
                  
               else
                  
                  lambda=coeff2*exp_aux*beta_aux/(1.d0-exp_aux)
c                  lambda=coeff2*beta_aux/(1.d0-exp_aux)
                  
               endif
               
               xtmp=gamma(coeff1*xsumj+
     &              xpoisson(lambda*ro_c/exp_aux))/lambda
               
c               ro(j)=xtmp/(1.d0+coeff3*xtmp)
               ro(j)=xtmp             
               IF (ro(j).gt.0d0) nactivs=nactivs+1

               ro_sum=ro_sum+ro(j)

               END IF

               E(j)=E_c+coeff4*Elap+coeff7-coeff5*E_c-coeff6*E_c*ro_c
               IF (E(j).lt.0d0) E(j)=0d0

               E_sum=E_sum+E(j)

            END DO  

            j=(N-1)*N+N

            ro_l=ro_c
            ro_c=ro_r
            ro_u=ro_aux(j2)
            ro_d=ro_aux2(j2)
            ro_r=ro_0
            
            ro_aux(j2)=ro(j)

            E_l=E_c
            E_c=E_r
            E_u=E_aux(j2)
            E_d=E_aux2(j2)
            E_r=E_0
            
            E_aux(j2)=E(j)

            xsumj=ro_l+ro_r+ro_u+ro_d
            Elap=E_l+E_r+E_u+E_d-4d0*E_c

            IF ((xsumj+ro_c).ne.0d0) THEN
            
            beta_aux=beta+w*E_c
            exp_aux=dexp(-beta_aux*dt)
            
            if (exp_aux.lt.1.d-5) then
               
               lambda=smallz
               
            else
               
               lambda=coeff2*exp_aux*beta_aux/(1.d0-exp_aux)
c               lambda=coeff2*beta_aux/(1.d0-exp_aux)
               
            endif
            
            xtmp=gamma(coeff1*xsumj+
     &           xpoisson(lambda*ro_c/exp_aux))/lambda
            
c            ro(j)=xtmp/(1.d0+coeff3*xtmp)
            ro(j)=xtmp
            IF (ro(j).gt.0d0) nactivs=nactivs+1
            
            ro_sum=ro_sum+ro(j)

            END IF

            E(j)=E_c+coeff4*Elap+coeff7-coeff5*E_c-coeff6*E_c*ro_c
            IF (E(j).lt.0d0) E(j)=0d0

            E_sum=E_sum+E(j)

            t=t+dt

         END DO

 10      CONTINUE
         
         open(33,FILE='shnerb_2D_R0.40825_L08',STATUS ='Unknown')
         
         tlogaux=0.d0
         
         do m=1,Nmax
            
            tlog=int(10**((m-1)/(2d2)))
            
            IF ((points(m).ge.1).and.(tlog.ne.tlogaux)) THEN
               
               x=ro_av(m)/points(m)
               x2=E_av(m)/points(m)
               
               Write(33,15) tlog,x,x2,k
               
            END IF
            
            tlogaux=tlog
            
         end do
         
         close(33)
         
      END DO

 15   format(3(E20.10,1x),I7)

      end

*******************************************************************
*******************************************************************

      function gauss(iset)
c     Numerical Recipes Box-Muller's algorithm (polar method)
      implicit real*8 (a-h,p-z)
      parameter(zsmall=1.d0/4294967296.d0,
     &     zshift=(1.d0+zsmall)/2.d0)

      save gset

      if (iset.eq.0) then
12345    u1=2.d0*(igrnd()*zsmall+zshift)-1.d0
         u2=2.d0*(igrnd()*zsmall+zshift)-1.d0
         rsq=u1*u1+u2*u2
         if (rsq.ge.1.d0.or.rsq.eq.0.d0) goto 12345
         fac=dsqrt(-2.d0*dlog(rsq)/rsq)
         gset=u1*fac
         gauss=u2*fac
         iset=1
         return
      else
         gauss=gset
         iset=0
         return
      endif

      return
      end

************************************************************************
      function gamma(a)
c     random number from a gamma law of parameter a>0:
c     (Proba(x \le gamma_a \le x+dx)= dx x^(a-1)e^(-x)/Gamma(a))
c     Remark: the absorbing state is dealt with through the 
c     natural convention gamma_0=delta(x)  
      
      implicit real*8 (a-h,p-z)
      parameter(zsmall=1.d0/4294967296.d0,
     &     zshift=(1.d0+zsmall)/2.d0)

c      if (a.lt.0.d0) then
c         write(*,*)'PARAMETER TURNED <0 !: a=',a
c         STOP
c      endif
c      write(*,*)a
      if (a.gt.1.d0) then
c     Best's algo for a>1
c     NB: this is this region of parameter space where the algo spends
c     most of its time for usual applications, since a # x(i)/dt # 1/dt
         gamma=gamma_gt1_best(a)
         return
      else if (a.eq.0.d0) then
c     Dirac peak at 0 if a=0.
         gamma=0.d0
         return
      else if (a.eq.1.d0) then
c     simple exponential variate if a=1.
         gamma=-dlog(igrnd()*zsmall+zshift)
         return
      else
c     Gamma_a=Gamma_{a+1}*U^{1/a} for 0<a<1 by Stuart's theorem
         gamma=gamma_gt1_best(a+1.d0)*
     &        (igrnd()*zsmall+zshift)**(1.d0/a)
         return
      endif

      return
      end
************************************************************************
      function gamma_gt1_best(a)
c     Best's method: comparison from a Student-t density with 2 degrees of
c     freedom, of which the cumulative pdf P(Y\le y)= 0.5*(1-y/sqrt(2+y**2)),
c     that is Y=sqrt(2)*(U-0.5)/sqrt(U*(1-U)) where  U is 
c     a uniform random number on ]0,1[.
c     If h_a(y) is the  pdf of (a-1)+Y*sqrt(1.5*(a-0.25)), then
c     Gamma_a(y) \le c_a h_a(y), with a rejection factor 
c     c_a=(2*sqrt(3*(a-0.25))/Gamma(a))*((a-1)/e)^(a-1) (which is \le
c     3 for a->1+, and  stricly nonincreasing with a 
c     towards the limit sqrt(6/Pi)=1.38... for a-> \infty (in practice 50)
c     Reference: Devroye's book, chapter IX, pp. 407-411
      implicit real*8 (a-h,p-z)
      parameter(zsmall=1.d0/4294967296.d0,
     &     zshift=(1.d0+zsmall)/2.d0)
      
      b=a-1.d0
      csq=dsqrt(3.d0*(a-0.25d0))
      iok=0
      do while(iok.eq.0)
         u=(igrnd()*zsmall+zshift)
         v=(igrnd()*zsmall+zshift)
         w=u*(1.d0-u)
         y=csq*(u-0.5d0)/dsqrt(w)
         x=b+y
         if (x.ge.0.d0) then
            z=64.d0*(w*w*w)*v*v
            if (z.le.1.d0-2.d0*y*y/x) then
               iok=1
            else
               if (dlog(z).le.2.d0*(b*dlog(x/b)-y)) iok=1
            endif
         endif
      enddo
      gamma_gt1_best=x
      
      return
      end
************************************************************************
      function xpoisson(xm)
c     Algorithm Numerical Recipes 7.3
c     The returned value xpoisson is a REAL number (to avoid overflow pbs)
c     with zero fractional part such that xpoisson is a random variate drawn
c      from a Poisson distribution with mean xm: 
c      Proba(xpoisson=n*1.d0)=(xm^n/n!)e^{-xm}
c     Note that with xm=0, P(xpoisson=n*1.)=delta_{n,0}, and again the
c     absorbing state is dealt with this way.
      real*8 pi
      parameter(pi=3.141592653589793116d0)
      real*8 zsmall,zshift
      parameter(zsmall=1.d0/4294967296.d0,
     &     zshift=(1.d0+zsmall)/2.d0)

      real*8 xm,alxm,xtmp,g,sq,t,y,gammln,xpoisson

c      write(*,*) xm
      if (xm.lt.12.d0) then
c     direct calculation
          g=dexp(-xm)
          xtmp=-1.d0
          t=1.d0
          do while(t.ge.g)
             xtmp=xtmp+1.d0
             t=t*(igrnd()*zsmall+zshift)
          enddo
          xpoisson=xtmp
          return
      else
c     rejection method from a Lorentzian shifted and scaled
         sq=dsqrt(2.d0*xm)
         alxm=dlog(xm)
         g=xm*alxm-gammln(xm+1.d0)
 1000    y=dtan(pi*(igrnd()*zsmall+zshift))
         xtmp=sq*y+xm
         if (xtmp.lt.0.d0) goto 1000
         xtmp=aint(xtmp)
         t=0.9d0*(1+y*y)*dexp(xtmp*alxm-gammln(xtmp+1.d0)-g)
         if ((igrnd()*zsmall+zshift).gt.t) goto 1000
      endif
      xpoisson=xtmp
      
      return
      end
************************************************************************
      function  gammln(xx)
c     NR 6.1
c     evaluation of ln\Gamma(x+1) (required for the generation of Poisson
c     variates) by Lanczos' algorithm: error smaller than 2.d-10 for all xx>0
      real*8 gammln
      real*8 xx,ser,stp,tmp,x,y,cof(6)
      save cof,stp
      data cof,stp/76.18009172947146d0, -86.50532032941677d0,
     &     24.01409824083091d0, -1.2317395724501155d0,
     &     0.1208650973866179d-2, -0.5395239384953d-5,
     &     2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*dlog(tmp)-tmp
      ser=1.000000000190015d0
      do j=1,6
         y=y+1.d0
         ser=ser+cof(j)/y
      enddo
      gammln=tmp+dlog(stp*ser/x)

      return
      end
************************************************************************
      subroutine sgrnd(iseed)
*
      implicit integer(a-z)
*
* Period parameters
      parameter(N     =  624)
*
      dimension mt(0:N-1)
*                     the array for the state vector
      common /block/mti,mt
      save   /block/
*
*      setting initial seeds to mt[N] using
*      the generator Line 25 of Table 1 in
*      [KNUTH 1981, The Art of Computer Programming
*         Vol. 2 (2nd Ed.), pp102]
*
      mt(0)= iand(iseed,-1)
      do 1000 mti=1,N-1
        mt(mti) = iand(69069 * mt(mti-1),-1)
 1000 continue
*
      return
      end
************************************************************************
      function igrnd()

      implicit integer(a-z)
c     Matsumoto's Mersenne twister generator: outputs 
c     2**32 integers, -2**31 \le igrnd() \le 2*31-1
c     equidistributed in 623 dimensions, and of period # 2^19639-1
c     (as far as I recall)
*
* Period parameters
      parameter(N     =  624)
      parameter(N1    =  N+1)
      parameter(M     =  397)
      parameter(MATA  = -1727483681)
*                                    constant vector a
      parameter(UMASK = -2147483648)
*                                    most significant w-r bits
      parameter(LMASK =  2147483647)
*                                    least significant r bits
* Tempering parameters
      parameter(TMASKB= -1658038656)
      parameter(TMASKC= -272236544)
*
      dimension mt(0:N-1)
*                     the array for the state vector
      common /block/mti,mt
      save   /block/
      data   mti/N1/
*                     mti==N+1 means mt[N] is not initialized
*
      dimension mag01(0:1)
      data mag01/0, MATA/
      save mag01
*                        mag01(x) = x * MATA for x=0,1
*
      TSHFTU(y)=ishft(y,-11)
      TSHFTS(y)=ishft(y,7)
      TSHFTT(y)=ishft(y,15)
      TSHFTL(y)=ishft(y,-18)
*
      if(mti.ge.N) then
*                       generate N words at one time
        if(mti.eq.N+1) then
*                            if sgrnd() has not been called,
          call sgrnd(4357)
*                              a default initial seed is used
        endif
*
        do 1000 kk=0,N-M-1
            y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
            mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
 1000   continue
        do 1100 kk=N-M,N-2
            y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
            mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
 1100   continue
        y=ior(iand(mt(N-1),UMASK),iand(mt(0),LMASK))
        mt(N-1)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
        mti = 0
      endif
*
      y=mt(mti)
      mti=mti+1
      y=ieor(y,TSHFTU(y))
      y=ieor(y,iand(TSHFTS(y),TMASKB))
      y=ieor(y,iand(TSHFTT(y),TMASKC))

      igrnd=ieor(y,TSHFTL(y))

      return
      end

     
 
*******************************************************************
*******************************************************************

      subroutine dran_ini(iseed0)
      implicit double precision(a-h,o-z)
      parameter(ip=1279)
      parameter(np=14)
      parameter(nbit=31)
      parameter(m=2**np,np1=nbit-np,nn=2**np1-1,nn1=nn+1)
      integer ix(ip)
      dimension g(0:m)
      common /ixx/ ix
      common /icc/ ic
      common /gg/ g

      data c0,c1,c2/2.515517,0.802853,0.010328/
      data d1,d2,d3/1.432788,0.189269,0.001308/
c
      dseed=iseed0
      do i=1,ip
          ix(i)=0
          do j=0,nbit-1
              if(rand_xx(dseed).lt.0.5) ix(i)=ibset(ix(i),j)
          enddo
      enddo
200   continue
      ic=0
c
      pi=4.0d0*datan(1.0d0)
      do 1 i=m/2,m
      p=1.0-real(i+1)/(m+2)
      t=sqrt(-2.0*log(p))
      x=t-(c0+t*(c1+c2*t))/(1.0+t*(d1+t*(d2+t*d3)))
      g(i)=x
      g(m-i)=-x
1     continue

      u2th=1.0-real(m+2)/m*sqrt(2.0/pi)*g(m)*exp(-g(m)*g(m)/2)
      u2th=nn1*sqrt(u2th)
      do 856 i=0,m
856   g(i)=g(i)/u2th

      return
      end

      function dran_u()
      implicit double precision(a-h,o-z)
      parameter(ip=1279)
      parameter(iq=418)
      parameter(is=ip-iq)
      parameter (rmax=2147483647.0)
      integer ix(ip)
      common /ixx/ ix
      common /icc/ ic
      ic=ic+1
        if(ic.gt.ip) ic=1
      if(ic.gt.iq) then
          ix(ic)=ieor(ix(ic),ix(ic-iq))
      else
            ix(ic)=ieor(ix(ic),ix(ic+is))
        endif
        dran_u=real(ix(ic))/rmax
      return
      end

      function i_dran(n)
      implicit double precision(a-h,o-z)
      parameter(ip=1279)
      parameter(iq=418)
      parameter(is=ip-iq)
      integer ix(ip)
      common /ixx/ ix
      common /icc/ ic
      ic=ic+1
      if(ic.gt.ip) ic=1
      if(ic.gt.iq) then
       ix(ic)=ieor(ix(ic),ix(ic-iq))
      else
       ix(ic)=ieor(ix(ic),ix(ic+is))
      endif
      i_ran=ix(ic)
      if (n.gt.0) i_dran=mod(i_ran,n)+1
      return
      end

      function rand_xx(dseed)
      double precision a,c,xm,rm,dseed,rand_xx
      parameter (xm=2.d0**32,rm=1.d0/xm,a=69069.d0,c=1.d0)
      dseed=mod(dseed*a+c,xm)
      rand_xx=dseed*rm
      return
      end

      function dran_g()
      implicit double precision(a-h,o-z)
      parameter(ip=1279)
      parameter(iq=418)
      parameter(np=14)
      parameter(nbit=31)
      parameter(m=2**np,np1=nbit-np,nn=2**np1-1,nn1=nn+1)
      parameter(is=ip-iq)
      integer ix(ip)
      dimension g(0:m)
      common /ixx/ ix
      common /icc/ ic
      common /gg/ g
      
      ic=ic+1
      if(ic.gt.ip) ic=1
      if(ic.gt.iq) then
         ix(ic)=ieor(ix(ic),ix(ic-iq))
      else
         ix(ic)=ieor(ix(ic),ix(ic+is))
      endif
      i=ishft(ix(ic),-np1)
      i2=iand(ix(ic),nn)
      dran_g=i2*g(i+1)+(nn1-i2)*g(i)
      return
      end



            
