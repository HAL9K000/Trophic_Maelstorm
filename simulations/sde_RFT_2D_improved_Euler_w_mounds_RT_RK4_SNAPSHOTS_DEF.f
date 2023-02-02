      Program sde_RFT_2D_with_Kadanoff_boxes
      Implicit none
      Double precision a,b,c,D,sigma,dt,dx,beta,beta_aux,coeff1,coeff2,
     &coeff3,smallz,tmax,rho_sum,rho,rho_av,points,t,rho_c,xsumj,exp_aux
     &,y,lambda,xtmp,x,tlog,tlogaux,dran_g,gamma,xpoisson,rho_l,rho_r,
     &rho_0,x2,rho_aux,rho_aux2,rho_u,rho_d,dran_u,coeff4,tsup,cont_sup,
     &xaux,gamma2,b_v,coeff3_v,delta_x,delta_y,width_x,width_y,r_change,
     &slope,rad,intercept,FL,cte1L,cte2L,dy,xaux2,Lmin,factorL,
     &rho_Nk,rho_Nk_av,rho_Nk_sum,b_Nk,b_Nk_av,b_Nk_sum,x3,x4,x5,RK1,RK2
     &,RK3,RK4,rho_av_k,points_k
      CHARACTER*7 ca
      CHARACTER*5 ck
      CHARACTER*3 cNk
      Integer N,L,itermax,Nmax,k,Kmax,i,j,iv,iter,m,nactives,Nbin,iv2,up
     &,down,left,right,j1,j2,iruns,m2,moccupied,nmounds,jran,ji,jj,jiaux
     &,jjaux,i_dran,encounters,mound_loc,ii,nmounds_max,itermaxMAX,
     &Nk,L2,lk,ll,i_old,iaux,nboxes,kk,k0,iterst
      PARAMETER (lk=08,N=2**lk,L=N*N,itermaxMAX=2.5*10**5,Nbin=10000
     &,Kmax=5000,k0=0,iterst=0.1*itermaxMAX)
      DIMENSION rho(L),rho_av(Nbin),points(Nbin),iv(0:L,4),iv2(0:N+1),
     &rho_aux(0:L+1),rho_aux2(0:L+1),iruns(Nbin),b_v(L),coeff3_v(L),
     &FL(L),moccupied(L),mound_loc(L),rho_av_k(Nbin),points_k(Nbin),
     &rho_Nk(L),rho_Nk_av(0:lk,Nbin),b_Nk(L),b_Nk_av(0:lk,Nbin)


      CALL dran_ini(91537)
      call sgrnd(3111)

      a=-0.991d0
      b=2d0
      c=1d0
      D=1d0
      sigma=1d0

      dt=0.1d0
      dx=1d0
      dy=dx

      beta=a-4d0*D/(dx*dx)

      coeff1=2d0*D/(sigma*dx)**2
      coeff2=2d0/(sigma*sigma)
cx      coeff3=b*dt
      coeff4=c*dt
      smallz=2d0/(sigma*sigma*dt)

      beta_aux=beta
      exp_aux=dexp(-beta_aux*dt)

      if (exp_aux.lt.1.d-5) then
                  
         lambda=smallz
                  
      else
                  
         lambda=coeff2*exp_aux*beta_aux/(1.d0-exp_aux)
cBAD         lambda=coeff2*beta_aux/(1.d0-exp_aux)
                  
      endif

      tmax=itermaxMAX*dt

      Nmax=int(dlog10(dble(tmax))*2d2)+1

      DO i=1,N

         iv2(i)=i

      END DO

      DO i=1,L

         rho_aux(i)=0d0
         rho_aux2(i)=0d0

      END DO

      rho_aux(0)=0d0
      rho_aux2(0)=0d0

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

c
c     TEST for the 2D Kadanoff box construction:
c     
cx      DO jj=0,lk
cx
cx         Nk=N/(2**jj)
cx
cx         write(*,*) 'BOX SIZE',Nk
cx         write(cNk,fmt='(i3.3)') Nk
cx       OPEN(33,FILE='checking_Kadanoff_box_size_'//cNk,STATUS='Unknown')
cx
cx         DO ii=1,L
cx
c
c     WARNING: I HAD TO DO THIS STUPID TRICK FOR i_old BECAUSE OTHERWISE
c     ifort WOULD LINK i_old and iaux IN A WAY THAT i_old WOULD
c     AUTOMATICALLY UPDATE IF iaux CHANGED!! YOU CAN SEE IT HAPPEN BY
c     JUST CHANGING THESE TWO i_old LINES FOR THE USUAL i_old=iaux AT
c     THE BEGINNING OF THE LOOP
c
cx            i_old=iaux+1
cx            iaux=(ii-1)/N+1
cx            i_old=i_old-1
cx
cx            
cx            i=((ii-1)/N)/Nk
cx            ll=mod(ii,N)
cx            IF (ll.eq.0) ll=N
cx            j=(ll-1)/Nk+1
cx
cx            write(*,*)
cx            write(*,*) 'Pixel',ii,'location',iaux,iv2(mod(ii,N))
cx
cx            kk=i*(N/Nk)+j
cx            
cxc            write(*,*) mod(ii,N),ll,(ll-1)/Nk,(ll-1)/Nk+1
cx            write(*,*) 'Box',k,'locationB',i+1,j
cx
cx            IF (i_old.ne.iaux) write(33,*)
cxc            write(33,*) iaux,iv2(mod(ii,N)),kk
cx            IF (mod((i+1)+j,2).eq.0) THEN
cx               write(33,*) iaux,iv2(mod(ii,N)),-1
cxc               write(33,*) iaux,iv2(mod(ii,N)),kk
cx            ELSE
cx               write(33,*) iaux,iv2(mod(ii,N)),+1
cxc               write(33,*) iaux,iv2(mod(ii,N)),-kk
cx            END IF
cx            
cx         END DO
cx         CLOSE(33)
cx      END DO
cx      STOP


      Do m=1,Nbin

         rho_av(m)=0d0
         points(m)=0d0
         iruns(m)=0d0

         DO jj=0,lk
            rho_Nk_av(jj,m)=0d0
            b_Nk_av(jj,m)=0d0
         END DO

      End Do

c
c     MOUND-RELATED PARAMETERS:
c
      nmounds_max=dx*dx*7/((2.4d-1)*2)**2 !We calculate the new number
                                          !of mounds to ensure a similar
                                          ! density of mounds as in
                                          ! Kenya
      nmounds_max=0

c      write(*,*) 'We need ',nmounds_max,' mounds!!'

      Lmin=1d0
      factorL=1.5d0

      delta_x=5d1               !distance between mounds, in m
      delta_y=delta_x
      
      width_x=5d0               !Radius of the mound, in m
      width_y=5d0

      r_change=width_x/2d0      !Radius from mound top at which the
                                !steep descent to zero starts for the L
                                ! and I functions.

c
c     Calculation of that steep descent function (we assume it is
c     linear):
c
      xaux=dexp(-((r_change/width_x)**2)*0.5d0) !WARNING: This assumes
                                                !width_x=width_y

      xaux=1d0
      slope=(0d0-xaux)/(width_x-r_change)
      intercept=0d0-slope*width_x


      write(ca,fmt='(f7.5)') abs(a)
c
c     Replicates:
c


      DO k=k0+1,Kmax
         
cx         write(*,*) 'Run',k,nmounds_max
         write(ck,fmt='(i5.5)') k

c
c     MOUND DEFINITION:
c

c
c     Definition of the improvement function, FL(x,y). Mounds are
c     separated a (REAL) distance delta_x and delta_y, or in terms of
c     the discretized system, a distance delta_x/dx and delta_y/dy.
c

      DO i=1,L

         moccupied(i)=0
         FL(i)=Lmin

      END DO

c
c     First we calculate the location of the mounds:
c
      nmounds=0
      DO WHILE (nmounds.lt.nmounds_max)
    
         jran=i_dran(L)

         ji=(jran-1)/N+1
         jj=iv2(mod(jran,N))       

c         write(*,*)
c         write(*,*) 'Trying at ',jran,ji,jj

         encounters=0
         DO i=1,L
     
c
c     In order to calculate distances between candidate mound and i, we
c     move everything to a system of reference in the center of the
c     system:
c     
            jiaux=(i-1)/N+1-ji
            jjaux=iv2(mod(i,N))-jj

            IF (jiaux.gt.N/2) THEN
               jiaux=jiaux-N
            ELSE
               IF (jiaux.lt.-N/2)  jiaux=jiaux+N
            END IF

            IF (jjaux.gt.N/2) THEN
               jjaux=jjaux-N
            ELSE
               IF (jjaux.lt.-N/2)  jjaux=jjaux+N
            END IF

            x=dx*jiaux
            y=dy*jjaux

            rad=dsqrt(x**2+y**2)

c
c     If any location within its radius had already been occupied, we
c     discard that j as a candidate. Otherwise, we mark all its covered
c     space as occupied.
c
            IF (rad.le.width_x) THEN

               IF (moccupied(i).eq.0) THEN
                                             
                  moccupied(i)=1
                
               ELSE

                  encounters=encounters+1
c                  write(*,*) i,' was occupied!!'
                  GO TO 22

               END IF

            END IF

         END DO

 22      CONTINUE

         IF (encounters.eq.0) THEN

            FL(i)=Lmin*(cte1L*xaux+cte2L)

            nmounds=nmounds+1
            mound_loc(nmounds)=jran

         END IF

      END DO

c      write(*,*) nmounds

      DO ii=1,nmounds

         j=mound_loc(ii)

         ji=(j-1)/N+1
         jj=iv2(mod(j,N))

c
c     Because we want our gaussian function below (which is into [0,1])
c     be into [1,factor], we need the change constants:
c

         xaux2=(factorL-1d0)*dran_u()+1d0 !number into factorL and 1d0

         cte1L=(xaux2-1d0)/(1d0-0d0)
         cte2L=xaux2-cte1L*1d0
         
         DO i=1,L

         jiaux=(j-1)/N+1
         jjaux=iv2(mod(j,N))

         ji=(i-1)/N+1-jiaux
         jj=iv2(mod(i,N))-jjaux

c
c     Periodic boundary conditions:
c

         IF (ji.gt.N/2) THEN
            ji=ji-N
         ELSE
            IF (ji.lt.-N/2)  ji=ji+N
         END IF
         
         IF (jj.gt.N/2) THEN
            jj=jj-N
         ELSE
            IF (jj.lt.-N/2)  jj=jj+N
         END IF

         x=dx*ji
         y=dy*jj

         rad=dsqrt(x**2+y**2)
         xaux=0d0
       
         IF (rad.le.width_x) THEN
            IF (rad.ge.r_change) THEN
cxcx               xaux=xaux*(-rad*0.4d0+1d0)
               xaux=slope*rad+intercept

            ELSE
               xaux=1d0
            END IF
            
            FL(i)=Lmin*(cte1L*xaux+cte2L)
         END IF

         END DO

      END DO

      IF (nmounds.eq.0) THEN

         DO i=1,L

cx            FL(i)=dran_u()
            FL(i)=Lmin

         END DO

      END IF


         rho_sum=0d0
         do i=1,L

            rho(i)=1d0

            rho_sum=rho_sum+rho(i)

            rho_Nk(i)=0d0       !This is just initialization because the
                                ! first thing we do at each time step is
                                ! update this array for each box
            b_Nk(i)=0d0

            b_v(i)=-b*FL(i)     !NOTE THE NEGATIVE SIGN FOR b!!!
            coeff3_v(i)=b_v(i)*dt

         end do
         nactives=L

cx      open(33,FILE='eluding_Lang_2D_improved_a-'//ca//
cx     &'_HEXA_L_function_matrix_'//ck,STATUS='Unknown')
cx      DO i=1,N
cx
cx         write(33,19) (b_v((i-1)*N+j),j=1,N)
cx                  
cx      END DO
cx      CLOSE(33)
cx      STOP

      Do m=1,Nbin
         rho_av_k(m)=0d0
         points_k(m)=0d0        
      End do

c
c     Time loop:
c

cx         itermax=iterst+i_dran(itermaxMAX-iterst)
         itermax=itermaxMAX
         t=1d0

         DO iter=1,itermax

c            IF (mod(iter,1000).eq.0) write(*,*) iter,rho_sum/L

            m=int(dlog10(t)*2d2)+1

            IF (rho_sum.eq.0d0) THEN
c
c     This is to calculate the survival time of a run:
c
               tsup=tsup+t
               cont_sup=cont_sup+1d0 
               GO TO 10

            END IF

c            IF (m.le.Nmax) THEN

c
c     Calculating the average rho within each box:
c
            DO jj=0,lk

               Nk=N/(2**jj)     !Box linear size, in pixels.
c
c     Because the linear size of the box is given by Nk, the number of
c     boxes is (N/Nk)**2, i.e. 2**(2*jj):
c
               nboxes=(N/Nk)*(N/Nk)

               DO ii=1,L
                  
                  i=((ii-1)/N)/Nk
                  ll=mod(ii,N)
                  IF (ll.eq.0) ll=N
                  j=(ll-1)/Nk+1
                  
                  kk=i*(N/Nk)+j

                  rho_Nk(kk)=rho_Nk(kk)+rho(ii)
                  b_Nk(kk)=b_Nk(kk)+b_v(ii)
                  
               END DO

c
c     If things work properly, boxes nboxes+1 to L should never be
c     used:
c
               rho_Nk_sum=0d0
               b_Nk_sum=0d0
               DO ii=1,nboxes
                  rho_Nk(ii)=rho_Nk(ii)/(Nk*Nk)
                  b_Nk(ii)=b_Nk(ii)/(Nk*Nk)

                  rho_Nk_sum=rho_Nk_sum+rho_Nk(ii)
                  b_Nk_sum=b_Nk_sum+b_Nk(ii)

                  rho_Nk(ii)=0d0
                  b_Nk(ii)=0d0
               END DO

               rho_Nk_av(jj,m)=rho_Nk_av(jj,m)+rho_Nk_sum/dble(nboxes)
               b_Nk_av(jj,m)=b_Nk_av(jj,m)+b_Nk_sum/dble(nboxes)
            END DO
            
            rho_av_k(m)=rho_av_k(m)+dble(rho_sum)/dble(L)
            points_k(m)=points_k(m)+1d0

            rho_av(m)=rho_av(m)+dble(rho_sum)/dble(L)
            points(m)=points(m)+1d0

c            END IF

            nactives=0
            rho_sum=0d0
  
c
c     The upper neighbours for the first row is just the last row
c     sites. The lower neighbour for the last row is just the first row.
c
            DO i=1,N
               rho_aux(i)=rho(N*(N-1)+i)
               rho_aux2(i)=rho(i)
            END DO
            
c
c     j1 is the number for the row and j2 for the column.
c
            DO j1=1,N-1

            rho_0=rho((j1-1)*N+1)
            rho_c=rho((j1-1)*N+N)
            rho_r=rho((j1-1)*N+1)
               
            DO j2=1,N-1

               j=(j1-1)*N+j2

               rho_l=rho_c
               rho_c=rho_r
               rho_u=rho_aux(j2)
               rho_d=rho(iv(j,3))
               rho_r=rho(iv(j,2))
               
               rho_aux(j2)=rho(j)

               xsumj=rho_l+rho_r+rho_u+rho_d

               IF ((xsumj+rho_c).ne.0) THEN

              
               xtmp=gamma2(coeff1*xsumj+
     &              xpoisson(lambda*rho_c/exp_aux))/lambda

               RK1=-(coeff3_v(j)+coeff4*xtmp)*xtmp*xtmp
               xaux=xtmp+0.5d0*RK1
               RK2=-(coeff3_v(j)+coeff4*xaux)*xaux*xaux
               xaux=xtmp+0.5d0*RK2
               RK3=-(coeff3_v(j)+coeff4*xaux)*xaux*xaux
               xaux=xtmp+RK3
               RK4=-(coeff3_v(j)+coeff4*xaux)*xaux*xaux
               rho(j)=xtmp+(RK1+2d0*(RK2+RK3)+RK4)/6d0

               IF (rho(j).gt.0d0) nactives=nactives+1

               rho_sum=rho_sum+rho(j)

               END IF

            END DO  


            j=(j1-1)*N+N

            rho_l=rho_c
            rho_c=rho_r
            rho_u=rho_aux(j2)
            rho_d=rho(iv(j,3))
            rho_r=rho_0
            
            rho_aux(j2)=rho(j)

            xsumj=rho_l+rho_r+rho_u+rho_d

            IF ((xsumj+rho_c).ne.0d0) THEN
            
            xtmp=gamma2(coeff1*xsumj+
     &           xpoisson(lambda*rho_c/exp_aux))/lambda
            
            RK1=-(coeff3_v(j)+coeff4*xtmp)*xtmp*xtmp
            xaux=xtmp+0.5d0*RK1
            RK2=-(coeff3_v(j)+coeff4*xaux)*xaux*xaux
            xaux=xtmp+0.5d0*RK2
            RK3=-(coeff3_v(j)+coeff4*xaux)*xaux*xaux
            xaux=xtmp+RK3
            RK4=-(coeff3_v(j)+coeff4*xaux)*xaux*xaux
            rho(j)=xtmp+(RK1+2d0*(RK2+RK3)+RK4)/6d0

            IF (rho(j).gt.0d0) nactives=nactives+1

            rho_sum=rho_sum+rho(j)

            END IF

            END DO

            rho_0=rho((N-1)*N+1)
            rho_c=rho((N-1)*N+N)
            rho_r=rho((N-1)*N+1)
               
            DO j2=1,N-1

               j=(N-1)*N+j2
                  
               rho_l=rho_c
               rho_c=rho_r
               rho_u=rho_aux(j2)
               rho_d=rho_aux2(j2)
               rho_r=rho(iv(j,2))

               rho_aux(j2)=rho(j)

               xsumj=rho_l+rho_r+rho_u+rho_d

               IF ((xsumj+rho_c).ne.0) THEN

               xtmp=gamma2(coeff1*xsumj+
     &              xpoisson(lambda*rho_c/exp_aux))/lambda

               RK1=-(coeff3_v(j)+coeff4*xtmp)*xtmp*xtmp
               xaux=xtmp+0.5d0*RK1
               RK2=-(coeff3_v(j)+coeff4*xaux)*xaux*xaux
               xaux=xtmp+0.5d0*RK2
               RK3=-(coeff3_v(j)+coeff4*xaux)*xaux*xaux
               xaux=xtmp+RK3
               RK4=-(coeff3_v(j)+coeff4*xaux)*xaux*xaux
               rho(j)=xtmp+(RK1+2d0*(RK2+RK3)+RK4)/6d0

               IF (rho(j).gt.0d0) nactives=nactives+1

               rho_sum=rho_sum+rho(j)

               END IF

            END DO  

            j=(N-1)*N+N

            rho_l=rho_c
            rho_c=rho_r
            rho_u=rho_aux(j2)
            rho_d=rho_aux2(j2)
            rho_r=rho_0
            
            rho_aux(j2)=rho(j)

            xsumj=rho_l+rho_r+rho_u+rho_d

            IF ((xsumj+rho_c).ne.0d0) THEN
            
            xtmp=gamma2(coeff1*xsumj+
     &           xpoisson(lambda*rho_c/exp_aux))/lambda

            RK1=-(coeff3_v(j)+coeff4*xtmp)*xtmp*xtmp
            xaux=xtmp+0.5d0*RK1
            RK2=-(coeff3_v(j)+coeff4*xaux)*xaux*xaux
            xaux=xtmp+0.5d0*RK2
            RK3=-(coeff3_v(j)+coeff4*xaux)*xaux*xaux
            xaux=xtmp+RK3
            RK4=-(coeff3_v(j)+coeff4*xaux)*xaux*xaux
            rho(j)=xtmp+(RK1+2d0*(RK2+RK3)+RK4)/6d0

            IF (rho(j).gt.0d0) nactives=nactives+1

            rho_sum=rho_sum+rho(j)

            END IF

            t=t+dt

         END DO
        
 10      CONTINUE
         
         Do m2=1,m
            iruns(m2)=iruns(m2)+1d0
         End do
         
         IF (rho_sum.ne.0d0) THEN
         open(32,FILE='eluding_Lang_2D_a-0.99100_L08_final_snap_'
     &   //ck,STATUS='Unknown')         
         do j=1,L
            
            IF (ji.ne.(j-1)/N+1) write(32,*)
            ji=(j-1)/N+1
            jj=iv2(mod(j,N))
            
            write(32,17) ji*dx,jj*dy,rho(j),b_v(j),itermax*dt
            
         end do
         close(32)
         END IF
         
         open(33,FILE='eluding_Lang_2D_improved_a-0.99100_L08_'//ck
     &   ,STATUS='Unknown')
        
         tlogaux=0.d0
         
         do m=1,Nmax
            
            tlog=int(10**((m-1)/(2d2)))
            
            IF ((points_k(m).ge.1).and.(tlog.ne.tlogaux)) THEN
               
               x2=rho_av_k(m)/points_k(m)
               
               Write(33,18) tlog,x2,points_k(m),iruns(m)+0d0,k-k0
               
            END IF
            
            tlogaux=tlog
            
         end do
         close(33)
         
         open(33,FILE='eluding_Lang_2D_improved_a-0.99100_L08',
     &        STATUS='Unknown')
        
         tlogaux=0.d0
         
         do m=1,Nmax
            
            tlog=int(10**((m-1)/(2d2)))
            
            IF ((points(m).ge.1).and.(tlog.ne.tlogaux)) THEN
               
               x2=rho_av(m)/points(m)
               x3=x2*(iruns(m)+0d0)/(k-k0)
               
               Write(33,15) tlog,x2,x3,tsup/cont_sup,points(m),iruns(m)
     &              +0d0,k-k0
               
            END IF
            
            tlogaux=tlog
            
         end do
         close(33)

         DO jj=0,lk

            Nk=N/(2**jj)

            write(cNk,fmt='(i3.3)') Nk
         OPEN(34,FILE='eluding_Lang_2D_improved_a-0.99100_L08_BOX_size_'
     &           //cNk,STATUS='Unknown')
            
            tlogaux=0.d0
         
            do m=1,Nmax
            
               tlog=int(10**((m-1)/(2d2)))
            
               IF ((points(m).ge.1).and.(tlog.ne.tlogaux)) THEN

                  x2=rho_Nk_av(jj,m)/points(m)
                  x3=x2*(iruns(m)+0d0)/(k-k0)
                  x4=b_Nk_av(jj,m)/points(m)
                  x5=x4*(iruns(m)+0d0)/(k-k0)

                  write(34,16) tlog,x2,x3,x4,x5,tsup/cont_sup,points(m),
     &              iruns(m)+0d0,k-k0

               END IF
            
               tlogaux=tlog
            end do
            CLOSE(34)

         END DO
         
      END DO

 15   format(6(E20.10,1x),I7)
 16   format(8(E20.10,1x),I7)
 17   format(5(E20.10,1x))
 18   format(4(E20.10,1x),I7) 
 19   format(<N>(E20.10,1x)) 

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
      function gamma2(a)
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
         gamma2=gamma_gt1_best(a)
         return
      else if (a.eq.0.d0) then
c     Dirac peak at 0 if a=0.
         gamma2=0.d0
         return
      else if (a.eq.1.d0) then
c     simple exponential variate if a=1.
         gamma2=-dlog(igrnd()*zsmall+zshift)
         return
      else
c     Gamma_a=Gamma_{a+1}*U^{1/a} for 0<a<1 by Stuart's theorem
         gamma2=gamma_gt1_best(a+1.d0)*
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



            
