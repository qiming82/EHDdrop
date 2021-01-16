      program eaxdrop
c------------------------------------------
c  viscous drop with insoluble surfactants
c------------------------------------------

c      implicit real*8 (a-h,o-y)
      implicit double precision (a-h,o-z)
	parameter(neq=1024)

	dimension z(neq),r(neq),z0(neq),r0(neq),s(neq),f(neq)
	dimension zold(neq),rold(neq),sold(neq),dsold(neq)
	dimension znew(neq),rnew(neq),snew(neq),v(neq)
	dimension dz(neq),ddz(neq),dr(neq),ddr(neq),ds(neq),indx(neq)
	dimension u1(neq),u2(neq),un(neq),fold(neq)
      dimension an1(neq),an2(neq),ut(neq)
      dimension anx(neq),anr(neq)
      dimension q1(neq),q2(neq),q1p(neq),q2p(neq),bb(neq),bbp(neq)
      dimension az(neq),bz(neq),cz(neq),ar(neq),br(neq),cr(neq)
      dimension d1(neq),b1(neq),c1(neq),d2(neq),b2(neq),c2(neq)
      dimension ag(neq),bg(neq),cg(neq)
c      dimension az(513),bz(513),cz(513),ar(513),br(513),cr(513)
c      dimension sumxx(neq,neq), sumxr(neq,neq), sumrx(neq,neq)
      dimension ml(neq,neq), mu(neq,neq),Am(neq,neq),hh(neq)
      dimension bbz(neq),bbr(neq),bbzp(neq),bbrp(neq),w_trp(neq)
	dimension Sxx(neq,neq),Sxr(neq,neq),Srx(neq,neq),Srr(neq,neq)
      dimension u_sxxtmp1(neq),u_sxrtmp1(neq),u_srxtmp1(neq)
	dimension u_srrtmp1(neq)
      dimension u_sxxtmp2(neq),u_sxrtmp2(neq),u_srxtmp2(neq)
	dimension u_srrtmp2(neq)
      Dimension gma(neq),sfct(neq),dsfct(neq),sftmp(neq),en(neq)
      dimension v0(neq),dgm(neq),dgm2(neq),gtmp2(neq)
      dimension pen(neq),pet(neq),pe(neq)
	integer ipiv(neq)


      pi=4*datan(1.d0)
	eps = 0.0000000001D0
	      
      open (22,file='axialv.dat')
      open (23,file='xyval.dat')
      open (24,file='fet.dat')
      open (25,file='fen.dat')
      open (26,file='fes.dat')
      open (27,file='tangv.dat')
      open (7,file='info.txt')
      open (11,file='hmin.dat')
      open (12,file='wmax.dat')
      open (10,file='unmax.dat')

c number of segments of the interface
      N=128
	Nm=2*N+2
      h=1.D0/N

      ram=1.0000D0
      beta=-(ram-1.0D0)/(ram+1.0D0)
	beta=beta/(4.0D0*pi)

c      eb=0.1800d0  ! electric bond number
c idef=1  desingularise
	idef=1

	Ca=0.0d0
	C=Ca*pi*4.0d0
c	coef1=-0.50d0/(2.0*pi*(ram+1.0d0))
      coef1=-0.5d0   ! after rescale the velocity

      aa1=1.6      !1.8   # 100
	bb1=1.2      !1.2   # 100
	cc1=0.2      !0.5   # 100
c initial settings...
	do j=1,N+1
	   s(j)=(j-1)*h
	   z0(j)=-dcos(pi*s(j))  !!
         z(j)=z0(j)
	   r0(j)=dsin(pi*s(j))
	   r(j)=r0(j)
	   write(23,*) z0(j), r0(j)
	end do    

c        cvol0=4.D0/3.D0*pi
        ssum0=1.0d0
         cvol0=0.0D0
c--------------- get the volume----------------------------------------
          do j=1,N
       	   hh(j)=z0(j+1)-z0(j)
	       cvol0=cvol0+hh(j)*(r0(j)**2+r0(j+1)**2)
	     end do
	     cvol0=dabs(cvol0)*pi
c----------------------------------------------------------------------
c       if(dabs(cvol0-4.0d0*pi/3.0d0).gt.0.00001d0)then
c          vscale=(cvol0/4.0d0/pi*3.0d0)**(1.D0/2.D0)
c           do j=1,N+1
c 	       r0(j)=r0(j)/vscale
c 	    end do	   
c	 endif
c	do j=1,N+1
c	   z(j)=z0(N+1-(j-1))
c	   r(j)=r0(N+1-(j-1))
c	   write(23,*) z(j), r(j)
c	end do
c-----------------------------------------------------------------------
	a0 = r(1)

      sum11=0.0D0
	sum12=0.0D0
	sum21=0.0D0
	sum22=0.0D0
      y11 = 0.0D0
	y12 = 0.0D0
	y21 = 0.0D0
	y22 = 0.0D0
c      nplt=1

c      itmax=100
c	reit=10
c      tol=0.000000001D0
      time=0.0D0
	dt=0.0005D0
c      dt=0.1d0*(1.0D0/N)**(1.5D0)
	dt2=dt/2.0d0
	dt3=dt2/2.0d0
	!
	tf=10.D0
c 	nplt=tf/dt/60
c	nplt2=nplt/10
c	tf=50.5D0
 	nplt=tf/dt/40
c      nplt=20
	nplt2=nplt/10
c-- temporary tf
c       tf=0.01d0

      jj=1
	jjj=1
	kkk=1

      eb = 0.32d0
      QQ = 1.370d0 
      RR = 0.100d0

c start evolution
      do while (time.lt.tf)
	   time = time + dt
c--------------------------------
c---- get derivatives first
c--
          call splc_clm
     +  (N
     +  ,s,z
     +  ,0.0D0
     +  ,0.0d0
     +  ,az,bz,cz
     +  )
          call splc_clm
     +  (N
     +  ,s,r
     +  ,pi
     +  ,-pi
     +  ,ar,br,cr
     +  )

	    cr(N+1)=-pi
	    cz(N+1)=0.0d0

	    do j=1,N+1
	       dr(j)=cr(j)
	       ddr(j)=2.0D0*br(j)
	       dz(j)=cz(j)
	       ddz(j)=2.0D0*bz(j)
	       
	       ds(j)=dsqrt(dz(j)**2+dr(j)**2)
               an1(j)=-dr(j)/ds(j)
	       an2(j)=dz(j)/ds(j)
	       anx(j)=-dr(j)
	       anr(j)=dz(j)

	       call df(r(j),dz(j),ddz(j),dr(j),ddr(j),ds(j),f(j))
	    end do
c--------
c instantaneous arclength 
c--------	    
       arcl=0.0d0
	   do k=1,N
	    arcl=arcl+(s(k+1)-s(k))/2.0d0*(r(k)*ds(k)+r(k+1)*ds(k+1))
	   enddo

c--- calculate the electric force
      if(eb.ne.0.0D0)then
c           call getef
c     +
c     +   (z,r
c     +   ,az,bz,cz
c     +   ,ar,br,cr
c     +   ,N,s
c     +   ,pen
c     +   )    ! for perfect conductor
           call getef2
     +
     +   (z,r,s,an1
     +   ,az,bz,cz
     +   ,ar,br,cr
     +   ,N
     +   ,RR
     +   ,pe,pet,pen
     +   )   ! for leaky dielectrics
      else
	  do kk=1,N+1
            pet(kk)=0.0d0
	    pen(kk)=0.0d0
	  end do
      end if

	    do j=1,N+1
               fen2 = pen(j)**2
               fet2 = (pet(j)/ds(j))**2
c-- perfect conductor
c               f(j) = f(j) - 0.5d0*eb*(-an1(j)+pen(j))**2
c	       q1(j)= coef1*f(j)*anx(j) 
c	       q2(j)= coef1*f(j)*anr(j) 

c-- leaky dielectrics
	       f(j) = f(j) - 0.5D0*eb*((1.0-QQ*RR*RR)*fen2-(1.0-QQ)*fet2)
	       q1(j)= coef1*(f(j)*anx(j) + 
     +                    -eb*(1.0-RR*QQ)*pen(j)*pet(j)*anr(j))
	       q2(j)= coef1*(f(j)*anr(j) +
     +                    -eb*(1.0-RR*QQ)*pen(j)*pet(j)*(-anx(j))) 
cc-- question: is the minus sign really there? seems so ! check
	    enddo

c---------------------------------------------------  
c prepare the integration of the kernels...
c--------------------------------------------------

c      step 1, single layers

	    do i=1,N+1
	        sum11=0.0D0
			sum12=0.0D0
			sum21=0.0D0
		    sum22=0.0D0
	       do k=1,N
	         call gsint
     +
     +   (z(k),r(k)
     +   ,s(k),s(k+1)
     +   ,q1(k),q1(k+1),q2(k),q2(k+1)
     +   ,s(i),z(i),r(i)
     +   ,az(k),bz(k),cz(k)
     +   ,ar(k),br(k),cr(k)
     +   ,y11,y12
     +   ,y21,y22
     +   )
	       sum11 = sum11 + y11
	       sum12 = sum12 + y12
	       sum21 = sum21 + y21
	       sum22 = sum22 + y22
	      end do
c  u1 = ux,  u2 = ur
	      bbz(i) = C*z(i) + sum11 + sum12
	      bbr(i) = -C/2.0*r(i) + sum21 + sum22
	    end do
          
	if(dabs(ram-1.0D0).lt.eps)then
         do k=1,N+1
	     u1(k)=bbz(k)
	     u2(k)=bbr(k)
	  end do
	else
c bb is the right hand side of the system
		do k=1,N+1
            bb(k)=bbz(k)
	      bb(k+N+1)=bbr(k)  
		end do
c       step 2, double layer potentials...

         do i=1,N+1
	   
	       do k=1,N
	if (idef.eq.0)then
          call strint
     +
     +   (z(k),r(k)
     +   ,s(k),s(k+1)
     +   ,anx(k),anx(k+1),anr(k),anr(k+1)
     +   ,z(i),r(i)
     +   ,1
     +   ,az(k),bz(k),cz(k)
     +   ,ar(k),br(k),cr(k)
     +   ,u_sxxtmp1(k),u_sxrtmp1(k)
     +   ,u_srxtmp1(k),u_srrtmp1(k)
     +   )

          call strint
     +
     +   (z(k),r(k)
     +   ,s(k),s(k+1)
     +   ,anx(k),anx(k+1),anr(k),anr(k+1)
     +   ,z(i),r(i)
     +   ,2
     +   ,az(k),bz(k),cz(k)
     +   ,ar(k),br(k),cr(k)
     +   ,u_sxxtmp2(k),u_sxrtmp2(k)
     +   ,u_srxtmp2(k),u_srrtmp2(k)
     +   )
	else
          call strint2
     +
     +   (z(k),r(k)
     +   ,s(k),s(k+1)
     +   ,anx(k),anx(k+1),anr(k),anr(k+1)
     +   ,z(i),r(i),an1(i),an2(i)
     +   ,1,arcL
     +   ,az(k),bz(k),cz(k)
     +   ,ar(k),br(k),cr(k)
     +   ,u_sxxtmp1(k),u_sxrtmp1(k)
     +   ,u_srxtmp1(k),u_srrtmp1(k)
     +   )

          call strint2
     +
     +   (z(k),r(k)
     +   ,s(k),s(k+1)
     +   ,anx(k),anx(k+1),anr(k),anr(k+1)
     +   ,z(i),r(i),an1(i),an2(i)
     +   ,2,arcL
     +   ,az(k),bz(k),cz(k)
     +   ,ar(k),br(k),cr(k)
     +   ,u_sxxtmp2(k),u_sxrtmp2(k)
     +   ,u_srxtmp2(k),u_srrtmp2(k)
     +   )
      endif
	   if (k.eq.1) then
	      Sxx(i,k)=u_sxxtmp1(k)
	      Sxr(i,k)=u_sxrtmp1(k)
	      Srx(i,k)=u_srxtmp1(k)
	      Srr(i,k)=u_srrtmp1(k)
	   else 
	      Sxx(i,k)=u_sxxtmp1(k)+u_sxxtmp2(k-1)
	      Sxr(i,k)=u_sxrtmp1(k)+u_sxrtmp2(k-1)
	      Srx(i,k)=u_srxtmp1(k)+u_srxtmp2(k-1)
	      Srr(i,k)=u_srrtmp1(k)+u_srrtmp2(k-1)
	   end if
	   end do

	    Sxx(i,N+1)=u_sxxtmp2(N)
	    Sxr(i,N+1)=u_sxrtmp2(N)
	    Srx(i,N+1)=u_srxtmp2(N)
	    Srr(i,N+1)=u_srrtmp2(N)
	   
	  end do

        Sxx=-beta*Sxx
	  Sxr=-beta*Sxr
	  Srx=-beta*Srx
	  Srr=-beta*Srr

	  do i=1,N+1
	    Sxx(i,i)=1.0D0+Sxx(i,i)
	    Srr(i,i)=1.0D0+Srr(i,i)
	  end do

c  Am=[Sxx -Sxr; -Srx Srr]
        do i=1,N+1
	    do k=1,N+1
	      Am(i,k)=Sxx(i,k)
	      Am(i,k+N+1)=Sxr(i,k)
	      Am(i+N+1,k)=Srx(i,k)
	      Am(i+N+1,k+N+1)=Srr(i,k)
	    end do
	  end do

        call dgesv(Nm, 1, Am, neq, ipiv, bb, neq, info)
	   v(1:Nm)=bb(1:Nm)       
	   u1(1:N+1)=v(1:N+1)
	   u2(1:N+1)=v(N+2:2*N+2)
	end if

c-----------
c  update marker points
c-----------
	    do k=1,N+1
	       un(k)=u1(k)*an1(k)+u2(k)*an2(k)
	       ut(k)=u1(k)*an2(k)-u2(k)*an1(k)
	       znew(k)=z(k)+dt*un(k)*an1(k)
	       rnew(k)=r(k)+dt*un(k)*an2(k)
c	       znew(k)=z(k)+dt*(un(k)*an1(k)+ut(k)*an2(k))
c   	       rnew(k)=r(k)+dt*(un(k)*an2(k)-ut(k)*an1(k))
	     end do

c+++++++++++++++++++++++++++++++++++++++++++++
           call regridf(znew,rnew,z,r,s,N,time,dt)
c------------------
         cvol=0.0D0
c       get the volume of the jet of one period,
          do j=1,N
       	   hh(j)=z(j+1)-z(j)
	       cvol=cvol+hh(j)*(r(j)**2+r(j+1)**2)
	     end do
	     cvol=cvol*pi
c	   write (*,*) 'ratio=', cvol/cvol0

c scale the jet
          vscale=(cvol/cvol0)**(1.D0/2.D0)
          do j=1,N+1
cc	       z(j)=z(j)/vscale
	       r(j)=r(j)/vscale
	    end do
c++++++++++++++++++++++++++++++++++++++++++++ END IF
c         rmin = minval(r(N/4:7*N/8))
         rmin = minval(r(5:N-4))
         vnmax = maxval(abs(un))
         v1max = maxval(abs(u1))

	   if (mod(jj,nplt2).eq.0) then
	      write (10,*) time, vnmax
	      write (11,*) time, rmin
	      write (12,*) time, v1max
         endif
	   if (mod(jj,nplt).eq.0) then
	      jjj=jjj+1
	     if(mod(jjj,2).eq.0)then
	        do j=N+1,1,-1
		      write (22,*) z(j), u1(j)
	          write(23,*) z(j), r(j)
	          write (24,*) z(j), pet(j)
	          write (25,*) z(j), pen(j)
	          write (26,*) s(j), pe(j)
	          write (27,*) s(j), ut(j)
	        end do
	     else
	        do j=1,N+1
		      write (22,*) z(j), u1(j)
	          write(23,*) z(j), r(j)
	          write (24,*) z(j), pet(j)
	          write (25,*) z(j), pen(j)
	          write (26,*) s(j), pe(j)
	          write (27,*) s(j), ut(j)
	        end do
	     end if
          end if	      
c------

c	   write(*,*) 'rmin=', rmin
c	   write(*,*) 'time=', time
c         write(*,*) 'unmax=', vnmax
c         if(rmin.lt.0.01.and.nplt.gt.40)then
c	       nplt=12
c	   end if
	   if(rmin.lt.0.02d0.and.dt.gt.dt2)then
	    dt=dt2
	   endif
	   if(rmin.lt.0.01d0.and.dt.gt.dt3)then
	    dt=dt3
	   endif
	   if(rmin.lt.0.005d0.and.dt.gt.dt3/2.0)then
	    dt=dt3/2.0
	   endif

         if (vnmax.gt.20.0D0) then
	       write(7,*) "velocity too big"
	       write(7,*) 'vmax=', vnmax
	       write(7,*) 'N=', N
	       exit
	    end if
         if (rmin.lt.0.001) then
	       write(7,*) "neck too small"
	       write(7,*) 'rmin=', rmin
	       exit
	    end if
         if (dabs(vnmax).lt.0.001D0) then
	       write(7,*) "velocity too small"
             write(7,*) 'vmax_min=', vnmax
	       exit
	    end if
c          if(dabs(vnmax).lt.0.0004d0)then
c	       Cs=Cs+0.005d0
c	       C=Cs*(1.0-E)*4.0*pi   
c             Ca=Ca+0.005d0
c	       C = Ca*4.0*pi
c		  write(*,*) 'Ca=', Ca   
c	    endif

         jj = jj + 1
c	   write(*,*) 'jj=', jj
	end do

      write(7,*) "stop time is"
	write(7,*) time
	RL=z(N+1)
	B=r(N/2+1)
	D=(RL-B)/(RL+B)
	write(7,*) 'ratio D=', D
      write(7,*) 'Ca=', Ca
      write(7,*) 'dt=', dt
      write(7,*) 'N=', N
      write(7,*) 'Eb=', eb, 'R=', RR, 'Q=', QQ
	
!
c output data
       open (21,file='initialc.dat')
       open (37,file='finalxy.dat')
c       open (38,file='finalgma.dat')
	 open (41,file='finalw.dat')

	      jjj=jjj+1
	     if(mod(jjj,2).eq.0)then
	        do j=N+1,1,-1
		      write (22,*) z(j), u1(j)
	          write(23,*) z(j), r(j)
	          write (24,*) z(j), pet(j)
	          write (25,*) z(j), pen(j)
	          write (26,*) s(j), pe(j)
	          write (27,*) s(j), ut(j)
	        end do
	      else
	        do j=1,N+1
		      write (22,*) z(j), u1(j)
	          write(23,*) z(j), r(j)
	          write (24,*) z(j), pet(j)
	          write (25,*) z(j), pen(j)
	          write (26,*) s(j), pe(j)
	          write (27,*) s(j), ut(j)
	        end do
	     end if

	 do j=1,N+1
	   write (21,*) z0(j), r0(j)
	   write (37,*) z(j), r(j)
	   write (41,*) z(j), u1(j)
	 end do

c close 23
      end program
