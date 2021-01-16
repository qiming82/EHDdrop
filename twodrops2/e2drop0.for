      program vdrop_axisym
c------------------------------------------

c      implicit real*8 (a-h,o-y)
      implicit double precision (a-h,o-z)
	parameter(neq=1024)

	dimension z1(neq),r1(neq),z10(neq),r10(neq),s1(neq),f1(neq)
	dimension z2(neq),r2(neq),z20(neq),r20(neq),s2(neq),f2(neq)
	dimension znew1(neq),rnew1(neq),snew1(neq),v1(neq)
	dimension znew2(neq),rnew2(neq),snew2(neq),v2(neq)
	dimension dz1(neq),ddz1(neq),dr1(neq),ddr1(neq),ds1(neq),indx(neq)
	dimension dz2(neq),ddz2(neq),dr2(neq),ddr2(neq),ds2(neq)
	dimension uz1(neq),uz2(neq),un1(neq),ut1(neq)
	dimension ur1(neq),ur2(neq),un2(neq),ut2(neq)
      dimension an11(neq),an12(neq),an21(neq),an22(neq)
      dimension anx1(neq),anr1(neq),bb1(neq),bb2(neq),bbr1(neq)
      dimension anx2(neq),anr2(neq),bbz1(neq)
      dimension qx1(neq),qr1(neq),qx2(neq),qr2(neq),bb(neq),bbp(neq)
      dimension az1(neq),bz1(neq),cz1(neq),ar1(neq),br1(neq),cr1(neq)
      dimension az2(neq),bz2(neq),cz2(neq),ar2(neq),br2(neq),cr2(neq)
      dimension d1(neq),b1(neq),c1(neq),d2(neq),b2(neq),c2(neq)
c      dimension az(513),bz(513),cz(513),ar(513),br(513),cr(513)
c      dimension sumxx(neq,neq), sumxr(neq,neq), sumrx(neq,neq)
      dimension ml(neq,neq), mu(neq,neq),hh(neq)
      dimension bbz(neq),bbr(neq),bbz2(neq),bbr2(neq),w_trp(neq)
	dimension Sxx1(neq,neq),Sxr1(neq,neq),Srx1(neq,neq),Srr1(neq,neq)
	dimension Sxx2(neq,neq),Sxr2(neq,neq),Srx2(neq,neq),Srr2(neq,neq)
	dimension Rxx1(neq,neq),Rxr1(neq,neq),Rrx1(neq,neq),Rrr1(neq,neq)
	dimension Rxx2(neq,neq),Rxr2(neq,neq),Rrx2(neq,neq),Rrr2(neq,neq)
      dimension u_sxxtmp1(neq),u_sxrtmp1(neq),u_srxtmp1(neq)
	dimension u_srrtmp1(neq)
      dimension u_sxxtmp2(neq),u_sxrtmp2(neq),u_srxtmp2(neq)
	dimension u_srrtmp2(neq)
      dimension v_sxxtmp1(neq),v_sxrtmp1(neq),v_srxtmp1(neq)
	dimension v_srrtmp1(neq)
      dimension v_sxxtmp2(neq),v_sxrtmp2(neq),v_srxtmp2(neq)
	dimension v_srrtmp2(neq)
	dimension ipiv(neq)
      dimension Am(neq,neq),Am11(neq,neq),Am12(neq,neq)
	dimension Am21(neq,neq),Am22(neq,neq)


      pi=4*datan(1.d0)
	pi4=pi*4.0d0
	eps = 0.0000000001D0
	      
      open (9,file='info.txt')
      open (18,file='initl1.dat')
      open (19,file='initl2.dat')
      open (20,file='xyval1.dat')
      open (21,file='xyval2.dat')
      open (22,file='axialv1.dat')
c      open (23,file='radialv1.dat')
      open (24,file='axialv2.dat')
c      open (25,file='radialv2.dat')
      open (11,file='hmin1.dat')
      open (12,file='wmax1.dat')
      open (13,file='hmin2.dat')
      open (14,file='wmax2.dat')

c number of segments of the interface
      N=32
	Nm=2*N+2
	Nm2=2*Nm
      h=1.D0/N

      a1=2.0d0  ! radius of drop 1
	a2=1.0d0  ! radius of drop 2
      ram1=1.000D0
	ram2=1.000d0
      beta1=-(ram1-1.0D0)/(ram1+1.0D0)
	beta1=beta1/(4.0D0*pi)
      beta2=-(ram2-1.0D0)/(ram1+1.0D0)
	beta2=beta2/(4.0D0*pi)
      beta3=-(ram1-1.0D0)/(ram2+1.0D0)
	beta3=beta3/(4.0D0*pi)
      beta4=-(ram2-1.0D0)/(ram2+1.0D0)
	beta4=beta4/(4.0D0*pi)

	Boi=0.10d0  ! Bond number
	idef1=0 ! idef=1 to deflate
	idef2=0
c	coef1=-1.0/(ram1+1.0d0)
c	coef2=-1.0/(ram2+1.0d0)
	coef1=-0.50d0/(2.0*pi*(ram1+1.0d0))
	coef2=-0.50d0/(2.0*pi*(ram2+1.0d0))
c      coef1=-0.5d0   ! after rescale the velocity
      ep=0.00d0
      dgap=0.50*a1
c initial settings...
	do j=1,N+1
	   s1(j)=(j-1)*h
	   s2(j)=s1(j)
	   P2  = 0.5d0*(3.0*dcos(pi*s1(j))*dcos(pi*s1(j))-1.0)
         z10(j)=dgap/2.0+a1-a1*dcos(pi*s1(j))*(1.0+ep*P2)  !!
         z1(j)=z10(j)
	   r10(j)=a1*dsin(pi*s1(j))*(1.0+ep*P2)
	   r1(j)=r10(j)
         z20(j)=-dgap/2.0-a2-a2*dcos(pi*s2(j))*(1.0-ep*P2)  !!
         z2(j)=z20(j)
	   r20(j)=a2*dsin(pi*s2(j))*(1.0-ep*P2)
	   r2(j)=r20(j)
	   write(18,*) z1(j), r10(j)
	   write(19,*) z2(j), r20(j)
	   write(20,*) z1(j), r10(j)
	   write(21,*) z2(j), r20(j)
	end do    

c        cvol0=4.D0/3.D0*pi
         cvol10=0.0D0
         cvol20=0.0D0
c--------------- get the volume----------------------------------------
          do j=1,N
       	   hh(j)=z10(j+1)-z10(j)
	       cvol10=cvol10+hh(j)*(r10(j)**2+r10(j+1)**2)
	     end do
	     cvol10=dabs(cvol10)*pi
          do j=1,N
       	   hh(j)=z20(j+1)-z20(j)
	       cvol20=cvol20+hh(j)*(r20(j)**2+r20(j+1)**2)
	     end do
	     cvol20=dabs(cvol20)*pi
c----------------------------------------------------------------------

      sum11=0.0D0
	sum12=0.0D0
	sum21=0.0D0
	sum22=0.0D0
      y11 = 0.0D0
	y12 = 0.0D0
	y21 = 0.0D0
	y22 = 0.0D0
      ssum11=0.0D0
	ssum12=0.0D0
	ssum21=0.0D0
	ssum22=0.0D0
      s11 = 0.0D0
	s12 = 0.0D0
	s21 = 0.0D0
	s22 = 0.0D0

      time=0.0D0
c	dt=0.0005D0
      dt=0.2d0*(1.0D0/N)**(1.5D0)
	dt2=dt/2.0d0
	dt3=dt2/2.0d0
	!
	tf=5.D0
c 	nplt=tf/dt/60
c	nplt2=nplt/10
c	tf=50.5D0
 	nplt=tf/dt/50
	nplt2=nplt/10

      jj=1
	jjj=1
	kkk=1

c start evolution
      do while (time.lt.tf)
	   time = time + dt
c--------------------------------
c---- get derivatives first
c--
          call splc_clm
     +  (N
     +  ,s1,z1
     +  ,0.0D0
     +  ,0.0d0
     +  ,az1,bz1,cz1
     +  )
          call splc_clm
     +  (N
     +  ,s1,r1
     +  ,a1*pi
     +  ,-a1*pi
     +  ,ar1,br1,cr1
     +  )

	    cr1(N+1)=-a1*pi
	    cz1(N+1)=0.0d0

          call splc_clm
     +  (N
     +  ,s2,z2
     +  ,0.0D0
     +  ,0.0d0
     +  ,az2,bz2,cz2
     +  )
          call splc_clm
     +  (N
     +  ,s2,r2
     +  ,a2*pi
     +  ,-a2*pi
     +  ,ar2,br2,cr2
     +  )

	    cr2(N+1)=-a2*pi
	    cz2(N+1)=0.0d0

	    do j=1,N+1
	       dr1(j)=cr1(j)
	       ddr1(j)=2.0D0*br1(j)
	       dz1(j)=cz1(j)
	       ddz1(j)=2.0D0*bz1(j)
		   ds1(j)=dsqrt(dz1(j)**2+dr1(j)**2)
             an11(j)=-dr1(j)/ds1(j)
	       an12(j)=dz1(j)/ds1(j)
	       anx1(j)=-dr1(j)
	       anr1(j)=dz1(j)
	       call df(r1(j),dz1(j),ddz1(j),dr1(j),ddr1(j),ds1(j),f1(j))
		   qx1(j)= (z1(j)+Boi*f1(j))*anx1(j)
	       qr1(j)= (z1(j)+Boi*f1(j))*anr1(j)

	       dr2(j)=cr2(j)
	       ddr2(j)=2.0D0*br2(j)
	       dz2(j)=cz2(j)
	       ddz2(j)=2.0D0*bz2(j)
		   ds2(j)=dsqrt(dz2(j)**2+dr2(j)**2)
             an21(j)=-dr2(j)/ds2(j)
	       an22(j)=dz2(j)/ds2(j)
	       anx2(j)=-dr2(j)
	       anr2(j)=dz2(j)
	       call df(r2(j),dz2(j),ddz2(j),dr2(j),ddr2(j),ds2(j),f2(j))
		   qx2(j)= (z2(j)+(a1/a2)**2*Boi*f2(j))*anx2(j)
	       qr2(j)= (z2(j)+(a1/a2)**2*Boi*f2(j))*anr2(j)
	    end do
c--------
c instantaneous arclength 
c--------	    
       arcl1=0.0d0
	   do k=1,N
	    arcl1=arcl1+
     +		(s1(k+1)-s1(k))/2.0d0*(r1(k)*ds1(k)+r1(k+1)*ds1(k+1))
	   enddo
       arcl2=0.0d0
	   do k=1,N
	    arcl2=arcl2+
     +		(s2(k+1)-s2(k))/2.0d0*(r2(k)*ds2(k)+r2(k+1)*ds2(k+1))
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
	        ssum11=0.0D0
			ssum12=0.0D0
			ssum21=0.0D0
		    ssum22=0.0D0
	       do k=1,N
c --- x0 on S1, over S1
	         call gsint
     +
     +   (z1(k),r1(k)
     +   ,s1(k),s1(k+1)
     +   ,qx1(k),qx1(k+1),qr1(k),qr1(k+1)
     +   ,s1(i),z1(i),r1(i)
     +   ,az1(k),bz1(k),cz1(k)
     +   ,ar1(k),br1(k),cr1(k)
     +   ,y11,y12
     +   ,y21,y22
     +   )
	       sum11 = sum11 + y11
	       sum12 = sum12 + y12
	       sum21 = sum21 + y21
	       sum22 = sum22 + y22
c --- x0 on S1, over S2 (regular integral..)
	         call gsint
     +
     +   (z2(k),r2(k)
     +   ,s2(k),s2(k+1)
     +   ,qx2(k),qx2(k+1),qr2(k),qr2(k+1)
     +   ,-10.0,z1(i),r1(i)
     +   ,az2(k),bz2(k),cz2(k)
     +   ,ar2(k),br2(k),cr2(k)
     +   ,y11,y12
     +   ,y21,y22
     +   )
	       sum11 = sum11 + y11
	       sum12 = sum12 + y12
	       sum21 = sum21 + y21
	       sum22 = sum22 + y22
c --- x0 on S2, over S1 (regular integral..)
	         call gsint
     +
     +   (z1(k),r1(k)
     +   ,s1(k),s1(k+1)
     +   ,qx1(k),qx1(k+1),qr1(k),qr1(k+1)
     +   ,-10.0,z2(i),r2(i)
     +   ,az1(k),bz1(k),cz1(k)
     +   ,ar1(k),br1(k),cr1(k)
     +   ,s11,s12
     +   ,s21,s22
     +   )
	       ssum11 = ssum11 + s11
	       ssum12 = ssum12 + s12
	       ssum21 = ssum21 + s21
	       ssum22 = ssum22 + s22
c --- x0 on S2, over S2
	         call gsint
     +
     +   (z2(k),r2(k)
     +   ,s2(k),s2(k+1)
     +   ,qx2(k),qx2(k+1),qr2(k),qr2(k+1)
     +   ,s2(i),z2(i),r2(i)
     +   ,az2(k),bz2(k),cz2(k)
     +   ,ar2(k),br2(k),cr2(k)
     +   ,s11,s12
     +   ,s21,s22
     +   )
	       ssum11 = ssum11 + s11
	       ssum12 = ssum12 + s12
	       ssum21 = ssum21 + s21
	       ssum22 = ssum22 + s22
	      end do
c  u1 = ux,  u2 = ur
	      bbz1(i) = coef1*(sum11 + sum12)
	      bbr1(i) = coef1*(sum21 + sum22)

	      bbz2(i) = coef2*(ssum11 + ssum12)
	      bbr2(i) = coef2*(ssum21 + ssum22)
	    end do
          
	if(dabs(ram1-1.0D0).lt.eps.and.
     +	dabs(ram2-1.0).lt.eps)then
         do k=1,N+1
	     uz1(k)=bbz1(k)
	     ur1(k)=bbr1(k)
	     uz2(k)=bbz2(k)
	     ur2(k)=bbr2(k)
	  end do
	else
c bb is the right hand side of the system
		do k=1,N+1
            bb1(k)=bbz1(k)
	      bb1(k+N+1)=bbr1(k)  
            bb2(k)=bbz2(k)
	      bb2(k+N+1)=bbr2(k)  
		end do
          do k=1,Nm
	      bb(k)=bb1(k)
	      bb(k+Nm)=bb2(k)
	    enddo

c       step 2, double layer potentials...

         do 9 i=1,N+1
	   
	       do k=1,N
	if (idef1.eq.0)then
c-- x0 on S1, over S1
          call strint
     +
     +   (z1(k),r1(k)
     +   ,s1(k),s1(k+1)
     +   ,anx1(k),anx1(k+1),anr1(k),anr1(k+1)
     +   ,z1(i),r1(i)
     +   ,1
     +   ,az1(k),bz1(k),cz1(k)
     +   ,ar1(k),br1(k),cr1(k)
     +   ,u_sxxtmp1(k),u_sxrtmp1(k)
     +   ,u_srxtmp1(k),u_srrtmp1(k)
     +   )

          call strint
     +
     +   (z1(k),r1(k)
     +   ,s1(k),s1(k+1)
     +   ,anx1(k),anx1(k+1),anr1(k),anr1(k+1)
     +   ,z1(i),r1(i)
     +   ,2
     +   ,az1(k),bz1(k),cz1(k)
     +   ,ar1(k),br1(k),cr1(k)
     +   ,u_sxxtmp2(k),u_sxrtmp2(k)
     +   ,u_srxtmp2(k),u_srrtmp2(k)
     +   )
c-- x0 on S2, over S1
          call strint
     +
     +   (z1(k),r1(k)
     +   ,s1(k),s1(k+1)
     +   ,anx1(k),anx1(k+1),anr1(k),anr1(k+1)
     +   ,z2(i),r2(i)
     +   ,1
     +   ,az1(k),bz1(k),cz1(k)
     +   ,ar1(k),br1(k),cr1(k)
     +   ,v_sxxtmp1(k),v_sxrtmp1(k)
     +   ,v_srxtmp1(k),v_srrtmp1(k)
     +   )

          call strint
     +
     +   (z1(k),r1(k)
     +   ,s1(k),s1(k+1)
     +   ,anx1(k),anx1(k+1),anr1(k),anr1(k+1)
     +   ,z2(i),r2(i)
     +   ,2
     +   ,az1(k),bz1(k),cz1(k)
     +   ,ar1(k),br1(k),cr1(k)
     +   ,v_sxxtmp2(k),v_sxrtmp2(k)
     +   ,v_srxtmp2(k),v_srrtmp2(k)
     +   )
	else
c-- x0 on S1, over S1
          call strint2
     +
     +   (z1(k),r1(k)
     +   ,s1(k),s1(k+1)
     +   ,anx1(k),anx1(k+1),anr1(k),anr1(k+1)
     +   ,z1(i),r1(i),an11(i),an12(i)
     +   ,1,arcL1
     +   ,az1(k),bz1(k),cz1(k)
     +   ,ar1(k),br1(k),cr1(k)
     +   ,u_sxxtmp1(k),u_sxrtmp1(k)
     +   ,u_srxtmp1(k),u_srrtmp1(k)
     +   )

          call strint2
     +
     +   (z1(k),r1(k)
     +   ,s1(k),s1(k+1)
     +   ,anx1(k),anx1(k+1),anr1(k),anr1(k+1)
     +   ,z1(i),r1(i),an11(i),an12(i)
     +   ,2,arcL1
     +   ,az1(k),bz1(k),cz1(k)
     +   ,ar1(k),br1(k),cr1(k)
     +   ,u_sxxtmp2(k),u_sxrtmp2(k)
     +   ,u_srxtmp2(k),u_srrtmp2(k)
     +   )
c-- x0 on S2, over S1
          call strint2
     +
     +   (z1(k),r1(k)
     +   ,s1(k),s1(k+1)
     +   ,anx1(k),anx1(k+1),anr1(k),anr1(k+1)
     +   ,z2(i),r2(i),an21(i),an22(i)
     +   ,1,arcL1
     +   ,az1(k),bz1(k),cz1(k)
     +   ,ar1(k),br1(k),cr1(k)
     +   ,u_sxxtmp1(k),u_sxrtmp1(k)
     +   ,u_srxtmp1(k),u_srrtmp1(k)
     +   )

          call strint2
     +
     +   (z1(k),r1(k)
     +   ,s1(k),s1(k+1)
     +   ,anx1(k),anx1(k+1),anr1(k),anr1(k+1)
     +   ,z2(i),r2(i),an21(i),an22(i)
     +   ,2,arcL1
     +   ,az1(k),bz1(k),cz1(k)
     +   ,ar1(k),br1(k),cr1(k)
     +   ,u_sxxtmp2(k),u_sxrtmp2(k)
     +   ,u_srxtmp2(k),u_srrtmp2(k)
     +   )
      endif
	   if (k.eq.1) then
	      Sxx1(i,k)=u_sxxtmp1(k)
	      Sxr1(i,k)=u_sxrtmp1(k)
	      Srx1(i,k)=u_srxtmp1(k)
	      Srr1(i,k)=u_srrtmp1(k)

	      Rxx1(i,k)=v_sxxtmp1(k)
	      Rxr1(i,k)=v_sxrtmp1(k)
	      Rrx1(i,k)=v_srxtmp1(k)
	      Rrr1(i,k)=v_srrtmp1(k)
	   else 
	      Sxx1(i,k)=u_sxxtmp1(k)+u_sxxtmp2(k-1)
	      Sxr1(i,k)=u_sxrtmp1(k)+u_sxrtmp2(k-1)
	      Srx1(i,k)=u_srxtmp1(k)+u_srxtmp2(k-1)
	      Srr1(i,k)=u_srrtmp1(k)+u_srrtmp2(k-1)

	      Rxx1(i,k)=v_sxxtmp1(k)+v_sxxtmp2(k-1)
	      Rxr1(i,k)=v_sxrtmp1(k)+v_sxrtmp2(k-1)
	      Rrx1(i,k)=v_srxtmp1(k)+v_srxtmp2(k-1)
	      Rrr1(i,k)=v_srrtmp1(k)+v_srrtmp2(k-1)
	   end if
	   end do

	    Sxx1(i,N+1)=u_sxxtmp2(N)
	    Sxr1(i,N+1)=u_sxrtmp2(N)
	    Srx1(i,N+1)=u_srxtmp2(N)
	    Srr1(i,N+1)=u_srrtmp2(N)
	   
	    Rxx1(i,N+1)=v_sxxtmp2(N)
	    Rxr1(i,N+1)=v_sxrtmp2(N)
	    Rrx1(i,N+1)=v_srxtmp2(N)
	    Rrr1(i,N+1)=v_srrtmp2(N)
9      continue

        Sxx1=-beta1*Sxx1
	  Sxr1=-beta1*Sxr1
	  Srx1=-beta1*Srx1
	  Srr1=-beta1*Srr1

        Rxx1=-beta3*Rxx1
	  Rxr1=-beta3*Rxr1
	  Rrx1=-beta3*Rrx1
	  Rrr1=-beta3*Rrr1

	  do i=1,N+1
	    Sxx1(i,i)=1.0D0+Sxx1(i,i)
	    Srr1(i,i)=1.0D0+Srr1(i,i)
	  end do
c******************++++++++++++++++++------------------00
         do 99 i=1,N+1
	   
	       do k=1,N
	if (idef2.eq.0)then
c-- x0 on S1, over S2
          call strint
     +
     +   (z2(k),r2(k)
     +   ,s2(k),s2(k+1)
     +   ,anx2(k),anx2(k+1),anr2(k),anr2(k+1)
     +   ,z1(i),r1(i)
     +   ,1
     +   ,az2(k),bz2(k),cz2(k)
     +   ,ar2(k),br2(k),cr2(k)
     +   ,u_sxxtmp1(k),u_sxrtmp1(k)
     +   ,u_srxtmp1(k),u_srrtmp1(k)
     +   )

          call strint
     +
     +   (z2(k),r2(k)
     +   ,s2(k),s2(k+1)
     +   ,anx2(k),anx2(k+1),anr2(k),anr2(k+1)
     +   ,z1(i),r1(i)
     +   ,2
     +   ,az2(k),bz2(k),cz2(k)
     +   ,ar2(k),br2(k),cr2(k)
     +   ,u_sxxtmp2(k),u_sxrtmp2(k)
     +   ,u_srxtmp2(k),u_srrtmp2(k)
     +   )
c-- x0 on S2, over S2
          call strint
     +
     +   (z2(k),r2(k)
     +   ,s2(k),s2(k+1)
     +   ,anx2(k),anx2(k+1),anr2(k),anr2(k+1)
     +   ,z2(i),r2(i)
     +   ,1
     +   ,az2(k),bz2(k),cz2(k)
     +   ,ar2(k),br2(k),cr2(k)
     +   ,v_sxxtmp1(k),v_sxrtmp1(k)
     +   ,v_srxtmp1(k),v_srrtmp1(k)
     +   )

          call strint
     +
     +   (z2(k),r2(k)
     +   ,s2(k),s2(k+1)
     +   ,anx2(k),anx2(k+1),anr2(k),anr2(k+1)
     +   ,z2(i),r2(i)
     +   ,2
     +   ,az2(k),bz2(k),cz2(k)
     +   ,ar2(k),br2(k),cr2(k)
     +   ,v_sxxtmp2(k),v_sxrtmp2(k)
     +   ,v_srxtmp2(k),v_srrtmp2(k)
     +   )
	else
c-- x0 on S1, over S2
          call strint2
     +
     +   (z2(k),r2(k)
     +   ,s2(k),s2(k+1)
     +   ,anx2(k),anx2(k+1),anr2(k),anr2(k+1)
     +   ,z1(i),r1(i),an11(i),an12(i)
     +   ,1,arcL2
     +   ,az2(k),bz2(k),cz2(k)
     +   ,ar2(k),br2(k),cr2(k)
     +   ,u_sxxtmp1(k),u_sxrtmp1(k)
     +   ,u_srxtmp1(k),u_srrtmp1(k)
     +   )

          call strint2
     +
     +   (z2(k),r2(k)
     +   ,s2(k),s2(k+1)
     +   ,anx2(k),anx2(k+1),anr2(k),anr2(k+1)
     +   ,z1(i),r1(i),an11(i),an12(i)
     +   ,2,arcL2
     +   ,az2(k),bz2(k),cz2(k)
     +   ,ar2(k),br2(k),cr2(k)
     +   ,u_sxxtmp2(k),u_sxrtmp2(k)
     +   ,u_srxtmp2(k),u_srrtmp2(k)
     +   )
c-- x0 on S2, over S2
          call strint2
     +
     +   (z2(k),r2(k)
     +   ,s2(k),s2(k+1)
     +   ,anx2(k),anx2(k+1),anr2(k),anr2(k+1)
     +   ,z2(i),r2(i),an21(i),an22(i)
     +   ,1,arcL2
     +   ,az2(k),bz2(k),cz2(k)
     +   ,ar2(k),br2(k),cr2(k)
     +   ,u_sxxtmp1(k),u_sxrtmp1(k)
     +   ,u_srxtmp1(k),u_srrtmp1(k)
     +   )

          call strint2
     +
     +   (z2(k),r2(k)
     +   ,s2(k),s2(k+1)
     +   ,anx2(k),anx2(k+1),anr2(k),anr2(k+1)
     +   ,z2(i),r2(i),an21(i),an22(i)
     +   ,2,arcL2
     +   ,az2(k),bz2(k),cz2(k)
     +   ,ar2(k),br2(k),cr2(k)
     +   ,u_sxxtmp2(k),u_sxrtmp2(k)
     +   ,u_srxtmp2(k),u_srrtmp2(k)
     +   )
      endif
	   if (k.eq.1) then
	      Sxx2(i,k)=u_sxxtmp1(k)
	      Sxr2(i,k)=u_sxrtmp1(k)
	      Srx2(i,k)=u_srxtmp1(k)
	      Srr2(i,k)=u_srrtmp1(k)

	      Rxx2(i,k)=v_sxxtmp1(k)
	      Rxr2(i,k)=v_sxrtmp1(k)
	      Rrx2(i,k)=v_srxtmp1(k)
	      Rrr2(i,k)=v_srrtmp1(k)
	   else 
	      Sxx2(i,k)=u_sxxtmp1(k)+u_sxxtmp2(k-1)
	      Sxr2(i,k)=u_sxrtmp1(k)+u_sxrtmp2(k-1)
	      Srx2(i,k)=u_srxtmp1(k)+u_srxtmp2(k-1)
	      Srr2(i,k)=u_srrtmp1(k)+u_srrtmp2(k-1)

	      Rxx2(i,k)=v_sxxtmp1(k)+v_sxxtmp2(k-1)
	      Rxr2(i,k)=v_sxrtmp1(k)+v_sxrtmp2(k-1)
	      Rrx2(i,k)=v_srxtmp1(k)+v_srxtmp2(k-1)
	      Rrr2(i,k)=v_srrtmp1(k)+v_srrtmp2(k-1)
	   end if
	   end do

	    Sxx2(i,N+1)=u_sxxtmp2(N)
	    Sxr2(i,N+1)=u_sxrtmp2(N)
	    Srx2(i,N+1)=u_srxtmp2(N)
	    Srr2(i,N+1)=u_srrtmp2(N)
	   
	    Rxx2(i,N+1)=v_sxxtmp2(N)
	    Rxr2(i,N+1)=v_sxrtmp2(N)
	    Rrx2(i,N+1)=v_srxtmp2(N)
	    Rrr2(i,N+1)=v_srrtmp2(N)
99      continue

        Sxx2=-beta2*Sxx1
	  Sxr2=-beta2*Sxr1
	  Srx2=-beta2*Srx1
	  Srr2=-beta2*Srr1

        Rxx2=-beta4*Rxx1
	  Rxr2=-beta4*Rxr1
	  Rrx2=-beta4*Rrx1
	  Rrr2=-beta4*Rrr1

	  do i=1,N+1
	    Rxx2(i,i)=1.0D0+Rxx2(i,i)
	    Rrr2(i,i)=1.0D0+Rrr2(i,i)
	  end do
c+++++++++++++++++++++++==*******************************
c  Am=[Sxx -Sxr; -Srx Srr]
        do i=1,N+1
	    do k=1,N+1
	      Am(i,k)=Sxx1(i,k)
	      Am(i,k+N+1)=Sxr1(i,k)
	      Am(i,k+Nm)=Sxx2(i,k)
	      Am(i,k+Nm+N+1)=Sxr2(i,k)

	      Am(i+N+1,k)=Srx1(i,k)
	      Am(i+N+1,k+N+1)=Srr1(i,k)
	      Am(i+N+1,k+Nm)=Srx2(i,k)
	      Am(i+N+1,k+Nm+N+1)=Srr2(i,k)

	      Am(i+Nm,k)=Rxx1(i,k)
	      Am(i+Nm,k+N+1)=Rxr1(i,k)
	      Am(i+Nm,k+Nm)=Rxx2(i,k)
	      Am(i+Nm,k+Nm+N+1)=Rxr2(i,k)

	      Am(i+Nm+N+1,k)=Rrx1(i,k)
	      Am(i+Nm+N+1,k+N+1)=Rrr1(i,k)
	      Am(i+Nm+N+1,k+Nm)=Rrx2(i,k)
	      Am(i+Nm+N+1,k+Nm+N+1)=Rrr2(i,k)
	    end do
	  end do

        call dgesv(Nm2, 1, Am, neq, ipiv, bb, neq, info)
	   uz1(1:N+1)=bb(1:N+1)       
	   ur1(1:N+1)=bb(N+2:Nm)       
	   uz2(1:N+1)=bb(Nm+1:Nm+N+1)
	   ur2(1:N+1)=bb(Nm+N+2:Nm2)
	end if

c-----------
c  update marker points
c-----------
	    do k=1,N+1
	       un1(k)=uz1(k)*an11(k)+ur1(k)*an12(k)
	       ut1(k)=uz1(k)*an12(k)-ur1(k)*an11(k)
	       znew1(k)=z1(k)+dt*un1(k)*an11(k)
	       rnew1(k)=r1(k)+dt*un1(k)*an12(k)

	       un2(k)=uz2(k)*an21(k)+ur2(k)*an22(k)
	       ut2(k)=uz2(k)*an22(k)-ur2(k)*an21(k)
	       znew2(k)=z2(k)+dt*un2(k)*an21(k)
	       rnew2(k)=r2(k)+dt*un2(k)*an22(k)
	     end do

c+++++++++++++++++++++++++++++++++++++++++++++
c         call regridf(znew,rnew,z,r,s,N,time,dt)
       call regridf(znew1,rnew1,z1,r1,a1,s1,N,time,dt)
       call regridf(znew2,rnew2,z2,r2,a2,s2,N,time,dt)
c------------------
         cvol1=0.0D0
          do j=1,N
       	   hh(j)=z1(j+1)-z1(j)
	       cvol1=cvol1+hh(j)*(r1(j)**2+r1(j+1)**2)
	     end do
	     cvol1=cvol1*pi
c	   write (*,*) 'ratio=', cvol1/cvol10

         cvol2=0.0D0
          do j=1,N
       	   hh(j)=z2(j+1)-z2(j)
	       cvol2=cvol2+hh(j)*(r2(j)**2+r2(j+1)**2)
	     end do
	     cvol2=cvol2*pi
c	   write (*,*) 'ratio=', cvol2/cvol20

c scale the drop
          vscale1=(cvol1/cvol10)**(1.D0/2.D0)
          vscale2=(cvol2/cvol20)**(1.D0/2.D0)
          do j=1,N+1
cc	       z(j)=z(j)/vscale
	       r1(j)=r1(j)/vscale1
	       r2(j)=r2(j)/vscale2
	    end do
c++++++++++++++++++++++++++++++++++++++++++++ END IF
c         rmin = minval(r(N/4:7*N/8))
         rmin1 = minval(r1(5:N-4))
         rmin2 = minval(r2(5:N-4))
         vnmax1 = maxval(abs(un1))
         vnmax2 = maxval(abs(un2))

	   if (mod(jj,nplt2).eq.0) then
	      write (10,*) time, vnmax
	      write (11,*) time, rmin
	      write (12,*) time, v1max
         endif
	   if (mod(jj,nplt).eq.0) then
	      jjj=jjj+1
	     if(mod(jjj,2).eq.0)then
	        do j=N+1,1,-1
		      write (20,*) z1(j), r1(j)
	          write(22,*) z1(j), uz1(j)
c	          write (23,*) z1(j), ur1(j)
		      write (21,*) z2(j), r2(j)
	          write(24,*) z2(j), uz2(j)
c	          write (25,*) z2(j), ur2(j)
	        end do
	     else
	        do j=1,N+1
		      write (20,*) z1(j), r1(j)
	          write(22,*) z1(j), uz1(j)
c	          write (23,*) z1(j), ur1(j)
		      write (21,*) z2(j), r2(j)
	          write(24,*) z2(j), uz2(j)
c	          write (25,*) z2(j), ur2(j)
	        end do
	     end if
          end if	      

c------

	   write(9,*) 'rmin1=', rmin1
	   write(9,*) 'rmin2=', rmin2
	   write(9,*) 'time=', time

c         if(rmin.lt.0.01.and.nplt.gt.40)then
c	       nplt=12
c	   end if
c	   if(rmin.lt.0.02d0.and.dt.gt.dt2)then
c	    dt=dt2
c	   endif
c	   if(rmin.lt.0.01d0.and.dt.gt.dt3)then
c	    dt=dt3
c	   endif
c	   if(rmin.lt.0.005d0.and.dt.gt.dt3/2.0)then
c	    dt=dt3/2.0
c	   endif

         if (vnmax1.gt.150.0D0) then
	       write(9,*) "velocity too big"
	       write(9,*) 'vmax=', vnmax
	       write(9,*) 'N=', N
	       exit
	    end if
	   da=minval(z1(1:N+1)-z2(1:N+1))
         if (da.lt.0.002D0) then
	       write(9,*) "touch"
	       write(9,*) 'da=', da
	       exit
	    end if
         if (dabs(vnmax1).lt.0.001D0.and.
     +	   dabs(vnmax2).lt.0.001D0) then
	       write(9,*) "velocity too small"
             write(9,*) 'vmax_min=', vnmax
	       exit
	    end if

         jj = jj + 1
	end do

      write(9,*) "stop time is"
	write(9,*) time
      write(9,*) 'dt=', dt
      write(9,*) 'N=', N
	write(9,*) 'inv(Bo)=', Boi
	
!
c output data
	open (41,file='xy_final1.dat')
	open (42,file='xy_final2.dat')

	   jjj=jjj+1
	     if(mod(jjj,2).eq.0)then
	        do j=N+1,1,-1
		      write (20,*) z1(j), r1(j)
	          write(22,*) z1(j), uz1(j)
c	          write (23,*) z1(j), ur1(j)
		      write (21,*) z2(j), r2(j)
	          write(24,*) z2(j), uz2(j)
c	          write (25,*) z2(j), ur2(j)
	        end do
	     else
	        do j=1,N+1
		      write (20,*) z1(j), r1(j)
	          write(22,*) z1(j), uz1(j)
c	          write (23,*) z1(j), ur1(j)
		      write (21,*) z2(j), r2(j)
	          write(24,*) z2(j), uz2(j)
c	          write (25,*) z2(j), ur2(j)
	        end do
	     end if

	 do j=1,N+1
	   write (41,*) z1(j), r1(j)
	   write (42,*) z2(j), r2(j)
	 end do
c close 23
      end program
