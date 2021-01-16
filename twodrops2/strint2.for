      subroutine strint2
     +
     +   (z1,r1
     +   ,s1,s2
     +   ,anx1,anx2,any1,any2
     +   ,z0,r0,anx0,any0
     +   ,n,arcl
     +   ,Az,Bz,Cz
     +   ,Ar,Br,Cr
     +   ,y1,y2
     +   ,y3,y4
     +   )

      Implicit Double Precision (a-h,o-z)

	dimension x1(10),w1(10)
      dimension sg(10),rg(10),zg(10),amg(10)
	dimension anzg(10),anrg(10)
      dimension Qxxx1(10),Qxxy1(10),Qxyx1(10),Qxyy1(10)
      dimension Qyxx1(10),Qyxy1(10),Qyyx1(10),Qyyy1(10)
      dimension Qxxx2(10),Qxxy2(10),Qxyx2(10),Qxyy2(10)
      dimension Qyxx2(10),Qyxy2(10),Qyyx2(10),Qyyy2(10)
      
	eps = 0.000000001D0
      pi  = 3.1415 92653 58979 32384 D0

	call Gauss_Legendre (8,x1,w1)

      if (n.eq.1) then
	  do i=1,8
	     amg(i)=(1.0D0-x1(i))/2.0D0
	   end do
	else if (n.eq.2) then
	   do i=1,8
	      amg(i)=(1.0D0+x1(i))/2.0D0
	    end do
	end if

	a = (s2-s1)/2.0D0
	b = (s2+s1)/2.0D0

	do j=1,8
	   sg(j) = a*x1(j) + b
	   anzg(j) = (anx2-anx1)/2.0D0*x1(j) + (anx2+anx1)/2.0D0
	   anrg(j) = (any2-any1)/2.0D0*x1(j) + (any2+any1)/2.0D0
	   xd = sg(j) - s1
	   rg(j) = Ar*xd**3 + Br*xd**2 + Cr*xd + r1
	   zg(j) = Az*xd**3 + Bz*xd**2 + Cz*xd + z1
	end do

      cst = -4.0*pi/arcl

c---   calculate the integral
      yxxtmp1=0.0D0
c	yxxtmp2=0.0D0
      yxrtmp1=0.0D0
c	yxrtmp2=0.0D0
      yrxtmp1=0.0D0
c	yrxtmp2=0.0D0
      yrrtmp1=0.0D0
c	yrrtmp2=0.0D0

	Iopt=2

      do j=1,8
	  call sgfax_str     
     + (Iopt,rg(j),zg(j),r0,z0
     + ,Qxxx1(j),Qxxy1(j),Qxyx1(j),Qxyy1(j)
     + ,Qyxx1(j),Qyxy1(j),Qyyx1(j),Qyyy1(j)
     + )


	  yxxtmp1=yxxtmp1+a*amg(j)*(Qxxx1(j)*anzg(j)+
     +	  Qxxy1(j)*anrg(j)+cst*anx0*anzg(j)*rg(j))*w1(j)

	  yxrtmp1=yxrtmp1+a*amg(j)*(Qxyx1(j)*anzg(j)+
     +	  Qxyy1(j)*anrg(j)+cst*anx0*anrg(j)*rg(j))*w1(j)

	  yrxtmp1=yrxtmp1+a*amg(j)*(Qyxx1(j)*anzg(j)+
     +	  Qyxy1(j)*anrg(j)+cst*any0*anzg(j)*rg(j))*w1(j)

	  yrrtmp1=yrrtmp1+a*amg(j)*(Qyyx1(j)*anzg(j)+
     +	  Qyyy1(j)*anrg(j)+cst*any0*anrg(j)*rg(j))*w1(j)
      end do

	y1 = yxxtmp1
	y2 = yxrtmp1
	y3 = yrxtmp1
	y4 = yrrtmp1

      return
	end

c----------------------------------------------------------




