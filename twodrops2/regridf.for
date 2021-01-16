      subroutine regridf(znew,rnew,z,r,aa,s,N,time,dt)
c
c
      implicit double precision (a-h,o-z)
	parameter(neq=1024)

	dimension z(neq),r(neq),znew(neq),rnew(neq),snew(neq),s(neq)
      dimension d1(neq),b1(neq),c1(neq),d2(neq),b2(neq),c2(neq)
	dimension coefz(neq,neq),coefr(neq,neq)

      call arc_evl(znew,rnew,snew,N,aa
     +  ,d2,c2,b2
     +  ,d1,c1,b1)
	RL=snew(N+1)
c-----	
c      rmin=minval(rnew(6:N-20))
c	if(time.gt.1.2D0 .and.N.lt.100)then
c	   N=128
c	   dt=0.3D0*(1.0D0/N)**(1.5D0)
c	endif
c      if(rmin.lt.0.025 .and. N.lt.160)then
c	   N=160
c	   dt=0.09D0*(1.0D0/N)**(1.5D0)
c	end if
c      if(rmin.lt.0.006 .and. N.lt.320)then
c	   N=320
c	   dt=0.12D0*(1.0D0/N)**(1.5D0)
c	end if

         h=RL/N

      	do j=1,N+1
	      s(j)=(j-1)*h
	    end do      

         do k=1,N+1
      	  z(k)=seval(N+1,s(k),snew,znew,b1,c1,d1)
	      r(k)=seval(N+1,s(k),snew,rnew,b2,c2,d2)
         end do

	r(1) = 0.0
	r(N+1)=0.0
      return
	end
c---------------------------------------------------------
      subroutine zabarankU(aa1,aa2,dd,RR1,RR2,eb,ram1,ram2,vel)
c
      implicit double precision (a-h,o-z)
	parameter(neq=1024)

c	dimension aa1(neq),aa2(neq)
c	dimension coefz(neq,neq),coefr(neq,neq)

c-----	
   	fac1 = 3.0*ram1+2.
	fac2 = 2.0*RR1+1.0
	fac22 = 2.0*RR2+1.0
	dd2 = dd**2
	dd3 = dd2*dd

	tmp1 = 1.0-3.*aa1**2*ram1/4.0/fac1/dd2 -aa2**2/4.0/dd2-
     +   (RR1-1.0)/2.0/fac2*aa1**3/dd3
	tmp2 = -3.*eb/(1.+ram2)/fac22*tmp1

	tmp3 = eb/2.0/fac1/(3.0*RR1+2.0)-
     +  (RR1-1.0)*((ram1+1.)/fac1-aa1/4./dd)

	tmp4 = 5.0*(RR2-1.)/fac2*aa1**2/dd2*tmp3

	tmp5 = 3.0*aa2**3/20./fac22/dd2*(tmp2+tmp4)
	
	tmp6 = 9.*(5.*ram2+2.0)*eb*aa1**2*aa2**3/160./(ram1+1.)/(ram2+1.)
	tmp6 = tmp6/fac2**2/dd2/dd3

	vel = tmp5 + tmp6

      return
	end