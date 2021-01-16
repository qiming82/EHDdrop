      subroutine regridfadp(znew,rnew,z,r,s,N,time,dt)
c
      implicit double precision (a-h,o-z)
	parameter(neq=1024)

	dimension z(neq),r(neq),znew(neq),rnew(neq),snew(neq),s(neq)
	dimension sold(neq),zold(neq),rold(neq)
      dimension d1(neq),b1(neq),c1(neq),d2(neq),b2(neq),c2(neq)
	dimension xtmp(neq),ytmp(neq),stmp(neq),sd(neq)
	dimension xtmp0(neq),ytmp0(neq),stmp0(neq)
      dimension dd1(neq),cc1(neq),bb1(neq)
      dimension dd2(neq),cc2(neq),bb2(neq)
      dimension dz(neq),dr(neq),ds(neq),f(neq),ddz(neq),ddr(neq)

      pi=4*datan(1.d0)

c first get the new arclength
      call arc_evl(znew,rnew,snew,N
     +  ,d2,c2,b2
     +  ,d1,c1,b1)
c   snew is the current arclength-spacing after advancing pts
      RL=snew(N+1) ! total lenght of the arclength
       rmin=minval(rnew(5:N/2+1)) ! rmin - before regridding
         h=RL/N
c-- calculate curvature first
	    do j=1,N+1
	       dr(j)=b2(j)
	       ddr(j)=2.0D0*c2(j)
	       dz(j)=b1(j)
	       ddz(j)=2.0D0*c1(j)
	       ds(j)=dsqrt(dz(j)**2+dr(j)**2)
	       call df(rnew(j),dz(j),ddz(j),dr(j),ddr(j),ds(j),f(j))
	    enddo
         fmax = maxval(f(1:N/2+1))
		 
       Ir=0  ! index of minimum
	 do 1 k=1,N/2+1
         if(dabs(fmax-f(k)).lt.0.0000001D0)then
	      Ir=k
	   endif
 1         continue		
 
c IF *****************	 
	  IF(fmax.le.9.0) then
c regular update
       write(*,*) 'regular', fmax, Ir, time
      	 do j=1,N+1
	    s(j)=(j-1)*h
	  enddo      
         do k=1,N+1
      	   z(k)=seval(N+1,s(k),snew,znew,b1,c1,d1)
	    r(k)=seval(N+1,s(k),snew,rnew,b2,c2,d2)
         enddo
c ELSE *****************	
      Elseif(fmax.gt.9.d0.and.Ir.eq.1)then
       write(*,*) 'adaptive', fmax, Ir, time

         do 21 j=1,N/2
            sd(j)=1.d0/(dabs(f(j))+1.0)
 21       continue

c-- change !
        RLH = snew(N/2+1)
        suml=0.d0
        do j=1,N/2
            suml = suml + sd(j)
        enddo
         hs1 = RLH/suml

         stmp(1)=0.0
          stmp(N/2+1)=snew(N/2+1)
         do j=1,N/2-1
             stmp(j+1)=stmp(j)+hs1*sd(j)
         enddo
		 do j= N/2+2,N
		   stmp(j) = stmp(j-1)+ hs1*sd(N+2-j)
		 enddo
		 stmp(N+1)=snew(N+1)
		 
        s = stmp   
         do k=1,N+1
      	  z(k)=seval(N+1,s(k),snew,znew,b1,c1,d1)
		r(k)=seval(N+1,s(k),snew,rnew,b2,c2,d2)
        end do
       r(1)=0.0d0
       r(N+1)=0.0d0
	  ELSE
c regular update
       write(*,*) 'regular', fmax, Ir
      	 do j=1,N+1
	    s(j)=(j-1)*h
	  enddo      
         do k=1,N+1
      	   z(k)=seval(N+1,s(k),snew,znew,b1,c1,d1)
	    r(k)=seval(N+1,s(k),snew,rnew,b2,c2,d2)
         enddo
		 ENDIF
c*************************** check max separation then
       r(1)=0.0d0
       r(N+1)=0.0d0
       dtnew=0.01d0/(5.0+fmax)
        if(dt.gt.dtnew) dt=dtnew
c DONE !
      return
	end
