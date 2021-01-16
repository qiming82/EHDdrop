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
c IF *****************	
      IF(rmin.gt.0.01d0)then ! equal-spacing regridding
      	 do j=1,N+1
	    s(j)=(j-1)*h
	  enddo      
         do k=1,N+1
      	   z(k)=seval(N+1,s(k),snew,znew,b1,c1,d1)
	    r(k)=seval(N+1,s(k),snew,rnew,b2,c2,d2)
         enddo
       dtnew=0.0025d0/fmax
          r(1)=0.d0
          r(N+1)=0.d0
        if(dt.gt.dtnew) dt=dtnew

c ELSE ***************
       Else
c       write(7,*) 'adaptive'
       Ir=0  ! index of minimum
	 do 1 k=5,N/2+1
         if(dabs(rmin-rnew(k)).lt.0.0000001D0)then
	      Ir=k
	   endif
 1         continue
         do 21 j=1,Ir-1
c            sd(j) = dsqrt((znew(j)-znew(Ir))**2+rnew(j)**2)  
c            sdtmp = 1.d0/dabs(f(j))
c            if(sd(j).lt.sdtmp) sd(j)=sdtmp
            sd(j)=1.d0/dabs(f(j))
 21       continue
         do 22 j=Ir+1,N/2+1
           sd(j-1)=dsqrt((znew(j)-znew(Ir))**2+rnew(j)**2)
           sdtmp=1.d0/dabs(f(j-1))
             if(sd(j-1).lt.sdtmp) sd(j-1)=sdtmp
 22         continue

c-- change !
        RLH = snew(Ir)
        suml=0.d0
        do j=1,Ir-1
            suml = suml + sd(j)
        enddo
         hs1 = RLH/suml
 
        suml=0.d0
        do j=Ir,N/2
            suml = suml + sd(j)
        enddo
         hs2 = (RL/2.d0-RLH)/suml

         stmp(1)=0.0
          stmp(Ir)=snew(Ir)
         stmp(N+1-(Ir-1))=RL-RLH
         stmp(N+1)=RL
         do j=Ir-1,2,-1
              stmp(j)=stmp(j+1)-hs1*sd(j)
             stmp(N+2-j)=stmp(N+1-j)+hs1*sd(j)
          enddo
         do j=Ir+1,N/2
             stmp(j)=stmp(j-1)+hs2*sd(j-1)
             stmp(N+2-j)=stmp(N+2-j+1)-hs2*sd(j-1)
         enddo
          stmp(N/2+1)=RL/2.d0
        s = stmp   
         do k=1,N+1
      	  z(k)=seval(N+1,s(k),snew,znew,b1,c1,d1)
	  r(k)=seval(N+1,s(k),snew,rnew,b2,c2,d2)
        end do
       r(1)=0.0d0
       r(N+1)=0.0d0

       dtnew=0.025d0/fmax
        if(dt.gt.dtnew) dt=dtnew

      ENDIF
c*************************** check max separation then
c DONE !
      return
	end
c-----------------------------------------------------
      subroutine regridfadp2(znew,rnew,z,r,s,N,umx,dt)
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
         fmax = maxval(dabs(f(1:N/2+1)))
       Ir=0  ! index of minimum
	 do 1 k=1,N/2+1
         if(dabs(fmax-f(k)).lt.0.000001D0)then
	      Ir=k
	   endif
 1         continue
c IF *****************	
      IF(fmax.le.6.d0)then ! equal-spacing regridding
      	 do j=1,N+1
	    s(j)=(j-1)*h
	  enddo      
         do k=1,N+1
      	   z(k)=seval(N+1,s(k),snew,znew,b1,c1,d1)
	    r(k)=seval(N+1,s(k),snew,rnew,b2,c2,d2)
         enddo
       dtnew=0.01d0/fmax/umx
          r(1)=0.d0
          r(N+1)=0.d0
        if(dt.gt.dtnew) dt=dtnew
c ELSE ***************
       Elseif (fmax.gt.6.d0.and.Ir.eq.1)then
         do 21 j=1,N/2
c            sd(j) = dsqrt((znew(j+1)-znew(1))**2+rnew(j+1)**2)
            sd(j) = 1.d0/dabs(f(j))
c            sd(N-j+1) = sd(j)
 21       continue
c-- change !
        RLH = snew(N/2+1)
        sum1=0.d0
c        sum2=0.d0
        do j=1,N/2
            sum1 = sum1 + sd(j)
        enddo
         hs1 = RLH/sum1
 
         stmp(1)=0.0
         stmp(N+1)=RL
         do j=1,N/2
             stmp(j+1)=stmp(j)+hs1*sd(j)
         enddo
          stmp(N/2+1)=RL/2.d0
         do j=N,N/2+2,-1
             stmp(j)=stmp(j+1)-hs1*sd(N+1-j)
         enddo
        s = stmp   
         do k=1,N+1
      	  z(k)=seval(N+1,s(k),snew,znew,b1,c1,d1)
	  r(k)=seval(N+1,s(k),snew,rnew,b2,c2,d2)
        end do
       r(1)=0.0d0
       r(N+1)=0.0d0
       dtnew=0.001d0/fmax/umx
        if(dt.gt.dtnew) dt=dtnew
      Else
         do j=1,Ir-1
           sd(j) = dsqrt((znew(j)-znew(Ir))**2+rnew(j)**2)
c           sd(j) = 0.5d0/(dabs(f(j))+1.d0)
         enddo
         do j=Ir,N/2
            sd(j) = dsqrt((znew(j+1)-znew(Ir))**2+rnew(j+1)**2)
          enddo

         stmp(1)=0.d0
          stmp(N+1)=RL
         stmp(N/2+1)=RL/2
         stmp(Ir)=snew(Ir)
        sum1=0.d0
        sum2=0.d0
        do j=1,Ir-1
            sum1 = sum1 + sd(j)
        enddo
        do j=Ir,N/2
            sum2 = sum2 + sd(j)
        enddo
         hs1 = snew(Ir)/sum1
         hs2 = (RLH-snew(Ir))/sum2
         do j=Ir-1,2,-1
             stmp(j)=stmp(j+1)-hs1*sd(j)
         enddo
         do j=Ir,N/2
            stmp(j+1)=stmp(j)+hs2*sd(j)
         enddo
         do j=2,N/2
            stmp(N/2+j)=stmp(N/2-j+2)
         enddo
c--
        s = stmp   
         do k=1,N+1
      	  z(k)=seval(N+1,s(k),snew,znew,b1,c1,d1)
	  r(k)=seval(N+1,s(k),snew,rnew,b2,c2,d2)
        end do
       r(1)=0.0d0
       r(N+1)=0.0d0
       dtnew=0.002d0/fmax/umx
        if(dt.gt.dtnew) dt=dtnew
      ENDIF
c*************************** check max separation then
c DONE !
      return
	end

c--------------------------------------------------------------
      subroutine regridfadp3(znew,rnew,z,r,s,N,umx,dt)
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
       dimension hh(neq)

      pi=4*datan(1.d0)

c   snew is the current arclength-spacing after advancing pts
      RL=snew(N+1) ! total lenght of the arclength
c       rmin=minval(rnew(5:N/2+1)) ! rmin - before regridding

        Nold = N

      call arc_evl(znew,rnew,snew,N
     +  ,d2,c2,b2
     +  ,d1,c1,b1)

       stmp(1:N+1)=snew(1:N+1)

 99      continue
         h=RL/N  ! h would be a benchmark 

         if(N.ge.300) go to 100
c-- get the grid size
         do j=1,N
             hh(j)= stmp(j+1) - stmp(j)
         enddo
c -- check if points are needed to be added or removed
          do j=1,N
            if(hh(j).gt.2.d0*h)then
               do k=j,N+1
                 s(k+1) = stmp(k)
               enddo
               do k=1,j-1
                  s(k) = stmp(k)
                enddo
                s(j) = (stmp(j)+stmp(j+1))/2.d0
               N = N + 1
                stmp(1:N+1)=s(1:N+1)
               go to 99
            endif
            if(hh(j).le.0.5d0*h)then
                 do k=1,j-1
                  s(k) = stmp(k)
                enddo
                 do k=j,N
                   s(k) = stmp(k+1)
                  enddo
                N = N - 1
                stmp(1:N+1)=s(1:N+1)
                go to 99
             endif
         enddo

 100     continue

          if(N.eq.Nold)then
      	 do j=1,N+1
	    s(j)=(j-1)*h
	  enddo      
            else
      	 do j=1,N+1
	    s(j)=stmp(j)
	  enddo      
         endif

         do k=1,N+1
      	   z(k)=seval(N+1,s(k),snew,znew,b1,c1,d1)
	    r(k)=seval(N+1,s(k),snew,rnew,b2,c2,d2)
         enddo
          r(1)=0.d0
          r(N+1)=0.d0

       dtnew=0.01d0/(umx+5.d0)
        if(dt.gt.dtnew) dt=dtnew

      return
	end

c-----------------------------------------------------
