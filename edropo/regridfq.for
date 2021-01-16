      subroutine regridf(znew,rnew,z,r,s,qs,N,time,dt)
c
c
      implicit double precision (a-h,o-z)
	parameter(neq=1024)

	dimension z(neq),r(neq),znew(neq),rnew(neq),snew(neq),s(neq)
      dimension d1(500),b1(500),c1(500),d2(500),b2(500),c2(500)
      dimension dq(neq),bq(neq),cq(neq)
      dimension qs(neq),qtmp(neq)

      call arc_evl(znew,rnew,snew,N
     +  ,d2,c2,b2
     +  ,d1,c1,b1)
c-----	
       qtmp=qs

         h=snew(N+1)/N
          call splc_clm
     +  (N
     +  ,snew,qtmp
     +  ,0.0d0
     +  ,0.0d0
     +  ,dq,cq,bq
     +  )
	 bq(N+1)=0.0d0

      	do j=1,N+1
	      s(j)=(j-1)*h
	    end do      

         do k=1,N+1
      	  z(k)=seval(N+1,s(k),snew,znew,b1,c1,d1)
	      r(k)=seval(N+1,s(k),snew,rnew,b2,c2,d2)
      	    qs(k)=seval(N+1,s(k),snew,qtmp,bq,cq,dq)
         end do
       qs(1)=qtmp(1)
       qs(N+1)=qtmp(N+1)

      return
	end
c--------
      subroutine regridfadp(znew,rnew,z,r,s,N,Umx,dt)
c
      implicit double precision (a-h,o-z)
	parameter(neq=1024)
	dimension z(neq),r(neq),znew(neq),rnew(neq),snew(neq),s(neq)
	dimension sold(neq),zold(neq),rold(neq)
      dimension d1(neq),b1(neq),c1(neq),d2(neq),b2(neq),c2(neq)
	dimension xtmp(neq),ytmp(neq),stmp(neq),sd(neq)
	dimension xtmp0(neq),ytmp0(neq),stmp0(neq),qs(neq)
      dimension dq(neq),cq(neq),bq(neq)
      dimension dd1(neq),cc1(neq),bb1(neq)
      dimension dd2(neq),cc2(neq),bb2(neq),qtmp(neq)
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
c--
c IF *****************	

      IF(fmax.le.8.d0)then ! equal-spacing regridding
      	 do j=1,N+1
	    s(j)=(j-1)*h
	  enddo      
         do k=1,N+1
      	   z(k)=seval(N+1,s(k),snew,znew,b1,c1,d1)
	    r(k)=seval(N+1,s(k),snew,rnew,b2,c2,d2)
c      	    qs(k)=seval(N+1,s(k),snew,qtmp,bq,cq,dq)
         enddo
       dtnew=0.01d0/fmax

c       qs(1)=qtmp(1)
c       qs(N+1)=qtmp(N+1)
          r(1)=0.d0
          r(N+1)=0.d0
        if(dt.gt.dtnew) dt=dtnew
c ELSE ***************
       Elseif (fmax.gt.8.d0.and.Ir.eq.1)then
         do 21 j=1,N/2
c            sd(j) = dsqrt((znew(j+1)-znew(1))**2+rnew(j+1)**2)
            sd(j) = 0.5d0/dabs(f(j)+1)
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
c      	    qs(k)=seval(N+1,s(k),snew,qtmp,bq,cq,dq)
        end do
c       qs(1)=qtmp(1)
c       qs(N+1)=qtmp(N+1)
       r(1)=0.0d0
       r(N+1)=0.0d0
c       dtnew=0.001d0/(fmax+2.d0)/(Umx+5.d0)
       dtnew=0.1d0/(fmax+10.d0)**2

        if(dt.gt.dtnew) dt=dtnew
c-------------------------------------------------------
      Else
c         do j=1,Ir-1
c           sd(j) = dsqrt((znew(j)-znew(Ir))**2+rnew(j)**2)
cc           sd(j) = 0.5d0/(dabs(f(j))+1.d0)
c         enddo
c         do j=Ir,N/2
c            sd(j) = dsqrt((znew(j+1)-znew(Ir))**2+rnew(j+1)**2)
c          enddo
         do j=1,N
           fk = dabs(f(j))
            if(fk.gt.1.0.and.fk.lt.100.0)then
               sd(j) = 1.0/fk
            elseif(fk.ge.100.0)then
               sd(j) = 1.2/fk
             else
               sd(j) = 0.1d0*f(j)**2
            endif
          enddo
         stmp(1)=0.d0
          stmp(N+1)=RL
         stmp(N/2+1)=RL/2
c         stmp(Ir)=snew(Ir)
        sum1=0.d0
        sum2=0.d0
c        do j=1,Ir-1
c            sum1 = sum1 + sd(j)
c        enddo
c        do j=Ir,N/2
c            sum2 = sum2 + sd(j)
c        enddo
c         hs1 = snew(Ir)/sum1
c         hs2 = (RLH-snew(Ir))/sum2
          do j=1,N
            sum1 = sum1 + sd(j)
          enddo
           hs = RL/sum1
          do j=2,N
             stmp(j) = stmp(j-1) + hs*sd(j-1)
          enddo
c         do j=Ir-1,2,-1
c             stmp(j)=stmp(j+1)-hs1*sd(j)
c         enddo
c         do j=Ir,N/2
c            stmp(j+1)=stmp(j)+hs2*sd(j)
c         enddo
         do j=2,N/2
            stmp(N/2+j)=stmp(N/2-j+2)
         enddo
c--
        s = stmp   
         do k=1,N+1
      	  z(k)=seval(N+1,s(k),snew,znew,b1,c1,d1)
	  r(k)=seval(N+1,s(k),snew,rnew,b2,c2,d2)
c      	    qs(k)=seval(N+1,s(k),snew,qtmp,bq,cq,dq)
        end do
c       qs(1)=qtmp(1)
c       qs(N+1)=qtmp(N+1)
       r(1)=0.0d0
       r(N+1)=0.0d0
c       dtnew=0.001d0/(fmax+10.d0)/(Umx+20.d0)
       dtnew=0.05d0/(Umx+20.d0)**2
        if(dt.gt.dtnew) dt=dtnew
      ENDIF
c*************************** check max separation then
c DONE !

      return
	end


c-----------------------------------------------------
