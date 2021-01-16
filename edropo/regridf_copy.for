      subroutine regridf(znew,rnew,z,r,s,N,time,dt)
c
c
      implicit double precision (a-h,o-z)
	parameter(neq=1024*2)

	dimension z(neq),r(neq),znew(neq),rnew(neq),snew(neq),s(neq)
      dimension d1(neq),b1(neq),c1(neq),d2(neq),b2(neq),c2(neq)
	dimension coefz(neq,neq),coefr(neq,neq)

      call arc_evl(znew,rnew,snew,N
     +  ,d2,c2,b2
     +  ,d1,c1,b1)
	RL=snew(N+1)
c-----	
         h=RL/N

      	do j=1,N+1
	      s(j)=(j-1)*h
	    end do      

         do k=1,N+1
      	  z(k)=seval(N+1,s(k),snew,znew,b1,c1,d1)
	      r(k)=seval(N+1,s(k),snew,rnew,b2,c2,d2)
         end do
c      call sm_ax(N,r)
       r(1)=0.d0
	 r(N+1)=0.d0
       dtnew = 0.05d0/(um+5.d0)   !(fmax+1.d0)
       if(dt.gt.dtnew) dt= dtnew
      return
	end
c----------------------------------------------
      subroutine regridf2(znew,rnew,z,r,s,N,f,dt,Um)
c
c
      implicit double precision (a-h,o-z)
	parameter(neq=1024*2)

	dimension z(neq),r(neq),znew(neq),rnew(neq),snew(neq),s(neq)
      dimension d1(neq),b1(neq),c1(neq),d2(neq),b2(neq),c2(neq)
c	dimension coefz(neq,neq),coefr(neq,neq)
      dimension stmp(neq),sd(neq),f(neq),sgm(neq),rtmp(neq),ztmp(neq)

      call arc_evl(znew,rnew,snew,N
     +  ,d2,c2,b2
     +  ,d1,c1,b1)
	RL=snew(N+1)
c-----	
        h=RL/N
        h3=2.d0*h
	  hh=h/10.d0
          Np = N  !old # of segments

         fmax = maxval(dabs(f(1:Np+1)))
        
       if(fmax.lt.5.d0)then
      	do j=1,Np+1
	      s(j)=(j-1)*h
	  enddo      

         do k=1,Np+1
      	  z(k)=seval(Np+1,s(k),snew,znew,b1,c1,d1)
	      r(k)=seval(Np+1,s(k),snew,rnew,b2,c2,d2)
         enddo
c--
        else  
 2        continue
c-- calculate each line segment
          do j=1,Np
            sgm(j)=snew(j+1)-snew(j)
          enddo
c           stmp= snew
c-- search for adding points
         do j=1,Np
          if(sgm(j).gt.h3.and.Np.lt.180) go to 1  ! regridding
c	    if(sgm(j).lt.hh.and.j.lt.Np.and.Np.gt.32) go to 4
         enddo

	   if(Np.eq.N)then
      	do j=1,Np+1
	      s(j)=(j-1)*h
	  enddo      
         do k=1,Np+1
      	  z(k)=seval(Np+1,s(k),snew,znew,b1,c1,d1)
	      r(k)=seval(Np+1,s(k),snew,rnew,b2,c2,d2)
         enddo
	   endif
c-- now arclength is updated
         go to 3

 1       continue
          write(7,*) 'regridding + pt'
         sint=0.5d0*(snew(j)+snew(j+1))
c         zint=0.5d0*(znew(j)+znew(j+1))
c         rint=0.5d0*(rnew(j)+rnew(j+1))
      if(j.lt.Np)then
      call INT_QUAD1
     +
     +  (snew(j),snew(j+1),snew(j+2),sint
     +  ,znew(j),znew(j+1),znew(j+2),zint
     +  ,rnew(j),rnew(j+1),rnew(j+2),rint
     +  )
	else
      call INT_QUAD1
     +
     +  (snew(j-1),snew(j),snew(j+1),sint
     +  ,znew(j-1),znew(j),znew(j+1),zint
     +  ,rnew(j-1),rnew(j),rnew(j+1),rint
     +  )
      endif
c
         do i=Np+1,j+1,-1
             i1=i+1
            znew(i1)=znew(i)
            rnew(i1)=rnew(i)
            snew(i1)=snew(i)
         enddo
          snew(j+1)=sint
          znew(j+1)=zint
          rnew(j+1)=rint
           Np = Np + 1
         go to 2
c--
4        continue
           write(7,*) 'remove a pt'
          do i=j+2,Np+1
	      i1=i-1
            znew(i1)=znew(i)
            rnew(i1)=rnew(i)
            snew(i1)=snew(i)
	   enddo
	      Np = Np - 1
         go to 2
c--
 3       continue
        N = Np
         s=snew
         z=znew
         r=rnew
c         write(7,*) 'Np=', Np, 'N=', N
c      call arc_evl(z,r,s,N
c     +  ,d2,c2,b2
c     +  ,d1,c1,b1)
c-- loop end
        endif
      call sm_ax(N,r)

        r(1)=0.d0
	  r(N+1)=0.d0
       dtnew = 0.1d0/(um+1.d0)   !(fmax+1.d0)
       if(dt.gt.dtnew) dt= dtnew
      return
	end
c-------------------------------------------------------
      subroutine INT_QUAD1
     +
     +  (x1,x2,x3,x
     +  ,u1,u2,u3,u
     +  ,v1,v2,v3,v
     +  )

c---------------------------------
c Quadratic Lagrange Interpolation
c---------------------------------

      Implicit Double Precision (A-H,O-Z)

      p1 = (x-x2)*(x-x3)/((x1-x2)*(x1-x3))
      p2 = (x-x1)*(x-x3)/((x2-x1)*(x2-x3))
      p3 = (x-x1)*(x-x2)/((x3-x1)*(x3-x2))

      u = u1*p1 + u2*p2 + u3*p3
      v = v1*p1 + v2*p2 + v3*p3

c-----
c Done
c-----

      Return
      End
c-----------------------------------
      subroutine sm_ax
     +
     +  (NSG
     +  ,P
     +  )

c-----------------------------
c Smoothing of property P
c using the 5-point formula of
c Longuett-Higgins and Cokelet
c-----------------------------

      Implicit Double Precision (a-h,o-z)

      Double Precision P(0:900),Ps(0:902)

c--------
c prepare
c--------

      NSGa = NSG-1
      NSG1 = NSG+1
      NSG2 = NSG+2
      NSG3 = NSG+3

c----------------
c Dump and extend
c----------------

      Do i=1,NSG1
       Ps(i) = P(i)
      End Do

      Ps(0)    = Ps(2)
      Ps(NSG2) = Ps(NSG)
      Ps(NSG3) = Ps(NSGa)

c----------
c smoothing
c----------

      P(1) = - Ps(3)+4.0*Ps(2)+10.0*Ps(1) 
     +              +4.0*Ps(2)     -Ps(3)

      Do i=2,NSG1
       P(i) = - Ps(i-2)+4.0*Ps(i-1)+10.0*Ps(i) 
     +                 +4.0*Ps(i+1)     -Ps(i+2)
      End Do

      Do i=1,NSG1
       P(i) = P(i)/16.00
      End Do
      Return
      End
