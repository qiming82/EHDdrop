      subroutine geten2
     +
     +   (z,r
     +   ,Az,Bz,Cz
     +   ,Ar,Br,Cr
     +   ,N
     +   ,s
     +   ,y
     +   )

      Implicit Double Precision (a-h,o-z)

c	dimension z(500),r(500),phi(500),s(500)
c	dimension F(500,500),RHS(500),y(500),H(500,500)
c	dimension h1(500),h2(500),f1(500),f2(500),ml(500,500),mu(500,500)
c	dimension Az(513),Bz(513),Cz(513),Ar(513),Br(513),Cr(513)
      parameter(neq=1024)
	dimension z(neq),r(neq),s(neq),phi(neq)
	dimension F(neq,neq),H(neq,neq),RHS(neq),x(neq),y(neq)
	dimension h1(neq),h2(neq),f1(neq),f2(neq),ml(neq,neq),mu(neq,neq)
	dimension Az(neq),Bz(neq),Cz(neq),Ar(neq),Br(neq),Cr(neq)
c	dimension s1(neq),s2(neq)
      dimension ipiv(neq)

      do i=1,N+1
	  do k=1,N
           call int_pot_s
     +
     +   (z(k),r(k)
     +   ,s(k),s(k+1)
     +   ,z(i),r(i),s(i)
     +   ,1
     +   ,Az(k),Bz(k),Cz(k)
     +   ,Ar(k),Br(k),Cr(k)
     +   ,f1(k)
     +   )
           call int_pot_s
     +
     +   (z(k),r(k)
     +   ,s(k),s(k+1)
     +   ,z(i),r(i),s(i)
     +   ,2
     +   ,Az(k),Bz(k),Cz(k)
     +   ,Ar(k),Br(k),Cr(k)
     +   ,f2(k)
     +   )
c	   if(k.eq.N)then
c	   write(*,*) f1(k)
c	   end if
           call int_pot_d2
     +
     +   (z(k),r(k)
     +   ,s(k),s(k+1)
     +   ,z(i),r(i),s(i)
     +   ,1
     +   ,Az(k),Bz(k),Cz(k)
     +   ,Ar(k),Br(k),Cr(k)
     +   ,h1(k)
     +   )
           call int_pot_d2
     +
     +   (z(k),r(k)
     +   ,s(k),s(k+1)
     +   ,z(i),r(i),s(i)
     +   ,2
     +   ,Az(k),Bz(k),Cz(k)
     +   ,Ar(k),Br(k),Cr(k)
     +   ,h2(k)
     +   )

	      if(k.eq.1) then
               F(i,k)=f1(k)
	         H(i,k)=h1(k)
	      else 
	         F(i,k)=f2(k-1)+f1(k)
	         H(i,k)=h2(k-1)+h1(k)
	       end if
	  end do
	   F(i,N+1)=f2(N)
	   H(i,N+1)=h2(N)

      end do
      
	do j=1,N+1
	   H(j,j)=H(j,j)-0.5D0
	 end do

	do j=1,N+1
	  tmp=0.0D0
	  do k=1,N+1
	    tmp =tmp + H(j,k)*z(k)
	  end do
         RHS(j)=tmp
	end do

         call dgesv(N+1, 1, F, neq, ipiv, RHS, neq, info)
      y(1:N+1)=RHS(1:N+1)

      return 
	end
c-----------------------------------------------------------------------------
      subroutine int_pot_s
     +
     +   (z1,r1
     +   ,s1,s2
     +   ,z0,r0,s0
     +   ,n
     +   ,Az,Bz,Cz
     +   ,Ar,Br,Cr
     +   ,y1
     +   )
      Implicit Double Precision (a-h,o-z)

	dimension x1(10),w1(10),x2(10),w2(10),amg(10),ams(10)
	dimension anzg(10),anrg(10),dzg(10),drg(10),dsg(10)
      dimension Bg(10),Bs(10)   
      dimension sg(10),ss(10),rg(10),rs(10),zg(10),zs(10)
      dimension dss(10),drs(10),dzs(10)
	dimension xh2(10),rh2(10),akg(10),aks(10)
	dimension pf1(10),pf2(10),pe1(10),pe2(10)
      dimension G(10),Fg(10),Eg(10)
	dimension Gz(10),Gr(10)
	   
	eps = 0.000000001D0
      
	pi  = 3.1415 92653 58979 32384 D0
      
	call Gauss_Legendre (8,x1,w1)
	call LogGL (8,x2,w2)

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
c	   anzg(j) = (anx2-anx1)/2.0D0*x1(j) + (anx2+anx1)/2.0D0
c	   anrg(j) = (any2-any1)/2.0D0*x1(j) + (any2+any1)/2.0D0
	   xd = sg(j) - s1
	   rg(j) = Ar*xd**3 + Br*xd**2 + Cr*xd + r1
	   drg(j) = 3.D0*Ar*xd**2 + 2.D0*Br*xd + Cr
	   zg(j) = Az*xd**3 + Bz*xd**2 + Cz*xd + z1
	   dzg(j) = 3.D0*Az*xd**2 + 2.D0*Bz*xd + Cz
	   dsg(j) = dsqrt(drg(j)**2+dzg(j)**2)
	   anzg(j) = -drg(j)/dsg(j)
	   anrg(j) = dzg(j)/dsg(j)
	   akg(j) = 4.0D0*r0*rg(j)/((zg(j)-z0)**2+(rg(j)+r0)**2)
	end do
 
      Iopt = 1  
c ipot ne 1    both green's function and its derivative

	do j=1,8
	    call B0(rg(j),zg(j),r0,z0,Bg(j))
	    call ell_int2(akg(j),Fg(j),Eg(j))
	    call polyn (akg(j),pf1(j),pf2(j),pe1(j),pe2(j))
          call lgf_ax_fs(Iopt,zg(j),rg(j),z0,r0,G(j),Gz(j),Gr(j))
      end do

c********
c first   single layer integrals
c -------
      if (dabs(s1-s0).lt.eps .and. r0.gt.eps) then 
c      if (s1.eq.s0) then 
	   xx_tmp1 = 0.0 D0
c	   xx_tmp2 = 0.0 D0
	   xx2 = 0.0D0
	   xx3 = 0.0D0
         
	  do j=1,8
	     Fg(j) = pf1(j) - dlog((1.D0-akg(j))/(sg(j)-s0)**2)*pf2(j)
        end do	   

c  regular part        
	  do j=1,8
	     xx_tmp1 = xx_tmp1 + a*amg(j)*Bg(j)*Fg(j)*dsg(j)*rg(j)*w1(j)
	     xx2 = xx2 + amg(j)*Bg(j)*dsg(j)*rg(j)*pf2(j)*w1(j)
	  end do

        xx2 = -xx2*dlog(s2-s1)*a

c  singular part
        do j=1,8
	     ss(j) = (s2-s1)*x2(j) + s1
	     if (n.eq.1)then
		      ams(j) = 1.D0 - x2(j)
	     else if (n.eq.2)then
	          ams(j) = x2(j)
	     end if
	     xd = ss(j)-s1
     	     rs(j) = Ar*xd**3 + Br*xd**2 + Cr*xd + r1
	     zs(j) = Az*xd**3 + Bz*xd**2 + Cz*xd + z1
	     drs(j) = 3.D0*Ar*xd**2 + 2.D0*Br*xd + Cr
	     dzs(j) = 3.D0*Az*xd**2 + 2.D0*Bz*xd + Cz
	     dss(j) = dsqrt(drs(j)**2+dzs(j)**2)
	     xh2(j) = (zs(j)-z0)**2
	     rh2(j) = xh2(j) + (rs(j)-r0)**2
	     Rd = dsqrt(xh2(j) + (rs(j)+r0)**2)
	     aks(j) = 4.0D0*r0*rs(j)/Rd**2
           call B0(rs(j),zs(j),r0,z0,Bs(j))
c	  end do
c	  do j=1,8
           call polyn (aks(j),Pff1,Pff2,Pee1,Pee2)
c	     R = sqrt(xh2(j) + (rs(j) + r0)**2)

	     xx3 = xx3 + ams(j)*Bs(j)*dss(j)*rs(j)*Pff2*w2(j)

	  end do
        xx3 = xx3*(s2-s1)

c       if (z0.lt.eps) then
c	     y1 = 0.0D0
c	   else
	     y1 = xx_tmp1 + 2.0D0*(xx2 + xx3)
c	  end if

	else if (dabs(s2-s0).lt.eps .and. r0.gt.eps) then
	   xx_tmp1 = 0.0 D0
	   xx2 = 0.0D0
	   xx3 = 0.0D0
         
	  do j=1,8
	     Fg(j) = pf1(j) - dlog((1.D0-akg(j))/(sg(j)-s0)**2)*pf2(j)
        end do	   

c  regular part        
	  do j=1,8
	     xx_tmp1 = xx_tmp1 + a*amg(j)*Bg(j)*Fg(j)*dsg(j)*rg(j)*w1(j)
	     xx2 = xx2 + amg(j)*Bg(j)*dsg(j)*rg(j)*pf2(j)*w1(j)
	  end do

        xx2 = -xx2*dlog(s2-s1)*a

c  singular part
        do j=1,8
	     ss(j) = (s1-s2)*x2(j) + s2
	     if (n.eq.1)then
	      ams(j) = x2(j)
	     else if (n.eq.2)then
	      ams(j) = 1.D0-x2(j)
	     end if
	     xd = ss(j)-s1
     	     rs(j) = Ar*xd**3 + Br*xd**2 + Cr*xd + r1
	     zs(j) = Az*xd**3 + Bz*xd**2 + Cz*xd + z1
	     drs(j) = 3.D0*Ar*xd**2 + 2.D0*Br*xd + Cr
	     dzs(j) = 3.D0*Az*xd**2 + 2.D0*Bz*xd + Cz
	     dss(j) = dsqrt(drs(j)**2+dzs(j)**2)
	     xh2(j) = (zs(j)-z0)**2
	     rh2(j) = xh2(j) + (rs(j)-r0)**2
	     Rd = dsqrt(xh2(j) + (rs(j)+r0)**2)
	     aks(j) = 4.0D0*r0*rs(j)/Rd**2
           call B0(rs(j),zs(j),r0,z0,Bs(j))
c	  end do
c	  do j=1,8
           call polyn (aks(j),Pff1,Pff2,Pee1,Pee2)
c	     R = sqrt(xh2(j) + (rs(j) + r0)**2)

	     xx3 = xx3 + ams(j)*Bs(j)*dss(j)*rs(j)*Pff2*w2(j)

	  end do
        xx3 = xx3*(s2-s1)
 
 	     y1 = xx_tmp1 + 2.0D0*(xx2 + xx3)
	else
	   xx_tmp1 = 0.0D0

	  do j=1,8
	     xx_tmp1=xx_tmp1+a*amg(j)*Bg(j)*Fg(j)*dsg(j)*rg(j)*w1(j)
	  end do

	  y1 = xx_tmp1

	end if
c
c ---- done
	     
	return
	end
c----------------------------------------------------------------------------
      subroutine int_pot_d2
     +
     +   (z1,r1
     +   ,s1,s2
     +   ,z0,r0,s0
     +   ,n
     +   ,Az,Bz,Cz
     +   ,Ar,Br,Cr
     +   ,y1
     +   )
      Implicit Double Precision (a-h,o-z)

	dimension x1(10),w1(10),x2(10),w2(10),amg(10),ams(10)
	dimension anzg(10),anrg(10),dzg(10),drg(10),dsg(10)
      dimension Bg(10),Bs(10),Cg(10),Cs(10)
      dimension sg(10),ss(10),rg(10),rs(10),zg(10),zs(10)
      dimension dss(10),drs(10),dzs(10)
	dimension xh2(10),rh2(10),akg(10),aks(10)
	dimension pf1(10),pf2(10),pe1(10),pe2(10)
      dimension G(10),Fg(10),Eg(10)
	dimension Gz(10),Gr(10)
	   
	eps = 0.000000001D0
      
	pi  = 3.1415 92653 58979 32384 D0
      
	call Gauss_Legendre (8,x1,w1)
	call LogGL (8,x2,w2)

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
c	   anzg(j) = (anx2-anx1)/2.0D0*x1(j) + (anx2+anx1)/2.0D0
c	   anrg(j) = (any2-any1)/2.0D0*x1(j) + (any2+any1)/2.0D0
	   xd = sg(j) - s1
	   rg(j) = Ar*xd**3 + Br*xd**2 + Cr*xd + r1
	   drg(j) = 3.D0*Ar*xd**2 + 2.D0*Br*xd + Cr
	   zg(j) = Az*xd**3 + Bz*xd**2 + Cz*xd + z1
	   dzg(j) = 3.D0*Az*xd**2 + 2.D0*Bz*xd + Cz
	   dsg(j) = dsqrt(drg(j)**2+dzg(j)**2)
	   anzg(j) = -drg(j)/dsg(j)
	   anrg(j) = dzg(j)/dsg(j)
	   akg(j) = 4.0D0*r0*rg(j)/((zg(j)-z0)**2+(rg(j)+r0)**2)
c	   akgi(j) = 4.0D0*r0*rg(j)/((zg(j)+z0)**2+(rg(j)+r0)**2)
	end do
 
      Iopt = 2  
c ipot ne 1    both green's function and its derivative

	do j=1,8
	    call B0(rg(j),zg(j),r0,z0,Bg(j))
	    call C1(rg(j),drg(j),zg(j),dzg(j),r0,z0,Cg(j))
	    call ell_int2(akg(j),Fg(j),Eg(j))
c	    call ell_int(akgi(j),Fgi(j),Egi(j))
	    call polyn (akg(j),pf1(j),pf2(j),pe1(j),pe2(j))
          call lgf_ax_fs(Iopt,zg(j),rg(j),z0,r0,G(j),Gz(j),Gr(j))
c          call lgf_ax_fs(Iopt,zg(j),rg(j),-z0,r0,Gi(j),Gzi(j),Gri(j))
      end do

c********
c first   single layer integrals
c -------
      if (dabs(s1-s0).lt.eps .and. r0.gt.eps) then 
c      if (s1.eq.s0) then 
	   xx_tmp1 = 0.0 D0
c	   xx_tmp2 = 0.0 D0
	   xx2 = 0.0D0
	   xx3 = 0.0D0
         
	  do j=1,8
	     Fg(j) = pf1(j) - dlog((1.D0-akg(j))/(sg(j)-s0)**2)*pf2(j)
	     Eg(j) = pe1(j) - dlog((1.D0-akg(j))/(sg(j)-s0)**2)*pe2(j)
        end do	   

c  regular part        
	  do j=1,8
	     xx_tmp1=xx_tmp1+a*amg(j)*Bg(j)*(-dzg(j)*Fg(j)+
     +	    Cg(j)*Eg(j))*w1(j)
c	     xx_tmp2 = xx_tmp2 + a*amg(j)*Bgi(j)*Fgi(j)*dsg(j)*rg(j)*w1(j)
	     xx2 = xx2 + amg(j)*Bg(j)*(-dzg(j)*pf2(j)+
     +		 Cg(j)*pe2(j))*w1(j)
	  end do

        xx2 = -xx2*dlog(s2-s1)*a

c  singular part
        do j=1,8
	     ss(j) = (s2-s1)*x2(j) + s1
	    if (n.eq.1)then
	     ams(j) = 1.D0 - x2(j)
	    else if (n.eq.2)then
	     ams(j) = x2(j)
	    end if
	     xd = ss(j)-s1
     	     rs(j) = Ar*xd**3 + Br*xd**2 + Cr*xd + r1
	     zs(j) = Az*xd**3 + Bz*xd**2 + Cz*xd + z1
	     drs(j) = 3.D0*Ar*xd**2 + 2.D0*Br*xd + Cr
	     dzs(j) = 3.D0*Az*xd**2 + 2.D0*Bz*xd + Cz
	     dss(j) = dsqrt(drs(j)**2+dzs(j)**2)
	     xh2(j) = (zs(j)-z0)**2
	     rh2(j) = xh2(j) + (rs(j)-r0)**2
	     Rd = dsqrt(xh2(j) + (rs(j)+r0)**2)
	     aks(j) = 4.0D0*r0*rs(j)/Rd**2
           call B0(rs(j),zs(j),r0,z0,Bs(j))
	     call C1(rs(j),drs(j),zs(j),dzs(j),r0,z0,Cs(j))
           call polyn (aks(j),Pff1,Pff2,Pee1,Pee2)
c	     R = sqrt(xh2(j) + (rs(j) + r0)**2)

	     xx3 = xx3 + ams(j)*Bs(j)*(-dzs(j)*Pff2+
     +        Cs(j)*Pee2)*w2(j)

	  end do
        xx3 = xx3*(s2-s1)

c       if (z0.lt.eps) then
c	     y1 = 0.0D0
c	   else
	     y1 = xx_tmp1 + 2.0D0*(xx2 + xx3)
c	  end if

	else if (dabs(s2-s0).lt.eps .and. r0.gt.eps) then
c	else if (s2.eq.s0) then
	   xx_tmp1 = 0.0 D0
c	   xx_tmp2 = 0.0D0
	   xx2 = 0.0D0
	   xx3 = 0.0D0
         
	  do j=1,8
	     Fg(j) = pf1(j) - dlog((1.D0-akg(j))/(sg(j)-s0)**2)*pf2(j)
	     Eg(j) = pe1(j) - dlog((1.D0-akg(j))/(sg(j)-s0)**2)*pe2(j)
        end do	   

c  regular part        
	  do j=1,8
	     xx_tmp1 = xx_tmp1 + a*amg(j)*Bg(j)*(-dzg(j)*Fg(j)+
     +       Cg(j)*Eg(j))*w1(j)
c	     xx_tmp2 = xx_tmp2 + a*amg(j)*Bgi(j)*Fgi(j)*dsg(j)*rg(j)*w1(j)
	     xx2 = xx2 + amg(j)*Bg(j)*(-dzg(j)*pf2(j)+
     +		 Cg(j)*pe2(j))*w1(j)
	  end do

        xx2 = -xx2*dlog(s2-s1)*a

c  singular part
        do j=1,8
	     ss(j) = (s1-s2)*x2(j) + s2
	     if (n.eq.1)then
	     ams(j) = x2(j)
	     else if (n.eq.2)then
	     ams(j) = 1.0D0-x2(j)
	     end if
	     xd = ss(j)-s1
     	     rs(j) = Ar*xd**3 + Br*xd**2 + Cr*xd + r1
	     zs(j) = Az*xd**3 + Bz*xd**2 + Cz*xd + z1
	     drs(j) = 3.D0*Ar*xd**2 + 2.D0*Br*xd + Cr
	     dzs(j) = 3.D0*Az*xd**2 + 2.D0*Bz*xd + Cz
	     dss(j) = dsqrt(drs(j)**2+dzs(j)**2)
	     xh2(j) = (zs(j)-z0)**2
	     rh2(j) = xh2(j) + (rs(j)-r0)**2
	     Rd = dsqrt(xh2(j) + (rs(j)+r0)**2)
	     aks(j) = 4.0D0*r0*rs(j)/Rd**2
           call B0(rs(j),zs(j),r0,z0,Bs(j))
	     call C1(rs(j),drs(j),zs(j),dzs(j),r0,z0,Cs(j))
           call polyn (aks(j),Pff1,Pff2,Pee1,Pee2)
c	     R = sqrt(xh2(j) + (rs(j) + r0)**2)

	     xx3 = xx3 + ams(j)*Bs(j)*(-dzs(j)*Pff2+
     +        Cs(j)*Pee2)*w2(j)

	  end do
        xx3 = xx3*(s2-s1)
 
 	     y1 = xx_tmp1 + 2.0D0*(xx2 + xx3)
c 	     y1 = xx_tmp1 + 2.0D0*(xx2 + xx3) + xx_tmp2
	else
	   xx_tmp1 = 0.0D0

	  do j=1,8
	     xx_tmp1=xx_tmp1+a*amg(j)*Bg(j)*(-dzg(j)*Fg(j)+
     +         Cg(j)*Eg(j))*w1(j)
	  end do

	  y1 = xx_tmp1

	end if
c
        y1 = 0.5D0*y1
c ---- done
	     
	return
	end
c------------------------------------------------------------------
      subroutine int_pot_d
     +
     +   (z1,r1
     +   ,s1,s2
     +   ,z0,r0,n
     +   ,Az,Bz,Cz
     +   ,Ar,Br,Cr
     +   ,y
     +   )

      Implicit Double Precision (a-h,o-z)

	dimension x1(10),w1(10)
	dimension anzg(10),anrg(10),dzg(10),drg(10),dsg(10)
      dimension sg(10),rg(10),zg(10),phig(10),amg(10)
      dimension Gi(10),G(10),Gz(10),Gr(10),Gzi(10),Gri(10)
      dimension Gci(10),Gc(10),Gzc(10),Grc(10),Gzci(10),Grci(10)
	   
	eps = 0.000000001D0
      pi  = 3.1415 92653 58979 32384 D0

      Iopt = 2  
c both green's function and its derivative

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
	   xd = sg(j) - s1
	   rg(j) = Ar*xd**3 + Br*xd**2 + Cr*xd + r1
	   drg(j) = 3.D0*Ar*xd**2 + 2.D0*Br*xd + Cr
	   zg(j) = Az*xd**3 + Bz*xd**2 + Cz*xd + z1
	   dzg(j) = 3.D0*Az*xd**2 + 2.D0*Bz*xd + Cz
	   dsg(j) = dsqrt(drg(j)**2+dzg(j)**2)
	   anzg(j) = -drg(j)/dsg(j)
	   anrg(j) = dzg(j)/dsg(j)
c	end do
 
c       Iopt = 2
c  periodic Green's function
           call lgf_ax_fs
     +   (Iopt
     +   ,zg(j),rg(j)
     +   ,z0,r0
     +   ,G(j)
     +   ,Gz(j),Gr(j)
     +   )
      end do


c  double layer integrals
c********
      yy_tmp1=0.0D0
c	yy_tmp2=0.0D0
c  nx = -dr/ds  nr = dz/ds	  
	do j=1,8
	    Gn=Gz(j)*anzg(j)+Gr(j)*anrg(j)
c	    Gni=Gzi(j)*anzg(j)+Gri(j)*anrg(j)
c	    Gn=(Gz(j)+Gzc(j))*anzg(j)+(Gr(j)+Grc(j))*anrg(j)
c	    Gni=(Gzi(j)+Gzci(j))*anzg(j)+(Gri(j)+Grci(j))*anrg(j)
	    yy_tmp1=yy_tmp1+a*amg(j)*Gn*dsg(j)*rg(j)*w1(j)
c	    yy_tmp1=yy_tmp1+a*(phig(j)-phi0)*(Gn+Gni)*dsg(j)*rg(j)*w1(j)
c	    yy_tmp2=yy_tmp2+a*(phig(j)-phi0)*Gni*dsg(j)*rg(j)*w1(j)
	end do
         
c	y = yy_tmp1 + yy_tmp2
      y = yy_tmp1
c ---- done
	     
	return
	end
c---------------------------------------------------------------
	subroutine C1(r,dr,z,dz,r0,z0,y)
      Implicit Double Precision (a-h,o-z)
      pi  = 3.1415 92653 58979 32384 D0
      xh = z-z0
	xh2 = xh*xh
	rh2 = xh2 + (r-r0)**2
	y = (dz*(xh2 + r0**2-r**2)-2.0D0*r*dr*(z0-z))/rh2

	return 
	end
c------------------------------------------------------------------
	subroutine B0(r,z,r0,z0,y)
      Implicit Double Precision (a-h,o-z)
      pi  = 3.1415 92653 58979 32384 D0
	rh = (z-z0)**2 + (r+r0)**2
	y = 1.D0/(pi*dsqrt(rh))
	return 
	end
