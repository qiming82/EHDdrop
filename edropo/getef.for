      subroutine getef
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
c	dimension h1(500),h2(500),f1(500),f2(500),ml(500,500), mu(500,500)
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
c           call int_pot_d
c     +
c     +   (z(k),r(k)
c     +   ,s(k),s(k+1)
c     +   ,z(i),r(i)
c     +   ,1
c     +   ,Az(k),Bz(k),Cz(k)
c     +   ,Ar(k),Br(k),Cr(k)
c     +   ,h1(k)
c     +   )
c           call int_pot_d
c     +
c     +   (z(k),r(k)
c     +   ,s(k),s(k+1)
c     +   ,z(i),r(i)
c     +   ,2
c     +   ,Az(k),Bz(k),Cz(k)
c     +   ,Ar(k),Br(k),Cr(k)
c     +   ,h2(k)
c     +   )

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
