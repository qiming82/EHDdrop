      subroutine getef22
     +
     +   (z1,r1,s1
     +   ,z2,r2,s2
     +   ,anz1,anz2
     +   ,Az1,Bz1,Cz1
     +   ,Ar1,Br1,Cr1
     +   ,Az2,Bz2,Cz2
     +   ,Ar2,Br2,Cr2
     +   ,N,RR
     +   ,phis1,dpdt1,dpdn1
     +   ,phis2,dpdt2,dpdn2
     +   )

      Implicit Double Precision (a-h,o-z)

      parameter (neq=1024)
	dimension z1(neq),r1(neq),s1(neq),phih1(neq),anz1(neq),anz2(neq)
        dimension z2(neq),r2(neq),s2(neq),phih2(neq)
	dimension RHS(neq),x(neq),y(neq),Am(neq,neq)
	dimension h1(neq),h2(neq),f1(neq),f2(neq),RHS1(neq),RHS2(neq)
	dimension hh1(neq),hh2(neq),ff1(neq),ff2(neq)
        dimension F11(neq,neq),H11(neq,neq),F12(neq,neq),H12(neq,neq)
        dimension F21(neq,neq),H21(neq,neq),F22(neq,neq),H22(neq,neq)
	dimension Az1(neq),Bz1(neq),Cz1(neq),Ar1(neq),Br1(neq),Cr1(neq)
	dimension Az2(neq),Bz2(neq),Cz2(neq),Ar2(neq),Br2(neq),Cr2(neq)
	dimension phis1(neq),dpdt1(neq),dpdn1(neq),fs1(neq),fs2(neq)       
	dimension phis2(neq),dpdt2(neq),dpdn2(neq)
      dimension ipiv(neq)


      do i=1,N+1
	  do k=1,N
c-- DOUBLE LAYER: x0 on S1, over S1       
            call int_pot_d2
     +
     +   (z1(k),r1(k)
     +   ,s1(k),s1(k+1)
     +   ,z1(i),r1(i),s1(i),1
     +   ,Az1(k),Bz1(k),Cz1(k)
     +   ,Ar1(k),Br1(k),Cr1(k)
     +   ,h1(k)
     +   )           
           call int_pot_d2
     +
     +   (z1(k),r1(k)
     +   ,s1(k),s1(k+1)
     +   ,z1(i),r1(i),s1(i),2
     +   ,Az1(k),Bz1(k),Cz1(k)
     +   ,Ar1(k),Br1(k),Cr1(k)
     +   ,h2(k)
     +   )           
c-- DOUBLE LAYER: x0 on S1, over S2       
           call int_pot_d2
     +
     +   (z2(k),r2(k)
     +   ,s2(k),s2(k+1)
     +   ,z1(i),r1(i),s1(i),1
     +   ,Az2(k),Bz2(k),Cz2(k)
     +   ,Ar2(k),Br2(k),Cr2(k)
     +   ,hh1(k)
     +   )           
           call int_pot_d2
     +
     +   (z2(k),r2(k)
     +   ,s2(k),s2(k+1)
     +   ,z1(i),r1(i),s1(i),2
     +   ,Az2(k),Bz2(k),Cz2(k)
     +   ,Ar2(k),Br2(k),Cr2(k)
     +   ,hh2(k)
     +   )           
c-
c-- SINGLE LAYER: x0 on S1, over S1
           call int_pot_s
     +
     +   (z1(k),r1(k)
     +   ,s1(k),s1(k+1)
     +   ,z1(i),r1(i),s1(i)
     +   ,1
     +   ,Az1(k),Bz1(k),Cz1(k)
     +   ,Ar1(k),Br1(k),Cr1(k)
     +   ,f1(k)
     +   )

           call int_pot_s
     +
     +   (z1(k),r1(k)
     +   ,s1(k),s1(k+1)
     +   ,z1(i),r1(i),s1(i)
     +   ,2
     +   ,Az1(k),Bz1(k),Cz1(k)
     +   ,Ar1(k),Br1(k),Cr1(k)
     +   ,f2(k)
     +   )
c-- SINGLE LAYER: x0 on S1, over S2 (regular...)
           call int_pot_s
     +
     +   (z2(k),r2(k)
     +   ,s2(k),s2(k+1)
     +   ,z1(i),r1(i),10.0d0   !s2(i)
     +   ,1
     +   ,Az2(k),Bz2(k),Cz2(k)
     +   ,Ar2(k),Br2(k),Cr2(k)
     +   ,ff1(k)
     +   )

           call int_pot_s
     +
     +   (z2(k),r2(k)
     +   ,s2(k),s2(k+1)
     +   ,z1(i),r1(i),10.0d0  !s1(i)
     +   ,2
     +   ,Az2(k),Bz2(k),Cz2(k)
     +   ,Ar2(k),Br2(k),Cr2(k)
     +   ,ff2(k)
     +   )

	      if(k.eq.1) then
               H11(i,k)=h1(k)
               H12(i,k)=hh1(k)
               F11(i,k)=f1(k)
               F12(i,k)=ff1(k)
	      else 
	         H11(i,k)=h2(k-1)+h1(k)    ! double layer
	         H12(i,k)=hh2(k-1)+hh1(k)    ! double layer
	         F11(i,k)=f2(k-1)+f1(k)  ! single layer
	         F12(i,k)=ff2(k-1)+ff1(k)  ! single layer
	       end if
	 enddo
	   
	   H11(i,N+1)=h2(N)
           H12(i,N+1)=hh2(N)
	   F11(i,N+1)=f2(N)
	   F12(i,N+1)=ff2(N)
      end do

c----
      do i=1,N+1
	  do k=1,N
c-- DOUBLE LAYER: x0 on S2, over S1       
           call int_pot_d2
     +
     +   (z1(k),r1(k)
     +   ,s1(k),s1(k+1)
     +   ,z2(i),r2(i),s2(i),1
     +   ,Az1(k),Bz1(k),Cz1(k)
     +   ,Ar1(k),Br1(k),Cr1(k)
     +   ,h1(k)
     +   )           
           call int_pot_d2
     +
     +   (z1(k),r1(k)
     +   ,s1(k),s1(k+1)
     +   ,z2(i),r2(i),s2(i),2
     +   ,Az1(k),Bz1(k),Cz1(k)
     +   ,Ar1(k),Br1(k),Cr1(k)
     +   ,h2(k)
     +   )           
c-- DOUBLE LAYER: x0 on S2, over S2       
           call int_pot_d2
     +
     +   (z2(k),r2(k)
     +   ,s2(k),s2(k+1)
     +   ,z2(i),r2(i),s2(i),1
     +   ,Az2(k),Bz2(k),Cz2(k)
     +   ,Ar2(k),Br2(k),Cr2(k)
     +   ,hh1(k)
     +   )           
           call int_pot_d2
     +
     +   (z2(k),r2(k)
     +   ,s2(k),s2(k+1)
     +   ,z2(i),r2(i),s2(i),2
     +   ,Az2(k),Bz2(k),Cz2(k)
     +   ,Ar2(k),Br2(k),Cr2(k)
     +   ,hh2(k)
     +   )           
c-
c-- SINGLE LAYER: x0 on S2, over S1 (regurlar...)
           call int_pot_s
     +
     +   (z1(k),r1(k)
     +   ,s1(k),s1(k+1)
     +   ,z2(i),r2(i),10.0d0 !s1(i)
     +   ,1
     +   ,Az1(k),Bz1(k),Cz1(k)
     +   ,Ar1(k),Br1(k),Cr1(k)
     +   ,f1(k)
     +   )

           call int_pot_s
     +
     +   (z1(k),r1(k)
     +   ,s1(k),s1(k+1)
     +   ,z2(i),r2(i),10.0d0 !s1(i)
     +   ,2
     +   ,Az1(k),Bz1(k),Cz1(k)
     +   ,Ar1(k),Br1(k),Cr1(k)
     +   ,f2(k)
     +   )
c-- SINGLE LAYER: x0 on S2, over S2
           call int_pot_s
     +
     +   (z2(k),r2(k)
     +   ,s2(k),s2(k+1)
     +   ,z2(i),r2(i),s2(i)
     +   ,1
     +   ,Az2(k),Bz2(k),Cz2(k)
     +   ,Ar2(k),Br2(k),Cr2(k)
     +   ,ff1(k)
     +   )

           call int_pot_s
     +
     +   (z2(k),r2(k)
     +   ,s2(k),s2(k+1)
     +   ,z2(i),r2(i),s2(i)
     +   ,2
     +   ,Az2(k),Bz2(k),Cz2(k)
     +   ,Ar2(k),Br2(k),Cr2(k)
     +   ,ff2(k)
     +   )
	      if(k.eq.1) then
               H21(i,k)=h1(k)
               H22(i,k)=hh1(k)
               F21(i,k)=f1(k)
               F22(i,k)=ff1(k)
	      else 
	         H21(i,k)=h2(k-1)+h1(k)    ! double layer
	         H22(i,k)=hh2(k-1)+hh1(k)    ! double layer
	         F21(i,k)=f2(k-1)+f1(k)  ! single layer
	         F22(i,k)=ff2(k-1)+ff1(k)  ! single layer
	       end if
	 enddo	   
	   H21(i,N+1)=h2(N)
           H22(i,N+1)=hh2(N)
	   F21(i,N+1)=f2(N)
	   F22(i,N+1)=ff2(N)
      end do

c----
	 H11 = H11*2.0d0*(1.0-RR)/(1.0+RR)
	 H12 = H12*2.0d0*(1.0-RR)/(1.0+RR)
	 H21 = H21*2.0d0*(1.0-RR)/(1.0+RR)
	 H22 = H22*2.0d0*(1.0-RR)/(1.0+RR)

	do i=1,N+1
	  H11(i,i)=1.0d0 + H11(i,i)
	  H22(i,i)=1.0d0 + H22(i,i)
	enddo
c--
       do i=1,N+1
         do k=1,N+1
           Am(i,k)=H11(i,k)
           Am(i,k+N+1)=H12(i,k)
           Am(i+N+1,k)=H21(i,k)
           Am(i+N+1,k+N+1)=H22(i,k)
         enddo
       enddo
c--
	do j=1,N+1
	  tmp1=0.0d0
          tmp2=0.0d0
	  do k=1,N+1
	    tmp1=tmp1+F11(j,k)*anz1(k)+F12(j,k)*anz2(k) 
	    tmp2=tmp2+F21(j,k)*anz1(k)+F22(j,k)*anz2(k) 
	  enddo
         RHS(j) = -2.0*(RR-1.0d0)/(RR+1.0d0)*tmp1
         RHS(j+N+1) = -2.0*(RR-1.0d0)/(RR+1.0d0)*tmp2
	enddo

         call dgesv(2*N+2, 1, Am, neq, ipiv, RHS, neq, info)
        phih1(1:N+1)=RHS(1:N+1)   ! surface potential phih = phis + z
        phih2(1:N+1)=RHS(N+2:2*N+2)
c-- phi_inf = -z

        do j=1,N+1
           phis1(j)=phih1(j)-z1(j)
           phis2(j)=phih2(j)-z2(j)
        enddo
c-- tangential derivative
       do j=2,N
	  dpdt1(j)=(phis1(j+1)-phis1(j-1))/(s1(j+1)-s1(j-1))
	  dpdt2(j)=(phis2(j+1)-phis2(j-1))/(s2(j+1)-s2(j-1))
	 enddo
	dpdt1(1)=0.0d0
	dpdt1(N+1)=0.0d0
	dpdt2(1)=0.0d0
	dpdt2(N+1)=0.0d0

c---- normal derivative
c
      do i=1,N+1
	  do k=1,N
c-- single layer: x0 on s1 over s1
           call int_pot_s
     +
     +   (z1(k),r1(k)
     +   ,s1(k),s1(k+1)
     +   ,z1(i),r1(i),s1(i)
     +   ,1
     +   ,Az1(k),Bz1(k),Cz1(k)
     +   ,Ar1(k),Br1(k),Cr1(k)
     +   ,f1(k)
     +   )
           call int_pot_s
     +
     +   (z1(k),r1(k)
     +   ,s1(k),s1(k+1)
     +   ,z1(i),r1(i),s1(i)
     +   ,2
     +   ,Az1(k),Bz1(k),Cz1(k)
     +   ,Ar1(k),Br1(k),Cr1(k)
     +   ,f2(k)
     +   )
c-- single layer: x0 on s1 over s2
           call int_pot_s
     +
     +   (z2(k),r2(k)
     +   ,s2(k),s2(k+1)
     +   ,z1(i),r1(i),10.0d0 !s1(i)
     +   ,1
     +   ,Az2(k),Bz2(k),Cz2(k)
     +   ,Ar2(k),Br2(k),Cr2(k)
     +   ,ff1(k)
     +   )
           call int_pot_s
     +
     +   (z2(k),r2(k)
     +   ,s2(k),s2(k+1)
     +   ,z1(i),r1(i),10.0d0 !s1(i)
     +   ,2
     +   ,Az2(k),Bz2(k),Cz2(k)
     +   ,Ar2(k),Br2(k),Cr2(k)
     +   ,ff2(k)
     +   )

	      if(k.eq.1) then
               F11(i,k)=f1(k)
	       F12(i,k)=ff1(k)
	      else 
	         F11(i,k)=f2(k-1)+f1(k)
	         F12(i,k)=ff2(k-1)+ff1(k)
	       end if
	  end do
	   F11(i,N+1)=f2(N)
	   F12(i,N+1)=ff2(N)
      enddo
c--      
         do i=1,N+1
	  do k=1,N
c-- single layer: x0 on s2 over s1
           call int_pot_s
     +
     +   (z1(k),r1(k)
     +   ,s1(k),s1(k+1)
     +   ,z2(i),r2(i),10.0d0 !s1(i)
     +   ,1
     +   ,Az1(k),Bz1(k),Cz1(k)
     +   ,Ar1(k),Br1(k),Cr1(k)
     +   ,f1(k)
     +   )
           call int_pot_s
     +
     +   (z1(k),r1(k)
     +   ,s1(k),s1(k+1)
     +   ,z2(i),r2(i),10.0d0 !s1(i)
     +   ,2
     +   ,Az1(k),Bz1(k),Cz1(k)
     +   ,Ar1(k),Br1(k),Cr1(k)
     +   ,f2(k)
     +   )
c-- single layer: x0 on s2 over s2
           call int_pot_s
     +
     +   (z2(k),r2(k)
     +   ,s2(k),s2(k+1)
     +   ,z2(i),r2(i),s2(i)
     +   ,1
     +   ,Az2(k),Bz2(k),Cz2(k)
     +   ,Ar2(k),Br2(k),Cr2(k)
     +   ,ff1(k)
     +   )
           call int_pot_s
     +
     +   (z2(k),r2(k)
     +   ,s2(k),s2(k+1)
     +   ,z2(i),r2(i),s2(i)
     +   ,2
     +   ,Az2(k),Bz2(k),Cz2(k)
     +   ,Ar2(k),Br2(k),Cr2(k)
     +   ,ff2(k)
     +   )

	      if(k.eq.1) then
               F21(i,k)=f1(k)
	       F22(i,k)=ff1(k)
	      else 
	         F21(i,k)=f2(k-1)+f1(k)
	         F22(i,k)=ff2(k-1)+ff1(k)
	       end if
	  end do
	   F21(i,N+1)=f2(N)
	   F22(i,N+1)=ff2(N)
      enddo
c--
       do i=1,N+1
         do k=1,N+1
           Am(i,k)=F11(i,k)
           Am(i,k+N+1)=F12(i,k)
           Am(i+N+1,k)=F21(i,k)
           Am(i+N+1,k+N+1)=F22(i,k)
         enddo
       enddo
      
	do j=1,N+1
         RHS(j) = phih1(j)/(RR-1.0d0)  ! RR ne 1
         RHS(j+N+1) = phih2(j)/(RR-1.0d0)  ! RR ne 1
	end do

         call dgesv(2*N+2, 1, Am, neq, ipiv, RHS, neq, info)
      dpdn1(1:N+1)=RHS(1:N+1)   !dphidn_2
      dpdn2(1:N+1)=RHS(N+1+1:2*N+2)   !dphidn_2
c--
      Return 
      END
c-----------------------------------------------------------------------------
