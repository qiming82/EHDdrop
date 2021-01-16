      subroutine polyn (ak2,pf1,pf2,pe1,pe2)

      Implicit Double Precision (a-h,o-z)

	dimension a(10), b(10), c(10), d(10)

      a(1)=0.09666344259 D0
      a(2)=0.03590092383 D0
      a(3)=0.03742563713 D0
      a(4)=0.01451196212 D0

      b(1)=0.12498593597 D0
      b(2)=0.06880248576 D0
      b(3)=0.03328355346 D0
      b(4)=0.00441787012 D0
	
	c(1)=0.44325141463 D0
      c(2)=0.06260601220 D0
      c(3)=0.04757383546 D0
      c(4)=0.01736506451 D0

      d(1)=0.24998368310 D0
      d(2)=0.09200180037 D0
      d(3)=0.04069697526 D0
      d(4)=0.00526449639 D0

	pf1 = 0.0 D0
	pf2 = 0.0 D0
	pe1 = 0.0 D0
	pe2 = 0.0 D0

	s1 = 0.0 D0
	s2 = 0.0 D0
	s3 = 0.0 D0
	s4 = 0.0 D0

      do i=1,4
	   s1 = s1 + a(i)*(1.0D0-ak2)**i
	   s2 = s2 + b(i)*(1.0D0-ak2)**i
         s3 = s3 + c(i)*(1.0D0-ak2)**i
	   s4 = s4 + d(i)*(1.0D0-ak2)**i
      end do
	pf1 = dlog(4.0D0) + s1
	pf2 = 0.5D0 + s2 
      pe1 = 1.0D0 + s3
	pe2 = s4
      
	return
	end
