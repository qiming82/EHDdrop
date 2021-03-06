      subroutine bess_I01K01
     +
     +   (x
     +   ,Iswitch
     +   ,BI0,BI1
     +   ,BK0,BK1
     +   )

c----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c----------------------------------------

c----------------------------------------
c  Modifed Bessel functions
c  of zeroth and first order
c  computed by approximation.
c
c  Iswitch = 1 to get BI0, BI1
c            2 to get BI0, BI1, BK0, BK1
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

c--------
c prepare
c--------

      xr = Dsqrt(x)
      ex =  Dexp(x)

      t  = x/3.75D0
      t2 = t*t
      ti = 1.0D0/t

c---------------------------
      If(abs(t).lt.1.0) then
c---------------------------

      BI0 = 1.0        + t2 * ( 3.5156229  +
     +                   t2 * ( 3.0899424  +
     +                   t2 * ( 1.2067492  +
     +                   t2 * ( 0.2659732  +
     +                   t2 * ( 0.0360768  +
     +                   t2 *   0.0045813 )))))
      BI1 = 0.5        + t2 * ( 0.87890594  +
     +                   t2 * ( 0.51498869  +
     +                   t2 * ( 0.15084934  +
     +                   t2 * ( 0.02658733  +
     +                   t2 * ( 0.00301532  +
     +                   t2 *   0.00032411  )))))
      BI1 = BI1 * x

c---------
      Else
c---------

      BI0 = 0.39894228 + ti * ( 0.01328592 +
     +                   ti * ( 0.00225319 +
     +                   ti * (-0.00157565 +
     +                   ti * ( 0.00916281 +
     +                   ti * (-0.02057706 +
     +                   ti * ( 0.02635537 +
     +                   ti * (-0.01647633 +
     +                   ti *   0.00392377 )))))))
      BI0 = BI0 * ex/xr
      BI1 = 0.39894228 + ti * (-0.03988024 +
     +                   ti * (-0.00362018 +
     +                   ti * ( 0.00163801 +
     +                   ti * (-0.01031555 +
     +                   ti * ( 0.02282967 +
     +                   ti * (-0.02895312 +
     +                   ti * ( 0.01787654 -
     +                   ti *   0.00420059 )))))))
      BI1 = BI1 * ex/xr

c-----------
      End If
c-----------

c------------------------------
      If(Iswitch.eq.1) Go to 99
c------------------------------

      z  = 0.50*x
      z2 = z*z
      zi = 1.0D0/z

c---
      If(abs(z).lt.1.0) then
c---

      BK0 =-0.57721566 + z2 * ( 0.42278420  +
     +                   z2 * ( 0.23069756  +
     +                   z2 * ( 0.03488590  +
     +                   z2 * ( 0.00262698  +
     +                   z2 * ( 0.00010750  +
     +                   z2 *   0.00000740 )))))
      BK0 = BK0 - Dlog(z)*BI0
      BK1 = 1.0        + z2 * ( 0.15443144  +
     +                   z2 * (-0.67278579  +
     +                   z2 * (-0.18156897  +
     +                   z2 * (-0.01919402  +
     +                   z2 * (-0.00110404  -
     +                   z2 *   0.00004686 )))))
      BK1 = ( x*log(z)*BI1+BK1 )/x

c---
      Else
c---

      BK0 = 1.25331414 + zi * (-0.07832358 +
     +                   zi * ( 0.02189568 +
     +                   zi * (-0.01062446 +
     +                   zi * ( 0.00587872 +
     +                   zi * (-0.00251540 +
     +                   zi *   0.00053208 )))))
      BK0 = BK0 / (xr*ex)
      BK1 = 1.25331414 + zi * ( 0.23498619 +
     +                   zi * (-0.03655620 +
     +                   zi * ( 0.01504268 +
     +                   zi * (-0.00780353 +
     +                   zi * ( 0.00325614 -
     +                   zi *   0.00068245 )))))
      BK1 = BK1 / (xr*ex)

c-----------
      End If
c-----------

  99  Continue

      Return
      End
