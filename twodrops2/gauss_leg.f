      subroutine Gauss_Legendre (NQ,Z,W)

c-----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c
c             C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c        Oxford University Press
c------------------------------------------------

c------------------------------------------------
c Abscissas and weights for the Gauss-Legendre
c quadrature with NQ points
c
c This table contains values for
c
c  NQ = 1,2,3,4,5,6,8,12,20
c
c  Default value is 20
c------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension Z(20),W(20)

c-----
c trap
c-----
 
      If(    NQ.ne. 1
     +  .and.NQ.ne. 2
     +  .and.NQ.ne. 3
     +  .and.NQ.ne. 4
     +  .and.NQ.ne. 5
     +  .and.NQ.ne. 6
     +  .and.NQ.ne. 8
     +  .and.NQ.ne.12
     +  .and.NQ.ne.20
     +  ) then

        write (6,*)
        write (6,*) ' Gauss_Legendre:'
        write (6,*)
        write (6,*) '   Chosen number of Gaussian points'
        write (6,*) '   is not available; Will take NQ=20'

        NQ = 20

      End If
 
c--------------------
      If(NQ.eq.1) then
c--------------------

      Z(1) = 0.0D0

      W(1) = 2.0D0

c-------------------------
      Else If(NQ.eq.2) then
c-------------------------

      Z(1) = -0.57735 02691 8962576450
      Z(2) = -Z(1)

      W(1) = 1.0D0
      W(2) = 1.0D0

c-------------------------
      Else If(NQ.eq.3) then
c-------------------------

      Z(1) = -0.77459 66692 4148337703
      Z(2) =  0.0D0
      Z(3) = -Z(1)

      W(1) = 0.55555 55555 55555555555
      W(2) = 0.88888 88888 88888888888
      W(3) = 0.55555 55555 55555555555

c-------------------------
      Else If(NQ.eq.4) then
c-------------------------

      Z(1) = -0.86113 63115 9405257522
      Z(2) = -0.33998 10435 8485626480 
      Z(3) = -Z(2)
      Z(4) = -Z(1)

      W(1) = 0.34785 48451 3745385737
      W(2) = 0.65214 51548 6254614262
      W(3) = W(2)
      W(4) = W(1)

c-------------------------
      Else If(NQ.eq.5) then
c-------------------------

      Z(1) = -0.90617 98459 3866399279
      Z(2) = -0.53846 93101 0568309103
      Z(3) =  0.0D0
      Z(4) = -Z(2)
      Z(5) = -Z(1)

      W(1) = 0.23692 68850 5618908751
      W(2) = 0.47862 86704 9936646804
      W(3) = 0.56888 88888 8888888889
      W(4) = W(2)
      W(5) = W(1)

c-------------------------
      Else If(NQ.eq.6) then
c-------------------------

      Z(1) = -0.93246 95142 03152
      Z(2) = -0.66120 93864 66265
      Z(3) = -0.23861 91860 83197

      Z(4) = -Z(3)
      Z(5) = -Z(2)
      Z(6) = -Z(1)

      W(1) = 0.17132 44923 79170
      W(2) = 0.36076 15730 48139
      W(3) = 0.46791 39345 72691

      W(4) = W(3)
      W(5) = W(2)
      W(6) = W(1)

c-------------------------
      Else If(NQ.eq.8) then
c-------------------------

      Z(1) = -0.96028 98564 9753623168
      Z(2) = -0.79666 64774 1362673959
      Z(3) = -0.52553 24099 1632898581
      Z(4) = -0.18343 46424 9564980493

      Z(5) = -Z(4)
      Z(6) = -Z(3)
      Z(7) = -Z(2)
      Z(8) = -Z(1)

      W(1) = 0.10122 85362 9037625915
      W(2) = 0.22238 10344 5337447054
      W(3) = 0.31370 66458 7788728733
      W(4) = 0.36268 37833 7836198296

      W(5) = W(4)
      W(6) = W(3)
      W(7) = W(2)
      W(8) = W(1)

c--------------------------
      Else If(NQ.eq.12) then
c--------------------------

      Z(1) = -0.98156 06342 46719
      Z(2) = -0.90411 72563 70475
      Z(3) = -0.76990 26741 94305
      Z(4) = -0.58731 79542 86617
      Z(5) = -0.36783 14989 98180
      Z(6) = -0.12523 34085 11469

      Z(7) = -Z(6)
      Z(8) = -Z(5)
      Z(9) = -Z(4)
      Z(10)= -Z(3)
      Z(11)= -Z(2)
      Z(12)= -Z(1)

      W(1) = 0.04717 53363 86511
      W(2) = 0.10693 93259 95318
      W(3) = 0.16007 83285 43346
      W(4) = 0.20316 74267 23066
      W(5) = 0.23349 25365 38355
      W(6) = 0.24914 70458 13403

      W(7) = W(6)
      W(8) = W(5)
      W(9) = W(4)
      W(10)= W(3)
      W(11)= W(2)
      W(12)= W(1)

c---------------------------
      Else If (NQ.eq.20) then
c---------------------------

      Z(1) = -0.99312 85991 85094 924786
      Z(2) = -0.96397 19272 77913 791268
      Z(3) = -0.91223 44282 51325 905868
      Z(4) = -0.83911 69718 22218 823395
      Z(5) = -0.74633 19064 60150 792614
      Z(6) = -0.63605 36807 26515 025453
      Z(7) = -0.51086 70019 50827 098004
      Z(8) = -0.37370 60887 15419 560673
      Z(9) = -0.22778 58511 41645 078080
      Z(10)= -0.07652 65211 33497 333755

      Z(11) = -Z(10)
      Z(12) = -Z(9)
      Z(13) = -Z(8)
      Z(14) = -Z(7)
      Z(15) = -Z(6)
      Z(16) = -Z(5)
      Z(17) = -Z(4)
      Z(18) = -Z(3)
      Z(19) = -Z(2)
      Z(20) = -Z(1)

      W(1) = 0.01761 40071 39152 118312
      W(2) = 0.04060 14298 00386 941331
      W(3) = 0.06267 20483 34109 063570
      W(4) = 0.08327 67415 76704 748725
      W(5) = 0.10193 01198 17240 435037
      W(6) = 0.11819 45319 61518 417312
      W(7) = 0.13168 86384 49176 626898
      W(8) = 0.14209 61093 18382 051329
      W(9) = 0.14917 29864 72603 746788
      W(10)= 0.15275 33871 30725 850698

      W(11) = W(10)
      W(12) = W(9)
      W(13) = W(8)
      W(14) = W(7)
      W(15) = W(6)
      W(16) = W(5)
      W(17) = W(4)
      W(18) = W(3)
      W(19) = W(2)
      W(20) = W(1)

c-----------
      End If
c-----------
 
      Return
      End
