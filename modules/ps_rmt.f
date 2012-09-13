      real function ZGOE(S)
      real S
c
      integer i
      real Z,tayW(0:43),tayZ(0:44)
      common /GOE/ tayW,tayZ
c
      if (S .lt. 2.2) then
         Z = tayZ(44)
         do i = 43,0,-1
            Z = Z*S + tayZ(i)
         end do
         ZGOE = Z
      else
         ZGOE = 1.04151836*exp(-(0.785398163*S)**2 - 0.785398163*S  
     .                         - 0.09375/(3.14159265*S + 1.0)**3)
     .        / (3.14159265*S + 1.0)**0.125
      end if
      return     
      end
c------------------------------
      real function WGOE(S)
      real S
c
      integer i
      real W,x,tayW(0:43),tayZ(0:44)
      common /GOE/ tayW,tayZ
c
      if (S .lt. 2.4) then
         W = tayW(43)
         do i = 42,0,-1
            W = W*S + tayW(i)
         end do
         WGOE = W
      else
         x = 3.14159265*S + 1.0
         WGOE = 1.0 + .409003303*exp(-(0.785398163*S)**2 - 0.785398163*S  
     .        - 0.09375/x**3)*(2.25/x**4 - 1./x - x - 1.)/x**0.125
      end if
      return     
      end
c------------------------
c      include 'graph.for'

      SUBROUTINE FIT(X,Y,SIG,N1,N2,A,B,SIGA,SIGB,CHI2)
      INTEGER N1,N2
      REAL*4 X(100000),Y(100000),SIG(100000),A,B,SIGA,SIGB,CHI2
C
      INTEGER I
      REAL*4 WT,SX,SY,SS,ST2,SXOSS
C
      SX = 0D0
      SY = 0D0
      ST2 = 0D0
      B = 0D0
      SS = 0D0
      DO I = N1,N2
         WT = 1D0/(SIG(I)**2)
         SS = SS + WT
         SX = SX + X(I)*WT
         SY = SY + Y(I)*WT
      END DO
      SXOSS = SX/SS
      DO I = N1,N2
         WT = (X(I) - SXOSS)/SIG(I)
         ST2 = ST2 + WT*WT
         B = B + WT*Y(I)/SIG(I)
      END DO 
      B = B/ST2
      A = (SY - SX*B)/SS
      SIGA = DSQRT((1D0 + SX*SX/(SS*ST2))/SS)
      SIGB = DSQRT(1D0/ST2)
      CHI2 = 0D0
      DO I = N1,N2
         CHI2 = CHI2 + ((Y(I) - A - B*X(I))/SIG(I))**2
      END DO
      RETURN
      END   
c----------------------------------------------------------------------
c
      SUBROUTINE SORT(N,RA)
      IMPLICIT REAL (A-H,O-Z)
      DIMENSION RA(N)
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
        ELSE
          RRA=RA(IR)
          RA(IR)=RA(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
      GO TO 10
      END
c-----------------------------------------------------------------------
c      
      subroutine BerryRobnik(Ni,Nf,S,W,dW,rho1,Chi2)
      integer Ni,Nf
      real S(100000),W(100000),dW(100000)                                    
      real rho1,Chi2
c
      integer maxcal,ifail
      real ax,bx,cx,tol,BRChi2
      integer N1,N2
      real S1(100000),W1(100000),dW1(100000)                           
      common /BRdata/ S1,W1,dW1,N1,N2
c
      external BRChi2
c
      N1 = Ni
      N2 = Nf
      do i = N1,N2
         S1(i) = S(i)
         W1(i) = W(i)
         dW1(i) = dW(i)
      end do
      ax =-0.5
      bx = 0.5
      cx = 1.5
      tol = 0.0001
      Chi2 = brent(ax,bx,cx,BRChi2,tol,rho1)
      return
      end
c-----------------------------------------------------------------------
c      
      real function BRChi2(rho1)
      real rho1
c
      integer N1,N2
      real S(100000),W(100000),dW(100000)                                    
      common /BRdata/ S,W,dW,N1,N2
c          
      integer i,j
      real rho2,a,Chi2,WBR,WGOE,ZGOE
c
      Chi2 = 0.0
      rho2 = 1.0 - rho1
      j = 0
      do i = N1,N2
         a = rho2*S(i)
         WBR = 1.0 - exp(-rho1*S(i))*(rho1*ZGOE(a)+rho2*(1.-WGOE(a)))
         Chi2 = Chi2 + ((WBR - W(i))/dW(i))**2
      end do
      write (6,*) 'rho1 =',rho1,', Chi2BR =',Chi2
      BRChi2 = Chi2
      return
      end
c-----------------------------------------------------------------------
c
      function brent(ax,bx,cx,f,tol,xmin)
      parameter (itmax = 100, cgold = .3819660, zeps = 1.0e-10)
c
      a = min(ax,cx)
      b = max(ax,cx)
      v = bx
      w = v
      x = v
      e = 0.0
      fx = f(x)
      fv = fx
      fw = fx
      do iter = 1,itmax
         xm = 0.5*(a + b)
         tol1 = tol*abs(x) + zeps
         tol2 = 2.0*tol1
         if (abs(x-xm) .le. (tol2 - 0.5*(b - a))) go to 3
         if (abs(e) .gt. tol1) then
            r = (x - w)*(fx - fv)
            q = (x - v)*(fx - fw)
            p = (x - v)*q - (x - w)*r
            q = 2.0*(q - r)
            if (q .gt. 0.0) p = -p
            q = abs(q)
            etemp = e
            e = d
            if (abs(p) .ge. abs(0.5*q*etemp) .or. p .le. q*(a - x) .or.
     .         p .ge. q*(b - x)) go to 1
            d = p/q
            u = x + d
            if (u - a .lt. tol2 .or. b - u .lt. tol2) d=sign(tol1,xm-x)
            go to 2
         end if
1        if (x .ge. xm) then
            e = a - x
         else
            e = b - x
         end if
         d = cgold*e
2        if (abs(d) .ge. tol1) then
            u = x + d
         else
            u = x + sign(tol1,d)
         end if
         fu = f(u)
         if (fu .le. fx) then
            if (u .ge. x) then
               a = x
            else
               b = x
            end if
            v = w
            fv = fw
            w = x
            fw = fx
            x = u
            fx = fu
         else
            if (u .lt. x) then
               a = u
            else
               b = u
            end if
            if (fu .le. fw .or. w .eq. x) then
               v = w
               fv = fw
               w = u
               fw = fu
            else if (fu .le. fv .or. v .eq. x .or. v. eq. w) then
               v = u
               fv = fu
            end if
         end if
      end do
      pause 'Brent exceed maximum iterations!'
3     xmin = x
      brent = fx
      return
      end
c---------------------------
      subroutine initGOE
      real tayP(0:42),tayW(0:43),tayZ(0:44)
      data tayP/
     . 0.00000000E+00,
     . 1.64493407E+00,
     . 0.00000000E+00,
     .-1.62348485E+00,
     . 0.36077441E+00,
     . 0.57225547E+00,
     .-0.20346861E+00,
     .-0.10459139E+00,
     . 0.51996863E-01,
     . 0.11730404E-01,
     .-0.84012045E-02,
     .-0.89057276E-03,
     . 0.10192294E-02,
     . 0.43354994E-04,
     .-0.10269378E-03,
     . 0.14456928E-06,
     . 0.89788109E-05,
     .-0.35907897E-06,
     .-0.68613431E-06,
     . 0.53373192E-07,
     . 0.45633878E-07,
     .-0.53051943E-08,
     .-0.26365437E-08,
     . 0.41366745E-09,
     . 0.13268190E-09,
     .-0.26832852E-10,
     .-0.58466725E-11,
     . 0.14984853E-11,
     . 0.22679440E-12,
     .-0.73871101E-13,
     .-0.77606149E-14,
     . 0.32799913E-14,
     . 0.23256905E-15,
     .-0.13336114E-15,
     .-0.59064987E-17,
     . 0.50305714E-17,
     . 0.11235768E-18,
     .-0.17770689E-18,
     .-0.55843340E-21,
     . 0.59124699E-20,
     .-0.86112533E-22,
     .-0.18575824E-21,
     . 0.58771977E-23/
      integer i
      common /GOE/ tayW,tayZ
c      
      tayW(0) = 0.0
      tayZ(0) = 1.0
      tayZ(1) =-1.0
      do i = 0,42
         tayW(i+1) = tayP(i)/float(i+1)
         tayZ(i+2) = tayW(i+1)/float(i+2)
      end do
      return
      end   
c------------------------------

