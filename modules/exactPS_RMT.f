      program b24br
c
c
      integer m,Ncl,Ninvnl,mi,mf,i,j,k,ifail,nspc,N,N0,N1,Nm
      real a,bup,bdown,eps,res,x,y,z,beta,dbeta,Chi2,gamma,dgamma
      real Smax,Swin,Wwin,W1,WBR,Uerr,r,rN,rho,UAM(100000)
      real minX,maxX,minY,maxY,rho1,Chi2BR,Uu(100000),Ul(100000)
      real lnS(100000),T(100000),dT(100000),U1(-2:2,100000)
      real S(100000),W(100000),dW(100000),U(100000),UB(100000)
      real Ep(100000),E(100000),Smin
      real*8 xx,E0,E1,c0,c1,c2,c3
      real WGOE,ZGOE
      character*20 fn
      character*50 comment
c
      call initGOE
      write (6,10) 'minimal spacing, Smin = '
      read (5,*) Smin
      write (6,10) 'maximal spacing, Smax = '
      read (5,*) Smax     
      write (6,10) ' small S window, Swin = '
      read (5,*) Swin
      write (6,10) 'rho = '
      read (5,*) rho
      write (6,10) '       input filename = '
      read (5,15) fn   
      write (6,10) 'Comment(none=data head) '
      read (5,15) comment
10    format (' ',A,$)
15    format (A)
c
      open (unit = 1,status = 'old',file = fn,err = 999)
      read (1,*)
      read (1,*) eps,nl,lmax
      c3 = 0d0
      c2 = 1d0/32d0
      c1 = 0d0
      c0 = 0d0
      N = 0
      N1 = 0
      N0 = 0
      do i = 1,nl
         read (1,*) xx
         E1 = c3*xx**3 + c2*xx**2 + c1*xx + c0
         if (i .eq. 1) E0 = E1
         E(i) = sngl(E1 - E0)
      end do
      call sort(nl,E)
      do i = 2,nl
         N = N + 1
         S(N) = E(i) - E(i-1)
         if (S(N) .le. 0.0) write (6,*) 'S(',N,') = ',S(N)
         if (S(N) .lt. Smax) N1 = N1 + 1
         if (S(N) .le. Smin) N0 = N0 + 1         
      end do
120   close (1)
      call sort(N,S)
      rN = float(N)
      do i = 1,N 
         lnS(i) = log(S(i)+1d-10)
         W(i) = (float(i) - 0.5)/rN
         dW(i) = sqrt(W(i)*(1. - W(i))/rN)
         T(i) = log(-log(1. - W(i)))
         dT(i) = -sqrt(W(i)/((1. - W(i))*rN))/log(1. - W(i))
      end do
      write (6,*) 'nl =',nl,', N1-N0=',N1-N0         
      write (6,*) 'Brody fit...'
      call fit(lnS,T,dT,N0+1,N1,gamma,beta,dgamma,dbeta,Chi2)
      write (6,*) 'beta =',beta-1.,', Chi2 =',Chi2
      write (6,*) 'Berry-Robnik fit...'
      call BerryRobnik(N0+1,N1,S,W,dW,rho1,Chi2BR)      
      write (6,*) 'rho1 =',rho1,', Chi2 =',Chi2BR
      minX = 0.0
      maxX = Smax
      minY = 0.0
      maxY = 1.0
      call InitGraph(2,0.,10.,0.,10.,'sso.W')
      if (comment(1:1) .ne. ' ') then
         call gotoXY(2,370,50)
         write (2,*) comment
         call gotoXY(2,470,300)
         write (2,*) '(a)'
         call gotoXY(2,470,2250)
         write (2,*) '(b)'
         go to 678
      end if
      call gotoXY(2,100,1)
      write (2,*) fn
      call gotoXY(2,100,61)
      write (2,810) gamma,100.*dgamma/abs(gamma),beta,100.*dbeta/beta,
     .nint(Chi2)
810   format ('T(lnS)=',F8.5,'(1+-',F5.2,'%)',
     .' +',F8.5,'(1+-',F5.2,'%)lnS, Chi2=',I9)
      call gotoXY(2,100,121)
      write (2,820) rho1,nint(Chi2BR)
820   format ('Berry-Robnik fit: rho1=',F7.4,', Chi2=',I8)
678   call scale(400,2100,200,1900,minX,maxX,minY,maxY)
      call frame(2,minX,maxX,minY,maxY,'S','W(S)',' ',' ')
      call SetStyle(1,2.1,1,1,1)
      call move(0.0,0.0)
      do j = 1,N1
         call draw(S(j),W(j))
      end do
      call SetStyle(1,2.3,11,1,1)
      call move(0.0,0.0)
      do j = 1,500
         x = float(j)*Smax/500.0
         call draw(x,1.0 - exp(-x))
      end do
      call move(0.0,0.0)
      do j = 1,500
         x = float(j)*Smax/500.0
         call draw(x,WGOE(x))
      end do
      call SetStyle(1,1.1,1,25,14)
      call move(0.0,0.0)
      do j = 1,500
         x = float(j)*Smax/500.0
         call draw(x,1.0 - exp(-exp(gamma)*x**beta))
      end do
      call SetStyle(1,1.1,6,24,18)
      call move(0.0,0.0)
      do j = 1,500
         x = float(j)*Smax/500.0         
         call draw(x,1.0-exp(-x*(qq + x*qq1)))
      end do  
      call SetStyle(1,1.1,1,1,1)
      call move(0.0,0.0)
      do j = 1,500
         x = float(j)*Smax/500.0
         y = (1-rho1)*x
         call 
     .   draw(x,1.0-exp(-rho1*x)*(rho1*ZGOE(y)+(1-rho1)*(1.-WGOE(y))))
      end do
      i = 1
      do while (S(i) .lt. Swin) 
         i = i + 1
      end do
      Wwin = (float(i) - 0.5)/float(N)
      call SetStyle(1,1.8,6,24,18)
      call move(Swin,0.0)
      call draw(Swin,Wwin)
      call draw(0.0,Wwin)
      call scale(1070,2070,870,1870,0.0,Swin,0.0,Wwin)
      call move(0.0,0.0)
      call draw(Swin,0.0)
      call draw(Swin,Wwin)
      call draw(0.0,Wwin)
      call draw(0.0,0.0)
      call SetStyle(1,2.1,1,1,1)
      call move(0.0,0.0)
      do j = 1,i-1
         call draw(S(j),W(j))
      end do
      call SetStyle(1,2.3,11,1,1)
      call move(0.0,0.0)
      do j = 1,500
         x = float(j)*Swin/500.0
         W1 = 1.0 - exp(-x)
         if (W1 .lt. Wwin) call draw(x,W1)
      end do
      call move(0.0,0.0)
      do j = 1,500
         x = float(j)*Swin/500.0
         W1 = WGOE(x)
         if (W1 .lt. Wwin) call draw(x,W1)
      end do
      call SetStyle(1,1.1,1,25,14)
      call move(0.0,0.0)
      do j = 1,500
         x = float(j)*Swin/500.0
         W1 = 1.0 - exp(-exp(gamma)*x**beta)
         if (W1 .lt. Wwin) call draw(x,W1)
      end do
      call SetStyle(1,1.1,1,1,1)
      call move(0.0,0.0)
      do j = 1,500
         x = float(j)*Swin/500.0
         y = (1-rho1)*x
         WBR = 1.0-exp(-rho1*x)*(rho1*ZGOE(y)+(1-rho1)*(1.-WGOE(y)))
         if (WBR .lt. Wwin) call draw(x,WBR)
      end do

c
c     U-function
c
      minY = 9E9
      maxY = -9E9
      Uerr = 1.0/(3.14159*sqrt(float(N)))
      mi1 = 1
      mi2 = 1
      ma1 = 1
      ma2 = 1
      do j = 1,N
         y = (1-rho1)*S(j)         
         WBR = 1.0-exp(-rho1*S(j))*(rho1*ZGOE(y)+(1-rho1)*(1.-WGOE(y)))
         x = acos(sqrt(1-WBR))
c        Wu = (W(j) + r)/(1. + r)
c        Wl = W(j)/(1. + r)
         Wu = W(j)
         Wl = Wu
c
         U(j) = 0.63661977*(acos(sqrt(1.-W(j)))-x)        
         if (U(j) .gt. U(ma1)) then
            ma2 = ma1
            ma1 = j
         else if (U(j) .gt. U(ma2)) then
            ma2 = j
         end if
         if (U(j) .lt. U(mi1)) then
            mi2 = mi1
            mi1 = j
         else if (U(j) .lt. U(mi2)) then
            mi2 = j
         end if
         Uu(j) = 0.63661977*(acos(sqrt(1.-Wu))-x)        
         Ul(j) = 0.63661977*(acos(sqrt(1.-Wl))-x)        
         do i = -2,2
            if (i .ne. 0) then
               r1 = rho1 + .01*float(i)
               y = (1-r1)*S(j)
               WBR = 1.0-exp(-r1*S(j))*(r1*ZGOE(y)+(1-r1)*(1.-WGOE(y)))
               U1(i,j) = 0.63661977*(acos(sqrt(1-WBR))-x)
               minY = min(minY,U1(i,j))
               maxY = max(maxY,U1(i,j))
            end if
         end do
         W1 = 1.0 - exp(-exp(gamma)*S(j)**beta)         
         UB(j) = 0.63661977*(acos(sqrt(1-W1))-x)

         minY = min(minY,min(U(j)-Uerr,min(Ul(j),UB(j))))
         maxY = max(maxY,max(U(j)+Uerr,max(Uu(j),UB(j))))
      end do
      do i = 2,nl
         S1 = E(i) - E(i-1)
         if (S1 .eq. S(ma1)) 
     .      write (6,*) 'max1,i=',i,' S=',S1,', E=',E(i-1),E(i)
         if (S1 .eq. S(ma2)) 
     .      write (6,*) 'max2,i=',i,' S=',S1,', E=',E(i-1),E(i)
         if (S1 .eq. S(mi1)) 
     .      write (6,*) 'min1,i=',i,' S=',S1,', E=',E(i-1),E(i)
         if (S1 .eq. S(mi2)) 
     .      write (6,*) 'min2,i=',i,' S=',S1,', E=',E(i-1),E(i)
      end do
      minX = 0.0
      maxX = 1.0
      call scale(400,2100,2150,3000,minX,maxX,minY,maxY)
      call SetStyle(-3,4.0,2,1,1)
      call move(W(1),U(1)-Uerr)
      do j = 2,N
         call draw(W(j),U(j)-Uerr)
      end do
      call move(W(1),U(1)+Uerr)
      do j = 2,N
         call draw(W(j),U(j)+Uerr)
      end do
      call SetStyle(1,2.3,2,1,1)
      call frame(2,minX,maxX,minY,maxY,'W','U-UBR',' ',' ')
      call move(minX,0.0)
      call draw(maxX,0.0)
      call SetStyle(1,2.8,2,1,1)
      call move(W(1),U(1))
      do j = 2,N
         call draw(W(j),U(j))
      end do
c     call SetStyle(1,1.5,1,1,1)
c     call move(W(1),Uu(1))
c     do j = 2,N
c        call draw(W(j),Uu(j))
c     end do
c     call move(W(1),Ul(1))
c     do j = 2,N
c        call draw(W(j),Ul(j))
c     end do
      call SetStyle(1,1.8,8,36,24)
      do i = -2,2
         if (i .ne. 0) then
            call move(W(1),U1(i,1))
            do j = 2,N,N/120
               call draw(W(j),U1(i,j))
            end do
         end if
      end do
      call SetStyle(1,1.8,1,25,14)
      call move(W(1),UB(1))    
      mm = N/500 + 1
      do j = 2,N,mm
         call draw(W(j),UB(j))
      end do
      call draw(W(N),UB(N))
      call print(2)
      call CloseGraph(2)
888   stop
999   write (6,*) 'File ',fn,' not found!'
      stop
      end
c----------------------------------------------------------------------
c
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
      include 'graph.for'

