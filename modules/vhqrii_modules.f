C----------------------------------------------------------
      SUBROUTINE VHQRII (N,IV,LV,IORD,AL_VEC,E,NVX,SCHMIDT,
     *                                      IWORK,WORK,IER)
C----------------------------------------------------------

C   Title: module VHQRII to diagonalize dense real-symmetric matrices
C          using vector processors.

C   Abstract: An input real-symmetric matrix AL_VEC is tridiagonalized
C             by the method of Householder.  The eigenvalues of the
C             tridiagonal matrix T are found by the QR method with 
C             origin shift.  Optionally, some or all eigenvectors of
C             matrix T are determined by inverse iteration.  The 
C             eigenvectors of matrix AL_VEC are found by 
C             left-multiplying the eigenvectors of matrix T times the
C             orthogonal matrix which brings AL_VEC to tridiagonal form.

C  Environment: Standard Fortran 77, 90 and 95 (with obsolescent GO TO).

C   Copyright by Manuel Berrondo, Annik Vivier Bunge, Carlos F. Bunge,
C                Yoshitaka Beppu and Ichizo Ninomiya, 1986, 2000.

C   Version  2-001: August, 2000.
C   Versions 1-001 through 1-005: March 1985 through December 1994.

C   Reference: M.Berrondo, A.V.Bunge and C.F.Bunge, Comput. Chem. 10,269
C              (1986).

C   Machine dependent parameter:

C   MACHEP      is the smallest number which makes 1D0 + MACHEP greater
C               than 1D0.  For an IBM computer MACHEP=1.11D-16.
C               For DEC Compaq and Vaxes MACHEP=1.4D-17.

C   Parameter values:

C   P5          is equal to 0.5D0.
C   BOUNDR=1D-6 is a small number involved in deciding whether two
C               eigenvalues are close enough to be considered degenerate
C               while computing eigenvectors of the tridiagonal matrix.
C   GOLDRT      GOLDRT=0.5D0*(SQRT(5D0)-1D0) is a golden random number
C               between 0 and 1.  Other choices may affect the accuracy
C               of nearly degenerate eigenvectors, for better or for
C               worse, depending on the particular case considered.

C   Argument values:

C   N           order of the matrix to be diagonalized in a given run.
C   IV          index of first wanted eigenvector.
C   LV          index of last wanted eigenvector.  If LV .LT. IV, no
C               eigenvectors are calculated.
C   IORD        controls the order of eigenvalues.  IORD .GE. 0 will
C               cause eigenvalues to be given in non-increasing order.
C               IORD .LT. 0 specifies non-decreasing order.
C   AL_VEC      On input:
C               real symmetric matrix of dimension at least N1=(N*N+N)/2
C               in rowwise lower triangular form, viz., AL_VEC(IJ)=
C               AA(I,J) where IJ=J+(I*I-I)/2 and AA is the usual square
C               symmetric form.
C               On output:
C               array of dimension NVX*N which contains the J-th eigen-
C               vector in positions AL_VEC(IJ) starting at 
C               IJ=(J-1)*NVX+1.
C   E           array of dimension at least N holding the eigenvalues
C               in the order determined by variable IORD.
C   NVX         column dimension of array AL_VEC in calling program
C               or within the subroutine in which the results are to be
C               used.
C               If AL_VEC is a one dimensional array in the calling
C               program, most likely NVX=N.
C   SCHMIDT     =.TRUE. for normal use.  If SCHMIDT=.FALSE. eigenvectors
C               of degenerate or nearly degenerate eigenvalues are not
C               orthogonal among themselves.
C   WORK        is an array of size at least 21*N + (N*N+N)/2.
C   IER         is a completion status indicator, see below.

C   Accuracy:
C            for the largest (in absolute value) eigenvalues E(I) which
C            are not highly degenerate, machine accuracy for MR(I),
C                               MR(I) = MODR(I)/E(I)
C            where MODR is the modulus of the residual vector R(I),
C                                R(I) = (A*V)(I) - E(I)*V(I),
C            and V(I) denotes the I-th eigenvector.  Up to M figures may
C            be lost for the smallest eigenvalues (and eigenvectors),
C            where 
C            M = ALOG10 (LARGEST/SMALLEST)
C            with LARGEST and SMALLEST being the largest and smallest
C            (in absolute value) eigenvalues, respectively.
C            In presence of highly degenerate eigenvalues, users may
C            enter SCHMIDT=.FALSE., which will produce very accurate
C            eigenvectors, but then it becomes their responsibility to
C            produce an orthogonal set. 

C   Completion status:

C   IER=0    normal successful completion.
C   IER=33   means that AL_VEC is exactly a null matrix.
C   IER=129  indicates a fatal error due to N, LV or NVX being outside
C            permissible bounds.

C   Caveats: matrices formed exclusively by very small or very large
C            matrix elements will give an arithmetic overflow.

C            In all processors tested so far, HQRII1 is more efficient
C            than VHQRII.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SCHMIDT
      DIMENSION AL_VEC(*),E(*),WORK(*),IWORK(*)


      CALL VHQRII_CODE (N,IV,LV,IORD,AL_VEC,E,NVX,SCHMIDT,
     *                 IWORK,WORK(N+1),WORK(2*N+1),WORK(3*N+1),
     *                  WORK(4*N+1),WORK(5*N+1),WORK(6*N+1),
     *                  WORK(7*N+1),WORK(8*N+1),WORK(9*N+1),
     *                                         WORK(10*N+1),IER)
      END
C---------------------------------------------------------------
      SUBROUTINE VHQRII_CODE (N_NL,IV_NL,LV_NL,IORD_NL,AL_VEC,E,
     *                        NVX_NL,SCHMIDT,IX,W1,W2,W1II,W7,
     *                        W5,W3,W4,W5I,W6,AU           ,IER)
C---------------------------------------------------------------

C   IX   is a working array of dimension N first defined by IX(I)=
C        (I*I-I)/2 so that the lower triangle matrix AL_VEC is 
C        indexed by IJ=IX(I)+J.  Further on, IX is defined by 
C        IX(1)=0 and IX(I)=IX(I-1)+N-I+1 for I.GT.1 so that the
C        upper triangle matrix AU is indexed by IJ=IX(I)+J.
C   W1-W7,W1II,W5I are nine working vectors each of dimension N.
C   AU   real symmetric matrix of dimension N1=(N*N+N)/2 given in
C        rowwise upper triangular form, viz., AU(IJ)= AA(I,J) where
C        IJ=IX(I)+J as explained below and AA is the usual square
C        symmetric form.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION MACHEP
      LOGICAL SCHMIDT
      PARAMETER (P5=0.5D0, BOUNDR=1D-6,
     *           MACHEP=1.4D-17, GOLDRT=0.618033988749894D0)
      DIMENSION AL_VEC(*),E(*),
     *          IX(*),W1(*),W2(*),W1II(*),W7(*),W5(*),
     *          W3(*),W4(*),W5I(*),W6(*),AU(*)

      N    = N_NL
      IV   = IV_NL
      LV   = LV_NL
      NVX  = NVX_NL
      IORD = IORD_NL
      IF (IV .LE. 0) IV = 1

      IER = 0
      IF (N .LT. 1      .OR.
     *   (N .GT. NVX  .AND. LV .GE. IV)  .OR.
     *                      LV .GT. N )  THEN
          PRINT *,' N OR LV OR NVX OUTSIDE PERMISSIBLE BOUNDS'
          PRINT *,' N=',N,' LV=',LV,' NVX=',NVX
          IER = 129
          RETURN
      END IF
      IF (N .EQ. 1) THEN
               E(1) = AL_VEC(1)
          AL_VEC(1) = 1D0
          RETURN
      END IF

      IX(1) = 0
      DO J=2,N
         IX(J) = IX(J-1) + J - 1
      END DO

      IER = 0
      NM1 = N - 1
      NM2 = N - 2
      NVF =(LV-1)*NVX
      IF (NVF .LT. 0) NVF = 0

      IF (N .GT. 2) THEN

C     Householder transformation.
          DO K=1,NM2
              KP1  = K + 1
             W2(K) = AL_VEC(K+IX(K))
             SCALE = 0D0
             DO J=KP1,N
                SCALE = SCALE + ABS (AL_VEC(IX(J)+K))
             END DO
             W1(K) = AL_VEC(IX(KP1)+K)
             IF (SCALE .GT. 0D0) THEN
                 SCALEI = 1D0/SCALE
                 DO J=KP1,N
                    W2(J) = AL_VEC(IX(J)+K)*SCALEI
                 END DO
                 DOTP = 0D0
                 DO J=KP1,N
                     DOTP = DOTP + W2(J)**2
                 END DO
                 S = SIGN (SQRT (DOTP),W2(KP1))
                 W1(K  ) = -S*SCALE
                 W2(KP1) = W2(KP1) + S
                 AL_VEC(IX(KP1)+K) = W2(KP1)*SCALE
                 H     = W2(KP1)*S
                 HUNS  = (H*SCALE)*SCALE
                 HI    = 1D0/H
                 IVSIZE= N - KP1 + 1

C*Vectorizable code, Section 1.
                 CALL ZRO1D (IVSIZE,W5(KP1))
C*End of Section 1.

C                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C                First numerically intensive kernel.
                 DO I=KP1,N
                    IM1   = I - 1
                    I0    = IX(I)
                    W2I   = W2(I)
                    IVSIZE= IM1 - KP1 + 1

C*Vectorizable code, Section 2.
                    CALL VDOTD_VSMA1D (IVSIZE,
     *                                 AL_VEC(I0+KP1),W2(KP1),DOTP,
     *                                                 W2I,W5(KP1))
                    W6(I) = W2I*AL_VEC(I0+I) + DOTP
                 END DO
C                End first numerically intensive kernel.
C                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

                 IVSIZE= N - KP1 + 1
                 CALL VSM2D (IVSIZE,HI,W5(KP1),W6(KP1),W1(KP1))
                 CALL VDOTD (IVSIZE,W1(KP1),W2(KP1),DOTP)
C*End of Section 2.

                 U = P5*DOTP*HI
C                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C                Second numerically intensive kernel.
                 DO I=KP1,N
                    I0    = IX(I)
                    W1(I) = W2(I)*U - W1(I)
                    W1I   = W1(I)
                    W2I   = W2(I)
                    IVSIZE= I- KP1 + 1

C*Vectorizable code, Section 3.
                    CALL VSMA2D (IVSIZE,W2I,W1I,W1(KP1),W2(KP1),
     *                                           AL_VEC(I0+KP1))
C*End of Section 3.
                 END DO
C                End second numerically intensive kernel.
C                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
             ELSE
                 HUNS = 0D0
             END IF
             AL_VEC(IX(K)+K) = HUNS
          END DO
      END IF

      NM1NM1  = IX(NM1) + NM1
      NM1N    = IX(N  ) + NM1
      NN      = NM1N    + 1
      W2(NM1) = AL_VEC(NM1NM1)
      W2(N  ) = AL_VEC(NN)
      W1(NM1) = AL_VEC(NM1N)
      W1(N  ) = 0D0

C     End Householder tridiagonalization.

      GERSCH  = ABS (W2(1)) + ABS (W1(1))
      DO I=1,NM1
         GERSCH = MAX (ABS(W2(I+1)) + ABS(W1(I)) + ABS(W1(I+1)), GERSCH)
      END DO

C     Trap null matrix before it is too late.

      IF (GERSCH .EQ. 0D0) THEN
          PRINT *,' NULL MATRIX IN SUBROUTINE VHQRII'
          IER = 33
          RETURN
      END IF

      SUMD  = 0D0
      SUMCOD= 0D0 
      DO I=1,N
         SUMCOD= ABS (W1(I)) + SUMCOD
          SUMD = ABS (W2(I)) + SUMD
      END DO
      SCALE = SUMD + SUMCOD
      SCALEI= 1D0/SCALE
      DO I=1,N
         W1(I) = W1(I)*SCALEI
         W2(I) = W2(I)*SCALEI
         W3(I) = W1(I)
          E(I) = W2(I)
      END DO
      EPS   = SQRT(MACHEP)
      GERSCH= GERSCH*SCALEI
      DEL   = GERSCH*EPS
      DELW5 = GERSCH*MACHEP

C     QR method with origin shift.

      DO K=N,2,-1
   10    L = K
   20    IF (ABS (W3(L-1)) .GT. DEL) THEN
             L = L - 1
             IF (L .GT. 1) GO TO 20
         END IF
         IF (L .NE. K) THEN
             WW = (E(K-1)+E(K))*P5
             R  =  E(K)-WW
             Z  = WW - SIGN (SQRT (W3(K-1)**2 + R*R),WW)
             EE = E(L) - Z
             FF = W3(L)
             R  = SQRT (EE*EE + FF*FF)
             RI = 1D0/R
            E(L)= EE
             C  =  E(L)*RI
             S  = W3(L)*RI
             WW =  E(L+1) - Z
            E(L)= (FF*C + WW*S)*S + EE + Z
          E(L+1)= C*WW - S*FF

C            //////////////// 
C            Start innerloop.
             DO J=L+1,K-1
                R  = SQRT (E(J)**2 + W3(J)**2)
                RI = 1D0/R
                EE =  E(J)*C + Z
                FF = W3(J)*C
                WW = E(J+1) - Z
                W3(J-1) = S*R
                C  =  E(J)*RI
                S  = W3(J)*RI
               E(J)= (FF*C + WW*S)*S + EE 
             E(J+1)= C*WW - S*FF
             END DO
C            End innerloop.
C            //////////////

             W3(K-1) = E(K)*S
              E(K  ) = E(K)*C + Z
             GO TO 10
         END IF
      END DO

C     Straight selection sort of eigenvalues.

      SORTER = 1D0
      IF (IORD .LT. 0) SORTER = -1D0
      J = N
  30  L = 1
      II= 1
      LL= 1            
      DO I=2,J
         IF ((E(I)-E(L))*SORTER .LE. 0D0) THEN
             L = I
         ELSE
             II= I
             LL= L
         END IF
      END DO
      IF (II .NE. LL) THEN
            WW = E(LL)
          E(LL)= E(II)
          E(II)= WW
      END IF
      J = II - 1
      IF (J .GT. 1) GO TO 30

      IF (LV .GE. IV) THEN

          CALL TRANSPOSE (N,AL_VEC,W4,W5,AU)
C         Inverse iteration for eigenvectors.

          DO J=1,N
             IF (W1(J) .EQ. 0D0) W1(J) = DEL
             W1II(J) = 1D0/W1(J)
          END DO


          FN  = DBLE(N)
          EPS1= SQRT(FN)*EPS
          SEPS= SQRT(EPS)
          EPS2=(GERSCH*BOUNDR)/(FN*SEPS)

C         Find if first wanted vector is too close to previous ones.
   40     CONTINUE
          IF (IV .EQ. 1) GO TO 50
          IF (ABS(E(IV-1)-E(IV)) .LT. EPS2) THEN
              IV = IV - 1
              GO TO 40
          END IF
   50     CONTINUE

          RN  = 0D0
          RA  = EPS*GOLDRT
          IF  = (IV-2)*NVX
          DO I=IV,LV
             IF  = IF + NVX 
             DO J=1,N
                W3(J) = 0D0
                W4(J) = W1(J)
                W5(J) = W2(J) - E(I)
                RN = RN + RA
                IF (RN .GE. EPS) RN = RN - EPS
                W6(J) = RN
             END DO
             DO J=1,NM1
                IF (ABS (W5(J)) .LE. ABS (W1(J))) THEN
                    W7(J) = -W5(J)*W1II(J)
                    W5(J) =  W1(J)
                   W5I(J) = W1II(J)
                    T = W5(J+1)
                   W5(J+1)=  W4(J)
                    W4(J) =  T
                    W3(J) =  W4(J+1)
                    IF (W3(J) .EQ. 0D0) W3(J) = DEL
                   W4(J+1)=  0D0
                ELSE
                    W5I(J) = 1D0/W5(J)
                    W7(J)  = -W1(J)*W5I(J)
                END IF
                W4(J+1) = W4(J+1) + W3(J)*W7(J)
                W5(J+1) = W5(J+1) + W4(J)*W7(J)
             END DO
             IF (W5(N) .EQ. 0D0) W5(N) = DELW5
             WNM15I = 1D0/W5(NM1)
             WN5I   = 1D0/W5(N  )

             W6(N  ) = W6(N  )*WN5I
             W6(NM1) =(W6(NM1)-W6(N )*W4(NM1))*WNM15I
                  VN = MAX (ABS (W6(N )), ABS (W6(NM1)))
             DO K = NM2,1,-1
                W6(K)=(W6(K)-W6(K+1)*W4(K)-W6(K+2)*W3(K))*W5I(K)
                  VN = MAX (ABS (W6(K)), VN)
             END DO
             S = EPS1/VN 
C*Vectorizable code, Section 4.
             CALL VSMD (N,S,W6,W6)
C*End of Section 4.

             DO J=1,NM1
                IF (W5(J) .EQ. W1(J)) THEN
                        T   = W6(J  )
                    W6(J  ) = W6(J+1)
                    W6(J+1) = T
                END IF
                W6(J+1) = W6(J+1) + W6(J)*W7(J)
             END DO
             W6(N ) = W6(N )*WN5I
             W6(NM1) =(W6(NM1)-W6(N)*W4(NM1))*WNM15I
                  VN = MAX (ABS (W6(N)), ABS (W6(NM1)))
             DO K = NM2,1,-1
                W6(K)=(W6(K)-W6(K+1)*W4(K)-W6(K+2)*W3(K))*W5I(K)
                  VN = MAX (ABS (W6(K)), VN)
             END DO
             S = EPS1/VN 
C*Vectorizable code, Section 4.
             CALL VSMD (N,S,W6,W6)
C*End of Section 4.

C*Vectorizable code, Section 5.
             CALL COPYV (N,W6,AL_VEC(IF+1))
C*End of Section 5.
          END DO

C         Use sequential row notation for upper triangular matrix.

          IX(1) = 0
          DO J=2,N
             IX(J) = IX(J-1) - J + 1 + N
          END DO

C         Back transformation of eigenvectors.

          IG = 1
          IF = (IV-2)*NVX
          DO I=IV,LV
             IF  = IF + NVX

C*Vectorizable code, Section 6.
             CALL COPYV (N,AL_VEC(IF+1),W6)
C*End of Section 6.

             IM1 = I - 1
             IF (N .GT. 2) THEN
                 DO J=1,NM2
                    K = N - J - 1
                    K0= IX(K)
                    IF (AU(K0+K) .NE. 0D0) THEN
                        KP1   = K + 1
                        IVSIZE= N - K

C*Vectorizable code, Section 7.
                        CALL VDOTD (IVSIZE,AU(K0+KP1),W6(KP1),DOTP)
                        S = -DOTP/AU(K0+K)
                        CALL VSMA1D (IVSIZE,S,AU(K0+KP1),W6(KP1))
C*End of Section 7.
                    END IF
                 END DO
             END IF

             DO J=IG,I
                JJ = J
                IF (ABS(E(J)-E(I)) .LT. EPS2) GO TO 60
             END DO
   60        IG = JJ

             IF (IG .NE. I  .AND.  SCHMIDT) THEN

C                Degenerate eigenvalues.  First, orthogonalize.

                 KF  = (IG-2)*NVX
                 DO K=IG,IM1
                    KF  = KF + NVX

C*Vectorizable code, Section 8.
                    CALL VDOTD (N,W6,W4,DOTP)
                    S = -DOTP
                    CALL VSMA1D (N,S,W4,W6)
C*End of Section 8.
                 END DO
             END IF

C            Normalization.

C*Vectorizable code, Section 9.
             CALL VDOTD (N,W6,W6,DOTP)
             S = 1D0/SQRT(DOTP)
             CALL VSMD (N,S,W6,AL_VEC(IF+1))
C*End of Section 9.
          END DO
      END IF

C*Vectorizable code, Section 10.
      CALL VSMD (N,SCALE,E,E)
C*End of Section 10.
      END
C---------------------------------
      SUBROUTINE VDOTD (N,V1,V2,S)
C---------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V1(*),V2(*)

      NL = N
      SL = 0D0
      NREST = MOD(NL,4)
      IF (NREST .EQ. 0) GO TO 10
      DO I=1,NREST
         SL =  SL + V1(I)*V2(I)
      END DO
 10   NRESTP1 = NREST + 1
      DO I=NRESTP1,NL,4
         SL =  SL + V1(I  )*V2(I  )
     *            + V1(I+1)*V2(I+1)
     *            + V1(I+2)*V2(I+2)
     *            + V1(I+3)*V2(I+3)
      END DO
      S = SL

      END
C-------------------------------
      SUBROUTINE VSMD (N,S,V1,V)
C-------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V1(*),V(*)

      NL = N
      SL = S
      NREST = MOD(NL,4)
      IF (NREST .EQ. 0) GO TO 10
      DO I=1,NREST
         V(I) = SL*V1(I)
      END DO
 10   NRESTP1 = NREST + 1
      DO I=NRESTP1,NL,4
         V(I  ) = SL*V1(I  )
         V(I+1) = SL*V1(I+1)
         V(I+2) = SL*V1(I+2)
         V(I+3) = SL*V1(I+3)
      END DO

      END
C-----------------------------------
      SUBROUTINE VSM2D (N,S,V1,V2,V)
C-----------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V1(*),V2(*),V(*)

      NL = N
      SL = S
      NREST = MOD(NL,4)
      IF (NREST .EQ. 0) GO TO 10
      DO I=1,NREST
         V(I) = SL*(V1(I)+V2(I))
      END DO
 10   NRESTP1 = NREST + 1
      DO I=NRESTP1,NL,4
         V(I  ) = SL*(V1(I  ) + V2(I  ))
         V(I+1) = SL*(V1(I+1) + V2(I+1))
         V(I+2) = SL*(V1(I+2) + V2(I+2))
         V(I+3) = SL*(V1(I+3) + V2(I+3))
      END DO
      END
C---------------------------
      SUBROUTINE ZRO1D (N,V)
C---------------------------
      IMPLICIT DOUBLE PRECISION (C,V)
      DIMENSION V(*)

      NL = N
      NREST = MOD(NL,4)
      IF (NREST .EQ. 0) GO TO 10
      DO I=1,NREST
         V(I) = 0D0
      END DO 
 10   NRESTP1 = NREST + 1
      DO I=NRESTP1,NL,4
         V(I  ) = 0D0
         V(I+1) = 0D0
         V(I+2) = 0D0
         V(I+3) = 0D0
      END DO
      END
C------------------------------
      SUBROUTINE COPYV (N,V1,V)
C------------------------------
      IMPLICIT DOUBLE PRECISION (C,V)
      DIMENSION V1(*),V(*)

      NL = N
      NREST = MOD(NL,4)
      IF (NREST .EQ. 0) GO TO 10
      DO I=1,NREST
         V(I) = V1(I)
      END DO 
 10   NRESTP1 = NREST + 1
      DO I=NRESTP1,NL,4
         V(I  ) = V1(I  )
         V(I+1) = V1(I+1)
         V(I+2) = V1(I+2)
         V(I+3) = V1(I+3)
      END DO
      END
C----------------------------------------
      SUBROUTINE VDOTD_VSMA1D (N,V1,V2,S,
     *                              S1,A)
C----------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V1(*),V2(*),A(*)

      NL = N
      SL = 0D0
      S1L= S1
      NREST = MOD(NL,4)
      IF (NREST .EQ. 0) GO TO 10
      DO I=1,NREST
           SL =  SL  + V2(I)*V1(I)
         A(I) = A(I) +   S1L*V1(I)
      END DO
 10   NRESTP1 = NREST + 1
      DO I=NRESTP1,NL,4
            SL =  SL + V2(I  )*V1(I  )
     *               + V2(I+1)*V1(I+1)
     *               + V2(I+2)*V1(I+2)
     *               + V2(I+3)*V1(I+3)
         A(I  ) = A(I  ) + S1L*V1(I  )
         A(I+1) = A(I+1) + S1L*V1(I+1)
         A(I+2) = A(I+2) + S1L*V1(I+2)
         A(I+3) = A(I+3) + S1L*V1(I+3)
      END DO
      S = SL

      END
C----------------------------------
      SUBROUTINE VSMA1D (N,S1,V1,A)
C----------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V1(*),A(*)
      NL = N
      S1L= S1
      NREST = MOD (NL,4)
      IF (NREST .EQ. 0) GO TO 10
      DO I=1,NREST
         A(I) = A(I) + S1L*V1(I)
      END DO
 10   NRESTP1 = NREST + 1
      DO I=NRESTP1,NL,4
         A(I  ) = A(I  ) + S1L*V1(I  )
         A(I+1) = A(I+1) + S1L*V1(I+1)
         A(I+2) = A(I+2) + S1L*V1(I+2)
         A(I+3) = A(I+3) + S1L*V1(I+3)
      END DO
      END
C----------------------------------------
      SUBROUTINE VSMA2D (N,S1,S2,V1,V2,A)
C----------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V1(*),V2(*),A(*)
      NL = N
      S1L= S1
      S2L= S2
      NREST = MOD (NL,4)
      IF (NREST .EQ. 0) GO TO 10
      DO I=1,NREST
         A(I) = A(I) + S1L*V1(I) + S2L*V2(I)
      END DO
 10   NRESTP1 = NREST + 1
      DO I=NRESTP1,NL,4
         AUX = A(I  ) + S1L*V1(I  )
         BUX = A(I+1) + S1L*V1(I+1)
         CUX = A(I+2) + S1L*V1(I+2)
         DUX = A(I+3) + S1L*V1(I+3)
         A(I  ) = AUX + S2L*V2(I  )
         A(I+1) = BUX + S2L*V2(I+1)
         A(I+2) = CUX + S2L*V2(I+2)
         A(I+3) = DUX + S2L*V2(I+3)
      END DO
      END
