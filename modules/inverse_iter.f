C---------------------------------------------------------------
      SUBROUTINE INVERSE_ITER (N,IV,LV,NVX,EPS,GERSCH,DEL,DELW5,
     *                                           E,V,W1,W2,EPS2,
     *                                   W1I,W7,W5,W3,W4,W5I,W6)
C---------------------------------------------------------------

C     Input:
C     N       order of matrix.
C     IV      first wanted eigenvector.
C     LV       last wanted eigenvector.
C     NVX     left dimension of vector V in subroutine where is
C             to be used as a two-dimensional array.
C     EPS     sqrt(MACHEP) where MACHEP is the machine epsilon
C             (smallest number which makes A+MACHEP different from A).
C     GERSCH  maximum tridiagonal element divided by SCALE, where
C             the latter is a scaling factor equal to the sum of the
C             absolute values of all the elements of the tridiagonal
C             matrix.
C     DEL     GERSCH*EPS
C     DELW5   GERSCH*MACHEP
C     E       array holding all eigenvalues.
C     W1      array holding codiagonal elements of tridiagonal matrix.
C     W2      array holding   diagonal elements of tridiagonal matrix.
C     Output:
C     EPS2    threshold to identify closely spaced eigenvalues.
C     V(*)    array holding eigenvectors of tridiagonal matrix.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (P5=0.5D0, BOUNDR=1D-6, 
     *           GOLDRT=0.618033988749894D0)
      DIMENSION E(*),V(*),W1(*),W2(*),
     *          W1I(*),W7(*),W5(*),W3(*),W4(*),W5I(*),W6(*)

      NM1 = N - 1
      NM2 = N - 2
C     Inverse iteration for eigenvectors.

      DO J=1,N
         IF (W1(J) .EQ. 0D0) W1(J) = DEL
         W1I(J) = 1D0/W1(J)
      END DO

      FN  = DBLE (N)
      EPS1= SQRT (FN)*EPS
      SEPS= SQRT (EPS)
      EPS2=(GERSCH*BOUNDR)/(FN*SEPS)

C     Find if first wanted vector is too close to previous ones.
   10 CONTINUE
      IF (IV .EQ. 1) GO TO 20
      IF (ABS(E(IV-1)-E(IV)) .LT. EPS2) THEN
          IV = IV - 1
          GO TO 10
      END IF
   20 CONTINUE

      RN  = 0D0
      RA  = EPS*GOLDRT
      DO I=IV,LV
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
                W7(J) = -W5(J)*W1I(J)
                W5(J) =  W1(J)
               W5I(J) =  W1I(J)
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
         WN5I   = 1D0/W5(N)

         W6(N  ) = W6(N  )*WN5I
         W6(NM1) =(W6(NM1)-W6(N)*W4(NM1))*WNM15I
              VN = MAX (ABS (W6(N)), ABS (W6(NM1)))
         DO K = NM2,1,-1
            W6(K)=(W6(K)-W6(K+1)*W4(K)-W6(K+2)*W3(K))*W5I(K)
              VN = MAX (ABS (W6(K)), VN)
         END DO
         S = EPS1/VN 
         DO J=1,N
            W6(J) = S*W6(J)
         END DO

         DO J=1,NM1
            IF (W5(J) .EQ. W1(J)) THEN
                    T   = W6(J  )
                W6(J  ) = W6(J+1)
                W6(J+1) = T
            END IF
            W6(J+1) = W6(J+1) + W6(J)*W7(J)
         END DO
         W6(N  ) = W6(N  )*WN5I
         W6(NM1) =(W6(NM1)-W6(N)*W4(NM1))*WNM15I
              VN = MAX (ABS (W6(N)), ABS (W6(NM1)))
         DO K = NM2,1,-1
            W6(K)=(W6(K)-W6(K+1)*W4(K)-W6(K+2)*W3(K))*W5I(K)
              VN = MAX (ABS (W6(K)), VN)
         END DO

         S = EPS1/VN 
         DO J=1,N
            V(J+(I-1)*NVX) = S*W6(J)
         END DO
      END DO
      END
