C-----------------------------------------------------------------
      SUBROUTINE QR (N,IORD,W1,W2,W3,E,EPS,GERSCH,DEL,DELW5,SCALE,
     *                                                        IER)
C-----------------------------------------------------------------

C     Input:
C     N      order of matrix
C     IORD   >0 for nonincreasing order of eigenvalues.
C     IORD   <0 for nondecreasing order of eigenvalues.
C     W1     codiagonal of the tridiagonal matrix.
C     W2     diagonal of the tridiagonal matrix.

C     Output:
C     EPS    sqrt(MACHEP) where MACHEP is the machine epsilon
C            (smallest number which makes A+MACHEP different from A).
C     GERSCH maximum tridiagonal element divided by SCALE.
C     DEL    GERSCH*EPS
C     DELW5  GERSCH*MACHEP
C     SCALE  overall scale factor, equal to the sum of all absolute
C            values of the elements of the tridiagonal matrix.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION MACHEP
      PARAMETER (P5=0.5D0, MACHEP=1.4D-17)
      DIMENSION W1(*),W2(*),W3(*),E(*)

      NM1 = N - 1
      GERSCH  = ABS (W2(1)) + ABS (W1(1))
      DO I=1,NM1
         GERSCH = MAX (ABS(W2(I+1)) + ABS(W1(I)) + ABS(W1(I+1)), GERSCH)
      END DO

C     Trap null matrix before it is too late.

      IF (GERSCH .EQ. 0D0) THEN
          PRINT *,' NULL MATRIX IN SUBROUTINE HQRII1'
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
      EPS   = SQRT (MACHEP)
      GERSCH= GERSCH*SCALEI
      DEL   = GERSCH*EPS
      DELW5 = GERSCH*MACHEP

      DO I=1,N
         W1(I) = W1(I)*SCALEI
         W2(I) = W2(I)*SCALEI
         W3(I) = W1(I)
          E(I) = W2(I)
      END DO

C     QR method with origin shift.

      DO K=N,2,-1
   10    L = K
   20    IF (ABS (W3(L-1)) .GT. DEL) THEN
             L = L - 1
             IF (L .GT. 1) GO TO 20
         END IF
         IF (L .NE. K) THEN
                WWW= (E(K-1)+E(K))*P5
                R  =  E(K)-WWW
                Z  = WWW - SIGN (SQRT (W3(K-1)**2 + R*R),WWW)
                EE = E(L) - Z
                FF = W3(L)
                R  = SQRT (EE*EE + FF*FF)
                RI = 1D0/R
               E(L)= EE
                C  =  E(L)*RI
                S  = W3(L)*RI
                WWW=  E(L+1) - Z
               E(L)= (FF*C + WWW*S)*S + EE + Z
             E(L+1)= C*WWW - S*FF

C            ////////////////
C            Start innerloop.
             DO J=L+1,K-1
                   R  = SQRT (E(J)**2 + W3(J)**2)
                   RI = 1D0/R
                   EE = E(J)*C + Z
                   FF = W3(J)*C
                   WWW= E(J+1) - Z
                   W3(J-1) = S*R
                   C  =  E(J)*RI
                   S  = W3(J)*RI
                  E(J)= (FF*C + WWW*S)*S + EE
                E(J+1)= C*WWW - S*FF
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
   30 L = 1
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
            WWW= E(LL)
          E(LL)= E(II)
          E(II)= WWW
      END IF
      J = II - 1
      IF (J .GT. 1) GO TO 30
      END
