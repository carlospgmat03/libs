C------------------------------------------------------------------
      SUBROUTINE HQRII1 (N,IV,LV,IORD,IUH,IUB,AL_VEC,E,NVX,SCHMIDT,
     *                   QRQL_LAPACK,               IWORK,WORK,IER)
C------------------------------------------------------------------

C  Title: module HQRII1 to diagonalize dense real-symmetric matrices.

C  Abstract: An input real-symmetric matrix AL_VEC is tridiagonalized
C            by the method of Householder.  The eigenvalues of the
C            tridiagonal matrix T are found by the QR method with
C            origin shift.  Optionally, some or all eigenvectors of
C            matrix T are determined by inverse iteration.  The
C            eigenvectors of matrix AL_VEC are found by 
C            left-multiplying the eigenvectors of matrix T times the
C            orthogonal matrix which brings AL_VEC to tridiagonal form,
C            being delivered in AL_VEC assuming a left dimension of NVX.

C  Environment: Standard Fortran 77, 90 and 95 (with obsolescent GO TO).

C  Copyright by Annik Vivier Bunge, Carlos F. Bunge, Yoshitaka Beppu,
C               Ichizo Ninomiya and Zdenko A. Tomasic, 1986, 2000.

C  Version  2-001: August, 2000.
C  Versions 1-001 through 1-005: March 1985 through December 1994.

C  References:  A.V.Bunge and C.F.Bunge, Comput. Chem. 10,259(1986).
C               C. F. Bunge, Comp. Phys. Commun. (2000).

C  Machine dependent parameter:  (in subroutine QR)

C  MACHEP       is the smallest number which makes 1. + MACHEP greater
C               than 1.  For DEC Compaq and Vaxes MACHEP=1.4D-17.
C               For IBM computers MACHEP=1.11D-16.

C  Parameter values:

C  P5           equal to 0.5D0.
C  BOUNDR=1D-6  is a small number involved in deciding whether two
C               eigenvalues are close enough to be considered degenerate
C               while computing eigenvectors of the tridiagonal matrix
C               (in subroutine INVERSE_ITER).
C  GOLDRT       GOLDRT=0.5D0*(SQRT(5D0)-1D0) is a random number between
C               0 and 1.  Other choices may affect the accuracy of 
C               nearly degenerate eigenvectors, for better or for worse,
C               depending on the particular case (in subroutine 
C               INVERSE_ITER).

C  Argument values:

C  N           order of the matrix to be diagonalized in a given run.
C  IV          index of first wanted eigenvector.  Default is IV=1.
C  LV          index of last wanted eigenvector.  If LV .LT. IV, no
C              eigenvectors are calculated.
C  IORD        controls the order of eigenvalues.  IORD .GE. 0 will
C              cause eigenvalues to be given in non-increasing order.
C              IORD .LT. 0 specifies non-decreasing order (used only
C              in subroutine QR).
C  **********  Below, optimal data for Alpha 21264 processors.
C  IUH         Order of hand-unrolling in subroutine HOUSEHOLDER_IUH.
C                  If IUH=0, program chooses best value.
C                     IUH=3 Compaq Fortran for Linux, N.GT.1000, 21264.
C                     IUH=6 Compaq Fortran for Linux, N.LE.1000, 21264.
C                     IUH=8 DEC Unix 4.0F Fortran, 21264A.
C  IUB         Order of hand-unrolling in subroutine 
C                                            BACK_TRANSFORMATION_IUB.
C                  If IUB=0, program chooses best value.
C                     IUB=4 Compaq Fortran for Linux, N.GT.1000, 21264.
C                     IUB=8 Compaq Fortran for Linux, N.LE.1000, 21264.
C                     IUB=4 DEC Unix 4.0F Fortran, N.GT.1000, 21264A.
C                     IUB=8 DEC Unix 4.0F Fortran, N.LE.1000, 21264A.
C  **********  Above, optimal data for Alpha 21264 processors.
C  AL_VEC      On input:
C              real symmetric matrix of dimension at least N1=(N*N+N)/2
C              given in rowwise lower triangular form, viz., AL_VEC(IJ)=
C              AA(J,I) where IJ=J+(I*I-I)/2 and AA is the usual square
C              symmetric form.
C              On output:
C              array of dimension NVX*N which contains the J-th eigen-
C              vector in positions AL_VEC(IJ) starting at 
C              IJ=(J-1)*NVX+1.
C   E          array of dimension at least N holding the eigenvalues
C              in the order determined by variable IORD.
C   NVX        column dimension of array AL_VEC in calling program.  If 
C              AL_VEC is a one dimensional array in the calling program,
C              most likely NVX=N.
C   SCHMIDT    =.TRUE. for normal use.  If SCHMIDT=.FALSE. eigenvectors
C              of degenerate or nearly degenerate eigenvalues are not
C              orthogonalized among themselves.
C   QRQL_LAPACK=.FALSE. for normal use.  If QRQL_LAPACK=.TRUE., the
C              program first attempts to evaluate the eigenvalues by
C              subroutine DSTERF from the LAPACK collection, which may
C              be faster for small matrices and certain compilers, like
C              Linux g77 for Alpha, but not for Pentium!.
C   WORK       working array of dimension 21*N + (N*N+N)/2.
C   IER        is a completion status indicator, see below.

C   Accuracy:
C       for the largest (in absolute value) eigenvalues E(I)
C       which are not highly degenerate, machine accuracy for MR(I),
C                            MR(I) = MODR(I)/E(I)
C       where MODR is the modulus of the residual vector R(I),
C                             R(I) = (A*V)(I) - E(I)*V(I),
C       and V(I) denotes the I-th eigenvector.  Up to M figures may be
C       lost for the smallest eigenvalues (and eigenvectors), where 
C                                M = ALOG10 (LARGEST/SMALLEST) 
C       with LARGEST and SMALLEST being the largest and smallest
C       (in absolute value) eigenvalues, respectively.
C       In presence of highly degenerate eigenvalues, users may enter
C       SCHMIDT=.FALSE., which will produce very accurate eigenvectors,
C       but then it becomes their responsibility to produce an
C       orthogonal set. 

C  Completion status:

C  IER=0     normal successful completion.
C  IER=33    means that input AL_VEC is exactly a null matrix.
C  IER=129   indicates a fatal error due to N, LV or NVX being outside
C            permissible bounds.

C  Caveat: matrices formed exclusively by very small or very large
C          matrix elements will give an arithmetic overflow.

C  Compilation suggestions (things will change with newer versions of
C                           compilers):

C  Compaq machines, DEC Unix 4.0F:
C    f90 or f95 -c -O5 -tune host -arch host -notransform_loops -inline all -unroll 6 $*.f
C    Replace f95 and f90 by fort, for Compaq Fortran for Linux.

C  Cray, UNICOS 8.0:
C     cf77 -Zv -c -Wf"-dp -a static -o inline"

C  Origin 2000 and Silicon Graphics:
C     f90 -c -O3 -mips4 -OPT:IEEE$\_$arithmetic=1:roundoff=0

C  g77 (Intel processors, viz., Pentium).
C     g77 -c -O3 -fforce-mem $*.f

C  g77 (DEC Alpha, etc, other than Intel, under Linux)
C     g77 -c -O3 -funroll-all-loops -fforce-mem $*.f

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION MACHEP
      REAL TIME,CPUT
      LOGICAL SCHMIDT,QRQL_LAPACK
      PARAMETER (P5=0.5D0, MACHEP=1.4D-17)
      DIMENSION AL_VEC(*),E(*),WORK(*),IWORK(*)
      EXTERNAL HOUSEHOLDER, BACK_TRANSFORMATION, TRANSPOSE


       IER = 0
      INFO = 0
C     Choice of IUH and IUB for DEC 21264A processors under DEC Unix 4.0F.
      IF (IUH .EQ. 0  .AND.  IUB .EQ. 0) THEN
           IF (N .GT. 1000) THEN
               IUH = 8
               IUB = 4
           ELSE
               IUH = 8
               IUB = 8
           END IF
      END IF

      IF (IV .LE. 0) IV = 1

      IF (N .LT. 1      .OR.
     *   (N .GT. NVX   .AND.  LV .GE. IV)  .OR.
     *                        LV .GT. N )  THEN
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
      CALL HOUSEHOLDER (IUH,N,AL_VEC,WORK,WORK(N+1),WORK(2*N+1),
     *                                 WORK(3*N+1),IWORK       )
          CALL QR (N,IORD,WORK,WORK(N+1),WORK(2*N+1),
     *              E,EPS,GERSCH,DEL,DELW5,SCALE,IER)
          IF (IER .NE. 0) RETURN

      IF (LV .LT. IV) THEN
          DO I=1,N
             E(I) = SCALE*E(I)
          END DO
C                     //////
                      RETURN  ! Return when eigenvectors not requested.
C                     //////
      END IF

!     CALL TIEMPO (.TRUE.,.TRUE.,' Begin TRANSPOSE.',6,TIME,CPUT)
      CALL TRANSPOSE (N,AL_VEC,IWORK,IWORK(N+1),
     *                                WORK(21*N+1))
!     CALL TIEMPO (.TRUE.,.FALSE.,' End TRANSPOSE.',6,TIME,CPUT)
!     CALL TIEMPO (.TRUE.,.TRUE.,' Begin INVERSE_ITER.',6,TIME,CPUT)
      CALL INVERSE_ITER (N,IV,LV,NVX,EPS,GERSCH,DEL,DELW5,
     *                   E,AL_VEC,WORK,WORK(N+1),EPS2,
     *                   WORK(2*N+1),WORK(3*N+1),WORK(4*N+1),
     *                   WORK(5*N+1),WORK(6*N+1),WORK(7*N+1),
     *                                           WORK(8*N+1))
!     CALL TIEMPO (.TRUE.,.FALSE.,' End INVERSE_ITER.',6,TIME,CPUT)

!     CALL TIEMPO (.TRUE.,.TRUE.,' Begin BACK_TRANSF.',6,TIME,CPUT)
      CALL BACK_TRANSFORMATION (IUB,N,IV,LV,SCHMIDT,EPS2,
     *                          WORK(21*N+1),E,NVX,AL_VEC,IWORK,WORK)
!     CALL TIEMPO (.TRUE.,.FALSE.,' End BACK_TRANSF.',6,TIME,CPUT)

      DO I=1,N
         E(I) = SCALE*E(I)
      END DO

      END
C-----------------------------------------------------
      SUBROUTINE HOUSEHOLDER (IUH,N,AL,W1,W2,W5,W6,IX)
C-----------------------------------------------------

C  Purpose: to select appropriate hand unrolling in tridiagonalization
C           subroutine HOUSEHOLDER_Un.
C  Working arrays:
C  Wn   n=1,2,5,6, are four working vectors each of dimension N_NL.
C  IX   is a working array of dimension N_NL, first defined by IX(I)=
C       (I*I-I)/2 so that the lower triangle matrix AL is indexed by
C       IJ=IX(I)+J.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL TIME,CPUT
      DIMENSION AL(*),W1(*),W2(*),W5(*),W6(*),IX(*)
      EXTERNAL HOUSEHOLDER_U1, HOUSEHOLDER_U2, HOUSEHOLDER_U3,
     *         HOUSEHOLDER_U4, HOUSEHOLDER_U6, HOUSEHOLDER_U8
 
C      CALL TIEMPO (.TRUE.,.TRUE.,'         Begin tridiag.',6,TIME,CPUT)
      IX(1) = 0
      IF (IUH .EQ. 8) THEN
          DO J=2,N
             IX(J) = IX(J-1) + N
          END DO
      ELSE
          DO J=2,N
             IX(J) = IX(J-1) + J - 1
          END DO
      END IF

      IF (IUH .LE. 0) IUH = 1
      IF (IUH .GT. 8) IUH = 8
      SELECT CASE (IUH)
             CASE (1)
             CALL HOUSEHOLDER_U1 (N,AL,W1,W2,W5,W6,IX)
             CASE (2)
             CALL HOUSEHOLDER_U2 (N,AL,W1,W2,W5,W6,IX)
             CASE (3)
             CALL HOUSEHOLDER_U3 (N,AL,W1,W2,W5,W6,IX)
             CASE (4)
             CALL HOUSEHOLDER_U4 (N,AL,W1,W2,W5,W6,IX)
             CASE (5:6)
             CALL HOUSEHOLDER_U6 (N,AL,W1,W2,W5,W6,IX)
             CASE (7:8)
             CALL HOUSEHOLDER_U8 (N,AL,W1,W2,W5,W6,IX)
      END SELECT
C      CALL TIEMPO (.TRUE.,.FALSE.,'           End tridiag.',6,TIME,CPUT)
      END
C-------------------------------------------------------
      SUBROUTINE HOUSEHOLDER_U1 (N_NL,AL,W1,W2,W5,W6,IX)
C-------------------------------------------------------

C  Purpose: to tridiagonalize matrix AL, hand unrolling=1.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (P5=0.5D0)
      DIMENSION AL(*),W1(*),W2(*),W5(*),W6(*),IX(*)

      N = N_NL
      NM1 = N - 1
      IF (N .GT. 2) THEN
C         Householder transformation.
          DO K=1,N-2
              KP1  = K + 1
             W2(K) = AL(K+IX(K))
             SCALE = 0D0
             DO J=KP1,N
                SCALE = SCALE + ABS (AL(IX(J)+K))
             END DO
             W1(K) = AL(IX(KP1)+K)
             IF (SCALE .GT. 0D0) THEN
                 SCALEI = 1D0/SCALE
                 DO J=KP1,N
                    W2(J) = AL(IX(J)+K)*SCALEI
                 END DO
                 SUM  = 0D0
                 DO J=KP1,N
                    SUM  = SUM + W2(J)**2
                 END DO
                 S = SIGN (SQRT (SUM),W2(KP1))
                 W1(K  ) = -S*SCALE
                 W2(KP1) = W2(KP1) + S
                 AL(IX(KP1)+K) = W2(KP1)*SCALE
                 H    = W2(KP1)*S
                 HUNS = (H*SCALE)*SCALE
                 HI   = 1D0/H
                 DO II=KP1,N
                    W5(II) = 0D0
                 END DO

C                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C                First numerically intensive kernel.
                 DO I=KP1,N
                    IM1  = I - 1
                    SUM  = 0D0
                    I0   = IX(I)
                    W2I  = W2(I)

C                   ////////////////////
C                   Start innermostloop.
!DIR$ IVDEP
!DIR$LOOP COUNT (4000)
                    DO J=KP1,IM1
                       SUM   =  SUM + W2(J)*AL(J+I0)
                       W5(J) =  W5(J) + W2I*AL(J+I0)
                    END DO
C                   End innermostloop.
C                   //////////////////

                    W6(I) = W2I*AL(I+I0) + SUM
                 END DO
C                End first numerically intensive kernel.
C                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

                 SUM = 0D0
                 DO I=KP1,N
                    W56 = W5(I) + W6(I)
                    W1(I) = HI*W56
                    SUM = SUM + HI*W56*W2(I)
                 END DO

                 U = P5*SUM*HI
                 DO I=KP1,N
                    W1(I) = W2(I)*U - W1(I)
                 END DO
C                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C                Second numerically intensive kernel.
                 DO I=N,KP1,-1
                    I0    = IX(I)
                    W1I   = W1(I)
                    W2I   = W2(I)

C                   ////////////////////
C                   Start innermostloop.
!DIR$ IVDEP
!DIR$LOOP COUNT (4000)
                    JREST = MOD(I-KP1+1,4)
                    DO J=KP1,KP1+JREST-1
                       AL(J+I0) = AL(J+I0) + W2I*W1(J) + W1I*W2(J)
                    END DO
                    DO J=KP1+JREST,I,4
                     AL(J  +I0) = AL(J  +I0) + W2I*W1(J  ) + W1I*W2(J  )
                     AL(J+1+I0) = AL(J+1+I0) + W2I*W1(J+1) + W1I*W2(J+1)
                     AL(J+2+I0) = AL(J+2+I0) + W2I*W1(J+2) + W1I*W2(J+2)
                     AL(J+3+I0) = AL(J+3+I0) + W2I*W1(J+3) + W1I*W2(J+3)
                    END DO
C                   End innermostloop.
C                   //////////////////

                 END DO
C                End second numerically intensive kernel.
C                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
             ELSE
                 HUNS = 0D0
             END IF
             AL(K+IX(K)) = HUNS
          END DO
      END IF
      NM1NM1  = IX(NM1) + NM1
      NM1N    = IX(N  ) + NM1
      NN      = NM1N    + 1
      W2(NM1) = AL(NM1NM1)
      W2(N  ) = AL(NN)
      W1(NM1) = AL(NM1N)
      W1(N  ) = 0D0
      END
C-------------------------------------------------------
      SUBROUTINE HOUSEHOLDER_U2 (N_NL,AL,W1,W2,W5,W6,IX)
C-------------------------------------------------------

C  Purpose: to tridiagonalize matrix AL, hand unrolling=2.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (P5=0.5D0)
      DIMENSION AL(*),W1(*),W2(*),W5(*),W6(*),IX(*)

      N = N_NL
      NM1 = N - 1
      IF (N .GT. 2) THEN
C         Householder transformation.
          DO K=1,N-2
              KP1  = K + 1
             W2(K) = AL(K+IX(K))
             SCALE = 0D0
             DO J=KP1,N
                SCALE = SCALE + ABS (AL(IX(J)+K))
             END DO
             W1(K) = AL(IX(KP1)+K)
             IF (SCALE .GT. 0D0) THEN
                 SCALEI = 1D0/SCALE
                 DO J=KP1,N
                    W2(J) = AL(IX(J)+K)*SCALEI
                 END DO
                 SUM  = 0D0
                 DO J=KP1,N
                    SUM  = SUM + W2(J)**2
                 END DO
                 S = SIGN (SQRT (SUM),W2(KP1))
                 W1(K  ) = -S*SCALE
                 W2(KP1) = W2(KP1) + S
                 AL(IX(KP1)+K) = W2(KP1)*SCALE
                 H    = W2(KP1)*S
                 HUNS = (H*SCALE)*SCALE
                 HI   = 1D0/H
                 DO II=KP1,N
                    W5(II) = 0D0
                 END DO
                 NREST = MOD(N-KP1+1,2)
                 DO I=KP1,KP1+NREST-1
                    IM1  = I - 1
                    SUM  = 0D0
                    I0   = IX(I)
                    W2I  = W2(I)
                    DO J=KP1,IM1
                       SUM   =  SUM + W2(J)*AL(J+I0)
                       W5(J) =  W5(J) + W2I*AL(J+I0)
                    END DO
                    W6(I) = W2I*AL(I+I0) + SUM
                 END DO

C                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C                First numerically intensive kernel.
                 DO I=KP1+NREST,N,2
                    IM1  = I - 1
                    I0   = IX(I  )
                    I1   = IX(I+1)
                    W2I  = W2(I  )
                    W2I1 = W2(I+1)
                    SUM0 = 0D0
                    SUM1 = 0D0

C                   ////////////////////
C                   Start innermostloop.
!DIR$ IVDEP
!DIR$LOOP COUNT (4000)
                    DO J=KP1,IM1
                        SUM0 = SUM0 + W2(J)*AL(J+I0)
                        SUM1 = SUM1 + W2(J)*AL(J+I1)
                       W5(J) = W5(J) + W2I *AL(J+I0)
     *                               + W2I1*AL(J+I1)
                    END DO
C                   End innermostloop.
C                   //////////////////

                    W5(I  ) = W5(I) + W2I1*AL(I+I1  )
                    W6(I  ) =  SUM0 + W2I *AL(I+I0  )
                    W6(I+1) =  SUM1 + W2I *AL(I+I1  ) + W2I1*AL(I+I1+1)
                 END DO
C                End first numerically intensive kernel.
C                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

                 SUM = 0D0
                 DO I=KP1,N
                    W56 = W5(I) + W6(I)
                    W1(I) = HI*W56
                    SUM = SUM + HI*W56*W2(I)
                 END DO

                 U = P5*SUM*HI
                 NREST = MOD(N-KP1+1,2)
                 DO I=KP1,KP1+NREST-1
                    I0    = IX(I)
                    W1(I) = W2(I)*U - W1(I)
                    W1I   = W1(I)
                    W2I   = W2(I)
                    DO J=KP1,I
                       AL(J+I0) = AL(J+I0) + W2I*W1(J) + W1I*W2(J)
                    END DO
                 END DO
                 DO I=KP1+NREST,N,2
                    W1(I  ) = W2(I  )*U - W1(I  )
                    W1(I+1) = W2(I+1)*U - W1(I+1)
                 END DO

C                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C                Second numerically intensive kernel.
                 DO I=N-1,KP1+NREST,-2
                    I0    = IX(I  )
                    I1    = IX(I+1)
                    W2I   = W2(I  )
                    W2I1  = W2(I+1)
                    W1I   = W1(I  )
                    W1I1  = W1(I+1)

C                   ////////////////////
C                   Start innermostloop.
!DIR$ IVDEP
!DIR$LOOP COUNT (4000)
                    DO J=KP1,I
!                      AL(J+I0) = AL(J+I0) + W2I *W1(J) + W1I *W2(J)
!                      AL(J+I1) = AL(J+I1) + W2I1*W1(J) + W1I1*W2(J)
                       TEMP     = AL(J+I0) + W2I *W1(J) + W1I *W2(J)
                       TEMP1    = AL(J+I1) + W2I1*W1(J) + W1I1*W2(J)
                       AL(J+I0) = TEMP
                       AL(J+I1) = TEMP1
                    END DO 
C                   End innermostloop.
C                   //////////////////

                   AL(I+I1+1) = AL(I+I1+1) + W2I1*W1(I+1) + W1I1*W2(I+1)
                 END DO
C                End second numerically intensive kernel.
C                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

             ELSE
                 HUNS = 0D0
             END IF
             AL(K+IX(K)) = HUNS
          END DO
      END IF
      NM1NM1  = IX(NM1) + NM1
      NM1N    = IX(N  ) + NM1
      NN      = NM1N    + 1
      W2(NM1) = AL(NM1NM1)
      W2(N  ) = AL(NN)
      W1(NM1) = AL(NM1N)
      W1(N  ) = 0D0
      END
C-------------------------------------------------------
      SUBROUTINE HOUSEHOLDER_U3 (N_NL,AL,W1,W2,W5,W6,IX)
C-------------------------------------------------------

C  Purpose: to tridiagonalize matrix AL, hand unrolling=3.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (P5=0.5D0)
      DIMENSION AL(*),W1(*),W2(*),W5(*),W6(*),IX(*)

      N = N_NL
      NM1 = N - 1
      IF (N .GT. 2) THEN
C         Householder transformation.
          DO K=1,N-2
              KP1  = K + 1
             W2(K) = AL(K+IX(K))
             SCALE = 0D0
             DO J=KP1,N
                SCALE = SCALE + ABS (AL(IX(J)+K))
             END DO
             W1(K) = AL(IX(KP1)+K)
             IF (SCALE .GT. 0D0) THEN
                 SCALEI = 1D0/SCALE
                 DO J=KP1,N
                    W2(J) = AL(IX(J)+K)*SCALEI
                 END DO
                 SUM  = 0D0
                 DO J=KP1,N
                    SUM  = SUM + W2(J)**2
                 END DO
                 S = SIGN (SQRT (SUM),W2(KP1))
                 W1(K  ) = -S*SCALE
                 W2(KP1) = W2(KP1) + S
                 AL(IX(KP1)+K) = W2(KP1)*SCALE
                 H    = W2(KP1)*S
                 HUNS = (H*SCALE)*SCALE
                 HI   = 1D0/H
                 DO II=KP1,N
                    W5(II) = 0D0
                 END DO
                 NREST = MOD(N-KP1+1,3)
                 DO I=KP1,KP1+NREST-1
                    IM1  = I - 1
                    SUM  = 0D0
                    I0   = IX(I)
                    W2I  = W2(I)
                    DO J=KP1,IM1
                       SUM   =  SUM + W2(J)*AL(J+I0)
                       W5(J) =  W5(J) + W2I*AL(J+I0)
                    END DO
                    W6(I) = W2I*AL(I+I0) + SUM
                 END DO

C                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C                First numerically intensive kernel.
                 DO I=KP1+NREST,N,3
                    IM1  = I - 1
                    SUM0 = 0D0
                    SUM1 = 0D0
                    SUM2 = 0D0
                    I0   = IX(I  )
                    I1   = IX(I+1)
                    I2   = IX(I+2)
                    W2I  = W2(I  )
                    W2I1 = W2(I+1)
                    W2I2 = W2(I+2)

C                   ////////////////////
C                   Start innermostloop.
!DIR$ IVDEP
!DIR$LOOP COUNT (4000)
                    DO J=KP1,IM1
                         AUX = W5(J) + W2I *AL(J+I0)
                        SUM0 = SUM0 + W2(J)*AL(J+I0)
                        SUM1 = SUM1 + W2(J)*AL(J+I1)
                        SUM2 = SUM2 + W2(J)*AL(J+I2)
                       W5(J) =  AUX  + W2I1*AL(J+I1)
     *                               + W2I2*AL(J+I2)
                    END DO
C                   End innermostloop.
C                   //////////////////

                       SUM1 =    SUM1 + W2I *AL(I+I1  )
                       SUM2 =    SUM2 + W2I *AL(I+I2  )
     *                                + W2I1*AL(I+I2+1)
                    W5(I  ) = W5(I  ) + W2I1*AL(I+I1  )
     *                                + W2I2*AL(I+I2  )
                    W5(I+1) = W5(I+1) + W2I2*AL(I+I2+1)
                    W6(I  ) =    SUM0 + W2I *AL(I+I0  )
                    W6(I+1) =    SUM1 + W2I1*AL(I+I1+1)
                    W6(I+2) =    SUM2 + W2I2*AL(I+I2+2)
                 END DO
C                End first numerically intensive kernel.
C                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

                 SUM = 0D0
                 DO I=KP1,N
                    W56 = W5(I) + W6(I)
                    W1(I) = HI*W56
                    SUM = SUM + HI*W56*W2(I)
                 END DO

                 U = P5*SUM*HI
                 NREST = MOD(N-KP1+1,3)
                 DO I=KP1,KP1+NREST-1
                    I0    = IX(I)
                    W1(I) = W2(I)*U - W1(I)
                    W1I   = W1(I)
                    W2I   = W2(I)
                    DO J=KP1,I
                       AL(J+I0) = AL(J+I0) + W2I*W1(J) + W1I*W2(J)
                    END DO
                 END DO
                 DO I=KP1+NREST,N,3
                    W1(I  ) = W2(I  )*U - W1(I  )
                    W1(I+1) = W2(I+1)*U - W1(I+1)
                    W1(I+2) = W2(I+2)*U - W1(I+2)
                 END DO

C                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C                Second numerically intensive kernel.
                 DO I=N-2,KP1+NREST,-3
                    I0    = IX(I  )
                    I1    = IX(I+1)
                    I2    = IX(I+2)
                    W2I   = W2(I  )
                    W2I1  = W2(I+1)
                    W2I2  = W2(I+2)
                    W1I   = W1(I  )
                    W1I1  = W1(I+1)
                    W1I2  = W1(I+2)

C                   ////////////////////
C                   Start innermostloop.
                    JREST = MOD(I-KP1+1,4)
                    DO J=KP1,KP1+JREST-1
                       AL(J+I0) = AL(J+I0) + W2I *W1(J) + W1I *W2(J)
                       AL(J+I1) = AL(J+I1) + W2I1*W1(J) + W1I1*W2(J)
                       AL(J+I2) = AL(J+I2) + W2I2*W1(J) + W1I2*W2(J)
                    END DO 
!DIR$ IVDEP
!DIR$LOOP COUNT (4000)
                    DO J=KP1+JREST,I,4
!!!                    TEMP  = AL(J+I0) + W2I *W1(J) + W1I *W2(J)
!!!                    TEMP1 = AL(J+I1) + W2I1*W1(J) + W1I1*W2(J)
!!!                    AL(J+I2) = AL(J+I2) + W2I2*W1(J) + W1I2*W2(J)
!!!                    AL(J+I0) = TEMP
!!!                    AL(J+I1) = TEMP1
                   TEMP  = AL(J  +I0) + W2I *W1(J  ) + W1I *W2(J  )
                   TEMP1 = AL(J+1+I0) + W2I *W1(J+1) + W1I *W2(J+1)
                   TEMP2 = AL(J+2+I0) + W2I *W1(J+2) + W1I *W2(J+2)
                   TEMP3 = AL(J+3+I0) + W2I *W1(J+3) + W1I *W2(J+3)
                   AL(J  +I0) = TEMP
                   AL(J+1+I0) = TEMP1
                   AL(J+2+I0) = TEMP2
                   AL(J+3+I0) = TEMP3
                   UEMP  = AL(J  +I1) + W2I1*W1(J  ) + W1I1*W2(J  )
                   UEMP1 = AL(J+1+I1) + W2I1*W1(J+1) + W1I1*W2(J+1)
                   UEMP2 = AL(J+2+I1) + W2I1*W1(J+2) + W1I1*W2(J+2)
                   UEMP3 = AL(J+3+I1) + W2I1*W1(J+3) + W1I1*W2(J+3)
                   AL(J  +I1) = UEMP
                   AL(J+1+I1) = UEMP1
                   AL(J+2+I1) = UEMP2
                   AL(J+3+I1) = UEMP3
                   VEMP  = AL(J  +I2) + W2I2*W1(J  ) + W1I2*W2(J  )
                   VEMP1 = AL(J+1+I2) + W2I2*W1(J+1) + W1I2*W2(J+1)
                   VEMP2 = AL(J+2+I2) + W2I2*W1(J+2) + W1I2*W2(J+2)
                   VEMP3 = AL(J+3+I2) + W2I2*W1(J+3) + W1I2*W2(J+3)
                   AL(J  +I2) = VEMP
                   AL(J+1+I2) = VEMP1
                   AL(J+2+I2) = VEMP2
                   AL(J+3+I2) = VEMP3
                    END DO 
C                   End innermostloop.
C                   //////////////////

                   AL(I+I1+1) = AL(I+I1+1) + W2I1*W1(I+1) + W1I1*W2(I+1)
                   AL(I+I2+1) = AL(I+I2+1) + W2I2*W1(I+1) + W1I2*W2(I+1)
                   AL(I+I2+2) = AL(I+I2+2) + W2I2*W1(I+2) + W1I2*W2(I+2)
                 END DO
C                End second numerically intensive kernel.
C                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
             ELSE
                 HUNS = 0D0
             END IF
             AL(K+IX(K)) = HUNS
          END DO
      END IF
      NM1NM1  = IX(NM1) + NM1
      NM1N    = IX(N  ) + NM1
      NN      = NM1N    + 1
      W2(NM1) = AL(NM1NM1)
      W2(N  ) = AL(NN)
      W1(NM1) = AL(NM1N)
      W1(N  ) = 0D0
      END
C-------------------------------------------------------
      SUBROUTINE HOUSEHOLDER_U4 (N_NL,AL_NL,
     *                       W1_NL,W2_NL,W5_NL,W6_NL,IX)
C-------------------------------------------------------

C  Purpose: to tridiagonalize matrix AL, hand unrolling=4.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (P5=0.5D0)
      DIMENSION AL_NL(*),W1_NL(*),W2_NL(*),W5_NL(*),W6_NL(*),IX(*)
      ALLOCATABLE :: W1(:), PAD1(:), W2(:), PADA(:), AL(:), PADB(:),
     *               W5(:), PAD2(:), W6(:)

      N = N_NL
      
      ALLOCATE (W1(N), PAD1(37), W2(N), 
     *          PADA(41), AL(N*N), PADB(47),
     *          W5(N), PAD2(53), W6(N))

      DO I=1,N*N
         AL(I) = AL_NL(I)
      END DO

      NM1 = N - 1
      IF (N .GT. 2) THEN
C         Householder transformation.
          DO K=1,N-2
              KP1  = K + 1
!!!!         JREST = MOD(K,4)
             W2(K) = AL(K+IX(K))
             SCALE = 0D0
             DO J=KP1,N
                SCALE = SCALE + ABS (AL(IX(J)+K))
             END DO
             W1(K) = AL(IX(KP1)+K)
             IF (SCALE .GT. 0D0) THEN
                 SCALEI = 1D0/SCALE
                 DO J=KP1,N
                    W2(J) = AL(IX(J)+K)*SCALEI
                 END DO
                 SUM  = 0D0
                 DO J=KP1,N
                    SUM  = SUM + W2(J)**2
                 END DO
                 S = SIGN (SQRT (SUM),W2(KP1))
                 W1(K  ) = -S*SCALE
                 W2(KP1) = W2(KP1) + S
                 AL(IX(KP1)+K) = W2(KP1)*SCALE
                 H    = W2(KP1)*S
                 HUNS = (H*SCALE)*SCALE
                 HI   = 1D0/H
                 DO II=KP1,N
                    W5(II) = 0D0
                 END DO
                 NREST = MOD(N-KP1+1,4)
                 DO I=KP1,KP1+NREST-1
                    IM1  = I - 1
                    SUM  = 0D0
                    I0   = IX(I)
                    W2I  = W2(I)
                    DO J=KP1,IM1
                       SUM   =  SUM + W2(J)*AL(J+I0)
                       W5(J) =  W5(J) + W2I*AL(J+I0)
                    END DO
                    W6(I) = W2I*AL(I+I0) + SUM
                 END DO

C                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C                First numerically intensive kernel.
!DIR$IVDEP:LOOP
                 DO I=KP1+NREST,N,4
                    IM1  = I - 1
                    SUM0 = 0D0
                    SUM1 = 0D0
                    SUM2 = 0D0
                    SUM3 = 0D0
                    I0   = IX(I  )
                    I1   = IX(I+1)
                    I2   = IX(I+2)
                    I3   = IX(I+3)
                    W2I  = W2(I  )
                    W2I1 = W2(I+1)
                    W2I2 = W2(I+2)
                    W2I3 = W2(I+3) 

C                   ////////////////////
C                   Start innermostloop.
!DIR$ IVDEP
!DIR$LOOP COUNT (4000)
                    DO J=KP1,IM1
                        SUM0 =  SUM0 + W2(J)*AL(J+I0)
                        SUM2 =  SUM2 + W2(J)*AL(J+I2)
                        SUM1 =  SUM1 + W2(J)*AL(J+I1)
                        SUM3 =  SUM3 + W2(J)*AL(J+I3)
                       W5(J) = W5(J) + W2I  *AL(J+I0)
     *                               + W2I2 *AL(J+I2)
     *                               + W2I1 *AL(J+I1)
     *                               + W2I3 *AL(J+I3)
                    END DO
C                   End innermostloop.
C                   //////////////////

                       SUM1 =    SUM1 + W2I *AL(I+I1  )
                       SUM2 =    SUM2 + W2I *AL(I+I2  )
     *                                + W2I1*AL(I+I2+1)
                       SUM3 =    SUM3 + W2I *AL(I+I3  )
     *                                + W2I1*AL(I+I3+1)
     *                                + W2I2*AL(I+I3+2)
                    W5(I  ) = W5(I  ) + W2I1*AL(I+I1  )
     *                                + W2I2*AL(I+I2  )
     *                                + W2I3*AL(I+I3  )
                    W5(I+1) = W5(I+1) + W2I2*AL(I+I2+1)
     *                                + W2I3*AL(I+I3+1)
                    W5(I+2) = W5(I+2) + W2I3*AL(I+I3+2)
                    W6(I  ) =    SUM0 + W2I *AL(I+I0  )
                    W6(I+1) =    SUM1 + W2I1*AL(I+I1+1)
                    W6(I+2) =    SUM2 + W2I2*AL(I+I2+2)
                    W6(I+3) =    SUM3 + W2I3*AL(I+I3+3)
                 END DO
C                End first numerically intensive kernel.
C                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

                 SUM = 0D0
                 DO I=KP1,N
                    W56 = W5(I) + W6(I)
                    W1(I) = HI*W56
                    SUM = SUM + HI*W56*W2(I)
                 END DO

                 U = P5*SUM*HI
                 NREST = MOD(N-KP1+1,4)
                 DO I=KP1,KP1+NREST-1
                    I0    = IX(I)
                    W1(I) = W2(I)*U - W1(I)
                    W1I   = W1(I)
                    W2I   = W2(I)
                    DO J=KP1,I
                       AL(J+I0) = AL(J+I0) + W2I*W1(J) + W1I*W2(J)
                    END DO
                 END DO
                 DO I=KP1+NREST,N,4
                    W1(I  ) = W2(I  )*U - W1(I  )
                    W1(I+1) = W2(I+1)*U - W1(I+1)
                    W1(I+2) = W2(I+2)*U - W1(I+2)
                    W1(I+3) = W2(I+3)*U - W1(I+3)
                 END DO

C                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C                Second numerically intensive kernel.
!!!!             DO I=KP1+NREST,N-3,4

!DIR$IVDEP:LOOP
CCC              DO I=N-3,KP1+NREST-2,-4
                 DO I=N-3,KP1+NREST,-4

CCCC             DO IZ=N-3,KP1+NREST,-4
CCCC                I_FIN = MAX(KP1+NREST,IZ-1)
CCCC             DO JZ=KP1,IZ,255
CCCC             DO I=IZ,I_FIN,-4

                    I0    = IX(I  )
                    I1    = IX(I+1)
                    I2    = IX(I+2)
                    I3    = IX(I+3)
                    W2I   = W2(I  )
                    W2I1  = W2(I+1)
                    W2I2  = W2(I+2)
                    W2I3  = W2(I+3)
                    W1I   = W1(I  )
                    W1I1  = W1(I+1)
                    W1I2  = W1(I+2)
                    W1I3  = W1(I+3)

C                   ////////////////////
C                   Start innermostloop.
CCCC                J_FIN = MIN(I,JZ+254)
CCCC                DO J=JZ,J_FIN

!!!!                JFIN = MIN(KP1+JREST-1,I)
!!!!                DO J=KP1,JFIN
!!!!                   AL(J+I0) = AL(J+I0) + W2I *W1(J) + W1I *W2(J)
!!!!                   AL(J+I1) = AL(J+I1) + W2I1*W1(J) + W1I1*W2(J)
!!!!                   AL(J+I2) = AL(J+I2) + W2I2*W1(J) + W1I2*W2(J)
!!!!                   AL(J+I3) = AL(J+I3) + W2I3*W1(J) + W1I3*W2(J)
!!!!                END DO
!DIR$IVDEP
!!!!                DO J=JFIN+1,I
!DIR$ IVDEP:LOOP
!DIR$LOOP COUNT (4000)
                    DO J=KP1,I
                       TEMP     = AL(J+I0) + W2I *W1(J) + W1I *W2(J)
!                      TEMP2    = AL(J+I2) + W2I2*W1(J) + W1I2*W2(J)
                       AL(J+I2) = AL(J+I2) + W2I2*W1(J) + W1I2*W2(J)
                       TEMP1    = AL(J+I1) + W2I1*W1(J) + W1I1*W2(J)
                       AL(J+I3) = AL(J+I3) + W2I3*W1(J) + W1I3*W2(J)
                       AL(J+I0) = TEMP
!                      AL(J+I2) = TEMP2
                       AL(J+I1) = TEMP1
!                      AL(J+I0) = AL(J+I0) + W2I *W1(J) + W1I *W2(J)
!                      AL(J+I1) = AL(J+I1) + W2I1*W1(J) + W1I1*W2(J)
!                      AL(J+I2) = AL(J+I2) + W2I2*W1(J) + W1I2*W2(J)
!                      AL(J+I3) = AL(J+I3) + W2I3*W1(J) + W1I3*W2(J)
CCCC                   AUX = AL(J+I0) + W2I *W1(J)
CCCC                   BUX = AL(J+I1) + W2I1*W1(J)
CCCC                   CUX = AL(J+I2) + W2I2*W1(J)
CCCC                   DUX = AL(J+I3) + W2I3*W1(J)
CCCC                   AL(J+I0) = AUX + W1I *W2(J)
CCCC                   AL(J+I1) = BUX + W1I1*W2(J)
CCCC                   AL(J+I2) = CUX + W1I2*W2(J)
CCCC                   AL(J+I3) = DUX + W1I3*W2(J)
                    END DO 
C                   End innermostloop.
C                   //////////////////

CCCC                IF (I .EQ. J_FIN) THEN
                   AL(I+I1+1) = AL(I+I1+1) + W2I1*W1(I+1) + W1I1*W2(I+1)
                   AL(I+I2+1) = AL(I+I2+1) + W2I2*W1(I+1) + W1I2*W2(I+1)
                   AL(I+I2+2) = AL(I+I2+2) + W2I2*W1(I+2) + W1I2*W2(I+2)
                   AL(I+I3+1) = AL(I+I3+1) + W2I3*W1(I+1) + W1I3*W2(I+1)
                   AL(I+I3+2) = AL(I+I3+2) + W2I3*W1(I+2) + W1I3*W2(I+2)
                   AL(I+I3+3) = AL(I+I3+3) + W2I3*W1(I+3) + W1I3*W2(I+3)
CCCC                END IF
CCCC             END DO
CCCC             END DO
                 END DO
C                End second numerically intensive kernel.
C                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
             ELSE
                 HUNS = 0D0
             END IF
             AL(K+IX(K)) = HUNS
          END DO
      END IF
      NM1NM1  = IX(NM1) + NM1
      NM1N    = IX(N  ) + NM1
      NN      = NM1N    + 1
      W2(NM1) = AL(NM1NM1)
      W2(N  ) = AL(NN)
      W1(NM1) = AL(NM1N)
      W1(N  ) = 0D0
      DO I=1,N
         W1_NL(I) = W1(I)
         W2_NL(I) = W2(I)
      END DO
      DO I=1,N*N
         AL_NL(I) = AL(I)
      END DO

      DEALLOCATE (W1,PAD1,W2,AL,W5,PAD2,W6)
      END
C-------------------------------------------------------
      SUBROUTINE HOUSEHOLDER_U6 (N_NL,AL,W1,W2,W5,W6,IX)
C-------------------------------------------------------

C  Purpose: to tridiagonalize matrix AL, hand unrolling=6.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (P5=0.5D0)
      DIMENSION AL(*),W1(*),W2(*),W5(*),W6(*),IX(*)

      N = N_NL
      NM1 = N - 1
      IF (N .GT. 2) THEN
C         Householder transformation.
          DO K=1,N-2
              KP1  = K + 1
             W2(K) = AL(K+IX(K))
             SCALE = 0D0
             DO J=KP1,N
                SCALE = SCALE + ABS (AL(IX(J)+K))
             END DO
             W1(K) = AL(IX(KP1)+K)
             IF (SCALE .GT. 0D0) THEN
                 SCALEI = 1D0/SCALE
                 DO J=KP1,N
                    W2(J) = AL(IX(J)+K)*SCALEI
                 END DO
                 SUM  = 0D0
                 DO J=KP1,N
                    SUM  = SUM + W2(J)**2
                 END DO
                 S = SIGN (SQRT (SUM),W2(KP1))
                 W1(K  ) = -S*SCALE
                 W2(KP1) = W2(KP1) + S
                 AL(IX(KP1)+K) = W2(KP1)*SCALE
                 H    = W2(KP1)*S
                 HUNS = (H*SCALE)*SCALE
                 HI   = 1D0/H
                 DO II=KP1,N
                    W5(II) = 0D0
                 END DO
                 NREST = MOD(N-KP1+1,6)
                 DO I=KP1,KP1+NREST-1
                    IM1  = I - 1
                    SUM  = 0D0
                    I0   = IX(I)
                    W2I  = W2(I)
                    DO J=KP1,IM1
                       SUM   =  SUM + W2(J)*AL(J+I0)
                       W5(J) =  W5(J) + W2I*AL(J+I0)
                    END DO
                    W6(I) = W2I*AL(I+I0) + SUM
                 END DO

C                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C                First numerically intensive kernel.
!DIR$IVDEP
                 DO I=KP1+NREST,N,6
                    IM1  = I - 1
                    SUM0 = 0D0
                    SUM1 = 0D0
                    SUM2 = 0D0
                    SUM3 = 0D0
                    SUM4 = 0D0
                    SUM5 = 0D0
                    I0   = IX(I  )
                    I1   = IX(I+1)
                    I2   = IX(I+2)
                    I3   = IX(I+3)
                    I4   = IX(I+4)
                    I5   = IX(I+5)
                    W2I  = W2(I  )
                    W2I1 = W2(I+1)
                    W2I2 = W2(I+2)
                    W2I3 = W2(I+3) 
                    W2I4 = W2(I+4) 
                    W2I5 = W2(I+5) 

C                   ////////////////////
C                   Start innermostloop.
!DIR$ IVDEP:LOOP
!DIR$LOOP COUNT (4000)
                    DO J=KP1,IM1
                       W5(J) = W5(J) + W2I  *AL(J+I0)
     *                               + W2I1 *AL(J+I1)
     *                               + W2I2 *AL(J+I2)
                        SUM0 =  SUM0 + W2(J)*AL(J+I0)
                        SUM1 =  SUM1 + W2(J)*AL(J+I1)
                        SUM2 =  SUM2 + W2(J)*AL(J+I2)
                    END DO
!DIR$ IVDEP:LOOP
!DIR$LOOP COUNT (4000)
                    DO J=KP1,IM1
                        SUM3 =  SUM3 + W2(J)*AL(J+I3)
                        SUM4 =  SUM4 + W2(J)*AL(J+I4)
                        SUM5 =  SUM5 + W2(J)*AL(J+I5)
                       W5(J) = W5(J) + W2I3 *AL(J+I3)
     *                               + W2I4 *AL(J+I4)
     *                               + W2I5 *AL(J+I5)
                    END DO
C                   End innermostloop.
C                   //////////////////

                       SUM1 =    SUM1 + W2I *AL(I+I1  )
                       SUM2 =    SUM2 + W2I *AL(I+I2  )
     *                                + W2I1*AL(I+I2+1)
                       SUM3 =    SUM3 + W2I *AL(I+I3  )
     *                                + W2I1*AL(I+I3+1)
     *                                + W2I2*AL(I+I3+2)
                       SUM4 =    SUM4 + W2I *AL(I+I4  )
     *                                + W2I1*AL(I+I4+1)
     *                                + W2I2*AL(I+I4+2)
     *                                + W2I3*AL(I+I4+3)
                       SUM5 =    SUM5 + W2I *AL(I+I5  )
     *                                + W2I1*AL(I+I5+1)
     *                                + W2I2*AL(I+I5+2)
     *                                + W2I3*AL(I+I5+3)
     *                                + W2I4*AL(I+I5+4)
                    W5(I  ) = W5(I  ) + W2I1*AL(I+I1  )
     *                                + W2I2*AL(I+I2  )
     *                                + W2I3*AL(I+I3  )
     *                                + W2I4*AL(I+I4  )
     *                                + W2I5*AL(I+I5  )
                    W5(I+1) = W5(I+1) + W2I2*AL(I+I2+1)
     *                                + W2I3*AL(I+I3+1)
     *                                + W2I4*AL(I+I4+1)
     *                                + W2I5*AL(I+I5+1)
                    W5(I+2) = W5(I+2) + W2I3*AL(I+I3+2)
     *                                + W2I4*AL(I+I4+2)
     *                                + W2I5*AL(I+I5+2)
                    W5(I+3) = W5(I+3) + W2I4*AL(I+I4+3)
     *                                + W2I5*AL(I+I5+3)
                    W5(I+4) = W5(I+4) + W2I5*AL(I+I5+4)
                    W6(I  ) =    SUM0 + W2I *AL(I+I0  )
                    W6(I+1) =    SUM1 + W2I1*AL(I+I1+1)
                    W6(I+2) =    SUM2 + W2I2*AL(I+I2+2)
                    W6(I+3) =    SUM3 + W2I3*AL(I+I3+3)
                    W6(I+4) =    SUM4 + W2I4*AL(I+I4+4)
                    W6(I+5) =    SUM5 + W2I5*AL(I+I5+5)
                 END DO
C                End first numerically intensive kernel.
C                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

                 SUM = 0D0
                 DO I=KP1,N
                    W56 = W5(I) + W6(I)
                    W1(I) = HI*W56
                    SUM = SUM + HI*W56*W2(I)
                 END DO

                 U = P5*SUM*HI
                 NREST = MOD(N-KP1+1,6)
                 DO I=KP1,KP1+NREST-1
                    I0    = IX(I)
                    W1(I) = W2(I)*U - W1(I)
                    W1I   = W1(I)
                    W2I   = W2(I)
                    DO J=KP1,I
                       AL(J+I0) = AL(J+I0) + W2I*W1(J) + W1I*W2(J)
                    END DO
                 END DO
                 DO I=KP1+NREST,N,6
                    W1(I  ) = W2(I  )*U - W1(I  )
                    W1(I+1) = W2(I+1)*U - W1(I+1)
                    W1(I+2) = W2(I+2)*U - W1(I+2)
                    W1(I+3) = W2(I+3)*U - W1(I+3)
                    W1(I+4) = W2(I+4)*U - W1(I+4)
                    W1(I+5) = W2(I+5)*U - W1(I+5)
                 END DO

C                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C                Second numerically intensive kernel.
CCCC             DO I=KP1+NREST,N,6
                 DO I=N-5,KP1+NREST,-6
                    I0    = IX(I  )
                    I1    = IX(I+1)
                    I2    = IX(I+2)
                    I3    = IX(I+3)
                    I4    = IX(I+4)
                    I5    = IX(I+5)
                    W2I   = W2(I  )
                    W2I1  = W2(I+1)
                    W2I2  = W2(I+2)
                    W2I3  = W2(I+3)
                    W2I4  = W2(I+4)
                    W2I5  = W2(I+5)
                    W1I   = W1(I  )
                    W1I1  = W1(I+1)
                    W1I2  = W1(I+2)
                    W1I3  = W1(I+3)
                    W1I4  = W1(I+4)
                    W1I5  = W1(I+5)

C                   ////////////////////
C                   Start innermostloop.
!DIR$ IVDEP:LOOP
!DIR$LOOP COUNT (4000)
                    DO J=KP1,I
                       AL(J+I0) = AL(J+I0) + W2I *W1(J) + W1I *W2(J)
                       AL(J+I1) = AL(J+I1) + W2I1*W1(J) + W1I1*W2(J)
                       AL(J+I2) = AL(J+I2) + W2I2*W1(J) + W1I2*W2(J)
                       AL(J+I3) = AL(J+I3) + W2I3*W1(J) + W1I3*W2(J)
                       AL(J+I4) = AL(J+I4) + W2I4*W1(J) + W1I4*W2(J)
                       AL(J+I5) = AL(J+I5) + W2I5*W1(J) + W1I5*W2(J)
!                      AUX = AL(J+I0) + W2I *W1(J)
!                      AL(J+I0) = AUX + W1I *W2(J)
!                      BUX = AL(J+I1) + W2I1*W1(J)
!                      AL(J+I1) = BUX + W1I1*W2(J)
!                      CUX = AL(J+I2) + W2I2*W1(J)
!                      AL(J+I2) = CUX + W1I2*W2(J)
!                      DUX = AL(J+I3) + W2I3*W1(J)
!                      AL(J+I3) = DUX + W1I3*W2(J)
!                      EUX = AL(J+I4) + W2I4*W1(J)
!                      AL(J+I4) = EUX + W1I4*W2(J)
!                      FUX = AL(J+I5) + W2I5*W1(J)
!                      AL(J+I5) = FUX + W1I5*W2(J)
                    END DO 
C                   End innermostloop.
C                   //////////////////

                   AL(I+I1+1) = AL(I+I1+1) + W2I1*W1(I+1) + W1I1*W2(I+1)
                   AL(I+I2+1) = AL(I+I2+1) + W2I2*W1(I+1) + W1I2*W2(I+1)
                   AL(I+I2+2) = AL(I+I2+2) + W2I2*W1(I+2) + W1I2*W2(I+2)
                   AL(I+I3+1) = AL(I+I3+1) + W2I3*W1(I+1) + W1I3*W2(I+1)
                   AL(I+I3+2) = AL(I+I3+2) + W2I3*W1(I+2) + W1I3*W2(I+2)
                   AL(I+I3+3) = AL(I+I3+3) + W2I3*W1(I+3) + W1I3*W2(I+3)
                   AL(I+I4+1) = AL(I+I4+1) + W2I4*W1(I+1) + W1I4*W2(I+1)
                   AL(I+I4+2) = AL(I+I4+2) + W2I4*W1(I+2) + W1I4*W2(I+2)
                   AL(I+I4+3) = AL(I+I4+3) + W2I4*W1(I+3) + W1I4*W2(I+3)
                   AL(I+I4+4) = AL(I+I4+4) + W2I4*W1(I+4) + W1I4*W2(I+4)
                   AL(I+I5+1) = AL(I+I5+1) + W2I5*W1(I+1) + W1I5*W2(I+1)
                   AL(I+I5+2) = AL(I+I5+2) + W2I5*W1(I+2) + W1I5*W2(I+2)
                   AL(I+I5+3) = AL(I+I5+3) + W2I5*W1(I+3) + W1I5*W2(I+3)
                   AL(I+I5+4) = AL(I+I5+4) + W2I5*W1(I+4) + W1I5*W2(I+4)
                   AL(I+I5+5) = AL(I+I5+5) + W2I5*W1(I+5) + W1I5*W2(I+5)
                 END DO
C                End second numerically intensive kernel.
C                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
             ELSE
                 HUNS = 0D0
             END IF
             AL(K+IX(K)) = HUNS
          END DO
      END IF
      NM1NM1  = IX(NM1) + NM1
      NM1N    = IX(N  ) + NM1
      NN      = NM1N    + 1
      W2(NM1) = AL(NM1NM1)
      W2(N  ) = AL(NN)
      W1(NM1) = AL(NM1N)
      W1(N  ) = 0D0
      END
C-------------------------------------------------------
      SUBROUTINE HOUSEHOLDER_U8 (N_NL,AL,W1_NL,W2_NL,W5_NL,W6_NL,IX)
C-------------------------------------------------------

C  Purpose: to tridiagonalize matrix AL, hand unrolling=8.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (P5=0.5D0)
      DIMENSION AL(*),W1_NL(*),W2_NL(*),W5_NL(*),W6_NL(*),IX(*)
      ALLOCATABLE :: W1(:), W2(:), W5(:), W6(:)

      N = N_NL
      ALLOCATE (W1(N), W2(N), W5(N), W6(N))

      NM1 = N - 1
      IF (N .GT. 2) THEN
C         Householder transformation.
          DO K=1,N-2
              KP1  = K + 1
             JREST = MOD(K,8)
             W2(K) = AL(K+IX(K))
             SCALE = 0D0
             DO J=KP1,N
                SCALE = SCALE + ABS (AL(IX(J)+K))
             END DO
             W1(K) = AL(IX(KP1)+K)
             IF (SCALE .GT. 0D0) THEN
                 SCALEI = 1D0/SCALE
                 DO J=KP1,N
                    W2(J) = AL(IX(J)+K)*SCALEI
                 END DO
                 SUM  = 0D0
                 DO J=KP1,N
                    SUM  = SUM + W2(J)**2
                 END DO
                 S = SIGN (SQRT (SUM),W2(KP1))
                 W1(K  ) = -S*SCALE
                 W2(KP1) = W2(KP1) + S
                 AL(IX(KP1)+K) = W2(KP1)*SCALE
                 H    = W2(KP1)*S
                 HUNS = (H*SCALE)*SCALE
                 HI   = 1D0/H
                 DO II=KP1,N
                    W5(II) = 0D0
                 END DO
                 NREST = MOD(N-KP1+1,8)
                 DO I=KP1,KP1+NREST-1
                    IM1  = I - 1
                    SUM  = 0D0
                    I0   = IX(I)
                    W2I  = W2(I)
                    DO J=KP1,IM1
                       SUM   =  SUM + W2(J)*AL(J+I0)
                       W5(J) =  W5(J) + W2I*AL(J+I0)
                    END DO
                    W6(I) = W2I*AL(I+I0) + SUM
                 END DO

C                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C                First numerically intensive kernel.
!DIR$IVDEP
                 DO I=KP1+NREST,N,8
                    IM1  = I - 1
                    SUM0 = 0D0
                    SUM1 = 0D0
                    SUM2 = 0D0
                    SUM3 = 0D0
                    SUM4 = 0D0
                    SUM5 = 0D0
                    SUM6 = 0D0
                    SUM7 = 0D0
                    I0   = IX(I  )
                    I1   = IX(I+1)
                    I2   = IX(I+2)
                    I3   = IX(I+3)
                    I4   = IX(I+4)
                    I5   = IX(I+5)
                    I6   = IX(I+6)
                    I7   = IX(I+7)
                    W2I  = W2(I  )
                    W2I1 = W2(I+1)
                    W2I2 = W2(I+2)
                    W2I3 = W2(I+3) 
                    W2I4 = W2(I+4) 
                    W2I5 = W2(I+5) 
                    W2I6 = W2(I+6) 
                    W2I7 = W2(I+7) 

C                   ////////////////////
C                   Start innermostloop.
!DIR$ IVDEP:LOOP
!DIR$LOOP COUNT (4000)
                    DO J=KP1,IM1
                       W5(J) = W5(J) + W2I  *AL(J+I0)
     *                               + W2I2 *AL(J+I2)
     *                               + W2I1 *AL(J+I1)
     *                               + W2I3 *AL(J+I3)
                        SUM0 =  SUM0 + W2(J)*AL(J+I0)
                        SUM2 =  SUM2 + W2(J)*AL(J+I2)
                        SUM1 =  SUM1 + W2(J)*AL(J+I1)
                        SUM3 =  SUM3 + W2(J)*AL(J+I3)
                    END DO
!DIR$IVDEP
                    DO J=KP1,IM1
                       W5(J) = W5(J) + W2I4 *AL(J+I4)
     *                               + W2I6 *AL(J+I6)
     *                               + W2I5 *AL(J+I5)
     *                               + W2I7 *AL(J+I7)
                        SUM4 =  SUM4 + W2(J)*AL(J+I4)
                        SUM6 =  SUM6 + W2(J)*AL(J+I6)
                        SUM5 =  SUM5 + W2(J)*AL(J+I5)
                        SUM7 =  SUM7 + W2(J)*AL(J+I7)
                    END DO
C                   End innermostloop.
C                   //////////////////

                       SUM1 =    SUM1 + W2I *AL(I+I1  )
                       SUM2 =    SUM2 + W2I *AL(I+I2  )
     *                                + W2I1*AL(I+I2+1)
                       SUM3 =    SUM3 + W2I *AL(I+I3  )
     *                                + W2I1*AL(I+I3+1)
     *                                + W2I2*AL(I+I3+2)
                       SUM4 =    SUM4 + W2I *AL(I+I4  )
     *                                + W2I1*AL(I+I4+1)
     *                                + W2I2*AL(I+I4+2)
     *                                + W2I3*AL(I+I4+3)
                       SUM5 =    SUM5 + W2I *AL(I+I5  )
     *                                + W2I1*AL(I+I5+1)
     *                                + W2I2*AL(I+I5+2)
     *                                + W2I3*AL(I+I5+3)
     *                                + W2I4*AL(I+I5+4)
                       SUM6 =    SUM6 + W2I *AL(I+I6  )
     *                                + W2I1*AL(I+I6+1)
     *                                + W2I2*AL(I+I6+2)
     *                                + W2I3*AL(I+I6+3)
     *                                + W2I4*AL(I+I6+4)
     *                                + W2I5*AL(I+I6+5)
                       SUM7 =    SUM7 + W2I *AL(I+I7  )
     *                                + W2I1*AL(I+I7+1)
     *                                + W2I2*AL(I+I7+2)
     *                                + W2I3*AL(I+I7+3)
     *                                + W2I4*AL(I+I7+4)
     *                                + W2I5*AL(I+I7+5)
     *                                + W2I6*AL(I+I7+6)
                    W5(I  ) = W5(I  ) + W2I1*AL(I+I1  )
     *                                + W2I2*AL(I+I2  )
     *                                + W2I3*AL(I+I3  )
     *                                + W2I4*AL(I+I4  )
     *                                + W2I5*AL(I+I5  )
     *                                + W2I6*AL(I+I6  )
     *                                + W2I7*AL(I+I7  )
                    W5(I+1) = W5(I+1) + W2I2*AL(I+I2+1)
     *                                + W2I3*AL(I+I3+1)
     *                                + W2I4*AL(I+I4+1)
     *                                + W2I5*AL(I+I5+1)
     *                                + W2I6*AL(I+I6+1)
     *                                + W2I7*AL(I+I7+1)
                    W5(I+2) = W5(I+2) + W2I3*AL(I+I3+2)
     *                                + W2I4*AL(I+I4+2)
     *                                + W2I5*AL(I+I5+2)
     *                                + W2I6*AL(I+I6+2)
     *                                + W2I7*AL(I+I7+2)
                    W5(I+3) = W5(I+3) + W2I4*AL(I+I4+3)
     *                                + W2I5*AL(I+I5+3)
     *                                + W2I6*AL(I+I6+3)
     *                                + W2I7*AL(I+I7+3)
                    W5(I+4) = W5(I+4) + W2I5*AL(I+I5+4)
     *                                + W2I6*AL(I+I6+4)
     *                                + W2I7*AL(I+I7+4)
                    W5(I+5) = W5(I+5) + W2I6*AL(I+I6+5)
     *                                + W2I7*AL(I+I7+5)
                    W5(I+6) = W5(I+6) + W2I7*AL(I+I7+6)
                    W6(I  ) =    SUM0 + W2I *AL(I+I0  )
                    W6(I+1) =    SUM1 + W2I1*AL(I+I1+1)
                    W6(I+2) =    SUM2 + W2I2*AL(I+I2+2)
                    W6(I+3) =    SUM3 + W2I3*AL(I+I3+3)
                    W6(I+4) =    SUM4 + W2I4*AL(I+I4+4)
                    W6(I+5) =    SUM5 + W2I5*AL(I+I5+5)
                    W6(I+6) =    SUM6 + W2I6*AL(I+I6+6)
                    W6(I+7) =    SUM7 + W2I7*AL(I+I7+7)
                 END DO
C                End first numerically intensive kernel.
C                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

                 SUM = 0D0
                 DO I=KP1,N
                    W56 = W5(I) + W6(I)
                    W1(I) = HI*W56
                    SUM = SUM + HI*W56*W2(I)
                 END DO

                 U = P5*SUM*HI
                 NREST = MOD(N-KP1+1,8)
                 DO I=KP1,KP1+NREST-1
                    I0    = IX(I)
                    W1(I) = W2(I)*U - W1(I)
                    W1I   = W1(I)
                    W2I   = W2(I)
                    DO J=KP1,I
                       AL(J+I0) = AL(J+I0) + W2I*W1(J) + W1I*W2(J)
                    END DO
                 END DO
                 DO I=KP1+NREST,N,8
                    W1(I  ) = W2(I  )*U - W1(I  )
                    W1(I+1) = W2(I+1)*U - W1(I+1)
                    W1(I+2) = W2(I+2)*U - W1(I+2)
                    W1(I+3) = W2(I+3)*U - W1(I+3)
                    W1(I+4) = W2(I+4)*U - W1(I+4)
                    W1(I+5) = W2(I+5)*U - W1(I+5)
                    W1(I+6) = W2(I+6)*U - W1(I+6)
                    W1(I+7) = W2(I+7)*U - W1(I+7)
                 END DO

C                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C                Second numerically intensive kernel.
!DIR$IVDEP
                 DO I=N-7,KP1+NREST,-8
                    I0    = IX(I  )
                    I1    = IX(I+1)
                    I2    = IX(I+2)
                    I3    = IX(I+3)
                    I4    = IX(I+4)
                    I5    = IX(I+5)
                    I6    = IX(I+6)
                    I7    = IX(I+7)
                    W2I   = W2(I  )
                    W2I1  = W2(I+1)
                    W2I2  = W2(I+2)
                    W2I3  = W2(I+3)
                    W2I4  = W2(I+4)
                    W2I5  = W2(I+5)
                    W2I6  = W2(I+6)
                    W2I7  = W2(I+7)
                    W1I   = W1(I  )
                    W1I1  = W1(I+1)
                    W1I2  = W1(I+2)
                    W1I3  = W1(I+3)
                    W1I4  = W1(I+4)
                    W1I5  = W1(I+5)
                    W1I6  = W1(I+6)
                    W1I7  = W1(I+7)

C                   ////////////////////
C                   Start innermostloop.
                    JFIN = MIN(KP1+JREST-1,I)
                    DO J=KP1,JFIN
                       AL(J+I0) = AL(J+I0) + W2I *W1(J) + W1I *W2(J)
                       AL(J+I1) = AL(J+I1) + W2I1*W1(J) + W1I1*W2(J)
                       AL(J+I2) = AL(J+I2) + W2I2*W1(J) + W1I2*W2(J)
                       AL(J+I3) = AL(J+I3) + W2I3*W1(J) + W1I3*W2(J)
                       AL(J+I4) = AL(J+I4) + W2I4*W1(J) + W1I4*W2(J)
                       AL(J+I5) = AL(J+I5) + W2I5*W1(J) + W1I5*W2(J)
                       AL(J+I6) = AL(J+I6) + W2I6*W1(J) + W1I6*W2(J)
                       AL(J+I7) = AL(J+I7) + W2I7*W1(J) + W1I7*W2(J)
                    END DO
!DIR$ IVDEP:LOOP
!DIR$LOOP COUNT (4000)
                    DO J=JFIN+1,I
!                      AL(J+I0) = AL(J+I0) + W2I *W1(J) + W1I *W2(J)
!                      AL(J+I1) = AL(J+I1) + W2I1*W1(J) + W1I1*W2(J)
!                      AL(J+I2) = AL(J+I2) + W2I2*W1(J) + W1I2*W2(J)
!                      AL(J+I3) = AL(J+I3) + W2I3*W1(J) + W1I3*W2(J)
!                      AL(J+I4) = AL(J+I4) + W2I4*W1(J) + W1I4*W2(J)
!                      AL(J+I5) = AL(J+I5) + W2I5*W1(J) + W1I5*W2(J)
!                      AL(J+I6) = AL(J+I6) + W2I6*W1(J) + W1I6*W2(J)

                       TEMP     = AL(J+I0) + W2I *W1(J) + W1I *W2(J)
                       TEMP2    = AL(J+I2) + W2I2*W1(J) + W1I2*W2(J)
                       TEMP1    = AL(J+I1) + W2I1*W1(J) + W1I1*W2(J)
                       TEMP4    = AL(J+I4) + W2I4*W1(J) + W1I4*W2(J)
                       TEMP3    = AL(J+I3) + W2I3*W1(J) + W1I3*W2(J)
                       TEMP6    = AL(J+I6) + W2I6*W1(J) + W1I6*W2(J)
                       TEMP5    = AL(J+I5) + W2I5*W1(J) + W1I5*W2(J)
                       AL(J+I7) = AL(J+I7) + W2I7*W1(J) + W1I7*W2(J)
                       AL(J+I0) = TEMP
                       AL(J+I2) = TEMP2
                       AL(J+I1) = TEMP1
                       AL(J+I4) = TEMP4
                       AL(J+I3) = TEMP3
                       AL(J+I6) = TEMP6
                       AL(J+I5) = TEMP5
!                      AUX = AL(J+I0) + W2I *W1(J)
!                      AL(J+I0) = AUX + W1I *W2(J)
!                      BUX = AL(J+I1) + W2I1*W1(J)
!                      AL(J+I1) = BUX + W1I1*W2(J)
!                      CUX = AL(J+I2) + W2I2*W1(J)
!                      AL(J+I2) = CUX + W1I2*W2(J)
!                      DUX = AL(J+I3) + W2I3*W1(J)
!                      AL(J+I3) = DUX + W1I3*W2(J)
!                      EUX = AL(J+I4) + W2I4*W1(J)
!                      AL(J+I4) = EUX + W1I4*W2(J)
!                      FUX = AL(J+I5) + W2I5*W1(J)
!                      AL(J+I5) = FUX + W1I5*W2(J)
!                      GUX = AL(J+I6) + W2I6*W1(J)
!                      AL(J+I6) = GUX + W1I6*W2(J)
!                      HUX = AL(J+I7) + W2I7*W1(J)
!                      AL(J+I7) = HUX + W1I7*W2(J)
                    END DO 
C                   End innermostloop.
C                   //////////////////

                   AL(I+I1+1) = AL(I+I1+1) + W2I1*W1(I+1) + W1I1*W2(I+1)
                   AL(I+I2+1) = AL(I+I2+1) + W2I2*W1(I+1) + W1I2*W2(I+1)
                   AL(I+I2+2) = AL(I+I2+2) + W2I2*W1(I+2) + W1I2*W2(I+2)
                   AL(I+I3+1) = AL(I+I3+1) + W2I3*W1(I+1) + W1I3*W2(I+1)
                   AL(I+I3+2) = AL(I+I3+2) + W2I3*W1(I+2) + W1I3*W2(I+2)
                   AL(I+I3+3) = AL(I+I3+3) + W2I3*W1(I+3) + W1I3*W2(I+3)
                   AL(I+I4+1) = AL(I+I4+1) + W2I4*W1(I+1) + W1I4*W2(I+1)
                   AL(I+I4+2) = AL(I+I4+2) + W2I4*W1(I+2) + W1I4*W2(I+2)
                   AL(I+I4+3) = AL(I+I4+3) + W2I4*W1(I+3) + W1I4*W2(I+3)
                   AL(I+I4+4) = AL(I+I4+4) + W2I4*W1(I+4) + W1I4*W2(I+4)
                   AL(I+I5+1) = AL(I+I5+1) + W2I5*W1(I+1) + W1I5*W2(I+1)
                   AL(I+I5+2) = AL(I+I5+2) + W2I5*W1(I+2) + W1I5*W2(I+2)
                   AL(I+I5+3) = AL(I+I5+3) + W2I5*W1(I+3) + W1I5*W2(I+3)
                   AL(I+I5+4) = AL(I+I5+4) + W2I5*W1(I+4) + W1I5*W2(I+4)
                   AL(I+I5+5) = AL(I+I5+5) + W2I5*W1(I+5) + W1I5*W2(I+5)
                   AL(I+I6+1) = AL(I+I6+1) + W2I6*W1(I+1) + W1I6*W2(I+1)
                   AL(I+I6+2) = AL(I+I6+2) + W2I6*W1(I+2) + W1I6*W2(I+2)
                   AL(I+I6+3) = AL(I+I6+3) + W2I6*W1(I+3) + W1I6*W2(I+3)
                   AL(I+I6+4) = AL(I+I6+4) + W2I6*W1(I+4) + W1I6*W2(I+4)
                   AL(I+I6+5) = AL(I+I6+5) + W2I6*W1(I+5) + W1I6*W2(I+5)
                   AL(I+I6+6) = AL(I+I6+6) + W2I6*W1(I+6) + W1I6*W2(I+6)
                   AL(I+I7+1) = AL(I+I7+1) + W2I7*W1(I+1) + W1I7*W2(I+1)
                   AL(I+I7+2) = AL(I+I7+2) + W2I7*W1(I+2) + W1I7*W2(I+2)
                   AL(I+I7+3) = AL(I+I7+3) + W2I7*W1(I+3) + W1I7*W2(I+3)
                   AL(I+I7+4) = AL(I+I7+4) + W2I7*W1(I+4) + W1I7*W2(I+4)
                   AL(I+I7+5) = AL(I+I7+5) + W2I7*W1(I+5) + W1I7*W2(I+5)
                   AL(I+I7+6) = AL(I+I7+6) + W2I7*W1(I+6) + W1I7*W2(I+6)
                   AL(I+I7+7) = AL(I+I7+7) + W2I7*W1(I+7) + W1I7*W2(I+7)
                 END DO
C                End second numerically intensive kernel.
C                $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
             ELSE
                 HUNS = 0D0
             END IF
             AL(K+IX(K)) = HUNS
          END DO
      END IF
      NM1NM1  = IX(NM1) + NM1
      NM1N    = IX(N  ) + NM1
      NN      = NM1N    + 1
      W2(NM1) = AL(NM1NM1)
      W2(N  ) = AL(NN)
      W1(NM1) = AL(NM1N)
      W1(N  ) = 0D0
      IJ = 0
      DO J=1,N
         JN = (J-1)*N
      DO I=1,J
         IJ = IJ + 1
         IJN = I + JN
         AL(IJ) = AL(IJN)
      END DO
      END DO
      DO JI=(N*N+N)/2+1,N*N
         AL(JI) = 0D0
      END DO
      DO I=1,N
         W1_NL(I) = W1(I)
         W2_NL(I) = W2(I)
      END DO
      DEALLOCATE (W1,W2,W5,W6)
      END
C--------------------------------------------
      SUBROUTINE TRANSPOSE (N_NL,AL,IX,JX,AU)
C--------------------------------------------

C     Purpose: to transpose matrix AL into matrix AU.

C     N_NL     is the nonlocal value of N, the matrix dimension.
C     AU       real symmetric matrix of dimension N1=(N*N+N)/2 given
C              in rowwise upper triangular form, viz., AU(IJ)= AA(I,J)
C              where IJ=IX(I)+J as explained below and AA is the usual
C              square symmetric form.
C     IX       Notation to access upper triangular part by rows.
C     JX       Notation to access upper triangular part by columns.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AL(*),AU(*),IX(*),JX(*)

      N = N_NL
C     Use sequential row notation for upper triangular matrix.
      IX(1) = 0
      JX(1) = 0
      DO J=2,N
         IX(J) = IX(J-1) - J + 1 + N
         JX(J) = JX(J-1) + J - 1
      END DO

C     Obtain transposed matrix AU.
      NREST_J = MOD(N,4)
      DO J=1,NREST_J
         JF = JX(J)
         DO I=1,J
            AU(IX(I)+J) = AL(JF+I)
         END DO
      END DO
      DO J=1+NREST_J,N,4
         JF = JX(J  )
         JF1= JX(J+1)
         JF2= JX(J+2)
         JF3= JX(J+3)
         NREST_I = MOD(J,4)
         DO I=1,J-NREST_I,4
            AU(IX(I  )+J  ) = AL(JF +I  )
            AU(IX(I  )+J+1) = AL(JF1+I  )
            AU(IX(I  )+J+2) = AL(JF2+I  )
            AU(IX(I  )+J+3) = AL(JF3+I  )
            AU(IX(I+1)+J  ) = AL(JF +I+1)
            AU(IX(I+1)+J+1) = AL(JF1+I+1)
            AU(IX(I+1)+J+2) = AL(JF2+I+1)
            AU(IX(I+1)+J+3) = AL(JF3+I+1)
            AU(IX(I+2)+J  ) = AL(JF +I+2)
            AU(IX(I+2)+J+1) = AL(JF1+I+2)
            AU(IX(I+2)+J+2) = AL(JF2+I+2)
            AU(IX(I+2)+J+3) = AL(JF3+I+2)
            AU(IX(I+3)+J  ) = AL(JF +I+3)
            AU(IX(I+3)+J+1) = AL(JF1+I+3)
            AU(IX(I+3)+J+2) = AL(JF2+I+3)
            AU(IX(I+3)+J+3) = AL(JF3+I+3)
         END DO
         DO I=J-NREST_I+1,J
            AU(IX(I)+J  ) = AL(JF +I)
            AU(IX(I)+J+1) = AL(JF1+I)
            AU(IX(I)+J+2) = AL(JF2+I)
            AU(IX(I)+J+3) = AL(JF3+I)
         END DO
            AU(IX(J+1)+J+1) = AL(JF1+J+1)
            AU(IX(J+1)+J+2) = AL(JF2+J+1)
            AU(IX(J+2)+J+2) = AL(JF2+J+2)
            AU(IX(J+1)+J+3) = AL(JF3+J+1)
            AU(IX(J+2)+J+3) = AL(JF3+J+2)
            AU(IX(J+3)+J+3) = AL(JF3+J+3)
      END DO
      END
