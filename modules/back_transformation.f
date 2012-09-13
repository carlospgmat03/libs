C-----------------------------------------------------------------------
      SUBROUTINE BACK_TRANSFORMATION (IUB,N,IV,LV,SCHMIDT,EPS2,AU,E,NVX,
     *                                                     V,IWORK,WORK)
C-----------------------------------------------------------------------

C     Purpose: to choose appropriate hand unrolling for the
C              back transformation step.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL TIME,CPUT
      LOGICAL SCHMIDT
      DIMENSION AU(*),E(*),V(*),WORK(*),IWORK(*)
      EXTERNAL BACK_TRANSFORMATION_U2, BACK_TRANSFORMATION_U3,
     *         BACK_TRANSFORMATION_U4, BACK_TRANSFORMATION_U5,
     *         BACK_TRANSFORMATION_U6, BACK_TRANSFORMATION_U8,
     *         BACK_TRANSFORMATION_U10,BACK_TRANSFORMATION_U16
 
C      CALL TIEMPO (.TRUE.,.TRUE.,'         Begin back_tr.',6,TIME,CPUT)
      IF (IUB .LE. 0 ) IUB = 1
      IF (IUB .GT. 16) IUB =16

      CALL DRIVE_BACK_TRANSFORMATION (N,IV,LV,SCHMIDT,EPS2,AU,E,NVX,
     *                                IUB,WORK,
     *                                WORK((IUB+1)*N+1),IWORK,
     *                                IWORK(N+1)       ,IWORK(2*N+1),
     *                                V,               IF,KRUN,NRVEC)

      SELECT CASE (IUB)
             CASE(2)
             CALL BACK_TRANSFORMATION_U2 (N,IV,LV,SCHMIDT,EPS2,AU,E,NVX,
     *                                    IF,KRUN,NRVEC,WORK,WORK(N+1),
     *                                   WORK((IUB+1)*N+1),IWORK,
     *                                   IWORK(N+1)       ,IWORK(2*N+1),
     *                                                                V)
             CASE(3)
             CALL BACK_TRANSFORMATION_U3 (N,IV,LV,SCHMIDT,EPS2,AU,E,NVX,
     *                                    IF,KRUN,NRVEC,WORK,WORK(N+1),
     *                                   WORK((IUB+1)*N+1),IWORK,
     *                                   IWORK(N+1)       ,IWORK(2*N+1),
     *                                                                V)
             CASE(4)
             CALL BACK_TRANSFORMATION_U4 (N,IV,LV,SCHMIDT,EPS2,AU,E,NVX,
     *                                    IF,KRUN,NRVEC,WORK,WORK(N+1),
     *                                   WORK((IUB+1)*N+1),IWORK,
     *                                   IWORK(N+1)       ,IWORK(2*N+1),
     *                                                                V)
             CASE(5)
             CALL BACK_TRANSFORMATION_U5 (N,IV,LV,SCHMIDT,EPS2,AU,E,NVX,
     *                                    IF,KRUN,NRVEC,WORK,WORK(N+1),
     *                                   WORK((IUB+1)*N+1),IWORK,
     *                                   IWORK(N+1)       ,IWORK(2*N+1),
     *                                                                V)
             CASE(6)
             CALL BACK_TRANSFORMATION_U6 (N,IV,LV,SCHMIDT,EPS2,AU,E,NVX,
     *                                    IF,KRUN,NRVEC,WORK,WORK(N+1),
     *                                   WORK((IUB+1)*N+1),IWORK,
     *                                   IWORK(N+1)       ,IWORK(2*N+1),
     *                                                                V)
             CASE(7:8)
             CALL BACK_TRANSFORMATION_U8 (N,IV,LV,SCHMIDT,EPS2,AU,E,NVX,
     *                                    IF,KRUN,NRVEC,WORK,WORK(N+1),
     *                                   WORK((IUB+1)*N+1),IWORK,
     *                                   IWORK(N+1)       ,IWORK(2*N+1),
     *                                                                V)
             CASE(9:10)
             CALL BACK_TRANSFORMATION_U10(N,IV,LV,SCHMIDT,EPS2,AU,E,NVX,
     *                                    IF,KRUN,NRVEC,WORK,WORK(N+1),
     *                                   WORK((IUB+1)*N+1),IWORK,
     *                                   IWORK(N+1)       ,IWORK(2*N+1),
     *                                                                V)
             CASE(11:12)
             CALL BACK_TRANSFORMATION_U12(N,IV,LV,SCHMIDT,EPS2,AU,E,NVX,
     *                                    IF,KRUN,NRVEC,WORK,WORK(N+1),
     *                                   WORK((IUB+1)*N+1),IWORK,
     *                                   IWORK(N+1)       ,IWORK(2*N+1),
     *                                                                V)
             CASE(13:16)
             CALL BACK_TRANSFORMATION_U16(N,IV,LV,SCHMIDT,EPS2,AU,E,NVX,
     *                                    IF,KRUN,NRVEC,WORK,WORK(N+1),
     *                                   WORK((IUB+1)*N+1),IWORK,
     *                                   IWORK(N+1)       ,IWORK(2*N+1),
     *                                                                V)
      END SELECT
C      CALL TIEMPO (.TRUE.,.FALSE.,'           End back_tr.',6,TIME,CPUT)
      END
C---------------------------------------------------------------------
      SUBROUTINE DRIVE_BACK_TRANSFORMATION (N,IV,LV,SCHMIDT,EPS2,AU,E,
     *                                      NVX,IUB,W6,AU_NEQ_0_I,IX,
     *                                      IG,K_TRUE,V,IF,KRUN,NRVEC)
C---------------------------------------------------------------------

C     AU       real symmetric matrix of dimension N1=(N*N+N)/2 given
C              in rowwise upper triangular form, viz., AU(IJ)= AA(I,J)
C              where IJ=IX(I)+J as explained below and AA is the usual
C              square symmetric form.
C     E        array holding eigenvalues.
C     W6       Working array  of dimension N.
C     AU_NEQ_0_I  Array holding nonzero diagonal elements of AU.
C     IX       Notation to access upper triangular part by rows.
C     IG       Working array.
C     K_TRUE   Array holding indices of nonzero diagonal elements.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SCHMIDT
      DIMENSION AU(*),E(*),V(*),
     *          W6(*),AU_NEQ_0_I(*),IX(*),IG(*),K_TRUE(*)

C     Use sequential row notation for upper triangular matrix.
      IX(1) = 0
      DO J=2,N
         IX(J) = IX(J-1) - J + 1 + N
      END DO

      DO I=1,IUB
         IG(I) = 1
      END DO

C     Collect inverses of nonzero diagonal matrix elements.
      KRUN = 0
       NM2 = N - 2
      DO K=NM2,1,-1
         IF (AU(K+IX(K)) .NE. 0D0) THEN
             KRUN = KRUN + 1
             K_TRUE(KRUN) = K
             AU_NEQ_0_I(KRUN) = 1D0/AU(K+IX(K))
         END IF
      END DO

C     ///////////////////////////////////////////////
C     Start short loops for end-effects of unrolling.
      IF = (IV-2)*NVX
      NRVEC = MOD(LV-IV+1,IUB)
      IF (IUB .EQ. 1) THEN
          IVFN = LV
      ELSE
          IVFN = IV + NRVEC - 1
      END IF
      DO I=IV,IVFN
         IF  = IF + NVX
         DO J=1,N
            W6(J) = V(J+IF)
         END DO
         IM1 = I  - 1
         IF (N .GT. 2) THEN
             DO KR=1,KRUN
                K = K_TRUE(KR)
                K0= IX(K)
                KP1   = K + 1
                SUM   = 0D0
                DO KK=N,KP1,-1
                   SUM = SUM + AU(KK+K0)*W6(KK)
                END DO
                S = -SUM*AU_NEQ_0_I(KR)
                DO KK=N,KP1,-1
                   W6(KK) = W6(KK) + S*AU(KK+K0)
                END DO
             END DO
         END IF

         DO J=IG(1),I-1
            IF (ABS(E(J)-E(I)) .LT. EPS2) THEN
                IG(1) = J
                GO TO 60
            END IF
         END DO
         IG(1) = I
   60    CONTINUE

         IF (IG(1) .NE. I  .AND.  SCHMIDT) THEN

C            Degenerate eigenvalues.  First, orthogonalize.
             KF  = (IG(1)-2)*NVX
             DO K=IG(1),IM1
                KF  = KF + NVX
                SUM = 0D0
                DO J=1,N
                   SUM = SUM + V(J+KF)*W6(J)
                END DO
                S = -SUM
                DO J=1,N
                   W6(J) = W6(J) + S*V(J+KF)
                END DO
             END DO

         END IF

C        Normalization.
         SUM = 0D0
         DO J=1,N
            SUM = SUM + W6(J)**2
         END DO
         S = 1D0/SQRT (SUM)
         DO J=1,N
            V(IF+J) = S*W6(J)
         END DO

      END DO
C     End short loops for end-effects of unrolling.
C     /////////////////////////////////////////////
      END
C----------------------------------------------------------------------
      SUBROUTINE BACK_TRANSFORMATION_U2 (N_NL,IV_NL,LV_NL,SCHMIDT,EPS2,
     *                                   AU,E,NVX,IF_NL,KRUN_NL,
     *                                   NRVEC_NL,
     *                                   W6,WW,AU_NEQ_0_I,IX,IG,K_TRUE,
     *                                                               V)
C----------------------------------------------------------------------

C     Array meanings explained in SUBROUTINE DRIVE_BACK_TRANSFORMATION.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SCHMIDT
      DIMENSION AU(*),E(*),V(*),
     *          W6(*),WW(*),AU_NEQ_0_I(*),IX(*),IG(*),K_TRUE(*)

      N = N_NL
      IV = IV_NL
      LV = LV_NL
      IF = IF_NL
      KRUN = KRUN_NL
      NRVEC = NRVEC_NL

C     //////////////////////////////////
C     Start full-fledged unrolled loops.
      IF1 = IF
      DO I=IV+NRVEC,LV,2

         IF  = IF1 + NVX
         IF1 = IF  + NVX
         JJ  = 0 
         DO J=1,2*N-1,2
            JJ = JJ + 1
            WW(J  ) = V(JJ+IF )
            WW(J+1) = V(JJ+IF1)
         END DO

         IM1 = I - 1
         IF (N .GT. 2) THEN
             DO KR=1,KRUN
                K = K_TRUE(KR)
                K0 = IX(K)
                KP1 = K + 1
                SUMA = 0D0
                SUMB = 0D0
                 KK0 = 2*N + 1

C               ////////////////////
C               Start innermostloop.
                 KKK = KK0
                DO KK=N,KP1,-1
                   KKK  = KKK - 2
                   SUMA = SUMA + AU(KK+K0)*WW(KKK  )
                   SUMB = SUMB + AU(KK+K0)*WW(KKK+1)
                END DO
C               End innermostloop.
C               //////////////////

                AUI = 1D0*AU_NEQ_0_I(KR)
                 SA = - SUMA*AUI
                 SB = - SUMB*AUI

C               ////////////////////
C               Start innermostloop.
                 KKK = KK0
                DO KK=N,KP1,-1
                         KKK = KKK - 2
                   WW(KKK  ) = WW(KKK  ) + SA*AU(KK+K0)
                   WW(KKK+1) = WW(KKK+1) + SB*AU(KK+K0)
                END DO
C               End innermostloop.
C               //////////////////

             END DO
         END IF

         DO II=1,2
            III = I + II - 1
            DO J=IG(II),III-1
               IF (ABS(E(J) - E(III)) .LT. EPS2) THEN
                   IG(II) = J
                   GO TO 70
               END IF
            END DO
            IG(II) = III
   70       CONTINUE

            IF (IG(II) .NE. III  .AND.  SCHMIDT) THEN

C               Degenerate eigenvalues.  First, orthogonalize.
                KF = (IG(II)-2)*NVX
                DO K=IG(II),III-1
                   KF  = KF + NVX
                   SUM = 0D0
                    JF = -2
                   DO J=1,N
                      JF  = JF + 2
                      SUM = SUM + V(J+KF)*WW(JF+II)
                   END DO
                    S = -SUM
                   JF = -2
                   DO J=1,N
                      JF = JF + 2
                      WW(JF+II) = WW(JF+II) + S*V(J+KF)
                   END DO
                END DO
            END IF
         END DO

C        Normalization.
         SUM0 = 0D0
         SUM1 = 0D0
         JF  = -2
         DO J=1,N
            JF = JF + 2
            SUM0 = SUM0 + WW(JF+1)**2
            SUM1 = SUM1 + WW(JF+2)**2
         END DO
         S0 = 1D0/SQRT(SUM0)
         S1 = 1D0/SQRT(SUM1)
         JF = -2
         DO J=1,N
            JF = JF + 2
            V(J+IF ) = S0*WW(JF+1)
            V(J+IF1) = S1*WW(JF+2)
         END DO
      END DO
      END
C----------------------------------------------------------------------
      SUBROUTINE BACK_TRANSFORMATION_U3 (N_NL,IV_NL,LV_NL,SCHMIDT,EPS2,
     *                                   AU,E,NVX,IF_NL,KRUN_NL,
     *                                   NRVEC_NL,
     *                                   W6,WW,AU_NEQ_0_I,IX,IG,K_TRUE,
     *                                                               V)
C----------------------------------------------------------------------

C     Array meanings explained in SUBROUTINE DRIVE_BACK_TRANSFORMATION.
C     AU       real symmetric matrix of dimension N1=(N*N+N)/2 given
C              in rowwise upper triangular form, viz., AU(IJ)= AA(I,J)
C              where IJ=IX(I)+J as explained below and AA is the usual
C              square symmetric form.
C     W6,WW    Working arrays of dimension N and 4*N, respectively.
C     AU_NEQ_0_I  Array holding nonzero diagonal elements of AU.
C     IX       Notation to access upper triangular part by rows.
C     IG       Working array.
C     K_TRUE   Array holding indices of nonzero diagonal elements.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SCHMIDT
      DIMENSION AU(*),E(*),V(*),
     *          W6(*),WW(*),AU_NEQ_0_I(*),IX(*),IG(*),K_TRUE(*)


      N = N_NL
      IV = IV_NL
      LV = LV_NL
      IF = IF_NL
      KRUN = KRUN_NL
      NRVEC = NRVEC_NL

C     //////////////////////////////////
C     Start full-fledged unrolled loops.
      IF2 = IF
      DO I=IV+NRVEC,LV,3

         IF  = IF2 + NVX
         IF1 = IF  + NVX
         IF2 = IF1 + NVX
         JJ  = 0 
         DO J=1,3*N-2,3
                 JJ = JJ + 1
            WW(J  ) = V(JJ+IF )
            WW(J+1) = V(JJ+IF1)
            WW(J+2) = V(JJ+IF2)
         END DO

         IM1 = I - 1
         IF (N .GT. 2) THEN
             DO KR=1,KRUN
                K = K_TRUE(KR)
                K0 = IX(K)
               KP1 = K + 1
              SUMA = 0D0
              SUMB = 0D0
              SUMC = 0D0
               KK0 = 3*N + 1

C               ////////////////////
C               Start innermostloop.
                 KKK = KK0
                DO KK=N,KP1,-1
                   KKK  = KKK - 3
                   SUMA = SUMA + AU(KK+K0)*WW(KKK  )
                   SUMB = SUMB + AU(KK+K0)*WW(KKK+1)
                   SUMC = SUMC + AU(KK+K0)*WW(KKK+2)
                END DO
C               End innermostloop.
C               //////////////////

                AUI = 1D0*AU_NEQ_0_I(KR)
                 SA = - SUMA*AUI
                 SB = - SUMB*AUI
                 SC = - SUMC*AUI

C               ////////////////////
C               Start innermostloop.
                 KKK = KK0
                DO KK=N,KP1,-1
                         KKK = KKK - 3
                   WW(KKK  ) = WW(KKK  ) + SA*AU(KK+K0)
                   WW(KKK+1) = WW(KKK+1) + SB*AU(KK+K0)
                   WW(KKK+2) = WW(KKK+2) + SC*AU(KK+K0)
                END DO
C               End innermostloop.
C               //////////////////

             END DO
         END IF

         DO II=1,3
            III = I + II - 1
            DO J=IG(II),III-1
               IF (ABS(E(J) - E(III)) .LT. EPS2) THEN
                   IG(II) = J
                   GO TO 70
               END IF
            END DO
            IG(II) = III
   70       CONTINUE

            IF (IG(II) .NE. III  .AND.  SCHMIDT) THEN

C               Degenerate eigenvalues.  First, orthogonalize.
                KF = (IG(II)-2)*NVX
                DO K=IG(II),III-1
                   KF  = KF + NVX
                   SUM = 0D0
                    JF = -3
                   DO J=1,N
                      JF  = JF + 3
                      SUM = SUM + V(J+KF)*WW(JF+II)
                   END DO
                    S = -SUM
                   JF = -3
                   DO J=1,N
                      JF = JF + 3
                      WW(JF+II) = WW(JF+II) + S*V(J+KF)
                   END DO
                END DO
            END IF
         END DO

C        Normalization.
         SUM0 = 0D0
         SUM1 = 0D0
         SUM2 = 0D0
         JF  = -3
         DO J=1,N
            JF = JF + 3
            SUM0 = SUM0 + WW(JF+1)**2
            SUM1 = SUM1 + WW(JF+2)**2
            SUM2 = SUM2 + WW(JF+3)**2
         END DO
         S0 = 1D0/SQRT(SUM0)
         S1 = 1D0/SQRT(SUM1)
         S2 = 1D0/SQRT(SUM2)
         JF = -3
         DO J=1,N
            JF = JF + 3
            V(J+IF ) = S0*WW(JF+1)
            V(J+IF1) = S1*WW(JF+2)
            V(J+IF2) = S2*WW(JF+3)
         END DO
      END DO
      END
C----------------------------------------------------------------------
      SUBROUTINE BACK_TRANSFORMATION_U4 (N_NL,IV_NL,LV_NL,SCHMIDT,EPS2,
     *                                   AU,E,NVX,IF_NL,KRUN_NL,
     *                                   NRVEC_NL,
     *                                   W6,WW,AU_NEQ_0_I,IX,IG,K_TRUE,
     *                                                               V)
C----------------------------------------------------------------------

C     Array meanings explained in SUBROUTINE DRIVE_BACK_TRANSFORMATION.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SCHMIDT
      DIMENSION AU(*),E(*),V(*),
     *          W6(*),WW(*),AU_NEQ_0_I(*),IX(*),IG(*),K_TRUE(*)

      N = N_NL
      IV = IV_NL
      LV = LV_NL
      IF = IF_NL
      KRUN = KRUN_NL
      NRVEC = NRVEC_NL

C     //////////////////////////////////
C     Start full-fledged unrolled loops.
      IF3 = IF
      DO I=IV+NRVEC,LV,4

         IF  = IF3 + NVX
         IF1 = IF  + NVX
         IF2 = IF1 + NVX
         IF3 = IF2 + NVX
         JJ  = 0 
         DO J=1,4*N-3,4
                 JJ = JJ + 1
            WW(J  ) = V(JJ+IF )
            WW(J+1) = V(JJ+IF1)
            WW(J+2) = V(JJ+IF2)
            WW(J+3) = V(JJ+IF3)
         END DO

         IM1 = I - 1
         IF (N .GT. 2) THEN
             DO KR=1,KRUN
                  K = K_TRUE(KR)
                 K0 = IX(K)
                KP1 = K + 1
               SUMA = 0D0
               SUMB = 0D0
               SUMC = 0D0
               SUMD = 0D0
                KK0 = 4*N + 1

C               ////////////////////
C               Start innermostloop.
                 KKK = KK0
                DO KK=N,KP1,-1
                   KKK  = KKK - 4
                   SUMA = SUMA + AU(KK+K0)*WW(KKK  )
                   SUMB = SUMB + AU(KK+K0)*WW(KKK+1)
                   SUMC = SUMC + AU(KK+K0)*WW(KKK+2)
                   SUMD = SUMD + AU(KK+K0)*WW(KKK+3)
                END DO
C               End innermostloop.
C               //////////////////

                AUI = 1D0*AU_NEQ_0_I(KR)
                 SA = - SUMA*AUI
                 SB = - SUMB*AUI
                 SC = - SUMC*AUI
                 SD = - SUMD*AUI 

C               ////////////////////
C               Start innermostloop.
                 KKK = KK0
                DO KK=N,KP1,-1
                         KKK = KKK - 4
                   WW(KKK  ) = WW(KKK  ) + SA*AU(KK+K0)
                   WW(KKK+1) = WW(KKK+1) + SB*AU(KK+K0)
                   WW(KKK+2) = WW(KKK+2) + SC*AU(KK+K0)
                   WW(KKK+3) = WW(KKK+3) + SD*AU(KK+K0)
                END DO
C               End innermostloop.
C               //////////////////

             END DO
         END IF

         DO II=1,4
            III = I + II - 1
            DO J=IG(II),III-1
               IF (ABS(E(J) - E(III)) .LT. EPS2) THEN
                   IG(II) = J
                   GO TO 70
               END IF
            END DO
            IG(II) = III
   70       CONTINUE

            IF (IG(II) .NE. III  .AND.  SCHMIDT) THEN

C               Degenerate eigenvalues.  First, orthogonalize.
                KF = (IG(II)-2)*NVX
                DO K=IG(II),III-1
                   KF  = KF + NVX
                   SUM = 0D0
                   JF = -4
                   DO J=1,N
                      JF  = JF + 4
                      SUM = SUM + V(J+KF)*WW(JF+II)
                   END DO
                    S = -SUM
                   JF = -4
                   DO J=1,N
                      JF = JF + 4
                      WW(JF+II) = WW(JF+II) + S*V(J+KF)
                   END DO
                END DO
            END IF
         END DO

C        Normalization.
         SUM0 = 0D0
         SUM1 = 0D0
         SUM2 = 0D0
         SUM3 = 0D0
         JF  = -4
         DO J=1,N
              JF = JF + 4
            SUM0 = SUM0 + WW(JF+1)**2
            SUM1 = SUM1 + WW(JF+2)**2
            SUM2 = SUM2 + WW(JF+3)**2
            SUM3 = SUM3 + WW(JF+4)**2
         END DO
         S0 = 1D0/SQRT(SUM0)
         S1 = 1D0/SQRT(SUM1)
         S2 = 1D0/SQRT(SUM2)
         S3 = 1D0/SQRT(SUM3)
         JF = -4
         DO J=1,N
                  JF = JF + 4
            V(J+IF ) = S0*WW(JF+1)
            V(J+IF1) = S1*WW(JF+2)
            V(J+IF2) = S2*WW(JF+3)
            V(J+IF3) = S3*WW(JF+4)
         END DO
      END DO
      END
C----------------------------------------------------------------------
      SUBROUTINE BACK_TRANSFORMATION_U5 (N_NL,IV_NL,LV_NL,SCHMIDT,EPS2,
     *                                   AU,E,NVX,IF_NL,KRUN_NL,
     *                                   NRVEC_NL,
     *                                   W6,WW,AU_NEQ_0_I,IX,IG,K_TRUE,
     *                                                               V)
C----------------------------------------------------------------------

C     Array meanings explained in SUBROUTINE DRIVE_BACK_TRANSFORMATION.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SCHMIDT
      DIMENSION AU(*),E(*),V(*),
     *          W6(*),WW(*),AU_NEQ_0_I(*),IX(*),IG(*),K_TRUE(*)
  
      N = N_NL
      IV = IV_NL
      LV = LV_NL
      IF = IF_NL
      KRUN = KRUN_NL
      NRVEC = NRVEC_NL

C     //////////////////////////////////
C     Start full-fledged unrolled loops.
      IF4 = IF
      DO I=IV+NRVEC,LV,5

         IF  = IF4 + NVX
         IF1 = IF  + NVX
         IF2 = IF1 + NVX
         IF3 = IF2 + NVX
         IF4 = IF3 + NVX
         JJ  = 0 
         DO J=1,5*N-4,5
                 JJ = JJ + 1
            WW(J  ) = V(JJ+IF )
            WW(J+1) = V(JJ+IF1)
            WW(J+2) = V(JJ+IF2)
            WW(J+3) = V(JJ+IF3)
            WW(J+4) = V(JJ+IF4)
         END DO

         IM1 = I - 1
         IF (N .GT. 2) THEN
             DO KR=1,KRUN
                   K = K_TRUE(KR)
                  K0 = IX(K)
                 KP1 = K + 1
                SUMA = 0D0
                SUMB = 0D0
                SUMC = 0D0
                SUMD = 0D0
                SUME = 0D0
                 KK0 = 5*N + 1

C               ////////////////////
C               Start innermostloop.
                 KKK = KK0
                DO KK=N,KP1,-1
                   KKK  = KKK - 5
                   SUMA = SUMA + AU(KK+K0)*WW(KKK  )
                   SUMB = SUMB + AU(KK+K0)*WW(KKK+1)
                   SUMC = SUMC + AU(KK+K0)*WW(KKK+2)
                   SUMD = SUMD + AU(KK+K0)*WW(KKK+3)
                   SUME = SUME + AU(KK+K0)*WW(KKK+4)
                END DO
C               End innermostloop.
C               //////////////////

                AUI = 1D0*AU_NEQ_0_I(KR)
                 SA = - SUMA*AUI
                 SB = - SUMB*AUI
                 SC = - SUMC*AUI
                 SD = - SUMD*AUI 
                 SE = - SUME*AUI 

C               ////////////////////
C               Start innermostloop.
                 KKK = KK0
                DO KK=N,KP1,-1
                         KKK = KKK - 5
                   WW(KKK  ) = WW(KKK  ) + SA*AU(KK+K0)
                   WW(KKK+1) = WW(KKK+1) + SB*AU(KK+K0)
                   WW(KKK+2) = WW(KKK+2) + SC*AU(KK+K0)
                   WW(KKK+3) = WW(KKK+3) + SD*AU(KK+K0)
                   WW(KKK+4) = WW(KKK+4) + SE*AU(KK+K0)
                END DO
C               End innermostloop.
C               //////////////////

             END DO
         END IF

         DO II=1,5
            III = I + II - 1
            DO J=IG(II),III-1
               IF (ABS(E(J) - E(III)) .LT. EPS2) THEN
                   IG(II) = J
                   GO TO 70
               END IF
            END DO
            IG(II) = III
   70       CONTINUE

            IF (IG(II) .NE. III  .AND.  SCHMIDT) THEN

C               Degenerate eigenvalues.  First, orthogonalize.
                KF = (IG(II)-2)*NVX
                DO K=IG(II),III-1
                   KF  = KF + NVX
                   SUM = 0D0
                    JF = -5
                   DO J=1,N
                      JF  = JF + 5
                      SUM = SUM + V(J+KF)*WW(JF+II)
                   END DO
                    S = -SUM
                   JF = -5
                   DO J=1,N
                      JF = JF + 5
                      WW(JF+II) = WW(JF+II) + S*V(J+KF)
                   END DO
                END DO
            END IF
         END DO

C        Normalization.
         SUM0 = 0D0
         SUM1 = 0D0
         SUM2 = 0D0
         SUM3 = 0D0
         SUM4 = 0D0
         JF  = -5
         DO J=1,N
              JF = JF + 5
            SUM0 = SUM0 + WW(JF+1)*WW(JF+1)
            SUM1 = SUM1 + WW(JF+2)*WW(JF+2)
            SUM2 = SUM2 + WW(JF+3)*WW(JF+3)
            SUM3 = SUM3 + WW(JF+4)*WW(JF+4)
            SUM4 = SUM4 + WW(JF+5)*WW(JF+5)
         END DO
         S0 = 1D0/SQRT(SUM0)
         S1 = 1D0/SQRT(SUM1)
         S2 = 1D0/SQRT(SUM2)
         S3 = 1D0/SQRT(SUM3)
         S4 = 1D0/SQRT(SUM4)
         JF = -5
         DO J=1,N
                  JF = JF + 5
            V(J+IF ) = S0*WW(JF+1)
            V(J+IF1) = S1*WW(JF+2)
            V(J+IF2) = S2*WW(JF+3)
            V(J+IF3) = S3*WW(JF+4)
            V(J+IF4) = S4*WW(JF+5)
         END DO
      END DO
      END
C----------------------------------------------------------------------
      SUBROUTINE BACK_TRANSFORMATION_U6 (N_NL,IV_NL,LV_NL,SCHMIDT,EPS2,
     *                                   AU,E,NVX,IF_NL,KRUN_NL,
     *                                   NRVEC_NL,
     *                                   W6,WW,AU_NEQ_0_I,IX,IG,K_TRUE,
     *                                                               V)
C----------------------------------------------------------------------

C     Array meanings explained in SUBROUTINE DRIVE_BACK_TRANSFORMATION.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SCHMIDT
      DIMENSION AU(*),E(*),V(*),
     *          W6(*),WW(*),AU_NEQ_0_I(*),IX(*),IG(*),K_TRUE(*)
  
      N = N_NL
      IV = IV_NL
      LV = LV_NL
      IF = IF_NL
      KRUN = KRUN_NL
      NRVEC = NRVEC_NL

C     //////////////////////////////////
C     Start full-fledged unrolled loops.
      IF5 = IF
      DO I=IV+NRVEC,LV,6

         IF  = IF5 + NVX
         IF1 = IF  + NVX
         IF2 = IF1 + NVX
         IF3 = IF2 + NVX
         IF4 = IF3 + NVX
         IF5 = IF4 + NVX
         JJ  = 0 
         DO J=1,6*N-5,6
                 JJ = JJ + 1
            WW(J  ) = V(JJ+IF )
            WW(J+1) = V(JJ+IF1)
            WW(J+2) = V(JJ+IF2)
            WW(J+3) = V(JJ+IF3)
            WW(J+4) = V(JJ+IF4)
            WW(J+5) = V(JJ+IF5)
         END DO

         IM1 = I - 1
         IF (N .GT. 2) THEN
             DO KR=1,KRUN
                   K = K_TRUE(KR)
                  K0 = IX(K)
                 KP1 = K + 1
                SUMA = 0D0
                SUMB = 0D0
                SUMC = 0D0
                SUMD = 0D0
                SUME = 0D0
                SUMF = 0D0
                 KK0 = 6*N + 1

C               ////////////////////
C               Start innermostloop.
                 KKK = KK0
                DO KK=N,KP1,-1
                   KKK  = KKK - 6
                   SUMA = SUMA + AU(KK+K0)*WW(KKK  )
                   SUMB = SUMB + AU(KK+K0)*WW(KKK+1)
                   SUMC = SUMC + AU(KK+K0)*WW(KKK+2)
                   SUMD = SUMD + AU(KK+K0)*WW(KKK+3)
                   SUME = SUME + AU(KK+K0)*WW(KKK+4)
                   SUMF = SUMF + AU(KK+K0)*WW(KKK+5)
                END DO
C               End innermostloop.
C               //////////////////

                AUI = 1D0*AU_NEQ_0_I(KR)
                 SA = - SUMA*AUI
                 SB = - SUMB*AUI
                 SC = - SUMC*AUI
                 SD = - SUMD*AUI 
                 SE = - SUME*AUI 
                 SF = - SUMF*AUI

C               ////////////////////
C               Start innermostloop.
                 KKK = KK0
                DO KK=N,KP1,-1
                         KKK = KKK - 6
                   WW(KKK  ) = WW(KKK  ) + SA*AU(KK+K0)
                   WW(KKK+1) = WW(KKK+1) + SB*AU(KK+K0)
                   WW(KKK+2) = WW(KKK+2) + SC*AU(KK+K0)
                   WW(KKK+3) = WW(KKK+3) + SD*AU(KK+K0)
                   WW(KKK+4) = WW(KKK+4) + SE*AU(KK+K0)
                   WW(KKK+5) = WW(KKK+5) + SF*AU(KK+K0)
                END DO
C               End innermostloop.
C               //////////////////

             END DO
         END IF

         DO II=1,6
            III = I + II - 1
            DO J=IG(II),III-1
               IF (ABS(E(J) - E(III)) .LT. EPS2) THEN
                   IG(II) = J
                   GO TO 70
               END IF
            END DO
            IG(II) = III
   70       CONTINUE

            IF (IG(II) .NE. III  .AND.  SCHMIDT) THEN

C               Degenerate eigenvalues.  First, orthogonalize.
                KF = (IG(II)-2)*NVX
                DO K=IG(II),III-1
                   KF  = KF + NVX
                   SUM = 0D0
                    JF = -6
                   DO J=1,N
                      JF  = JF + 6
                      SUM = SUM + V(J+KF)*WW(JF+II)
                   END DO
                    S = -SUM
                   JF = -6
                   DO J=1,N
                      JF = JF + 6
                      WW(JF+II) = WW(JF+II) + S*V(J+KF)
                   END DO
                END DO
            END IF
         END DO

C        Normalization.
         SUM0 = 0D0
         SUM1 = 0D0
         SUM2 = 0D0
         SUM3 = 0D0
         SUM4 = 0D0
         SUM5 = 0D0
         JF  = -6
         DO J=1,N
              JF = JF + 6
            SUM0 = SUM0 + WW(JF+1)*WW(JF+1)
            SUM1 = SUM1 + WW(JF+2)*WW(JF+2)
            SUM2 = SUM2 + WW(JF+3)*WW(JF+3)
            SUM3 = SUM3 + WW(JF+4)*WW(JF+4)
            SUM4 = SUM4 + WW(JF+5)*WW(JF+5)
            SUM5 = SUM5 + WW(JF+6)*WW(JF+6)
         END DO
         S0 = 1D0/SQRT(SUM0)
         S1 = 1D0/SQRT(SUM1)
         S2 = 1D0/SQRT(SUM2)
         S3 = 1D0/SQRT(SUM3)
         S4 = 1D0/SQRT(SUM4)
         S5 = 1D0/SQRT(SUM5)
         JF = -6
         DO J=1,N
                  JF = JF + 6
            V(J+IF ) = S0*WW(JF+1)
            V(J+IF1) = S1*WW(JF+2)
            V(J+IF2) = S2*WW(JF+3)
            V(J+IF3) = S3*WW(JF+4)
            V(J+IF4) = S4*WW(JF+5)
            V(J+IF5) = S5*WW(JF+6)
         END DO
      END DO
      END
C----------------------------------------------------------------------
      SUBROUTINE BACK_TRANSFORMATION_U8 (N_NL,IV_NL,LV_NL,SCHMIDT,EPS2,
     *                                   AU,E,NVX,IF_NL,KRUN_NL,
     *                                   NRVEC_NL,
     *                                   W6,WW,AU_NEQ_0_I,IX,IG,K_TRUE,
     *                                                               V)
C----------------------------------------------------------------------

C     Array meanings explained in SUBROUTINE DRIVE_BACK_TRANSFORMATION.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SCHMIDT
      DIMENSION AU(*),E(*),V(*),
     *          W6(*),WW(*),AU_NEQ_0_I(*),IX(*),IG(*),K_TRUE(*)
  
      N = N_NL
      IV = IV_NL
      LV = LV_NL
      IF = IF_NL
      KRUN = KRUN_NL
      NRVEC = NRVEC_NL

C     //////////////////////////////////
C     Start full-fledged unrolled loops.
      IF7 = IF
      DO I=IV+NRVEC,LV,8

         IF  = IF7 + NVX
         IF1 = IF  + NVX
         IF2 = IF1 + NVX
         IF3 = IF2 + NVX
         IF4 = IF3 + NVX
         IF5 = IF4 + NVX
         IF6 = IF5 + NVX
         IF7 = IF6 + NVX
         JJ  = 0 
         DO J=1,8*N-7,8
                 JJ = JJ + 1
            WW(J  ) = V(JJ+IF )
            WW(J+1) = V(JJ+IF1)
            WW(J+2) = V(JJ+IF2)
            WW(J+3) = V(JJ+IF3)
            WW(J+4) = V(JJ+IF4)
            WW(J+5) = V(JJ+IF5)
            WW(J+6) = V(JJ+IF6)
            WW(J+7) = V(JJ+IF7)
         END DO

         IM1 = I - 1
         IF (N .GT. 2) THEN
             DO KR=1,KRUN
                K = K_TRUE(KR)
                K0 = IX(K)
                KP1 = K + 1
                SUMA = 0D0
                SUMB = 0D0
                SUMC = 0D0
                SUMD = 0D0
                SUME = 0D0
                SUMF = 0D0
                SUMG = 0D0
                SUMH = 0D0
                 KK0 = 8*N + 1

C               ////////////////////
C               Start innermostloop.
                 KKK = KK0
                DO KK=N,KP1,-1
                   KKK  = KKK - 8
                   SUMA = SUMA + AU(KK+K0)*WW(KKK  )
                   SUMB = SUMB + AU(KK+K0)*WW(KKK+1)
                   SUMC = SUMC + AU(KK+K0)*WW(KKK+2)
                   SUMD = SUMD + AU(KK+K0)*WW(KKK+3)
                   SUME = SUME + AU(KK+K0)*WW(KKK+4)
                   SUMF = SUMF + AU(KK+K0)*WW(KKK+5)
                   SUMG = SUMG + AU(KK+K0)*WW(KKK+6)
                   SUMH = SUMH + AU(KK+K0)*WW(KKK+7)
                END DO
C               End innermostloop.
C               //////////////////

                AUI = 1D0*AU_NEQ_0_I(KR)
                SA = - SUMA*AUI
                SB = - SUMB*AUI
                SC = - SUMC*AUI
                SD = - SUMD*AUI
                SE = - SUME*AUI
                SF = - SUMF*AUI
                SG = - SUMG*AUI
                SH = - SUMH*AUI

C               ////////////////////
C               Start innermostloop.
                 KKK = KK0
                DO KK=N,KP1,-1
                         KKK = KKK - 8
                   WW(KKK  ) = WW(KKK  ) + SA*AU(KK+K0)
                   WW(KKK+1) = WW(KKK+1) + SB*AU(KK+K0)
                   WW(KKK+2) = WW(KKK+2) + SC*AU(KK+K0)
                   WW(KKK+3) = WW(KKK+3) + SD*AU(KK+K0)
                   WW(KKK+4) = WW(KKK+4) + SE*AU(KK+K0)
                   WW(KKK+5) = WW(KKK+5) + SF*AU(KK+K0)
                   WW(KKK+6) = WW(KKK+6) + SG*AU(KK+K0)
                   WW(KKK+7) = WW(KKK+7) + SH*AU(KK+K0)
                END DO
C               End innermostloop.
C               //////////////////

             END DO
         END IF

         DO II=1,8
            III = I + II - 1
            DO J=IG(II),III-1
               IF (ABS(E(J) - E(III)) .LT. EPS2) THEN
                   IG(II) = J
                   GO TO 70
               END IF
            END DO
            IG(II) = III
   70       CONTINUE

            IF (IG(II) .NE. III  .AND.  SCHMIDT) THEN

C               Degenerate eigenvalues.  First, orthogonalize.
                KF = (IG(II)-2)*NVX
                DO K=IG(II),III-1
                   KF  = KF + NVX
                   SUM = 0D0
                    JF = -8
                   DO J=1,N
                      JF  = JF + 8
                      SUM = SUM + V(J+KF)*WW(JF+II)
                   END DO
                    S = -SUM
                   JF = -8
                   DO J=1,N
                             JF = JF + 8
                      WW(JF+II) = WW(JF+II) + S*V(J+KF)
                   END DO
                END DO
            END IF
         END DO

C        Normalization.
         SUM0 = 0D0
         SUM1 = 0D0
         SUM2 = 0D0
         SUM3 = 0D0
         SUM4 = 0D0
         SUM5 = 0D0
         SUM6 = 0D0
         SUM7 = 0D0
          JF  = -8
         DO J=1,N
            JF = JF + 8
            SUM0 = SUM0 + WW(JF+1)*WW(JF+1)
            SUM1 = SUM1 + WW(JF+2)*WW(JF+2)
            SUM2 = SUM2 + WW(JF+3)*WW(JF+3)
            SUM3 = SUM3 + WW(JF+4)*WW(JF+4)
            SUM4 = SUM4 + WW(JF+5)*WW(JF+5)
            SUM5 = SUM5 + WW(JF+6)*WW(JF+6)
            SUM6 = SUM6 + WW(JF+7)*WW(JF+7)
            SUM7 = SUM7 + WW(JF+8)*WW(JF+8)
         END DO
         S0 = 1D0/SQRT(SUM0)
         S1 = 1D0/SQRT(SUM1)
         S2 = 1D0/SQRT(SUM2)
         S3 = 1D0/SQRT(SUM3)
         S4 = 1D0/SQRT(SUM4)
         S5 = 1D0/SQRT(SUM5)
         S6 = 1D0/SQRT(SUM6)
         S7 = 1D0/SQRT(SUM7)
         JF = -8
         DO J=1,N
            JF = JF + 8
            V(J+IF ) = S0*WW(JF+1)
            V(J+IF1) = S1*WW(JF+2)
            V(J+IF2) = S2*WW(JF+3)
            V(J+IF3) = S3*WW(JF+4)
            V(J+IF4) = S4*WW(JF+5)
            V(J+IF5) = S5*WW(JF+6)
            V(J+IF6) = S6*WW(JF+7)
            V(J+IF7) = S7*WW(JF+8)
         END DO
      END DO
      END
C----------------------------------------------------------------------
      SUBROUTINE BACK_TRANSFORMATION_U10(N_NL,IV_NL,LV_NL,SCHMIDT,EPS2,
     *                                   AU,E,NVX,IF_NL,KRUN_NL,
     *                                   NRVEC_NL,
     *                                   W6,WW,AU_NEQ_0_I,IX,IG,K_TRUE,
     *                                                               V)
C----------------------------------------------------------------------

C     Array meanings explained in SUBROUTINE DRIVE_BACK_TRANSFORMATION.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SCHMIDT
      DIMENSION AU(*),E(*),V(*),
     *          W6(*),WW(*),AU_NEQ_0_I(*),IX(*),IG(*),K_TRUE(*)
  
      N = N_NL
      IV = IV_NL
      LV = LV_NL
      IF = IF_NL
      KRUN = KRUN_NL
      NRVEC = NRVEC_NL

C     //////////////////////////////////
C     Start full-fledged unrolled loops.
      IF9 = IF
      DO I=IV+NRVEC,LV,10

         IF   = IF9  + NVX
         IF1  = IF   + NVX
         IF2  = IF1  + NVX
         IF3  = IF2  + NVX
         IF4  = IF3  + NVX
         IF5  = IF4  + NVX
         IF6  = IF5  + NVX
         IF7  = IF6  + NVX
         IF8  = IF7  + NVX
         IF9  = IF8  + NVX
          JJ  = 0 
         DO J=1,10*N-9,10
                  JJ = JJ + 1
            WW(J   ) = V(JJ+IF  )
            WW(J+1 ) = V(JJ+IF1 )
            WW(J+2 ) = V(JJ+IF2 )
            WW(J+3 ) = V(JJ+IF3 )
            WW(J+4 ) = V(JJ+IF4 )
            WW(J+5 ) = V(JJ+IF5 )
            WW(J+6 ) = V(JJ+IF6 )
            WW(J+7 ) = V(JJ+IF7 )
            WW(J+8 ) = V(JJ+IF8 )
            WW(J+9 ) = V(JJ+IF9 )
         END DO

         IM1 = I - 1
         IF (N .GT. 2) THEN
             DO KR=1,KRUN
                K = K_TRUE(KR)
                K0 = IX(K)
                KP1 = K + 1
                SUMA = 0D0
                SUMB = 0D0
                SUMC = 0D0
                SUMD = 0D0
                SUME = 0D0
                SUMF = 0D0
                SUMG = 0D0
                SUMH = 0D0
                SUMI = 0D0
                SUMJ = 0D0
                 KK0 = 10*N + 1

C               ////////////////////
C               Start innermostloop.
                 KKK = KK0
                DO KK=N,KP1,-1
                   KKK = KKK - 10
                  SUMA = SUMA + AU(KK+K0)*WW(KKK  )
                  SUMB = SUMB + AU(KK+K0)*WW(KKK+1)
                  SUMC = SUMC + AU(KK+K0)*WW(KKK+2)
                  SUMD = SUMD + AU(KK+K0)*WW(KKK+3)
                  SUME = SUME + AU(KK+K0)*WW(KKK+4)
                  SUMF = SUMF + AU(KK+K0)*WW(KKK+5)
                  SUMG = SUMG + AU(KK+K0)*WW(KKK+6)
                  SUMH = SUMH + AU(KK+K0)*WW(KKK+7)
                  SUMI = SUMI + AU(KK+K0)*WW(KKK+8)
                  SUMJ = SUMJ + AU(KK+K0)*WW(KKK+9)
                END DO
C               End innermostloop.
C               //////////////////

                AUI = 1D0*AU_NEQ_0_I(KR)
                SA = - SUMA*AUI
                SB = - SUMB*AUI
                SC = - SUMC*AUI
                SD = - SUMD*AUI
                SE = - SUME*AUI
                SF = - SUMF*AUI
                SG = - SUMG*AUI
                SH = - SUMH*AUI
                SI = - SUMI*AUI
                SJ = - SUMJ*AUI

C               ////////////////////
C               Start innermostloop.
                 KKK = KK0
                DO KK=N,KP1,-1
                          KKK = KKK - 10
                   WW(KKK   ) = WW(KKK   ) + SA*AU(KK+K0)
                   WW(KKK+1 ) = WW(KKK+1 ) + SB*AU(KK+K0)
                   WW(KKK+2 ) = WW(KKK+2 ) + SC*AU(KK+K0)
                   WW(KKK+3 ) = WW(KKK+3 ) + SD*AU(KK+K0)
                   WW(KKK+4 ) = WW(KKK+4 ) + SE*AU(KK+K0)
                   WW(KKK+5 ) = WW(KKK+5 ) + SF*AU(KK+K0)
                   WW(KKK+6 ) = WW(KKK+6 ) + SG*AU(KK+K0)
                   WW(KKK+7 ) = WW(KKK+7 ) + SH*AU(KK+K0)
                   WW(KKK+8 ) = WW(KKK+8 ) + SI*AU(KK+K0)
                   WW(KKK+9 ) = WW(KKK+9 ) + SJ*AU(KK+K0)
                END DO
C               End innermostloop.
C               //////////////////

             END DO
         END IF

         DO II=1,10
            III = I + II - 1
            DO J=IG(II),III-1
               IF (ABS(E(J) - E(III)) .LT. EPS2) THEN
                   IG(II) = J
                   GO TO 70
               END IF
            END DO
            IG(II) = III
   70       CONTINUE

            IF (IG(II) .NE. III  .AND.  SCHMIDT) THEN

C               Degenerate eigenvalues.  First, orthogonalize.
                KF = (IG(II)-2)*NVX
                DO K=IG(II),III-1
                   KF  = KF + NVX
                   SUM = 0D0
                   JF = -10
                   DO J=1,N
                      JF  = JF + 10
                      SUM = SUM + V(J+KF)*WW(JF+II)
                   END DO
                    S = -SUM
                   JF = -10
                   DO J=1,N
                             JF = JF + 10
                      WW(JF+II) = WW(JF+II) + S*V(J+KF)
                   END DO
                END DO
            END IF
         END DO

C        Normalization.
         SUM0 = 0D0
         SUM1 = 0D0
         SUM2 = 0D0
         SUM3 = 0D0
         SUM4 = 0D0
         SUM5 = 0D0
         SUM6 = 0D0
         SUM7 = 0D0
         SUM8 = 0D0
         SUM9 = 0D0
         JF  = -10
         DO J=1,N
            JF = JF + 10
            SUM0 = SUM0 + WW(JF+1) *WW(JF+1)
            SUM1 = SUM1 + WW(JF+2) *WW(JF+2)
            SUM2 = SUM2 + WW(JF+3) *WW(JF+3)
            SUM3 = SUM3 + WW(JF+4) *WW(JF+4)
            SUM4 = SUM4 + WW(JF+5) *WW(JF+5)
            SUM5 = SUM5 + WW(JF+6) *WW(JF+6)
            SUM6 = SUM6 + WW(JF+7) *WW(JF+7)
            SUM7 = SUM7 + WW(JF+8) *WW(JF+8)
            SUM8 = SUM8 + WW(JF+9) *WW(JF+9)
            SUM9 = SUM9 + WW(JF+10)*WW(JF+10)
         END DO
         S0 = 1D0/SQRT(SUM0)
         S1 = 1D0/SQRT(SUM1)
         S2 = 1D0/SQRT(SUM2)
         S3 = 1D0/SQRT(SUM3)
         S4 = 1D0/SQRT(SUM4)
         S5 = 1D0/SQRT(SUM5)
         S6 = 1D0/SQRT(SUM6)
         S7 = 1D0/SQRT(SUM7)
         S8 = 1D0/SQRT(SUM8)
         S9 = 1D0/SQRT(SUM9)
         JF = -10
         DO J=1,N
            JF = JF + 10
            V(J+IF ) = S0*WW(JF+1)
            V(J+IF1) = S1*WW(JF+2)
            V(J+IF2) = S2*WW(JF+3)
            V(J+IF3) = S3*WW(JF+4)
            V(J+IF4) = S4*WW(JF+5)
            V(J+IF5) = S5*WW(JF+6)
            V(J+IF6) = S6*WW(JF+7)
            V(J+IF7) = S7*WW(JF+8)
            V(J+IF8) = S8*WW(JF+9)
            V(J+IF9) = S9*WW(JF+10)
         END DO
      END DO
      END
C----------------------------------------------------------------------
      SUBROUTINE BACK_TRANSFORMATION_U12(N_NL,IV_NL,LV_NL,SCHMIDT,EPS2,
     *                                   AU,E,NVX,IF_NL,KRUN_NL,
     *                                   NRVEC_NL,
     *                                   W6,WW,AU_NEQ_0_I,IX,IG,K_TRUE,
     *                                                               V)
C----------------------------------------------------------------------

C     Array meanings explained in SUBROUTINE DRIVE_BACK_TRANSFORMATION.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SCHMIDT
      DIMENSION AU(*),E(*),V(*),
     *          W6(*),WW(*),AU_NEQ_0_I(*),IX(*),IG(*),K_TRUE(*)
  
      N = N_NL
      IV = IV_NL
      LV = LV_NL
      IF = IF_NL
      KRUN = KRUN_NL
      NRVEC = NRVEC_NL

C     //////////////////////////////////
C     Start full-fledged unrolled loops.
      IF11 = IF
      DO I=IV+NRVEC,LV,12

         IF   = IF11 + NVX
         IF1  = IF   + NVX
         IF2  = IF1  + NVX
         IF3  = IF2  + NVX
         IF4  = IF3  + NVX
         IF5  = IF4  + NVX
         IF6  = IF5  + NVX
         IF7  = IF6  + NVX
         IF8  = IF7  + NVX
         IF9  = IF8  + NVX
         IF10 = IF9  + NVX
         IF11 = IF10 + NVX
           JJ = 0 
         DO J=1,12*N-11,12
            JJ = JJ + 1
            WW(J   ) = V(JJ+IF  )
            WW(J+1 ) = V(JJ+IF1 )
            WW(J+2 ) = V(JJ+IF2 )
            WW(J+3 ) = V(JJ+IF3 )
            WW(J+4 ) = V(JJ+IF4 )
            WW(J+5 ) = V(JJ+IF5 )
            WW(J+6 ) = V(JJ+IF6 )
            WW(J+7 ) = V(JJ+IF7 )
            WW(J+8 ) = V(JJ+IF8 )
            WW(J+9 ) = V(JJ+IF9 )
            WW(J+10) = V(JJ+IF10)
            WW(J+11) = V(JJ+IF11)
         END DO

         IM1 = I - 1
         IF (N .GT. 2) THEN
             DO KR=1,KRUN
                K = K_TRUE(KR)
                K0 = IX(K)
                KP1 = K + 1
                SUMA = 0D0
                SUMB = 0D0
                SUMC = 0D0
                SUMD = 0D0
                SUME = 0D0
                SUMF = 0D0
                SUMG = 0D0
                SUMH = 0D0
                SUMI = 0D0
                SUMJ = 0D0
                SUMK = 0D0
                SUML = 0D0
                 KK0 = 12*N + 1

C               ////////////////////
C               Start innermostloop.
                 KKK = KK0
                DO KK=N,KP1,-1
                   KKK  = KKK - 12
                   SUMA = SUMA + AU(KK+K0)*WW(KKK   )
                   SUMB = SUMB + AU(KK+K0)*WW(KKK+1 )
                   SUMC = SUMC + AU(KK+K0)*WW(KKK+2 )
                   SUMD = SUMD + AU(KK+K0)*WW(KKK+3 )
                   SUME = SUME + AU(KK+K0)*WW(KKK+4 )
                   SUMF = SUMF + AU(KK+K0)*WW(KKK+5 )
                   SUMG = SUMG + AU(KK+K0)*WW(KKK+6 )
                   SUMH = SUMH + AU(KK+K0)*WW(KKK+7 )
                   SUMI = SUMI + AU(KK+K0)*WW(KKK+8 )
                   SUMJ = SUMJ + AU(KK+K0)*WW(KKK+9 )
                   SUMK = SUMK + AU(KK+K0)*WW(KKK+10)
                   SUML = SUML + AU(KK+K0)*WW(KKK+11)
                END DO
C               End innermostloop.
C               //////////////////

                AUI = 1D0*AU_NEQ_0_I(KR)
                 SA = - SUMA*AUI
                 SB = - SUMB*AUI
                 SC = - SUMC*AUI
                 SD = - SUMD*AUI
                 SE = - SUME*AUI
                 SF = - SUMF*AUI
                 SG = - SUMG*AUI
                 SH = - SUMH*AUI
                 SI = - SUMI*AUI
                 SJ = - SUMJ*AUI
                 SK = - SUMK*AUI
                 SL = - SUML*AUI

C               ////////////////////
C               Start innermostloop.
                 KKK = KK0
                DO KK=N,KP1,-1
                          KKK = KKK - 12
                   WW(KKK   ) = WW(KKK   ) + SA*AU(KK+K0)
                   WW(KKK+1 ) = WW(KKK+1 ) + SB*AU(KK+K0)
                   WW(KKK+2 ) = WW(KKK+2 ) + SC*AU(KK+K0)
                   WW(KKK+3 ) = WW(KKK+3 ) + SD*AU(KK+K0)
                   WW(KKK+4 ) = WW(KKK+4 ) + SE*AU(KK+K0)
                   WW(KKK+5 ) = WW(KKK+5 ) + SF*AU(KK+K0)
                   WW(KKK+6 ) = WW(KKK+6 ) + SG*AU(KK+K0)
                   WW(KKK+7 ) = WW(KKK+7 ) + SH*AU(KK+K0)
                   WW(KKK+8 ) = WW(KKK+8 ) + SI*AU(KK+K0)
                   WW(KKK+9 ) = WW(KKK+9 ) + SJ*AU(KK+K0)
                   WW(KKK+10) = WW(KKK+10) + SK*AU(KK+K0)
                   WW(KKK+11) = WW(KKK+11) + SL*AU(KK+K0)
                END DO
C               End innermostloop.
C               //////////////////

             END DO
         END IF

         DO II=1,12
            III = I + II - 1
            DO J=IG(II),III-1
               IF (ABS(E(J) - E(III)) .LT. EPS2) THEN
                   IG(II) = J
                   GO TO 70
               END IF
            END DO
            IG(II) = III
   70       CONTINUE

            IF (IG(II) .NE. III  .AND.  SCHMIDT) THEN

C               Degenerate eigenvalues.  First, orthogonalize.
                KF = (IG(II)-2)*NVX
                DO K=IG(II),III-1
                   KF  = KF + NVX
                   SUM = 0D0
                    JF = -12
                   DO J=1,N
                      JF  = JF + 12
                      SUM = SUM + V(J+KF)*WW(JF+II)
                   END DO
                    S = -SUM
                   JF = -12
                   DO J=1,N
                      JF = JF + 12
                      WW(JF+II) = WW(JF+II) + S*V(J+KF)
                   END DO
                END DO
            END IF
         END DO

C        Normalization.
         SUM0 = 0D0
         SUM1 = 0D0
         SUM2 = 0D0
         SUM3 = 0D0
         SUM4 = 0D0
         SUM5 = 0D0
         SUM6 = 0D0
         SUM7 = 0D0
         SUM8 = 0D0
         SUM9 = 0D0
         SUM10= 0D0
         SUM11= 0D0
         JF  = -12
         DO J=1,N
            JF = JF + 12
            SUM0 = SUM0 + WW(JF+1) *WW(JF+1)
            SUM1 = SUM1 + WW(JF+2) *WW(JF+2)
            SUM2 = SUM2 + WW(JF+3) *WW(JF+3)
            SUM3 = SUM3 + WW(JF+4) *WW(JF+4)
            SUM4 = SUM4 + WW(JF+5) *WW(JF+5)
            SUM5 = SUM5 + WW(JF+6) *WW(JF+6)
            SUM6 = SUM6 + WW(JF+7) *WW(JF+7)
            SUM7 = SUM7 + WW(JF+8) *WW(JF+8)
            SUM8 = SUM8 + WW(JF+9) *WW(JF+9)
            SUM9 = SUM9 + WW(JF+10)*WW(JF+10)
            SUM10= SUM10+ WW(JF+11)*WW(JF+11)
            SUM11= SUM11+ WW(JF+12)*WW(JF+12)
         END DO
         S0 = 1D0/SQRT(SUM0)
         S1 = 1D0/SQRT(SUM1)
         S2 = 1D0/SQRT(SUM2)
         S3 = 1D0/SQRT(SUM3)
         S4 = 1D0/SQRT(SUM4)
         S5 = 1D0/SQRT(SUM5)
         S6 = 1D0/SQRT(SUM6)
         S7 = 1D0/SQRT(SUM7)
         S8 = 1D0/SQRT(SUM8)
         S9 = 1D0/SQRT(SUM9)
         S10= 1D0/SQRT(SUM10)
         S11= 1D0/SQRT(SUM11)
         JF = -12
         DO J=1,N
            JF = JF + 12
            V(J+IF ) = S0*WW(JF+1)
            V(J+IF1) = S1*WW(JF+2)
            V(J+IF2) = S2*WW(JF+3)
            V(J+IF3) = S3*WW(JF+4)
            V(J+IF4) = S4*WW(JF+5)
            V(J+IF5) = S5*WW(JF+6)
            V(J+IF6) = S6*WW(JF+7)
            V(J+IF7) = S7*WW(JF+8)
            V(J+IF8) = S8*WW(JF+9)
            V(J+IF9) = S9*WW(JF+10)
            V(J+IF10)=S10*WW(JF+11)
            V(J+IF11)=S11*WW(JF+12)
         END DO
      END DO
      END
C----------------------------------------------------------------------
      SUBROUTINE BACK_TRANSFORMATION_U16(N_NL,IV_NL,LV_NL,SCHMIDT,EPS2,
     *                                   AU,E,NVX,IF_NL,KRUN_NL,
     *                                   NRVEC_NL,
     *                                   W6,WW,AU_NEQ_0_I,IX,IG,K_TRUE,
     *                                                               V)
C----------------------------------------------------------------------

C     Array meanings explained in SUBROUTINE DRIVE_BACK_TRANSFORMATION.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SCHMIDT
      DIMENSION AU(*),E(*),V(*),
     *          W6(*),WW(*),AU_NEQ_0_I(*),IX(*),IG(*),K_TRUE(*)

      N = N_NL
      IV = IV_NL
      LV = LV_NL
      IF = IF_NL
      KRUN = KRUN_NL
      NRVEC = NRVEC_NL

C     //////////////////////////////////
C     Start full-fledged unrolled loops.
      IF15 = IF
      DO I=IV+NRVEC,LV,16

         IF   = IF15 + NVX
         IF1  = IF   + NVX
         IF2  = IF1  + NVX
         IF3  = IF2  + NVX
         IF4  = IF3  + NVX
         IF5  = IF4  + NVX
         IF6  = IF5  + NVX
         IF7  = IF6  + NVX
         IF8  = IF7  + NVX
         IF9  = IF8  + NVX
         IF10 = IF9  + NVX
         IF11 = IF10 + NVX
         IF12 = IF11 + NVX
         IF13 = IF12 + NVX
         IF14 = IF13 + NVX
         IF15 = IF14 + NVX
           JJ = 0 
         DO J=1,16*N-15,16
            JJ = JJ + 1
            WW(J   ) = V(JJ+IF  )
            WW(J+1 ) = V(JJ+IF1 )
            WW(J+2 ) = V(JJ+IF2 )
            WW(J+3 ) = V(JJ+IF3 )
            WW(J+4 ) = V(JJ+IF4 )
            WW(J+5 ) = V(JJ+IF5 )
            WW(J+6 ) = V(JJ+IF6 )
            WW(J+7 ) = V(JJ+IF7 )
            WW(J+8 ) = V(JJ+IF8 )
            WW(J+9 ) = V(JJ+IF9 )
            WW(J+10) = V(JJ+IF10)
            WW(J+11) = V(JJ+IF11)
            WW(J+12) = V(JJ+IF12)
            WW(J+13) = V(JJ+IF13)
            WW(J+14) = V(JJ+IF14)
            WW(J+15) = V(JJ+IF15)
         END DO

         IM1 = I - 1
         IF (N .GT. 2) THEN
             DO KR=1,KRUN
                K = K_TRUE(KR)
                K0 = IX(K)
                KP1 = K + 1
                SUMA = 0D0
                SUMB = 0D0
                SUMC = 0D0
                SUMD = 0D0
                SUME = 0D0
                SUMF = 0D0
                SUMG = 0D0
                SUMH = 0D0
                SUMI = 0D0
                SUMJ = 0D0
                SUMK = 0D0
                SUML = 0D0
                SUMM = 0D0
                SUMN = 0D0
                SUMO = 0D0
                SUMP = 0D0
                 KK0 = 16*N + 1

C               ////////////////////
C               Start innermostloop.
                 KKK = KK0
                DO KK=N,KP1,-1
                   KKK  = KKK - 16
                   SUMA = SUMA + AU(KK+K0)*WW(KKK   )
                   SUMB = SUMB + AU(KK+K0)*WW(KKK+1 )
                   SUMC = SUMC + AU(KK+K0)*WW(KKK+2 )
                   SUMD = SUMD + AU(KK+K0)*WW(KKK+3 )
                   SUME = SUME + AU(KK+K0)*WW(KKK+4 )
                   SUMF = SUMF + AU(KK+K0)*WW(KKK+5 )
                   SUMG = SUMG + AU(KK+K0)*WW(KKK+6 )
                   SUMH = SUMH + AU(KK+K0)*WW(KKK+7 )
                   SUMI = SUMI + AU(KK+K0)*WW(KKK+8 )
                   SUMJ = SUMJ + AU(KK+K0)*WW(KKK+9 )
                   SUMK = SUMK + AU(KK+K0)*WW(KKK+10)
                   SUML = SUML + AU(KK+K0)*WW(KKK+11)
                   SUMM = SUMM + AU(KK+K0)*WW(KKK+12)
                   SUMN = SUMN + AU(KK+K0)*WW(KKK+13)
                   SUMO = SUMO + AU(KK+K0)*WW(KKK+14)
                   SUMP = SUMP + AU(KK+K0)*WW(KKK+15)
                END DO
C               End innermostloop.
C               //////////////////

                AUI = 1D0*AU_NEQ_0_I(KR)
                 SA = - SUMA*AUI
                 SB = - SUMB*AUI
                 SC = - SUMC*AUI
                 SD = - SUMD*AUI
                 SE = - SUME*AUI
                 SF = - SUMF*AUI
                 SG = - SUMG*AUI
                 SH = - SUMH*AUI
                 SI = - SUMI*AUI
                 SJ = - SUMJ*AUI
                 SK = - SUMK*AUI
                 SL = - SUML*AUI
                 SM = - SUMM*AUI
                 SN = - SUMN*AUI
                 SO = - SUMO*AUI
                 SP = - SUMP*AUI

C               ////////////////////
C               Start innermostloop.
                 KKK = KK0
                DO KK=N,KP1,-1
                          KKK = KKK - 16
                   WW(KKK   ) = WW(KKK   ) + SA*AU(KK+K0)
                   WW(KKK+1 ) = WW(KKK+1 ) + SB*AU(KK+K0)
                   WW(KKK+2 ) = WW(KKK+2 ) + SC*AU(KK+K0)
                   WW(KKK+3 ) = WW(KKK+3 ) + SD*AU(KK+K0)
                   WW(KKK+4 ) = WW(KKK+4 ) + SE*AU(KK+K0)
                   WW(KKK+5 ) = WW(KKK+5 ) + SF*AU(KK+K0)
                   WW(KKK+6 ) = WW(KKK+6 ) + SG*AU(KK+K0)
                   WW(KKK+7 ) = WW(KKK+7 ) + SH*AU(KK+K0)
                   WW(KKK+8 ) = WW(KKK+8 ) + SI*AU(KK+K0)
                   WW(KKK+9 ) = WW(KKK+9 ) + SJ*AU(KK+K0)
                   WW(KKK+10) = WW(KKK+10) + SK*AU(KK+K0)
                   WW(KKK+11) = WW(KKK+11) + SL*AU(KK+K0)
                   WW(KKK+12) = WW(KKK+12) + SM*AU(KK+K0)
                   WW(KKK+13) = WW(KKK+13) + SN*AU(KK+K0)
                   WW(KKK+14) = WW(KKK+14) + SO*AU(KK+K0)
                   WW(KKK+15) = WW(KKK+15) + SP*AU(KK+K0)
                END DO
C               End innermostloop.
C               //////////////////

             END DO
         END IF

         DO II=1,16
            III = I + II - 1
            DO J=IG(II),III-1
               IF (ABS(E(J) - E(III)) .LT. EPS2) THEN
                   IG(II) = J
                   GO TO 70
               END IF
            END DO
            IG(II) = III
   70       CONTINUE

            IF (IG(II) .NE. III  .AND.  SCHMIDT) THEN

C               Degenerate eigenvalues.  First, orthogonalize.
                KF = (IG(II)-2)*NVX
                DO K=IG(II),III-1
                   KF  = KF + NVX
                   SUM = 0D0
                   JF = -16
                   DO J=1,N
                      JF  = JF + 16
                      SUM = SUM + V(J+KF)*WW(JF+II)
                   END DO
                    S = -SUM
                   JF = -16
                   DO J=1,N
                      JF = JF + 16
                      WW(JF+II) = WW(JF+II) + S*V(J+KF)
                   END DO
                END DO
            END IF
         END DO

C        Normalization.
         SUM0 = 0D0
         SUM1 = 0D0
         SUM2 = 0D0
         SUM3 = 0D0
         SUM4 = 0D0
         SUM5 = 0D0
         SUM6 = 0D0
         SUM7 = 0D0
         SUM8 = 0D0
         SUM9 = 0D0
         SUM10= 0D0
         SUM11= 0D0
         SUM12= 0D0
         SUM13= 0D0
         SUM14= 0D0
         SUM15= 0D0
         JF  = -16
         DO J=1,N
            JF = JF + 16
            SUM0 = SUM0 + WW(JF+1) *WW(JF+1)
            SUM1 = SUM1 + WW(JF+2) *WW(JF+2)
            SUM2 = SUM2 + WW(JF+3) *WW(JF+3)
            SUM3 = SUM3 + WW(JF+4) *WW(JF+4)
            SUM4 = SUM4 + WW(JF+5) *WW(JF+5)
            SUM5 = SUM5 + WW(JF+6) *WW(JF+6)
            SUM6 = SUM6 + WW(JF+7) *WW(JF+7)
            SUM7 = SUM7 + WW(JF+8) *WW(JF+8)
            SUM8 = SUM8 + WW(JF+9) *WW(JF+9)
            SUM9 = SUM9 + WW(JF+10)*WW(JF+10)
            SUM10= SUM10+ WW(JF+11)*WW(JF+11)
            SUM11= SUM11+ WW(JF+12)*WW(JF+12)
            SUM12= SUM12+ WW(JF+13)*WW(JF+13)
            SUM13= SUM13+ WW(JF+14)*WW(JF+14)
            SUM14= SUM14+ WW(JF+15)*WW(JF+15)
            SUM15= SUM15+ WW(JF+16)*WW(JF+16)
         END DO
         S0 = 1D0/SQRT(SUM0)
         S1 = 1D0/SQRT(SUM1)
         S2 = 1D0/SQRT(SUM2)
         S3 = 1D0/SQRT(SUM3)
         S4 = 1D0/SQRT(SUM4)
         S5 = 1D0/SQRT(SUM5)
         S6 = 1D0/SQRT(SUM6)
         S7 = 1D0/SQRT(SUM7)
         S8 = 1D0/SQRT(SUM8)
         S9 = 1D0/SQRT(SUM9)
         S10= 1D0/SQRT(SUM10)
         S11= 1D0/SQRT(SUM11)
         S12= 1D0/SQRT(SUM12)
         S13= 1D0/SQRT(SUM13)
         S14= 1D0/SQRT(SUM14)
         S15= 1D0/SQRT(SUM15)
         JF = -16
         DO J=1,N
            JF = JF + 16
            V(J+IF ) = S0*WW(JF+1)
            V(J+IF1) = S1*WW(JF+2)
            V(J+IF2) = S2*WW(JF+3)
            V(J+IF3) = S3*WW(JF+4)
            V(J+IF4) = S4*WW(JF+5)
            V(J+IF5) = S5*WW(JF+6)
            V(J+IF6) = S6*WW(JF+7)
            V(J+IF7) = S7*WW(JF+8)
            V(J+IF8) = S8*WW(JF+9)
            V(J+IF9) = S9*WW(JF+10)
            V(J+IF10)=S10*WW(JF+11)
            V(J+IF11)=S11*WW(JF+12)
            V(J+IF12)=S12*WW(JF+13)
            V(J+IF13)=S13*WW(JF+14)
            V(J+IF14)=S14*WW(JF+15)
            V(J+IF15)=S15*WW(JF+16)
         END DO
      END DO
      END
