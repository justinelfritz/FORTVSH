c     Functions designed to test numerical results
c     of VSH module functions

      MODULE TESTS
      USE KINDS
      USE GLOBALS
      USE VSH
      IMPLICIT NONE

      CONTAINS

c  Analytic result for Test 1
        FUNCTION TEST1_ANALYTIC(NTH)
        IMPLICIT NONE
        INTEGER(KIND=i4) :: NTH,I
        REAL(KIND=dp) :: TH
        COMPLEX(KIND=dp), DIMENSION(NTH) :: TEST1_ANALYTIC

        DO I=1,NTH
          TH = (I-1)*pi/(NTH-1)
          TEST1_ANALYTIC(I)=3.d0*DSQRT(5.d0)/64.d0/pi**2*
     &      (DSIN(TH)*DCOS(TH))**2
        ENDDO
        RETURN
        END FUNCTION TEST1_ANALYTIC

c  Numerical result for Test 1 from VSH inner products
        FUNCTION TEST1_VSH(NTH)
        IMPLICIT NONE
        INTEGER(KIND=i4) :: NTH,I
        REAL(KIND=dp) :: TH,PH
        COMPLEX(KIND=dp), DIMENSION(NTH) :: TEST1_VSH

        PH = 0.123d0   !- arbitrary choice
        DO I=1,NTH
          TH = (I-1)*pi/(NTH-1)
          TEST1_VSH(I) = DOT(VSH_TOR(2,0,TH,PH),VSH_TOR(1,1,TH,PH))*
     &         DOT(VSH_POL_DN(1,-1,TH,PH),CONJG(VSH_POL_UP(0,0,TH,PH)))
        ENDDO
        RETURN
        END FUNCTION TEST1_VSH

c  Numerical result for Test 1 from Geppert & Wiebicke forms
        FUNCTION TEST1_GW(NTH)
        IMPLICIT NONE
        INTEGER(KIND=i4) :: NTH,I,I1,I2
        REAL(KIND=dp) :: TH,PH
        COMPLEX(KIND=dp) :: TEMP1, TEMP2
        COMPLEX(KIND=dp), DIMENSION(NTH) :: TEST1_GW

        PH = 0.123d0   !- arbitrary choice
        DO I=1,NTH
          TH = (I-1)*pi/(NTH-1)
          TEMP1 = DCMPLX(0.d0,0.d0)
          TEMP2 = DCMPLX(0.d0,0.d0)
          DO I1=1,3,2
            TEMP1 = TEMP1-((8.d0-I1*(I1+1.d0))/2.d0/DSQRT(12.d0))*
     &        GWI(2,0,1,1,I1,1)*SSH(I1,1,TH,PH)
          ENDDO
          I2 = 1
          TEMP2 = -(1.d0/DSQRT(3.d0))*GWI(1,-1,1,1,0,0)*
     &        CONJG(SSH(1,1,TH,PH))

          TEST1_GW(I) = TEMP1*TEMP2
        ENDDO
        RETURN
        END FUNCTION TEST1_GW

c  Analytic result for Test 2
        FUNCTION TEST2_ANALYTIC(NTH)
        IMPLICIT NONE
        INTEGER(KIND=i4) :: NTH,I
        REAL(KIND=dp) :: TH,PH
        COMPLEX(KIND=dp), DIMENSION(NTH) :: TEST2_ANALYTIC
        PH = 1.006d0   !- arbitrary choice
        DO I=1,NTH
          TH = (I-1)*pi/(NTH-1)
          TEST2_ANALYTIC(I)=15.d0*DSQRT(15.d0*7.d0)/128.d0/pi**2*
     &      (DSIN(TH)*DCOS(TH))**3*DSIN(TH)/DSQRT(2.d0)*
     &      exp(4.d0*j*PH)
        ENDDO
        RETURN
        END FUNCTION TEST2_ANALYTIC

c  Numerical result for Test 2 from VSH inner products
        FUNCTION TEST2_VSH(NTH)
        IMPLICIT NONE
        INTEGER(KIND=i4) :: NTH,I
        REAL(KIND=dp) :: TH,PH
        COMPLEX(KIND=dp), DIMENSION(NTH) :: TEST2_VSH

        PH = 1.006d0   !- arbitrary choice
        DO I=1,NTH
          TH = (I-1)*pi/(NTH-1)
          TEST2_VSH(I) = DOT(PVSH_POL(2,0,TH,PH),PVSH_TOR(3,2,TH,PH))*
     &         DOT(PVSH_RAD(1,0,TH,PH),CONJG(PVSH_RAD(2,-2,TH,PH)))
        ENDDO
        RETURN
        END FUNCTION TEST2_VSH

c  Numerical result for Test 2 from Geppert & Wiebicke forms
        FUNCTION TEST2_GW(NTH)
        IMPLICIT NONE
        INTEGER(KIND=i4) :: NTH,I,I1,I2
        REAL(KIND=dp) :: TH,PH
        COMPLEX(KIND=dp) :: TEMP1, TEMP2
        COMPLEX(KIND=dp), DIMENSION(NTH) :: TEST2_GW

        PH = 1.006d0   !- arbitrary choice
        DO I=1,NTH
          TH = (I-1)*pi/(NTH-1)
          TEMP1 = DCMPLX(0.d0,0.d0)
          TEMP2 = DCMPLX(0.d0,0.d0)
          DO I1=2,4,2
            TEMP1 = TEMP1+GWJ(3,2,2,0,I1,2)*SSH(I1,2,TH,PH)
          ENDDO
          TEMP1=-TEMP1*j/DSQRT(6.d0*12.d0)
          DO I2=1,3,2
            TEMP2 = TEMP2+GWI(1,0,I2,-2,2,-2)*
     &        CONJG(SSH(I2,-2,TH,PH))
          ENDDO
          TEST2_GW(I) = TEMP1*TEMP2
        ENDDO
        RETURN
        END FUNCTION TEST2_GW


C  Test ASSOC_LEGENDRE_ALL against individual ASSOC_LEGENDRE calls
C  Columns: L  M  X  P_batch  P_single  abs_diff
      SUBROUTINE BATCH_ALM_CONS(LMAX, NX, OUTUNIT)
      IMPLICIT NONE
      INTEGER(KIND=i4), INTENT(IN) :: LMAX, NX, OUTUNIT
      INTEGER(KIND=i4) :: L, M, IX, PSIZE
      REAL(KIND=dp) :: X, DX, P_SINGLE
      REAL(KIND=dp), ALLOCATABLE :: P_BATCH(:)

      PSIZE = (LMAX+1)*(LMAX+2)/2
      ALLOCATE(P_BATCH(PSIZE))
      DX = 1.8_dp / (NX - 1)

      WRITE(OUTUNIT, '(A)') '# L  M  X  P_batch  P_single  abs_diff'
      DO IX = 0, NX-1
        X = -0.9_dp + IX * DX
        CALL ASSOC_LEGENDRE_ALL(P_BATCH, LMAX, X)
        DO L = 0, LMAX
          DO M = 0, L
            P_SINGLE = ASSOC_LEGENDRE(L, M, X)
            WRITE(OUTUNIT, *) L, M, X,
     &          P_BATCH(PLM_INDEX(L, M)), P_SINGLE,
     &          ABS(P_BATCH(PLM_INDEX(L, M)) - P_SINGLE)
          END DO
        END DO
      END DO

      DEALLOCATE(P_BATCH)
      END SUBROUTINE BATCH_ALM_CONS


C  Test DDX_ASSOC_LEGENDRE_ALL against individual DDX_ASSOC_LEGENDRE calls
C  Columns: L  M  X  dP_batch  dP_single  abs_diff
      SUBROUTINE BATCH_DALM_CONS(LMAX, NX, OUTUNIT)
      IMPLICIT NONE
      INTEGER(KIND=i4), INTENT(IN) :: LMAX, NX, OUTUNIT
      INTEGER(KIND=i4) :: L, M, IX, PSIZE
      REAL(KIND=dp) :: X, DX, DP_SINGLE
      REAL(KIND=dp), ALLOCATABLE :: P_BATCH(:), DP_BATCH(:)

      PSIZE = (LMAX+1)*(LMAX+2)/2
      ALLOCATE(P_BATCH(PSIZE), DP_BATCH(PSIZE))
      DX = 1.8_dp / (NX - 1)

      WRITE(OUTUNIT, '(A)') '# L  M  X  dP_batch  dP_single  abs_diff'
      DO IX = 0, NX-1
        X = -0.9_dp + IX * DX
        CALL ASSOC_LEGENDRE_ALL(P_BATCH, LMAX, X)
        CALL DDX_ASSOC_LEGENDRE_ALL(DP_BATCH, P_BATCH, LMAX, X)
        DO L = 0, LMAX
          DO M = 0, L
            DP_SINGLE = DDX_ASSOC_LEGENDRE(L, M, X)
            WRITE(OUTUNIT, *) L, M, X,
     &          DP_BATCH(PLM_INDEX(L, M)), DP_SINGLE,
     &          ABS(DP_BATCH(PLM_INDEX(L, M)) - DP_SINGLE)
          END DO
        END DO
      END DO

      DEALLOCATE(P_BATCH, DP_BATCH)
      END SUBROUTINE BATCH_DALM_CONS


C  Test SSH_ALL against individual SSH calls
C  Columns: L  M  theta  phi  re(batch)  im(batch)  re(single)  im(single)
      SUBROUTINE BATCH_SSH_CONS(LMAX, NTH, OUTUNIT)
      IMPLICIT NONE
      INTEGER(KIND=i4), INTENT(IN) :: LMAX, NTH, OUTUNIT
      INTEGER(KIND=i4) :: L, M, ITH, NYLM
      REAL(KIND=dp) :: THETA, PHI, DTH
      COMPLEX(KIND=dp) :: Y_SINGLE
      COMPLEX(KIND=dp), ALLOCATABLE :: YLM(:)

      NYLM = (LMAX+1)**2
      ALLOCATE(YLM(NYLM))
      DTH = pi / (NTH + 1)
      PHI = pi / 4.d0

      WRITE(OUTUNIT, '(A)')
     &    '# L  M  theta  phi  re(batch)  im(batch)'//
     &    '  re(single)  im(single)'
      DO ITH = 1, NTH
        THETA = ITH * DTH
        CALL SSH_ALL(YLM, LMAX, THETA, PHI)
        DO L = 0, LMAX
          DO M = -L, L
            Y_SINGLE = SSH(L, M, THETA, PHI)
            WRITE(OUTUNIT, *) L, M, THETA, PHI,
     &          DREAL(YLM(YLM_INDEX(L, M))),
     &          DIMAG(YLM(YLM_INDEX(L, M))),
     &          DREAL(Y_SINGLE), DIMAG(Y_SINGLE)
          END DO
        END DO
      END DO

      DEALLOCATE(YLM)
      END SUBROUTINE BATCH_SSH_CONS


C  Verify orthonormality of SSH_ALL using midpoint quadrature in u=cos(theta)
C  Phi integral is done analytically (= 2*pi for matching m, 0 otherwise)
C  Writes upper triangle of inner product matrix: L1  M  L2  integral  expected
      SUBROUTINE SSH_ORTHO(LMAX, NQUAD, OUTUNIT)
      IMPLICIT NONE
      INTEGER(KIND=i4), INTENT(IN) :: LMAX, NQUAD, OUTUNIT
      INTEGER(KIND=i4) :: L1, L2, M, IQ, PSIZE, EXPECTED
      REAL(KIND=dp) :: U, DU, N_L1M, N_L2M, FULL_INT
      REAL(KIND=dp), ALLOCATABLE :: P(:), ORTHO(:,:,:)

      PSIZE = (LMAX+1)*(LMAX+2)/2
      ALLOCATE(P(PSIZE))
      ALLOCATE(ORTHO(0:LMAX, 0:LMAX, 0:LMAX))
      ORTHO = 0.0_dp
      DU = 2.0_dp / NQUAD

C     Accumulate int_{-1}^{1} P_l1^m(u) * P_l2^m(u) du via midpoint rule
      DO IQ = 1, NQUAD
        U = -1.0_dp + (IQ - 0.5_dp) * DU
        CALL ASSOC_LEGENDRE_ALL(P, LMAX, U)
        DO M = 0, LMAX
          DO L1 = M, LMAX
            DO L2 = L1, LMAX
              ORTHO(L1, L2, M) = ORTHO(L1, L2, M) +
     &            P(PLM_INDEX(L1, M)) * P(PLM_INDEX(L2, M))
            END DO
          END DO
        END DO
      END DO
      ORTHO = ORTHO * DU

      WRITE(OUTUNIT, '(A)')
     &    '# L1  M  L2  inner_product  expected(0_or_1)'
      DO M = 0, LMAX
        DO L1 = M, LMAX
          DO L2 = L1, LMAX
            N_L1M = DSQRT((2.d0*L1+1.d0)/(4.d0*pi)) *
     &          EXP(0.5_dp*(LOG_FACT(L1-M) - LOG_FACT(L1+M)))
            N_L2M = DSQRT((2.d0*L2+1.d0)/(4.d0*pi)) *
     &          EXP(0.5_dp*(LOG_FACT(L2-M) - LOG_FACT(L2+M)))
            FULL_INT = ORTHO(L1, L2, M) * 2.d0*pi * N_L1M * N_L2M
            IF (L1 .EQ. L2) THEN
              EXPECTED = 1
            ELSE
              EXPECTED = 0
            END IF
            WRITE(OUTUNIT, *) L1, M, L2, FULL_INT, EXPECTED
          END DO
        END DO
      END DO

      DEALLOCATE(P, ORTHO)
      END SUBROUTINE SSH_ORTHO


C  Test GRAD_SSH_ALL against GRAD_SSH for each (l,m,theta)
C  Columns: L  M  theta  phi  absdiff_r  absdiff_th  absdiff_ph
      SUBROUTINE BATCH_GRAD_SSH_CONS(LMAX, NTH, OUTUNIT)
      IMPLICIT NONE
      INTEGER(KIND=i4), INTENT(IN) :: LMAX, NTH, OUTUNIT
      INTEGER(KIND=i4) :: L, M, ITH, NYLM
      REAL(KIND=dp) :: THETA, PHI, DTH
      COMPLEX(KIND=dp), DIMENSION(3) :: VS
      COMPLEX(KIND=dp), ALLOCATABLE :: OUT(:,:)
      NYLM = (LMAX+1)**2
      ALLOCATE(OUT(3, NYLM))
      DTH = pi / (NTH + 1)
      PHI = pi / 4.d0
      WRITE(OUTUNIT,'(A)')
     &    '# L  M  theta  phi  absdiff_r  absdiff_th  absdiff_ph'
      DO ITH = 1, NTH
        THETA = ITH * DTH
        CALL GRAD_SSH_ALL(OUT, LMAX, THETA, PHI)
        DO L = 0, LMAX
          DO M = -L, L
            VS = GRAD_SSH(L, M, THETA, PHI)
            WRITE(OUTUNIT, *) L, M, THETA, PHI,
     &          ABS(OUT(1,YLM_INDEX(L,M))-VS(1)),
     &          ABS(OUT(2,YLM_INDEX(L,M))-VS(2)),
     &          ABS(OUT(3,YLM_INDEX(L,M))-VS(3))
          END DO
        END DO
      END DO
      DEALLOCATE(OUT)
      END SUBROUTINE BATCH_GRAD_SSH_CONS


C  Test L_SSH_ALL against L_SSH for each (l,m,theta)
C  Columns: L  M  theta  phi  absdiff_r  absdiff_th  absdiff_ph
      SUBROUTINE BATCH_L_SSH_CONS(LMAX, NTH, OUTUNIT)
      IMPLICIT NONE
      INTEGER(KIND=i4), INTENT(IN) :: LMAX, NTH, OUTUNIT
      INTEGER(KIND=i4) :: L, M, ITH, NYLM
      REAL(KIND=dp) :: THETA, PHI, DTH
      COMPLEX(KIND=dp), DIMENSION(3) :: VS
      COMPLEX(KIND=dp), ALLOCATABLE :: OUT(:,:)
      NYLM = (LMAX+1)**2
      ALLOCATE(OUT(3, NYLM))
      DTH = pi / (NTH + 1)
      PHI = pi / 4.d0
      WRITE(OUTUNIT,'(A)')
     &    '# L  M  theta  phi  absdiff_r  absdiff_th  absdiff_ph'
      DO ITH = 1, NTH
        THETA = ITH * DTH
        CALL L_SSH_ALL(OUT, LMAX, THETA, PHI)
        DO L = 0, LMAX
          DO M = -L, L
            VS = L_SSH(L, M, THETA, PHI)
            WRITE(OUTUNIT, *) L, M, THETA, PHI,
     &          ABS(OUT(1,YLM_INDEX(L,M))-VS(1)),
     &          ABS(OUT(2,YLM_INDEX(L,M))-VS(2)),
     &          ABS(OUT(3,YLM_INDEX(L,M))-VS(3))
          END DO
        END DO
      END DO
      DEALLOCATE(OUT)
      END SUBROUTINE BATCH_L_SSH_CONS


C  Test PVSH_RAD_ALL against PVSH_RAD for each (l,m,theta)
C  Columns: L  M  theta  phi  absdiff_r  absdiff_th  absdiff_ph
      SUBROUTINE BATCH_PVSH_RAD_CONS(LMAX, NTH, OUTUNIT)
      IMPLICIT NONE
      INTEGER(KIND=i4), INTENT(IN) :: LMAX, NTH, OUTUNIT
      INTEGER(KIND=i4) :: L, M, ITH, NYLM
      REAL(KIND=dp) :: THETA, PHI, DTH
      COMPLEX(KIND=dp), DIMENSION(3) :: VS
      COMPLEX(KIND=dp), ALLOCATABLE :: OUT(:,:)
      NYLM = (LMAX+1)**2
      ALLOCATE(OUT(3, NYLM))
      DTH = pi / (NTH + 1)
      PHI = pi / 4.d0
      WRITE(OUTUNIT,'(A)')
     &    '# L  M  theta  phi  absdiff_r  absdiff_th  absdiff_ph'
      DO ITH = 1, NTH
        THETA = ITH * DTH
        CALL PVSH_RAD_ALL(OUT, LMAX, THETA, PHI)
        DO L = 0, LMAX
          DO M = -L, L
            VS = PVSH_RAD(L, M, THETA, PHI)
            WRITE(OUTUNIT, *) L, M, THETA, PHI,
     &          ABS(OUT(1,YLM_INDEX(L,M))-VS(1)),
     &          ABS(OUT(2,YLM_INDEX(L,M))-VS(2)),
     &          ABS(OUT(3,YLM_INDEX(L,M))-VS(3))
          END DO
        END DO
      END DO
      DEALLOCATE(OUT)
      END SUBROUTINE BATCH_PVSH_RAD_CONS


C  Test PVSH_POL_ALL against PVSH_POL for each (l,m,theta)
C  Columns: L  M  theta  phi  absdiff_r  absdiff_th  absdiff_ph
      SUBROUTINE BATCH_PVSH_POL_CONS(LMAX, NTH, OUTUNIT)
      IMPLICIT NONE
      INTEGER(KIND=i4), INTENT(IN) :: LMAX, NTH, OUTUNIT
      INTEGER(KIND=i4) :: L, M, ITH, NYLM
      REAL(KIND=dp) :: THETA, PHI, DTH
      COMPLEX(KIND=dp), DIMENSION(3) :: VS
      COMPLEX(KIND=dp), ALLOCATABLE :: OUT(:,:)
      NYLM = (LMAX+1)**2
      ALLOCATE(OUT(3, NYLM))
      DTH = pi / (NTH + 1)
      PHI = pi / 4.d0
      WRITE(OUTUNIT,'(A)')
     &    '# L  M  theta  phi  absdiff_r  absdiff_th  absdiff_ph'
      DO ITH = 1, NTH
        THETA = ITH * DTH
        CALL PVSH_POL_ALL(OUT, LMAX, THETA, PHI)
        DO L = 0, LMAX
          DO M = -L, L
            VS = PVSH_POL(L, M, THETA, PHI)
            WRITE(OUTUNIT, *) L, M, THETA, PHI,
     &          ABS(OUT(1,YLM_INDEX(L,M))-VS(1)),
     &          ABS(OUT(2,YLM_INDEX(L,M))-VS(2)),
     &          ABS(OUT(3,YLM_INDEX(L,M))-VS(3))
          END DO
        END DO
      END DO
      DEALLOCATE(OUT)
      END SUBROUTINE BATCH_PVSH_POL_CONS


C  Test PVSH_TOR_ALL against PVSH_TOR for l>=1 (single-mode PVSH_TOR
C  has a 0/0 form at l=0 that the batch routine resolves to zero)
C  Columns: L  M  theta  phi  absdiff_r  absdiff_th  absdiff_ph
      SUBROUTINE BATCH_PVSH_TOR_CONS(LMAX, NTH, OUTUNIT)
      IMPLICIT NONE
      INTEGER(KIND=i4), INTENT(IN) :: LMAX, NTH, OUTUNIT
      INTEGER(KIND=i4) :: L, M, ITH, NYLM
      REAL(KIND=dp) :: THETA, PHI, DTH
      COMPLEX(KIND=dp), DIMENSION(3) :: VS
      COMPLEX(KIND=dp), ALLOCATABLE :: OUT(:,:)
      NYLM = (LMAX+1)**2
      ALLOCATE(OUT(3, NYLM))
      DTH = pi / (NTH + 1)
      PHI = pi / 4.d0
      WRITE(OUTUNIT,'(A)')
     &    '# L  M  theta  phi  absdiff_r  absdiff_th  absdiff_ph'
      DO ITH = 1, NTH
        THETA = ITH * DTH
        CALL PVSH_TOR_ALL(OUT, LMAX, THETA, PHI)
        DO L = 0, LMAX
          DO M = -L, L
            VS = PVSH_TOR(L, M, THETA, PHI)
            WRITE(OUTUNIT, *) L, M, THETA, PHI,
     &          ABS(OUT(1,YLM_INDEX(L,M))-VS(1)),
     &          ABS(OUT(2,YLM_INDEX(L,M))-VS(2)),
     &          ABS(OUT(3,YLM_INDEX(L,M))-VS(3))
          END DO
        END DO
      END DO
      DEALLOCATE(OUT)
      END SUBROUTINE BATCH_PVSH_TOR_CONS


C  Test VSH_TOR_ALL against VSH_TOR for l>=1
C  Columns: L  M  theta  phi  absdiff_r  absdiff_th  absdiff_ph
      SUBROUTINE BATCH_VSH_TOR_CONS(LMAX, NTH, OUTUNIT)
      IMPLICIT NONE
      INTEGER(KIND=i4), INTENT(IN) :: LMAX, NTH, OUTUNIT
      INTEGER(KIND=i4) :: L, M, ITH, NYLM
      REAL(KIND=dp) :: THETA, PHI, DTH
      COMPLEX(KIND=dp), DIMENSION(3) :: VS
      COMPLEX(KIND=dp), ALLOCATABLE :: OUT(:,:)
      NYLM = (LMAX+1)**2
      ALLOCATE(OUT(3, NYLM))
      DTH = pi / (NTH + 1)
      PHI = pi / 4.d0
      WRITE(OUTUNIT,'(A)')
     &    '# L  M  theta  phi  absdiff_r  absdiff_th  absdiff_ph'
      DO ITH = 1, NTH
        THETA = ITH * DTH
        CALL VSH_TOR_ALL(OUT, LMAX, THETA, PHI)
        DO L = 1, LMAX
          DO M = -L, L
            VS = VSH_TOR(L, M, THETA, PHI)
            WRITE(OUTUNIT, *) L, M, THETA, PHI,
     &          ABS(OUT(1,YLM_INDEX(L,M))-VS(1)),
     &          ABS(OUT(2,YLM_INDEX(L,M))-VS(2)),
     &          ABS(OUT(3,YLM_INDEX(L,M))-VS(3))
          END DO
        END DO
      END DO
      DEALLOCATE(OUT)
      END SUBROUTINE BATCH_VSH_TOR_CONS


C  Test VSH_POL_UP_ALL against VSH_POL_UP for each (l,m,theta)
C  Columns: L  M  theta  phi  absdiff_r  absdiff_th  absdiff_ph
      SUBROUTINE BATCH_VSH_POL_UP_CONS(LMAX, NTH, OUTUNIT)
      IMPLICIT NONE
      INTEGER(KIND=i4), INTENT(IN) :: LMAX, NTH, OUTUNIT
      INTEGER(KIND=i4) :: L, M, ITH, NYLM
      REAL(KIND=dp) :: THETA, PHI, DTH
      COMPLEX(KIND=dp), DIMENSION(3) :: VS
      COMPLEX(KIND=dp), ALLOCATABLE :: OUT(:,:)
      NYLM = (LMAX+1)**2
      ALLOCATE(OUT(3, NYLM))
      DTH = pi / (NTH + 1)
      PHI = pi / 4.d0
      WRITE(OUTUNIT,'(A)')
     &    '# L  M  theta  phi  absdiff_r  absdiff_th  absdiff_ph'
      DO ITH = 1, NTH
        THETA = ITH * DTH
        CALL VSH_POL_UP_ALL(OUT, LMAX, THETA, PHI)
        DO L = 0, LMAX
          DO M = -L, L
            VS = VSH_POL_UP(L, M, THETA, PHI)
            WRITE(OUTUNIT, *) L, M, THETA, PHI,
     &          ABS(OUT(1,YLM_INDEX(L,M))-VS(1)),
     &          ABS(OUT(2,YLM_INDEX(L,M))-VS(2)),
     &          ABS(OUT(3,YLM_INDEX(L,M))-VS(3))
          END DO
        END DO
      END DO
      DEALLOCATE(OUT)
      END SUBROUTINE BATCH_VSH_POL_UP_CONS


C  Test VSH_POL_DN_ALL against VSH_POL_DN for each (l,m,theta)
C  Columns: L  M  theta  phi  absdiff_r  absdiff_th  absdiff_ph
      SUBROUTINE BATCH_VSH_POL_DN_CONS(LMAX, NTH, OUTUNIT)
      IMPLICIT NONE
      INTEGER(KIND=i4), INTENT(IN) :: LMAX, NTH, OUTUNIT
      INTEGER(KIND=i4) :: L, M, ITH, NYLM
      REAL(KIND=dp) :: THETA, PHI, DTH
      COMPLEX(KIND=dp), DIMENSION(3) :: VS
      COMPLEX(KIND=dp), ALLOCATABLE :: OUT(:,:)
      NYLM = (LMAX+1)**2
      ALLOCATE(OUT(3, NYLM))
      DTH = pi / (NTH + 1)
      PHI = pi / 4.d0
      WRITE(OUTUNIT,'(A)')
     &    '# L  M  theta  phi  absdiff_r  absdiff_th  absdiff_ph'
      DO ITH = 1, NTH
        THETA = ITH * DTH
        CALL VSH_POL_DN_ALL(OUT, LMAX, THETA, PHI)
        DO L = 0, LMAX
          DO M = -L, L
            VS = VSH_POL_DN(L, M, THETA, PHI)
            WRITE(OUTUNIT, *) L, M, THETA, PHI,
     &          ABS(OUT(1,YLM_INDEX(L,M))-VS(1)),
     &          ABS(OUT(2,YLM_INDEX(L,M))-VS(2)),
     &          ABS(OUT(3,YLM_INDEX(L,M))-VS(3))
          END DO
        END DO
      END DO
      DEALLOCATE(OUT)
      END SUBROUTINE BATCH_VSH_POL_DN_CONS


C  Test pointwise bilinear orthogonality of PVSH_POL and PVSH_TOR:
C  DOT(PVSH_POL(l,m), PVSH_TOR(l,m)) = 0 for all l>=1, m, theta, phi.
C  Verified analytically; any nonzero value indicates numerical error.
C  Columns: L  M  theta  phi  |DOT(PVSH_POL, PVSH_TOR)|
      SUBROUTINE PVSH_POL_TOR_ORTHO(LMAX, NTH, OUTUNIT)
      IMPLICIT NONE
      INTEGER(KIND=i4), INTENT(IN) :: LMAX, NTH, OUTUNIT
      INTEGER(KIND=i4) :: L, M, ITH, NYLM, IDX
      REAL(KIND=dp) :: THETA, PHI, DTH
      COMPLEX(KIND=dp) :: DOTVAL
      COMPLEX(KIND=dp), ALLOCATABLE :: POL(:,:), TOR(:,:)
      NYLM = (LMAX+1)**2
      ALLOCATE(POL(3, NYLM), TOR(3, NYLM))
      DTH = pi / (NTH + 1)
      PHI = pi / 4.d0
      WRITE(OUTUNIT,'(A)')
     &    '# L  M  theta  phi  |DOT(PVSH_POL,PVSH_TOR)|'
      DO ITH = 1, NTH
        THETA = ITH * DTH
        CALL PVSH_POL_ALL(POL, LMAX, THETA, PHI)
        CALL PVSH_TOR_ALL(TOR, LMAX, THETA, PHI)
        DO L = 1, LMAX
          DO M = -L, L
            IDX = YLM_INDEX(L, M)
            DOTVAL = DOT(POL(:,IDX), TOR(:,IDX))
            WRITE(OUTUNIT, *) L, M, THETA, PHI, ABS(DOTVAL)
          END DO
        END DO
      END DO
      DEALLOCATE(POL, TOR)
      END SUBROUTINE PVSH_POL_TOR_ORTHO


C  Test VSH_POL_UP/DN inversion back to PVSH_POL and PVSH_RAD:
C    sqrt(l/(2l+1))*VSH_POL_UP + sqrt((l+1)/(2l+1))*VSH_POL_DN = PVSH_POL
C   -sqrt((l+1)/(2l+1))*VSH_POL_UP + sqrt(l/(2l+1))*VSH_POL_DN = PVSH_RAD
C  Columns: L  M  theta  phi  max_absdiff_pol  max_absdiff_rad
      SUBROUTINE VSH_POL_INVERSION(LMAX, NTH, OUTUNIT)
      IMPLICIT NONE
      INTEGER(KIND=i4), INTENT(IN) :: LMAX, NTH, OUTUNIT
      INTEGER(KIND=i4) :: L, M, C, ITH, NYLM, IDX
      REAL(KIND=dp) :: THETA, PHI, DTH, SC_POL, SC_RAD
      REAL(KIND=dp) :: DIFF_POL, DIFF_RAD
      COMPLEX(KIND=dp), ALLOCATABLE :: UP(:,:), DN(:,:)
      COMPLEX(KIND=dp), ALLOCATABLE :: POL(:,:), RAD(:,:)
      NYLM = (LMAX+1)**2
      ALLOCATE(UP(3,NYLM), DN(3,NYLM), POL(3,NYLM), RAD(3,NYLM))
      DTH = pi / (NTH + 1)
      PHI = pi / 4.d0
      WRITE(OUTUNIT,'(A)')
     &    '# L  M  theta  phi  max_absdiff_pol  max_absdiff_rad'
      DO ITH = 1, NTH
        THETA = ITH * DTH
        CALL VSH_POL_UP_ALL(UP,  LMAX, THETA, PHI)
        CALL VSH_POL_DN_ALL(DN,  LMAX, THETA, PHI)
        CALL PVSH_POL_ALL(POL, LMAX, THETA, PHI)
        CALL PVSH_RAD_ALL(RAD, LMAX, THETA, PHI)
        DO L = 0, LMAX
          SC_POL = DSQRT(DBLE(L)/(2*L+1.d0))
          SC_RAD = DSQRT((L+1.d0)/(2*L+1.d0))
          DO M = -L, L
            IDX = YLM_INDEX(L, M)
            DIFF_POL = 0.d0
            DIFF_RAD = 0.d0
            DO C = 1, 3
              DIFF_POL = MAX(DIFF_POL, ABS(
     &            SC_POL*UP(C,IDX) + SC_RAD*DN(C,IDX) - POL(C,IDX)))
              DIFF_RAD = MAX(DIFF_RAD, ABS(
     &           -SC_RAD*UP(C,IDX) + SC_POL*DN(C,IDX) - RAD(C,IDX)))
            END DO
            WRITE(OUTUNIT, *) L, M, THETA, PHI, DIFF_POL, DIFF_RAD
          END DO
        END DO
      END DO
      DEALLOCATE(UP, DN, POL, RAD)
      END SUBROUTINE VSH_POL_INVERSION


      END MODULE TESTS