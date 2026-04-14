c     Function definitions for scalar spherical harmonics (SSH),
c     vector spherical harmonics (VSH), and related
c     angular functions

      MODULE VSH
      USE KINDS
      USE GLOBALS
      IMPLICIT NONE

      CONTAINS
c     Compute Factorial k! for k>=0
        FUNCTION FACTORIAL(K)
        IMPLICIT NONE
        REAL(KIND=dp) :: FACTORIAL
        INTEGER(KIND=i4) :: I,K
        FACTORIAL = 1.0d0
        DO I = 1,K
            FACTORIAL = I * FACTORIAL
        END DO
        RETURN
        END FUNCTION FACTORIAL

c     Compute Log of Gamma function for (L+K)!/(L-K)! SSH pre-factor
        FUNCTION SSH_LG(L,K)
        IMPLICIT NONE
        REAL(KIND=dp) :: SSH_LG, TEMP1, TEMP2
        INTEGER(KIND=i4) :: L,K
        TEMP1 = L+K+1.d0
        TEMP2 = L-K+1.d0
        SSH_LG = EXP(LOG_GAMMA(TEMP1)-LOG_GAMMA(TEMP2))
        RETURN
        END FUNCTION SSH_LG


c     Compute Legendre polynomial P_l(x) order l at -1<=x<=+1
        FUNCTION LEGENDRE(L,X)       
        IMPLICIT NONE
        REAL(KIND=dp) :: P0,PIM1,PIM2       
        REAL(KIND=dp) :: FI,X,LEGENDRE  
        INTEGER(KIND=I4) :: I,L
C     CHECK FOR VALID VALUES OF N AND X HERE BEFORE PROCEEDING C 
C     PROCEEDING WITH CALCULATIONS       
        IF (L.EQ.0) THEN         
            LEGENDRE=1.d0       
        ELSEIF (L.EQ.1) THEN         
            LEGENDRE=X       
        ELSE 
C     USE RECURRENCE RELATIONS       
            PIM1=1.d0         
            P0=X         
            DO I=2,L    
                FI=I        
                PIM2=PIM1         
                PIM1=P0         
                P0=((2.d0*I-1.d0)*X*PIM1-(I-1.d0)*PIM2)/FI    
            ENDDO
            LEGENDRE=P0
        ENDIF       
        RETURN
        END FUNCTION LEGENDRE 

C     Compute d/dx of Legendre polynomial P_l(x) order l at -1<=x<=+1
        FUNCTION DDX_LEGENDRE(L,X)
        IMPLICIT NONE
        REAL(KIND=dp) :: X,DDX_LEGENDRE
        INTEGER(KIND=i4) :: L
        IF (L.EQ.0) THEN
            DDX_LEGENDRE=0.d0
        ELSEIF (L.EQ.1) THEN
            DDX_LEGENDRE=1.d0
        ELSEIF (L.GT.1 .AND. ABS(X).LT.1.d0) THEN
            DDX_LEGENDRE=L*(X*LEGENDRE(L,X)-
     &      LEGENDRE(L-1,X))/(X**2-1.d0)
        ELSEIF (L.GT.1 .AND. X.EQ.1.d0) THEN
            DDX_LEGENDRE=L*(L+1)/2.d0
        ELSEIF (L.GT.1 .AND. X.EQ.-1.d0) THEN
            DDX_LEGENDRE=(L*(L+1)/2.d0)*(-1)**L
        ENDIF
        RETURN
        END FUNCTION DDX_LEGENDRE

c     Compute d/dx of Associated Legendre polynomial P_l^k(x)
        FUNCTION DDX_ASSOC_LEGENDRE_DEPRC(L,K,X)
        IMPLICIT NONE
        INTEGER(KIND=i4) :: L,K,ABSK
        REAL(KIND=dp) :: X, DDX_ASSOC_LEGENDRE_DEPRC,F0,F1,F2
        ABSK = ABS(K)
        IF (K.EQ.0) THEN
            DDX_ASSOC_LEGENDRE_DEPRC = DDX_LEGENDRE(L, X)
        ELSEIF (K.NE.0 .AND. ABSK.LE.L .AND. ABS(X).LT.1.d0) THEN
            F0 = ASSOC_LEGENDRE(L-1,K,X)
            F1 = ASSOC_LEGENDRE(L,K,X)
            F2 = (L*X*F1-(L+K)*F0)/(-1.d0+X**2)

            DDX_ASSOC_LEGENDRE_DEPRC = F2  
        ELSE 
            DDX_ASSOC_LEGENDRE_DEPRC = 0.d0          
        ENDIF
        RETURN
        END FUNCTION DDX_ASSOC_LEGENDRE_DEPRC

C     Compute Associated Legendre polynomial P_l^k(x) 
        FUNCTION ASSOC_LEGENDRE_DEPRC(L,K,X)
        IMPLICIT NONE
        INTEGER(KIND=i4) :: L,K,I,ABSK
        REAL(KIND=dp) :: X,ASSOC_LEGENDRE_DEPRC,F0, F1, F2
c        REAL(KIND=dp) FACTORIAL
        ABSK = ABS(K)
        IF (K .EQ. 0) THEN
            ASSOC_LEGENDRE_DEPRC = LEGENDRE(L,X)
        ELSE IF (K .NE. 0 .AND. ABSK .LE. L) THEN
            IF (ABS(X) .LT. 1.d0) THEN
                F0 = LEGENDRE(L-1,X)
                F1 = LEGENDRE(L,X)
                F2 = (L*X*F1-L*F0)/DSQRT(1.d0-X**2)
                I = 2
                DO WHILE (I .LE. ABSK)
                    F0 = F1
                    F1 = F2
                    F2 = -2.d0*X*(I-1.d0)*F1/DSQRT(1.d0-X**2)-
     &             (L+I-1.d0)*(L-I+2.d0)*F0
                    I = I+1
                END DO
                IF (K .LT. 0) THEN
                    F2=(-1)**ABSK*FACTORIAL(L-ABSK)*F2/FACTORIAL(L+ABSK)
                ENDIF
                ASSOC_LEGENDRE_DEPRC = F2  
            ELSE
                ASSOC_LEGENDRE_DEPRC = 0.d0
            ENDIF
        ELSE 
            ASSOC_LEGENDRE_DEPRC = 0.d0          
        ENDIF
        RETURN 
        END FUNCTION ASSOC_LEGENDRE_DEPRC

        FUNCTION ASSOC_LEGENDRE(L, K, X)
            IMPLICIT NONE
            INTEGER :: L, K, ABSK, I
            REAL(KIND=dp) :: X, ASSOC_LEGENDRE
            REAL(KIND=dp) :: Pmm, Pmm1, Plm, somx2
            
            ABSK = ABS(K)
            
            IF (ABSK > L .OR. ABS(X) > 1.0_dp) THEN
                ASSOC_LEGENDRE = 0.0_dp
                RETURN
            END IF
        
            Pmm = 1.0_dp
            IF (ABSK > 0) THEN
                somx2 = sqrt((1.0_dp - X) * (1.0_dp + X))
                DO I = 1, ABSK
                    Pmm = -Pmm * (2*I - 1) * somx2
                END DO
            END IF
        
            IF (L == ABSK) THEN
                Plm = Pmm
            ELSE
                Pmm1 = X * (2 * ABSK + 1) * Pmm
                IF (L == ABSK + 1) THEN
                    Plm = Pmm1
                ELSE
                    DO I = ABSK + 2, L
                        Plm = (X*(2*I-1)*Pmm1-(I+ABSK-1)*Pmm) / (I-ABSK)
                        Pmm = Pmm1
                        Pmm1 = Plm
                    END DO
                END IF
            END IF
        
            IF (K < 0) THEN
                Plm = ((-1.0_dp)**ABSK)*
     &                (FACTORIAL(L-ABSK)/FACTORIAL(L+ABSK))*Plm
            END IF
        
            ASSOC_LEGENDRE = Plm
        END FUNCTION ASSOC_LEGENDRE

        FUNCTION DDX_ASSOC_LEGENDRE(L, K, X) RESULT(RES)
            IMPLICIT NONE
            INTEGER(KIND=i4), INTENT(IN) :: L, K
            REAL(KIND=dp), INTENT(IN)    :: X
            REAL(KIND=dp)                :: RES
            
            INTEGER(KIND=i4) :: ABSK
            REAL(KIND=dp)    :: F0, F1, DOM
            
            ABSK = ABS(K)
            DOM = 1.0_dp - X**2
        
            IF (ABS(X) >= 1.0_dp) THEN
                ! Mathematically, the derivative is often undefined or 
                ! requires a limit-based approach at the poles.
                RES = 0.0_dp ! Or handle as a specific boundary case
                RETURN
            END IF
        
            IF (ABSK > L) THEN
                RES = 0.0_dp
            ELSE IF (K == 0) THEN
                RES = DDX_LEGENDRE(L, X)
            ELSE
                ! Ensure L-1 is valid
                IF (L > 0) THEN
                    F0 = ASSOC_LEGENDRE(L-1, K, X)
                    F1 = ASSOC_LEGENDRE(L, K, X)
                    ! Using the standard recurrence: ( (l+k)P_{l-1} - l*x*P_l ) / (1 - x^2)
                    RES = ((L + K) * F0 - L * X * F1) / DOM
                ELSE
                    RES = 0.0_dp
                END IF
            END IF
        
        END FUNCTION DDX_ASSOC_LEGENDRE

C     Compute scalar spherical harmonic Y_l^k(theta,phi) at 0<=theta<=pi, 0<=phi<=2pi
C     Ylk = sqrt{(2l+1)/4pi}sqrt{(l-k)!/(l+k)!}P_l^k*exp(i k phi)
        FUNCTION SSH(L,K,THETA,PHI)
        IMPLICIT NONE
        INTEGER(KIND=i4) :: L,K
        REAL(KIND=dp) :: THETA,PHI
        COMPLEX(KIND=dp) :: SSH
        SSH = DSQRT((2.d0*L+1.d0)*FACTORIAL(L-K))/
     &      DSQRT((4.d0*pi)*FACTORIAL(L+K))*
     &      ASSOC_LEGENDRE(L,K,DCOS(THETA))*EXP(j*K*PHI)
        RETURN
        END FUNCTION SSH

c        FUNCTION SCALAR_SH(L,K,THETA,PHI)
c        IMPLICIT NONE
c        INTEGER(KIND=i4) :: L,K
c        REAL(KIND=dp) :: THETA,PHI
c        COMPLEX(KIND=dp) :: SCALAR_SH
c        SSH = DSQRT((2.d0*L+1.d0)*FACTORIAL(L-K))/
c     &      DSQRT((4.d0*pi)*FACTORIAL(L+K))*
c     &      ASSOC_LEGENDRE(L,K,DCOS(THETA))*EXP(j*K*PHI)

c        SCALAR_SH = 0.5d0*DSQRT((2.d0*L+1.d0)/pi) * 
c     &              SSH_LG(L,K) * 
c        RETURN
c        END FUNCTION SCALAR_SH    

C     Compute angular gradient of scalar spherical harmonic Ylk at 0<=theta<=pi, 0<=phi<=2pi
        FUNCTION GRAD_SSH(L,K,THETA,PHI)
        IMPLICIT NONE
        INTEGER(KIND=i4) :: L,K
        REAL(KIND=dp) :: THETA,PHI
        COMPLEX(KIND=dp) :: YLM
        COMPLEX(KIND=dp), DIMENSION(3) :: GRAD_SSH
        YLM = SSH(L,K,THETA,PHI)
        GRAD_SSH(1) = CMPLX(0.d0, 0.d0)
        GRAD_SSH(2) =-DSIN(THETA)*DSQRT((2.d0*L+1.d0)*FACTORIAL(L-K))/
     &       DSQRT((4.d0*pi)*FACTORIAL(L+K))*
     &       DDX_ASSOC_LEGENDRE(L,K,DCOS(THETA))*EXP(j*K*PHI)
        IF(ABS(K) .EQ. 1) THEN
          GRAD_SSH(3) = j*K*YLM/DSIN(THETA)
        ELSE
          GRAD_SSH(3) = j*K*YLM/DSIN(THETA)
        ENDIF
        RETURN
        END FUNCTION GRAD_SSH 


C     Compute rhat X angular gradient of scalar spherical harmonic Ylk at 0<=theta<=pi, 0<=phi<=2pi
        FUNCTION L_SSH(L,K,THETA,PHI)
        IMPLICIT NONE
        INTEGER(KIND=i4) :: L,K
        REAL(KIND=dp) :: THETA,PHI
        COMPLEX(KIND=dp) :: YLM
        COMPLEX(KIND=dp), DIMENSION(3) :: L_SSH
        YLM = SSH(L,K,THETA,PHI)
        L_SSH(1) = DCMPLX(0.d0, 0.d0)
        L_SSH(2) = -j*K*YLM/DSIN(THETA)
        L_SSH(3) = -DSIN(THETA)*DSQRT((2.d0*L+1.d0)*FACTORIAL(L-K))/
     &       DSQRT((4.d0*pi)*FACTORIAL(L+K))*
     &       DDX_ASSOC_LEGENDRE(L,K,DCOS(THETA))*EXP(j*K*PHI)
        RETURN
        END FUNCTION L_SSH


C     Compute longitudinal/radial component of polar vector spherical harmonic Y_(L,K)^(-1) = rhat.Ylm
        FUNCTION PVSH_RAD(L,K,THETA,PHI)
        IMPLICIT NONE
        INTEGER(KIND=i4) :: L,K
        REAL(KIND=dp) :: THETA,PHI
        COMPLEX(KIND=dp) :: YLM
        COMPLEX(KIND=dp), DIMENSION(3) :: PVSH_RAD
        YLM = SSH(L,K,THETA,PHI)
        PVSH_RAD(1) = YLM
        PVSH_RAD(2) = DCMPLX(0.d0, 0.d0)
        PVSH_RAD(3) = DCMPLX(0.d0, 0.d0)
        RETURN
        END FUNCTION PVSH_RAD


C     Compute transverse/toroidal component of polar vector spherical harmonic Y_(L,K)^(0) ~ L(SSH)
        FUNCTION PVSH_TOR(L,K,THETA,PHI)
        IMPLICIT NONE
        INTEGER(KIND=i4) :: L,K
        REAL(KIND=dp) :: THETA,PHI
        COMPLEX(KIND=dp), DIMENSION(3) :: LHATSSH
        COMPLEX(KIND=dp), DIMENSION(3) :: PVSH_TOR
        LHATSSH = L_SSH(L,K,THETA,PHI)
        PVSH_TOR(1) = DCMPLX(0.d0, 0.d0)
        PVSH_TOR(2) = -j*LHATSSH(2)/DSQRT(L*(L+1.d0))
        PVSH_TOR(3) = -j*LHATSSH(3)/DSQRT(L*(L+1.d0))
        RETURN
        END FUNCTION PVSH_TOR


C     Compute transverse/polar component of polar vector spherical harmonic Y_(L,K)^(+1) ~ GRADSSH
        FUNCTION PVSH_POL(L,K,THETA,PHI)
        IMPLICIT NONE
        INTEGER(KIND=i4), INTENT (IN) :: L,K
        REAL(KIND=dp), INTENT (IN) :: THETA,PHI
        COMPLEX(KIND=dp), DIMENSION(3) :: GRDSSH
        COMPLEX(KIND=dp), DIMENSION(3) :: PVSH_POL
        GRDSSH = GRAD_SSH(L,K,THETA,PHI)
        PVSH_POL(1) = DCMPLX(0.d0, 0.d0)
        IF (L .GT. 0) THEN
          PVSH_POL(2) = GRDSSH(2)/DSQRT(L*(L+1.d0))
          PVSH_POL(3) = GRDSSH(3)/DSQRT(L*(L+1.d0))
        ELSE
          PVSH_POL(2) = DCMPLX(0.d0, 0.d0)
          PVSH_POL(3) = DCMPLX(0.d0, 0.d0)
        ENDIF
        RETURN
        END FUNCTION PVSH_POL


C     Compute transverse/toroidal component of standard vector spherical harmonic Y_(L,K)^L
        FUNCTION VSH_TOR(L,K,THETA,PHI)
        IMPLICIT NONE
        INTEGER(KIND=i4) :: L,K
        REAL(KIND=dp) :: THETA,PHI
        COMPLEX(KIND=dp) :: YLM
        COMPLEX(KIND=dp), DIMENSION(3) :: VSH_TOR
        VSH_TOR = PVSH_TOR(L,K,THETA,PHI)
        RETURN
        END FUNCTION VSH_TOR


C     Compute transverse component of standard vector spherical harmonic Y_(L,K)^L-1
        FUNCTION VSH_POL_DN(L,K,THETA,PHI)
        IMPLICIT NONE
        INTEGER(KIND=i4) :: L,K
        REAL(KIND=dp) :: THETA,PHI 
        COMPLEX(KIND=dp), DIMENSION(3) :: UR, UGRAD
        COMPLEX(KIND=dp), DIMENSION(3) :: VSH_POL_DN
        UR = PVSH_RAD(L,K,THETA,PHI)
        UGRAD = PVSH_POL(L,K,THETA,PHI)
        VSH_POL_DN = DSQRT((L+1.d0)/(2.d0*L+1.d0))*UGRAD + 
     &          DSQRT((L)/(2.d0*L+1.d0))*UR
        RETURN
        END FUNCTION VSH_POL_DN


C     Compute transverse component of polar vector spherical harmonic Y_(L,K)^L+1
        FUNCTION VSH_POL_UP(L,K,THETA,PHI)
        IMPLICIT NONE
        INTEGER(KIND=i4) :: L,K
        REAL(KIND=dp) :: THETA,PHI
        COMPLEX(KIND=dp), DIMENSION(3) :: UR, UGRAD
        COMPLEX(KIND=dp), DIMENSION(3) :: VSH_POL_UP
        UR = PVSH_RAD(L,K,THETA,PHI)
        UGRAD = PVSH_POL(L,K,THETA,PHI)
        VSH_POL_UP = DSQRT((L)/(2.d0*L+1.d0))*UGRAD - 
     &          DSQRT((L+1.d0)/(2.d0*L+1.d0))*UR
        RETURN
        END FUNCTION VSH_POL_UP

C     Compute Clebsch-Gordan coefficient
C     Adapted from David Simpson (NASA GSFC)
        FUNCTION CGCOEFF(J1,M1,J2,M2,J3,M3)
        IMPLICIT NONE
        INTEGER(KIND=i4) :: J1,M1,J2,M2,J3,M3,K,KMIN,KMAX
        REAL(KIND=dp) :: SUMK,TERM,CGCOEFF

        IF ((J3 .LT. ABS(J1-J2)) .OR.  
     &      (J3 .GT. (J1+J2))    .OR.  
     &      (ABS(M1) .GT. J1)    .OR.  
     &      (ABS(M2) .GT. J2)    .OR.  
     &      (ABS(M3) .GT. J3)  .OR. 
     &      ((M1+M2) .NE. M3))  THEN
            CGCOEFF = 0.0d0
        ELSE
            CGCOEFF = DSQRT(2.0_dp*J3+1.0_dp)*
     &       EXP(-0.5_dp*LOG_FACT(J1+J2+J3+1))
            CGCOEFF = CGCOEFF * EXP(0.5_dp*(LOG_FACT(J1+J2-J3)+
     &                LOG_FACT(J2+J3-J1)+LOG_FACT(J3+J1-J2)))

            CGCOEFF = CGCOEFF * EXP(0.5_dp*(LOG_FACT(J1+M1)+
     &                LOG_FACT(J1-M1)+LOG_FACT(J2+M2)+LOG_FACT(J2-M2)+
     &                LOG_FACT(J3+M3)+LOG_FACT(J3-M3)))

            SUMK = 0.0_dp
            KMIN = MAX(0, J1-J3+M2, J2-J3-M1)
            KMAX = MIN(J1+J2-J3, J1-M1, J2+M2)
            DO K = KMIN,KMAX
                TERM = LOG_FACT(J1+J2-J3-K)+LOG_FACT(J3-J1-M2+K)+ 
     &                 LOG_FACT(J3-J2+M1+K)+LOG_FACT(J1-M1-K)+  
     &                 LOG_FACT(J2+M2-K)+LOG_FACT(K)
                IF (MOD(K,2) == 1) THEN
                    SUMK = SUMK - EXP(-TERM)
                ELSE
                    SUMK = SUMK + EXP(-TERM)
                END IF
            END DO
            CGCOEFF = CGCOEFF*SUMK
        ENDIF
        RETURN
        END FUNCTION CGCOEFF

C     Compute log(factorial()) for stability  
        FUNCTION LOG_FACT(N)
            INTEGER, INTENT(IN) :: N
            REAL(KIND=dp) :: LOG_FACT
            IF (N <= 1) THEN
                LOG_FACT = 0.0_dp
            ELSE
                ! ln(n!) = ln(Gamma(n+1))
                LOG_FACT = LOG_GAMMA(REAL(N + 1, KIND=dp))
            END IF
        END FUNCTION LOG_FACT

C     Compute Wigner 3j symbol
        FUNCTION SYMBOL3J(J1,M1,J2,M2,J3,M3)
        IMPLICIT NONE
        INTEGER(KIND=i4) :: J1,M1,J2,M2,J3,M3,K
        REAL(KIND=dp) :: TERM,CG,SYMBOL3J
        TERM = (-1)**(J3+M3+NINT(2.d0*J1))/DSQRT(2.d0*J3+1.d0)
        CG = CGCOEFF(J1,-M1,J2,-M2,J3,M3)
        SYMBOL3J = TERM*CG
        RETURN
        END FUNCTION SYMBOL3J

C     Compute Wigner 6j symbol
C-----Placeholder for future extension

C     Compute Wigner 9j symbol
C-----Placeholder for future extension

C     Compute dot product of two VSH evaluated at theta,phi coordinates
        FUNCTION DOT(VSH1, VSH2)
        IMPLICIT NONE
        COMPLEX(KIND=dp), DIMENSION(3) :: VSH1, VSH2
        COMPLEX(KIND=dp) :: DOT
        DOT = VSH1(1)*VSH2(1) + VSH1(2)*VSH2(2) + VSH1(3)*VSH2(3)
        RETURN
        END FUNCTION DOT
 

        FUNCTION GWI(J1,M1,J2,M2,L,M)
        IMPLICIT NONE
        INTEGER(KIND=i4) :: J1,M1,J2,M2,L,M
        COMPLEX(KIND=dp) :: GWI
        GWI=DSQRT((2.d0*J1+1.d0)*(2.d0*J2+1.d0)/(2.d0*L+1.d0)/4.d0/pi)*
     &    CGCOEFF(J1,0,J2,0,L,0)*CGCOEFF(J1,M1,J2,M2,L,M)
        RETURN
        END FUNCTION GWI

        
        FUNCTION GWJ(J1,M1,J2,M2,L,M)
        IMPLICIT NONE
        INTEGER(KIND=i4) :: J1,M1,J2,M2,L,M
        COMPLEX(KIND=dp) :: GWJ
        GWJ=-DSQRT((2.d0*J1+1.d0)*(2.d0*J2+1.d0)/(2.d0*L+1.d0)/4.d0/pi)*
     &    CGCOEFF(J1+1,0,J2,0,L,0)*CGCOEFF(J1,M1,J2,M2,L,M)*
     &    DCMPLX(0.d0, 1.d0)*DSQRT((J1+J2+L+2.d0)*(J2+L-J1)*
     &    (J1+J2-L+1.d0)*(J1-J2+L+1.d0))/2.d0
        RETURN
        END FUNCTION GWJ

      END MODULE VSH
