


      PROGRAM MAIN
      USE KINDS
      USE GLOBALS
      USE VSH
      USE TESTS
      IMPLICIT NONE

      INTEGER(KIND=i4) :: NTHETA, I
      PARAMETER (NTHETA=100)
      COMPLEX(KIND=dp), DIMENSION(NTHETA) :: TEST1A,TEST1B,TEST1C
      COMPLEX(KIND=dp), DIMENSION(NTHETA) :: TEST2A,TEST2B,TEST2C

      TEST1A = TEST1_ANALYTIC(NTHETA)
      TEST1B = TEST1_VSH(NTHETA)
      TEST1C = TEST1_GW(NTHETA)

      TEST2A = TEST2_ANALYTIC(NTHETA)
      TEST2B = TEST2_VSH(NTHETA)
      TEST2C = TEST2_GW(NTHETA)

      OPEN(UNIT=11, FILE="./validation/TEST1.DAT")
      OPEN(UNIT=12, FILE="./validation/TEST2.DAT")

      DO I=1,NTHETA
        WRITE(11,*) (I-1)*pi/(NTHETA-1),TEST1A(I),TEST1B(I),TEST1C(I)
      ENDDO

      DO I=1,NTHETA
        WRITE(12,*) (I-1)*pi/(NTHETA-1),TEST2A(I),TEST2B(I),TEST2C(I)
      ENDDO

      CLOSE(11)
      CLOSE(12)

      OPEN(UNIT=21, FILE="./out/batch_alm_cons.dat")
      CALL BATCH_ALM_CONS(8, 20, 21)
      CLOSE(21)

      OPEN(UNIT=22, FILE="./out/batch_dalm_cons.dat")
      CALL BATCH_DALM_CONS(8, 20, 22)
      CLOSE(22)

      OPEN(UNIT=23, FILE="./out/batch_ssh_cons.dat")
      CALL BATCH_SSH_CONS(5, 15, 23)
      CLOSE(23)

      OPEN(UNIT=24, FILE="./out/ssh_ortho.dat")
      CALL SSH_ORTHO(6, 10000, 24)
      CLOSE(24)

      OPEN(UNIT=31, FILE="./out/batch_grad_ssh_cons.dat")
      CALL BATCH_GRAD_SSH_CONS(5, 15, 31)
      CLOSE(31)

      OPEN(UNIT=32, FILE="./out/batch_l_ssh_cons.dat")
      CALL BATCH_L_SSH_CONS(5, 15, 32)
      CLOSE(32)

      OPEN(UNIT=33, FILE="./out/batch_pvsh_rad_cons.dat")
      CALL BATCH_PVSH_RAD_CONS(5, 15, 33)
      CLOSE(33)

      OPEN(UNIT=34, FILE="./out/batch_pvsh_pol_cons.dat")
      CALL BATCH_PVSH_POL_CONS(5, 15, 34)
      CLOSE(34)

      OPEN(UNIT=35, FILE="./out/batch_pvsh_tor_cons.dat")
      CALL BATCH_PVSH_TOR_CONS(5, 15, 35)
      CLOSE(35)

      OPEN(UNIT=36, FILE="./out/batch_vsh_tor_cons.dat")
      CALL BATCH_VSH_TOR_CONS(5, 15, 36)
      CLOSE(36)

      OPEN(UNIT=37, FILE="./out/batch_vsh_pol_up_cons.dat")
      CALL BATCH_VSH_POL_UP_CONS(5, 15, 37)
      CLOSE(37)

      OPEN(UNIT=38, FILE="./out/batch_vsh_pol_dn_cons.dat")
      CALL BATCH_VSH_POL_DN_CONS(5, 15, 38)
      CLOSE(38)

      OPEN(UNIT=39, FILE="./out/pvsh_pol_tor_ortho.dat")
      CALL PVSH_POL_TOR_ORTHO(5, 15, 39)
      CLOSE(39)

      OPEN(UNIT=40, FILE="./out/vsh_pol_inversion.dat")
      CALL VSH_POL_INVERSION(5, 15, 40)
      CLOSE(40)

      END
