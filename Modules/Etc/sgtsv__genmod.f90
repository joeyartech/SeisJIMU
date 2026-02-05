        !COMPILER-GENERATED INTERFACE MODULE: Tue Feb  3 23:23:19 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SGTSV__genmod
          INTERFACE 
            SUBROUTINE SGTSV(N,NRHS,DL,D,DU,B,LDB,INFO)
              INTEGER(KIND=4) :: LDB
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: NRHS
              REAL(KIND=4) :: DL(*)
              REAL(KIND=4) :: D(*)
              REAL(KIND=4) :: DU(*)
              REAL(KIND=4) :: B(LDB,*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE SGTSV
          END INTERFACE 
        END MODULE SGTSV__genmod
