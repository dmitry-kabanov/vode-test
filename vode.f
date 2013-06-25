      MODULE VODE_INT
      IMPLICIT NONE
      PRIVATE

!     Fortran 90 Interface
      INTERFACE VODE
        MODULE PROCEDURE D_VODE
      END INTERFACE

      PUBLIC :: VODE

      CONTAINS

        SUBROUTINE D_VODE(ISTATE, F, JAC, Y, T, TOUT, TOL, PARAM)
        INTEGER, INTENT(INOUT) :: ISTATE
        DOUBLE PRECISION, INTENT(INOUT) :: Y(:)
        DOUBLE PRECISION, INTENT(INOUT) :: T
        DOUBLE PRECISION, INTENT(IN) :: TOUT
        DOUBLE PRECISION, INTENT(IN), OPTIONAL :: TOL
        DOUBLE PRECISION, INTENT(INOUT), OPTIONAL :: PARAM(50)

        INTERFACE
          SUBROUTINE F(NEQ, T, Y, YDOT, RPAR, IPAR)
          INTEGER, INTENT(IN) :: NEQ
          DOUBLE PRECISION, INTENT(IN) :: T
          DOUBLE PRECISION, INTENT(IN) :: Y(NEQ)
          DOUBLE PRECISION, INTENT(OUT) :: YDOT(NEQ)
          DOUBLE PRECISION, INTENT(INOUT) :: RPAR(*)
          INTEGER, INTENT(INOUT) :: IPAR(*)
          END SUBROUTINE

          SUBROUTINE JAC(NEQ, T, Y, ML, MU, PD, NROWPD, RPAR, IPAR)
          INTEGER, INTENT(IN) :: NEQ
          DOUBLE PRECISION, INTENT(IN) :: T
          DOUBLE PRECISION, INTENT(IN) :: Y(NEQ)
          INTEGER, INTENT(IN) :: ML
          INTEGER, INTENT(IN) :: MU
          INTEGER, INTENT(IN) :: NROWPD
          DOUBLE PRECISION, INTENT(INOUT) :: PD(NROWPD,NEQ)
          DOUBLE PRECISION, INTENT(INOUT) :: RPAR(*)
          INTEGER, INTENT(INOUT) :: IPAR(*)
          END SUBROUTINE
        END INTERFACE

      END SUBROUTINE

      END MODULE
