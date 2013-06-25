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

      INTEGER, SAVE :: NEQ
      INTEGER, PARAMETER :: ITOL = 1 ! Always use scalar tolerances
      DOUBLE PRECISION, SAVE :: RTOL(1)
      DOUBLE PRECISION, SAVE :: ATOL(1)
      INTEGER, SAVE :: MF
      INTEGER :: MITER, METH, LWM, MAXORD
      INTEGER, SAVE :: ITASK
      INTEGER, SAVE :: IOPT
      INTEGER, SAVE :: LIW, LRW
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: RWORK(:)
      INTEGER, ALLOCATABLE, SAVE :: IWORK(:)
      DOUBLE PRECISION :: RPAR(1)
      INTEGER :: IPAR(1)
      INTEGER :: ERROR, I

      ! If this is the first call, process input and allocate storage
      IF ( ISTATE .EQ. 1 ) THEN
        NEQ = SIZE(Y,1)

        ! Process tolerances
        RTOL(1) = 1.D-3
        ATOL(1) = 1.D-3
        IF ( PRESENT(TOL) ) THEN
          RTOL(1) = TOL
          ATOL(1) = TOL
        END IF

        ! Process MF related data
        METH = 1 ! Use implicit Adams method by default
        MITER = 0 ! Use functional iteration by default
        IF ( PRESENT(PARAM) ) THEN
          IF ( PARAM(12) .NE. 0.D0 ) METH = INT(PARAM(12))
          IF ( PARAM(13) .NE. 0.D0 ) MITER = INT(PARAM(13))
        END IF
        IF ( METH .EQ. 1 ) MAXORD = 12
        IF ( METH .EQ. 2 ) MAXORD = 5
        MF = 10*METH + MITER

        ! Calculate work array sizes
        LIW = 30 + NEQ
        IF ( ( MITER .EQ. 0 ) .OR. ( MITER .EQ. 3 ) ) THEN
          LIW = 30
        END IF
        IF ( PRESENT(PARAM) ) THEN
          IF ( PARAM(6) .NE. 0.D0 ) MAXORD = INT(PARAM(6))
        END IF
        SELECT CASE ( MITER )
        CASE (0)
          LWM = 0
        CASE (1:2)
          LWM = 2*NEQ**2 + 2
        CASE (3)
          LWM = NEQ + 2
        CASE (4:5)
          WRITE (*,*) 'MITER = 4 OR 5 NOT IMPLEMENTED'
          STOP
        END SELECT
        LRW = 20 + NEQ*(MAXORD + 1) + 3*NEQ + LWM

        ! Allocate work arrays
        IF ( ALLOCATED(RWORK) ) THEN
          IF ( SIZE(RWORK) .GE. LRW ) THEN
            LRW = SIZE(RWORK)
          ELSE
            DEALLOCATE(RWORK)
            ALLOCATE(RWORK(LRW),STAT=ERROR)
            IF ( ERROR. NE. 0 ) THEN
              WRITE (*,*) 'Failed to allocate storage.  STAT=', ERROR
              STOP
            END IF
          END IF
        ELSE
          ALLOCATE(RWORK(LRW),STAT=ERROR)
          IF ( ERROR. NE. 0 ) THEN
            WRITE (*,*) 'Failed to allocate storage.  STAT=', ERROR
            STOP
          END IF
        END IF
        IF ( ALLOCATED(IWORK) ) THEN
          IF ( SIZE(IWORK) .GE. LIW ) THEN
            LIW = SIZE(IWORK)
          ELSE
            DEALLOCATE(IWORK)
            ALLOCATE(IWORK(LIW),STAT=ERROR)
            IF ( ERROR. NE. 0 ) THEN
              WRITE (*,*) 'Failed to allocate storage.  STAT=', ERROR
              STOP
            END IF
          END IF
        ELSE
          ALLOCATE(IWORK(LIW),STAT=ERROR)
          IF ( ERROR. NE. 0 ) THEN
            WRITE (*,*) 'Failed to allocate storage.  STAT=', ERROR
            STOP
          END IF
        END IF

        ! Figure out if we have optional input
        IOPT = 0
        IF ( PRESENT(PARAM) ) THEN
          DO I = 1, 6 ! The first 6 PARAM entries are optional input to VODE
            IF ( PARAM(I) .NE. 0.D0 ) THEN
              IOPT = 1
              RWORK(5:10) = 0.D0
              IWORK(5:10) = 0.D0
              EXIT
            END IF
          END DO
        END IF

        ! Process optional input other than MAXORD (see lrw section above)
        ITASK = 1 ! Integrate to TOUT with interpolation past TCRIT
        IF ( PRESENT(PARAM) ) THEN
          IF ( PARAM(1) .NE. 0.D0 ) RWORK(5) = PARAM(1) ! Initial step size
          IF ( PARAM(2) .NE. 0.D0 ) RWORK(7) = PARAM(2) ! Min step size
          IF ( PARAM(3) .NE. 0.D0 ) RWORK(6) = PARAM(3) ! Max step size
          IF ( PARAM(4) .NE. 0.D0 ) IWORK(6) = INT(PARAM(4)) ! Max # of steps
          IF ( PARAM(5) .NE. 0.D0 ) THEN
            ITASK = INT(PARAM(5)) ! Mode of operation
            IF ( ( ITASK .EQ. 4 ) .OR. ( ITASK .EQ. 5 ) ) THEN
              RWORK(1) = PARAM(7) ! Set TCRIT
            END IF
          END IF
        END IF
      END IF

      CALL DVODE(F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK, ISTATE,
     &   IOPT, RWORK, LRW, IWORK, LIW, DDUMMY_JAC, MF, RPAR, IPAR)
      END SUBROUTINE

      SUBROUTINE DDUMMY_JAC(NEQ, T, Y, ML, MU, PD, NROWPD, RPAR, IPAR)
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


      SUBROUTINE SDUMMY_JAC(NEQ, T, Y, ML, MU, PD, NROWPD, RPAR, IPAR)
      INTEGER, INTENT(IN) :: NEQ
      REAL, INTENT(IN) :: T
      REAL, INTENT(IN) :: Y(NEQ)
      INTEGER, INTENT(IN) :: ML
      INTEGER, INTENT(IN) :: MU
      INTEGER, INTENT(IN) :: NROWPD
      REAL, INTENT(INOUT) :: PD(NROWPD,NEQ)
      REAL, INTENT(INOUT) :: RPAR(*)
      INTEGER, INTENT(INOUT) :: IPAR(*)
      END SUBROUTINE

      END MODULE
