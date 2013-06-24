module fcns
contains
subroutine lambda_fcn(n, x, lambda, rhs, rp, ip)
    ! Computes the RHS of the ODE: dl/dx = k*(1-lambda)*exp(-e/(p*v))/u
    integer, intent(in) :: n
    double precision, intent(in) :: x, lambda(1)
    double precision, intent(out) :: rhs(1)
    double precision, intent(inout) :: rp(1)
    integer, intent(inout) :: ip(1)
end subroutine lambda_fcn

subroutine dummy_jac(NEQ, T, Y, ML, MU, PD, NROWPD, RPAR, IPAR)
    integer, intent(in) :: NEQ
    double precision, intent(in) :: T
    double precision, intent(in) :: Y(NEQ)
    integer, intent(in) :: ML
    integer, intent(in) :: MU
    integer, intent(in) :: NROWPD
    double precision, intent(inout) :: PD(NROWPD,NEQ)
    double precision, intent(inout) :: RPAR(:)
    integer, intent(inout) :: IPAR(:)
end subroutine dummy_jac
end module

program dummy
    use vode_int
    use fcns
    implicit none
    integer :: istate
    double precision :: x_tmp, x_end, lambda(1), tol, pm(50)

    call vode(istate, lambda_fcn, dummy_jac, lambda, x_tmp, x_end, tol, pm)
!    call vode(istate)
end program dummy
