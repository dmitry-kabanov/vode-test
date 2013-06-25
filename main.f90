program dummy
    use vode_int
    use fcns
    implicit none
    integer :: istate
    double precision :: x_tmp, x_end, lambda(1), tol, pm(50)

    call vode(istate, lambda_fcn, dummy_jac, lambda, x_tmp, x_end, tol, pm)
!    call vode(istate)
end program dummy
