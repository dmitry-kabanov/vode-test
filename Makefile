all:
	gfortran -c vode.f
	gfortran -c fcns.f90
	gfortran -c dummy.f90