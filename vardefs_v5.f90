! Module to define to define variables used in driver_mpimc_1D.f90
!
! E. Mucciolo & C. Chamon (February 22, 2018)
!
 module vardefs_v5
!
   use double_precision_mod
   use matrix_mod
!
   implicit none
!
! number of spins
   integer :: Nspin
!
! lattice dimensions
   integer :: Nx
!
! number of spin pairs
   integer :: Npairs
!
! number of steps
   integer :: nsteps
!
! maximum bond dimension allowed
   integer :: Dmax, Dcap
!
! array containing bond dimensions
   integer, allocatable, dimension(:) :: bondsize
!
! array of allocatable matrices for the bit string
   type(matrix), allocatable, dimension(:) :: string
!
! array with one-body matrix elements
   real (kind=double), dimension(0:1,0:1) :: matelem1x, matelem1z
!
! array with two-body matrix elements
   real (kind=double), dimension(0:1,0:1,0:1,0:1) :: matelem2zz
!
! array with all spin pairs
   integer, allocatable, dimension(:,:) :: pairs
!
! auxiliary arrays
   real (kind=double), allocatable, dimension(:,:,:) :: A, B, Atil, Btil
!
! initial spin projection amplitudes
   real (kind=double) :: ampl_0, ampl_1
!
! moments of the magnetization
   real (kind=double) :: moment_z, moment_x
!
! ground state energy
   real (kind=double) :: energy1, energy2
!
! other real variables
   real (kind=double) :: Bx, Jzz, tau, Z, aux, alphaB, alphaJ, eps, alphaB2
!
! pi
   real (kind=double), parameter :: pi = 3.141592653590d0
!
! indices
   integer :: k, i, j, n, nn, kk
   integer :: MA, NA, NB, R, Rtil, Rnew
!
! truncation choice
   character (len=1) :: trmode
!
! boundary conditions
   character (len=1) :: bcond
!files storing final matrices
character(*),parameter :: fileplace = "/mnt/d/mat/"
character*13 :: filename1 , filename2, filename3!
! error flags
   integer :: ifail, ifail_pair, ifail_swap, ifail_zz, ifail_expval
!
! logical variables
   logical :: screen_count
!
! external subroutines passed as arguments 
   external :: swap_v3, twobody_zz
!
 end module vardefs_v5
