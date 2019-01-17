! Subroutine to read the input parameters used in the program
! driver_mpimc_1D.f90
!
! E. Mucciolo & C. Chamon (February 22, 2018)
!
subroutine readdata_v3 (Nx, Bx, Jzz, tau, nsteps, Dmax, Dcap, eps, trmode, &
     bcond, ampl_0, ampl_1, screen_count)
!
 use double_precision_mod
!
 implicit none
!
! external variables
 integer, intent(out) :: Nx, nsteps, Dmax, Dcap
 real (kind=double), intent(out) :: Bx, Jzz, tau, eps, ampl_0, ampl_1
 logical :: screen_count
 character (len=1) :: trmode, bcond
!
! read input data
 open (unit=15,file='input3',status='unknown')
 read (15,*) Nx             ! chain length
 read (15,*) nsteps         ! number of Trotter-Suzuki steps
 read (15,*) Bx             ! transverse field
 read (15,*) Jzz            ! Ising exchange interaction strength
 read (15,*) tau            ! imaginary time interval
 read (15,*) Dmax           ! maximum bond dimension allowed
 read (15,*) Dcap           ! cap for the bond dimension truncation
 read (15,*) eps            ! smallest error
 read (15,*) trmode         ! choice of truncation ('t') or not ('f')
 read (15,*) bcond          ! choice of periodic (1) or open (2) boundaries
 read (15,*) ampl_0, ampl_1 ! initial spin state amplitudes
 read (15,*) screen_count   ! print counting on screen or not
 close (unit=15)
!
 if ((trmode /= 'f').and.(trmode /= 't')) then
    write (55,'("error in reading truncation mode")')
    trmode = 'f'
 end if
!
 end subroutine readdata_v3
