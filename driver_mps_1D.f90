! Driver program to compute the expectation value of one- and two-body
! observables using a matrix product state representation of a
! one-dimensional (chain) spin model with nearest-neighbor two-body
! interactions and a uniform external magnetic field.
!
! E. Mucciolo & C. Chamon (February 22, 2018)
!
 program driver_mps_1D
!
! load module with variable definitions
 use vardefs_v5
!
 implicit none
!
! read input data
 call readdata_v3 (Nx, Bx, Jzz, tau, nsteps, Dmax, Dcap, eps, trmode, &
      bcond, ampl_0, ampl_1, screen_count)
!
! number of spins
 Nspin = Nx
!
! define the boundary conditions
 if (bcond == 'p') then
!
! periodic boundary conditions
    if (Nx < 2) then
       Npairs = 0
    else if (Nx > 2) then
       Npairs = Nx
    end if
!
! open boundary conditions
 else
    Npairs = Nx - 1
!
 end if
!
 if (Npairs == 0) then
    write (20,'("No spin pairs")')
 end if
!
! the calculations needs more than one site and a proper initial state:
 if ((Npairs > 0).and.(abs(ampl_0) + abs(ampl_1) /= 0)) then
!
! allocate array with spin pairs
    allocate (pairs(Npairs,2))
!
! determine all interacting spin pairs
    call findpairs3 (Nx, bcond, Npairs, pairs, ifail_pair)
    if (ifail_pair /= 0) then
       write (20,'("error in counting pairs")')
    end if
!
! allocate bond dimension size array
    allocate (bondsize(0:Nspin))
!
! initialize bond dimension array
    bondsize = 1
!
! allocate bit string
    allocate (string(Nspin))
!
! allocate bit matrices in the string and initialize (product state)
    do k = 1, Nspin
       allocate (string(k)%M(bondsize(k),bondsize(k-1),0:1))
       string(k)%M = 0.d0
       string(k)%M(1,1,0) = ampl_0
       string(k)%M(1,1,1) = ampl_1
    end do
!
! define one-body matrix element of the magnetization (sigma_z):
! matelem1z(a,b) = <a|sigma_z|b>
    matelem1z = 0
    matelem1z(0,0) = 1.d0
    matelem1z(1,1) = -1.d0
!
! define one-body matrix element of the magnetization (sigma_x)
! matelem1x(a,b) = <a|sigma_x|b>
    matelem1x = 0.d0
    matelem1x(0,1) = 1.d0
    matelem1x(1,0) = 1.d0
!
! find normalization constant: Z = <psi|psi>
    call normal (Nspin, string, bondsize, Z)
!
! evaluate the expectation value of the magnetization along z
    moment_z = 0.d0
    LATTICE1: do k = 1, Nspin
       call expval_onebody (Nspin, string, bondsize, matelem1z, k, Z, aux)
       moment_z = moment_z + aux
    end do LATTICE1
!
! evaluate the expectation value of the magnetization along x
    moment_x = 0.d0
    LATTICE2: do k = 1, Nspin
       call expval_onebody (Nspin, string, bondsize, matelem1x, k, Z, aux)
       moment_x = moment_x + aux
    end do LATTICE2
!
! define two-body matrix element of the interaction (sigma_z*sigma_z)
! matelem2zz(d,a,b,c) = <c,b|sigma_z(2) sigma_z(1)|d,a>
    matelem2zz = 0.d0
    matelem2zz(0,0,0,0) = 1.d0
    matelem2zz(0,1,1,0) = -1.d0
    matelem2zz(1,0,0,1) = -1.d0
    matelem2zz(1,1,1,1) = 1.d0
!
! set one-body operation angle
    alphaB = tau*Bx/nsteps
    alphaB2 = alphaB/2.d0
!
! set two-body operation angle
    alphaJ = tau*Jzz/nsteps
!
! open output file with expectation values and write first line
    open (unit=46,file='out',status='unknown')
    energy1 = 0.d0
    energy2 = 0.d0
    write (46,'(5(1x,e12.5))') moment_z/Nspin, moment_x/Nspin, &
            energy1, energy2, energy1 + energy2
!
!****************************************************************************
!** Block for evolving matrices
!
!============================================================================
! Initial stage: one-body operations
    do k = 1, Nspin
       MA = bondsize(k)
       NA = bondsize(k-1)
       allocate (A(MA,NA,0:1),Atil(MA,NA,0:1))
       A(1:MA,1:NA,0:1) = string(k)%M(1:MA,1:NA,0:1)
       call onebody_x (MA, NA, A, Atil, alphaB2)
       string(k)%M(1:MA,1:NA,0:1) = Atil(1:MA,1:NA,0:1)
       deallocate (A,Atil)
    end do
!
! find normalization constant
    call normal (Nspin, string, bondsize, Z)
!
! rescale the matrices in the string
    do k = 1, Nspin
       string(k)%M = string(k)%M/Z**(.5d0/Nspin)
    end do
!
!============================================================================
! Initial stage: two-body operations
!
! >> loop over pairs
    PAIRS1: do k = 1, Npairs
!
       i = pairs(k,1)
       j = pairs(k,2)       
!
!----------------------------------------------------------------------------
! Apply the two-body operator
!----------------------------------------------------------------------------
!
! get bond dimensions
       MA = bondsize(j)
       R = bondsize(j-1)
       NB = bondsize(j-2)
       Rtil = 2*min(MA,NB)
!
! allocate auxiliary matrices - input
       allocate (A(MA,R,0:1))
       allocate (B(R,NB,0:1))
!
! transfer neighboring bit matrices to auxiliary matrices
       B(1:R,1:NB,0:1) = string(j-1)%M(1:R,1:NB,0:1)
       A(1:MA,1:R,0:1) = string(j)%M(1:MA,1:R,0:1)
!
! allocate auxilirary matrices - output
       allocate (Btil(Rtil,NB,0:1))
       allocate (Atil(MA,Rtil,0:1))
!
! implement the two-body operation
       if (alphaJ /= 0.d0) then
          call twobody_op4 (A, B, MA, R, NB, Atil, Btil, Dmax, Dcap, & 
               Rtil, Rnew, twobody_zz, alphaJ, eps, trmode, ifail_zz)
!
       else
          Rnew = R
          Atil = 0.d0
          Atil(1:MA,1:Rnew,0:1) = A(1:MA,1:Rnew,0:1)
          Btil = 0.d0
          Btil(1:Rnew,1:NB,0:1) = B(1:Rnew,1:NB,0:1)
       end if
!
! check maximum bond dimension allowed and check for errors in operation
       if (Rnew > Dmax) then
          ifail = 4
       else if (ifail_zz /= 0) then
          ifail = 5
       else
!
! change bond dimension in global array
          bondsize(j-1) = Rnew
!
! deallocate/allocate bit matrices
          if (Rnew /= R) then
             deallocate (string(j-1)%M)
             deallocate (string(j)%M)
             allocate (string(j-1)%M(Rnew,NB,0:1))
             allocate (string(j)%M(MA,Rnew,0:1))
          end if
!
! transfer the resulting matrices back to the bit string
          string(j-1)%M(1:Rnew,1:NB,0:1) = Btil(1:Rnew,1:NB,0:1)
          string(j)%M(1:MA,1:Rnew,0:1) = Atil(1:MA,1:Rnew,0:1)
!
       end if
!
! deallocate auxiliary arrays 
       deallocate (A,B)
       deallocate (Atil,Btil)
!
! find normalization constant
       call normal (Nspin, string, bondsize, Z)
!
! rescale the matrices in the string
       do kk = 1, Nspin
          string(kk)%M = string(kk)%M/Z**(.5d0/Nspin)
       end do
!
!----------------------------------------------------------------------------
!
    end do PAIRS1
!
!============================================================================
! >> Loop over intermediate steps
    STEPS: do nn = 1, nsteps-1

!!!       print *,'nn=',nn
       
!
!============================================================================
! Intermediate stage: one-body operations
       do k = 1, Nspin
          MA = bondsize(k)
          NA = bondsize(k-1)
          allocate (A(MA,NA,0:1),Atil(MA,NA,0:1))
          A(1:MA,1:NA,0:1) = string(k)%M(1:MA,1:NA,0:1)
          call onebody_x (MA, NA, A, Atil, alphaB)
          string(k)%M(1:MA,1:NA,0:1) = Atil(1:MA,1:NA,0:1)
          deallocate (A,Atil)
       end do
!
! find normalization constant
       call normal (Nspin, string, bondsize, Z)
!
! rescale the matrices in the string
       do k = 1, Nspin
          string(k)%M = string(k)%M/Z**(.5d0/Nspin)
       end do
!
!============================================================================
! Intermediate stage: two-body operations
!
! >> loop over pairs
       PAIRS2: do k = 1, Npairs
          i = pairs(k,1)
          j = pairs(k,2)
!          
!----------------------------------------------------------------------------
! Apply the two-body operator
!----------------------------------------------------------------------------
! get bond dimensions
          MA = bondsize(j)
          R = bondsize(j-1)
          NB = bondsize(j-2)
          Rtil = 2*min(MA,NB)
!
! allocate auxiliary matrices - input
          allocate (A(MA,R,0:1))
          allocate (B(R,NB,0:1))
!
! transfer neighboring bit matrices to auxiliary matrices
          B(1:R,1:NB,0:1) = string(j-1)%M(1:R,1:NB,0:1)
          A(1:MA,1:R,0:1) = string(j)%M(1:MA,1:R,0:1)
!
! allocate auxilirary matrices - output
          allocate (Btil(Rtil,NB,0:1))
          allocate (Atil(MA,Rtil,0:1))
!
! implement the two-body operation
          if (alphaJ /= 0.d0) then
             call twobody_op4 (A, B, MA, R, NB, Atil, Btil, Dmax, Dcap, &
                  Rtil, Rnew, twobody_zz, alphaJ, eps, trmode, ifail_zz)
!
          else
             Rnew = R
             Atil = 0.d0
             Atil(1:MA,1:Rnew,0:1) = A(1:MA,1:Rnew,0:1)
             Btil = 0.d0
             Btil(1:Rnew,1:NB,0:1) = B(1:Rnew,1:NB,0:1)
          end if
!
! check maximum bond dimension allowed and check for errors in operation
          if (Rnew > Dmax) then
             ifail = 4
          else if (ifail_zz /= 0) then
             ifail = 5
          else
!
! change bond dimension in global array
             bondsize(j-1) = Rnew
!
! deallocate/allocate bit matrices
             if (Rnew /= R) then
                deallocate (string(j-1)%M)
                deallocate (string(j)%M)
                allocate (string(j-1)%M(Rnew,NB,0:1))
                allocate (string(j)%M(MA,Rnew,0:1))
             end if
!
! transfer the resulting matrices back to the bit string
             string(j-1)%M(1:Rnew,1:NB,0:1) = Btil(1:Rnew,1:NB,0:1)
             string(j)%M(1:MA,1:Rnew,0:1) = Atil(1:MA,1:Rnew,0:1)
!
          end if
!
! deallocate auxiliary arrays
          deallocate (A,B)
          deallocate (Atil,Btil)
!
! find normalization constant
          call normal (Nspin, string, bondsize, Z)
!
! rescale the matrices in the string
          do kk = 1, Nspin
             string(kk)%M = string(kk)%M/Z**(.5d0/Nspin)
          end do
!
!----------------------------------------------------------------------------
!
       end do PAIRS2
!
!============================================================================
!
! >> close loop over intermediate steps
    end do STEPS
!
!============================================================================
! Final stage: one-body operations
!
    do k = 1, Nspin
       MA = bondsize(k)
       NA = bondsize(k-1)
       allocate (A(MA,NA,0:1),Atil(MA,NA,0:1))
       A(1:MA,1:NA,0:1) = string(k)%M(1:MA,1:NA,0:1)
       call onebody_x (MA, NA, A, Atil, alphaB2)
       string(k)%M(1:MA,1:NA,0:1) = Atil(1:MA,1:NA,0:1)
       deallocate (A,Atil)
    end do
!
! find normalization constant
    call normal (Nspin, string, bondsize, Z)
!
! rescale the matrices in the string
    do k = 1, Nspin
       write(filename1,'("matrixu",I2,".txt")')k
       write(filename2,'("matrixd",I2,".txt")')k
       open(unit = k+49, file = fileplace//filename1,status="replace")
       open(unit = k+109, file =fileplace//filename2,status="replace")
       string(k)%M = string(k)%M/Z**(.5d0/Nspin)
       MA = bondsize(k)
       NA = bondsize(k-1)
       do i = 1,MA
       do j = 1,NA
         write(k+49,*) CMPLX(string(k)%M(i,j,0),0.d0)
         write(k+109,*) CMPLX(string(k)%M(i,j,1),0.d0)
       enddo
       enddo 
     close(k+49)
     close(k+109)
    end do

open (unit=1009 ,file = fileplace//"bonddim.txt",status ="replace")
    do k=0,Nspin
     write(1009,*) bondsize(k)
    enddo
close(1009)      
   
!
!============================================================================
!
!****************************************************************************
!
! find normalization constant
    call normal (Nspin, string, bondsize, Z)
!
! evaluate the expectation value of the magnetization along z
    moment_z = 0.d0
    LATTICE3: do k = 1, Nspin
       call expval_onebody (Nspin, string, bondsize, matelem1z, k, Z, aux)
       moment_z = moment_z + aux
    end do LATTICE3
!
! evaluate the expectation value of the magnetization along x
    moment_x = 0.d0
    LATTICE4: do k = 1, Nspin
       call expval_onebody (Nspin, string, bondsize, matelem1x, k, Z, aux)
       moment_x = moment_x + aux
    end do LATTICE4
!
! evaluate the expectation value of the total energy over horizontal segments
! (open boundaries)
!
    energy1 = 0.d0
    LATTICE5: do k = 1, Nspin
       call expval_onebody (Nspin, string, bondsize, matelem1x, k, Z, aux)
       energy1 = energy1 - aux*Bx
    end do LATTICE5
!
    energy2 = 0.d0
    PAIRS3: do k = 1, Npairs
       i = pairs(k,1)
       j = pairs(k,2)
       call expval_twobody_v2 (Nspin, string, bondsize, matelem2zz, &
            i, j, Z, aux, ifail_expval)
       if (ifail_expval /= 0) then
          ifail = 6
       end if
       energy2 = energy2 - aux*Jzz
    end  do PAIRS3
!
! save results
    write (46,'(6(1x,e22.15))') moment_z/Nspin, moment_x/Nspin, &
              energy1/Nspin, energy2/Nspin, (energy1 + energy2)/Nspin
!
! screen counting?
    if (screen_count) then
       write (*,'(" Z = ",e12.5," Dmax = ",i3," <E_0> = ",f15.12)') &
            Z, maxval(bondsize), (energy1+energy2)/Nspin
    end if
!
! close output file with expectation values
    close (unit=46)
!
! deallocate bit string matrices
    do k = 1, Nspin
       deallocate (string(k)%M)
    end do
!
! print final results on the screen, if requested
    if (screen_count) then
       write (*,'("Nx = ",i3)') Nx
       if (bcond == 'o') then
          write (*,'("boundary conditions: open")')
       else if (bcond == 'p') then
          write (*,'("boundary conditions: periodic")')
       end if
       if (trmode == 'f') then
          write (*,'("full matrices")')
       else if (trmode == 't') then
          write (*,'("truncation at D = ",i3)') Dcap
       else if (trmode == 's') then
          write (*,'("sampling up to D = ",i3)') Dcap
       end if
       write (*,'("Jzz = ",f11.8)') Jzz
       write (*,'("B_x = ",f11.8)') Bx
       write (*,'("m_z = ",f11.8)') moment_z/Nspin
       write (*,'("m_x = ",f11.8)') moment_x/Nspin
       write (*,'("E_0 = ",f12.8, &
            " = ",f15.11," (2-body) + ",f15.11," (1-body)")') &
            (energy1+energy2)/Nspin, energy2/Nspin, energy1/Nspin
    end if
!
 end if
!
! deallocate arrays
 deallocate (string)
 deallocate (bondsize)
 deallocate (pairs)
!
 end program driver_mps_1D
!
 include "readdata_v3.f90"
 include "findpairs3.f90"
 include "singvaldec.f90"
 include "factorize_seq.f90"
 include "twobody_op4.f90"
 include "swap_v3.f90"
 include "twobody_zz.f90"
 include "normal.f90"
 include "expval_onebody.f90"
 include "onebody_x.f90"
 include "expval_twobody_v2.f90"
