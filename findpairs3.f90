! Subroutine to determine all pairs of nearest neighbor spins in a
! one-dimensional lattice (chain).
!
! INPUT: Nx, Npairs
! OUTPUT: pairs, ifail
!
! The chain sites are disposed in the following order (for either
! periodic or open boundary conditions):
!
!  -1--2--3--4-
!
! The array pairs(k,j) contains the neighboring site numbers
! j1=pairs(k,1) and j2 = pairs(k,2) of the k-th pair. Each pair is
! counted only once.
!
! E. Mucciolo & C. Chamon (Feb 21, 2018)
!
 subroutine findpairs3 (Nx, bcond, Npairs, pairs, ifail)
!
 implicit none
!
! external variables
 integer, intent(in) :: Nx, Npairs
 integer, intent(out) :: ifail
 integer, intent(out), dimension(Npairs,2) :: pairs
 character (len=1) :: bcond
!
! internal variables
 integer :: i, k
!
 ifail = 0
 k = 0
!
 if (bcond == 'p') then
    k = 1
    pairs(k,1) = 1
    pairs(k,2) = Nx
 end if
!
 do i = 1, Nx-1
    k = k + 1
    pairs(k,1) = i
    pairs(k,2) = i + 1
 end do
!
 if (k /= Npairs) then
    ifail = 1
 end if
!
 end subroutine findpairs3
