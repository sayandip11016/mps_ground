! Subroutine to compute the expectation value of a one-body operator
! with predefined matrix elements. The operator acts on spin i.
!
! INPUT: Nspin, string, bondsize, matelem, i, Z
! OUTPUT: value
!
! Notice: string defined as (Nspin,Nspin-1,...,3,2,1)
!
! E. Mucciolo & C. Chamon (Mar 04, 2013)
!
 subroutine expval_onebody (Nspin, string, bondsize, matelem, i, Z, value)
!
 use double_precision_mod
 use matrix_mod
!
 implicit none
!
! external variables
 integer, intent(in) :: Nspin, i
 type(matrix), intent(in), dimension(Nspin) :: string
 integer, intent(in), dimension(0:Nspin) :: bondsize
 real (kind=double), intent(in) :: Z
 real (kind=double), intent(in), dimension(0:1,0:1) :: matelem
 real (kind=double), intent(out) :: value
!
! internal variables
 integer :: n, k, MA, NA, MD
 real (kind=double), allocatable, dimension(:,:) :: A0, A1, B0, B1
 real (kind=double), allocatable, dimension(:,:) :: P
!
! allocate and initialize auxiliary matrix
 MD = bondsize(0)
 allocate (P(MD,MD))
 P = 1.
!
! recursive loop over the string (first i-1 spins)
 do n = 1, i-1
    NA = bondsize(n-1)
    MA = bondsize(n)
    allocate (A0(MA,NA),A1(MA,NA))
    A0(:,:) = string(n)%M(1:MA,1:NA,0)
    A1(:,:) = string(n)%M(1:MA,1:NA,1)
    allocate (B0(NA,MA),B1(NA,MA))
    B0 = matmul(P,transpose(A0))
    B1 = matmul(P,transpose(A1))
    deallocate (P)
    MD = MA
    allocate (P(MD,MD))
    P = matmul(A0,B0) + matmul(A1,B1)
    deallocate (A0,A1)
    deallocate (B0,B1)
 end do
!
! spin i
    NA = bondsize(i-1)
    MA = bondsize(i)
    allocate (A0(MA,NA),A1(MA,NA))
    A0(:,:) = string(i)%M(1:MA,1:NA,0)
    A1(:,:) = string(i)%M(1:MA,1:NA,1)
    allocate (B0(NA,MA),B1(NA,MA))
    B0 = matmul(P,transpose(A0))
    B1 = matmul(P,transpose(A1))
    deallocate (P)
    MD = MA
    allocate (P(MD,MD))
    P = matmul(A0,B0)*matelem(0,0)
    P = P + matmul(A0,B1)*matelem(0,1)
    P = P + matmul(A1,B0)*matelem(1,0)
    P = P + matmul(A1,B1)*matelem(1,1)
    deallocate (A0,A1)
    deallocate (B0,B1)
!
! recursive loop over the string (last Nspin-i spins)
 do n = i+1, Nspin
    NA = bondsize(n-1)
    MA = bondsize(n)
    allocate (A0(MA,NA),A1(MA,NA))
    A0(:,:) = string(n)%M(1:MA,1:NA,0)
    A1(:,:) = string(n)%M(1:MA,1:NA,1)
    allocate (B0(NA,MA),B1(NA,MA))
    B0 = matmul(P,transpose(A0))
    B1 = matmul(P,transpose(A1))
    deallocate (P)
    MD = MA
    allocate (P(MD,MD))
    P = matmul(A0,B0) + matmul(A1,B1)
    deallocate (A0,A1)
    deallocate (B0,B1)
 end do
!
! compute normalization factor
 value = 0.d0
 do k = 1, MD
    value = value + P(k,k)
 end do
 value = value/Z
!
! deallocate auxiliary matrix
 deallocate (P)
!
 end subroutine expval_onebody
