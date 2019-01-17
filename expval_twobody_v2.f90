! Subroutine to compute the expectation value of a two-body operator
! with predefined matrix elements. The operator acts on arbitrary
! spins i and j (j /= i).
!
! INPUT: Nspin, string, bondsize, matelem, i, j, Z
! OUTPUT: val
!
! Notice: string defined as (Nspin,Nspin-1,...,3,2,1)
!
! matelem(d,a,b,c) = <c,b|sigma_z(j) sigma_z(i)|d,a>
!                  = <c|sigma_z(j)|d> * <b|sigma_z(i)|a>
!
! E. Mucciolo & C. Chamon (July 10, 2013)
!
 subroutine expval_twobody_v2 (Nspin, string, bondsize, matelem, i, j, Z, &
      val, ifail)
!
 use double_precision_mod
 use matrix_mod
!
 implicit none
!
! external variables
 integer, intent(in) :: Nspin, i, j
 integer, intent(out) :: ifail
 type(matrix), intent(in), dimension(Nspin) :: string
 integer, intent(in), dimension(0:Nspin) :: bondsize
 real (kind=double), intent(in) :: Z
 real (kind=double), intent(in), dimension(0:1,0:1,0:1,0:1) :: matelem
 real (kind=double), intent(out) :: val
!
! internal variables
 integer :: n, k, MA, NA, MD
 integer :: ii, jj
 real (kind=double), allocatable, dimension(:,:) :: A0, A1, B0, B1
 real (kind=double), allocatable, dimension(:,:) :: C0, C1, D0, D1
 real (kind=double), allocatable, dimension(:,:) :: P, P00, P01, P10, P11
!
! check relative position
 ifail = 0
 if (j > i) then
    jj = j
    ii = i
 else if (j < i) then
    jj = i
    ii = j
 else
    ifail = 1
 end if
!
! allocate and initialize auxiliary matrix
 MD = bondsize(0)
 allocate (P(MD,MD))
 P = 1.
!
! recursive loop over the string (first ii-1 spins)
 do n = 1, ii-1
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
! spin ii
 NA = bondsize(ii-1)
 MA = bondsize(ii)
 allocate (A0(MA,NA),A1(MA,NA))
 A0(:,:) = string(ii)%M(1:MA,1:NA,0)
 A1(:,:) = string(ii)%M(1:MA,1:NA,1)
 allocate (B0(NA,MA),B1(NA,MA))
 B0 = matmul(P,transpose(A0))
 B1 = matmul(P,transpose(A1))
 deallocate (P)
 MD = MA
 allocate (P00(MD,MD),P01(MD,MD),P10(MD,MD),P11(MD,MD))
 P00 = 0.
 P00 = P00 + matmul(A0,B0)*matelem(0,0,0,0)
 P00 = P00 + matmul(A0,B1)*matelem(0,0,1,0)
 P00 = P00 + matmul(A1,B0)*matelem(0,1,0,0)
 P00 = P00 + matmul(A1,B1)*matelem(0,1,1,0)
 P01 = 0.
 P01 = P01 + matmul(A0,B0)*matelem(0,0,0,1)
 P01 = P01 + matmul(A0,B1)*matelem(0,0,1,1)
 P01 = P01 + matmul(A1,B0)*matelem(0,1,0,1)
 P01 = P01 + matmul(A1,B1)*matelem(0,1,1,1)
 P10 = 0.
 P10 = P10 + matmul(A0,B0)*matelem(1,0,0,0)
 P10 = P10 + matmul(A0,B1)*matelem(1,0,1,0)
 P10 = P10 + matmul(A1,B0)*matelem(1,1,0,0)
 P10 = P10 + matmul(A1,B1)*matelem(1,1,1,0)
 P11 = 0.
 P11 = P11 + matmul(A0,B0)*matelem(1,0,0,1)
 P11 = P11 + matmul(A0,B1)*matelem(1,0,1,1)
 P11 = P11 + matmul(A1,B0)*matelem(1,1,0,1)
 P11 = P11 + matmul(A1,B1)*matelem(1,1,1,1)
 deallocate (A0,A1)
 deallocate (B0,B1)
!
! recursive loop over the string (spin ii+1 until spin jj-1)
 do n = ii+1, jj-1
    NA = bondsize(n-1)
    MA = bondsize(n)
    allocate (A0(MA,NA),A1(MA,NA))
    A0(:,:) = string(n)%M(1:MA,1:NA,0)
    A1(:,:) = string(n)%M(1:MA,1:NA,1)
    allocate (B0(NA,MA),B1(NA,MA))
    MD = MA
    B0 = matmul(P00,transpose(A0))
    B1 = matmul(P00,transpose(A1))
    deallocate (P00)
    allocate (P00(MD,MD))
    P00 = matmul(A0,B0) + matmul(A1,B1)
    B0 = matmul(P01,transpose(A0))
    B1 = matmul(P01,transpose(A1))
    deallocate (P01)
    allocate (P01(MD,MD))
    P01 = matmul(A0,B0) + matmul(A1,B1)
    B0 = matmul(P10,transpose(A0))
    B1 = matmul(P10,transpose(A1))
    deallocate (P10)
    allocate (P10(MD,MD))
    P10 = matmul(A0,B0) + matmul(A1,B1)
    B0 = matmul(P11,transpose(A0))
    B1 = matmul(P11,transpose(A1))
    deallocate (P11)
    allocate (P11(MD,MD))
    P11 = matmul(A0,B0) + matmul(A1,B1)
    deallocate (A0,A1)
    deallocate (B0,B1)
 end do
!
! spin jj
 NA = bondsize(jj-1)
 MA = bondsize(jj)
 allocate (A0(MA,NA),A1(MA,NA))
 A0(:,:) = string(jj)%M(1:MA,1:NA,0)
 A1(:,:) = string(jj)%M(1:MA,1:NA,1)
 allocate (B0(NA,MA),B1(NA,MA))
 MD = MA
 B0 = matmul(P00,transpose(A0)) + matmul(P01,transpose(A1))
 B1 = matmul(P10,transpose(A0)) + matmul(P11,transpose(A1))
 deallocate (P00,P01,P10,P11)
 allocate (P(MD,MD))
 P = matmul(A0,B0) + matmul(A1,B1)
 deallocate (A0,A1)
 deallocate (B0,B1)
!
! recursive loop over the string (last Nspin-jj spins)
 do n = jj+1, Nspin
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
 val = 0.d0
 do k = 1, MD
    val = val + P(k,k)
 end do
 val = val/Z
!
! deallocate auxiliary matrix
 deallocate (P)
!
 end subroutine expval_twobody_v2
