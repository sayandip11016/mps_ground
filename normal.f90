! Subroutine to compute the normalization constant for the matrix
! product state of a spin string.
!
! Notice: string defined as (Nbit,Nbit-1,...,3,2,1)
!
! E. Mucciolo & C. Chamon (Feb 26, 2013)
!
 subroutine normal (Nbit, string, bondsize, Z)
!
 use double_precision_mod
 use matrix_mod
!
 implicit none
!
! external variables
 integer, intent(in) :: Nbit
 type(matrix), intent(in), dimension(Nbit) :: string
 integer, intent(in), dimension(0:Nbit) :: bondsize
 real (kind=double), intent(out) :: Z
!
! internal variables
 integer :: n, i, MA, NA
 real (kind=double), allocatable, dimension(:,:) :: A0, A1, B0, B1, P
!
! allocate and initialize auxiliary matrix
 NA = bondsize(0)
 allocate (P(NA,NA))
 P = 1.
!
! recursive loop over the string
 do n = 1, Nbit
    NA = bondsize(n-1)
    MA = bondsize(n)
    allocate (A0(MA,NA),A1(MA,NA))
    A0(:,:) = string(n)%M(1:MA,1:NA,0)
    A1(:,:) = string(n)%M(1:MA,1:NA,1)
    allocate (B0(NA,MA),B1(NA,MA))
    B0 = matmul(P,transpose(A0))
    B1 = matmul(P,transpose(A1))
    deallocate (P)
    allocate (P(MA,MA))
    P = matmul(A0,B0) + matmul(A1,B1)
    deallocate (A0,A1)
    deallocate (B0,B1)
 end do
!
! compute normalization factor
 Z = 0.d0
 do i = 1, MA 
    Z = Z + P(i,i)
 end do
!! Z = sqrt(Z)
!
! deallocate auxiliary matrix
 deallocate (P)
!
 end subroutine normal
