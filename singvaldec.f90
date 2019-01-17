! PURPOSE: 
!   To find the singular value decomposition of a real rectangular
!   matrix A of dimensions MxN,
!
!   A = U * SIGMA * V
!
!   where SIGMA is a semi-positive diagonal matrix and U and V are
!   orthogonal matrices. It calls the standard LAPACK subroutine
!   DGESVD (fortran 77).
!
!   INPUT: A, N, M
!   OUTPUT: S, U, V, INFO
!
!   E. Mucciolo & C. Chamon (May/25/2012)
!
 subroutine singvaldec (A, M, N, S, U, V, INFO)
!
 use double_precision_mod
!
 implicit none
!
! external variables and arrays
 integer, intent(in) :: N, M
 integer, intent(out) :: INFO
 real (kind=double), intent(inout), dimension(M,N) :: A
 real (kind=double), intent(inout), dimension(M:M) :: U
 real (kind=double), intent(inout), dimension(N:N) :: V
! S contains the singular values of A sorted so that S(i)>=S(i+1)
 real (kind=double), intent(inout), dimension(min(N,M)) :: S  
!
! internal variables and arrays
 integer :: LWORK
 character(len=1) :: JOBU, JOBV
 real (kind=double), allocatable :: WORK(:)
!
! if LWORK=-1 on entry, WORK(1) returns optimal LWORK
! LWORK >= max(1,3*min(M,N)+max(M,N),5*min(M,N))
! For good performance, LWORK should generally be larger.
 LWORK = 2*max(1,3*min(M,N)+max(M,N),5*min(M,N))
!
! memory allocation for auxiliary array
 allocate (WORK(LWORK))
!
! option to return all columns of U in array U
 JOBU = 'A'
!
! option to return all rows of V in array V
 JOBV = 'A'
!
! call LAPACK subroutine (written in FORTRAN 77)
 call DGESVD( JOBU, JOBV, M, N, A, M, S, U, M, V, N, WORK, & 
              LWORK, INFO )
!
! Note: successful exit only if INFO=0 (test it outside subroutine)
!
 end subroutine singvaldec
!
