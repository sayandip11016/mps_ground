! Purpose: 

!    To assemble the matrices resulting from the action of a one-body
!    operator of the form exp(alpha*sigma_x).
!
!    INPUT: A, MA, NA, alpha
!    OUTPUT: Atilde
!
! E. Mucciolo & C. Chamon (Jan 29, 2013)
!
 subroutine onebody_x (MA, NA, A, Atilde, alpha)
!
 use double_precision_mod
!
 implicit none
!
! external variables and arrays
 integer, intent(in) :: MA, NA
 real (kind=double), intent(in) :: alpha
 real (kind=double), intent(in), dimension(MA,NA,0:1) :: A
 real (kind=double), intent(out), dimension(MA,NA,0:1) :: Atilde
!
! initialization
 Atilde = 0.
!
! implement the transformation
 Atilde(:,:,0) = cosh(alpha)*A(1:MA,1:NA,0) + sinh(alpha)*A(1:MA,1:NA,1)
 Atilde(:,:,1) = sinh(alpha)*A(1:MA,1:NA,0) + cosh(alpha)*A(1:MA,1:NA,1)
!
 end subroutine onebody_x
