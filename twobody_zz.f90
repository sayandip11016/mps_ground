! Purpose:
!    To assemble the matrix resulting from the action of a two-body
!    operator of the form exp(alpha*sigma_z*sigma_z) on two
!    neighboring spins described by the rectangular matrices A and
!    B. The result is a larger matrix.
!
!    INPUT: A0, A1, B0, B1, MA, R, NB, alpha
!    OUTPUT: matrix
!
! E. Mucciolo & C. Chamon (Jan 29, 2013)
!
 subroutine twobody_zz (MA, R, NB, A, B, alpha, matrix)
!
 use double_precision_mod
!
 implicit none
!
! external variables and arrays
 integer, intent(in) :: MA, R, NB
 real (kind=double), intent(in) :: alpha
 real (kind=double), intent(in), dimension(MA,R,0:1) :: A
 real (kind=double), intent(in), dimension(R,NB,0:1) :: B
 real (kind=double), intent(out), dimension(2*MA,2*NB) :: matrix
!
! internal variables
 integer :: i, j, k, i_shift, j_shift
 real (kind=double) :: aux
!
! initialization
 matrix = 0.
!
! upper left block (00)
 do j = 1, NB
    do k = 1, R
       aux = B(k,j,0)*exp(alpha)
       if (aux /= 0.) then
          do i = 1, MA
             matrix(i,j) = matrix(i,j) + A(i,k,0)*aux
          end do
       end if
    end do
 end do
!
! upper right block (01)
 do j = 1, NB
    j_shift = j + NB
    do k = 1, R
       aux = B(k,j,1)*exp(-alpha)
       if (aux /= 0.) then
          do i = 1, MA
             matrix(i,j_shift) = matrix(i,j_shift) + A(i,k,0)*aux
          end do
       end if
    end do
 end do
!
! lower left block (10)
 do j = 1, NB
    do k = 1, R
       aux = B(k,j,0)*exp(-alpha)
       if (aux /= 0.) then
          do i = 1, MA
             i_shift = i + MA
             matrix(i_shift,j) = matrix(i_shift,j) + A(i,k,1)*aux
          end do
       end if
    end do
 end do
!
! lower right block (11)
 do j = 1, NB
    j_shift = j + NB
    do k = 1, R
       aux = B(k,j,1)*exp(alpha)
       if (aux /= 0.) then
          do i = 1, MA
             i_shift = i + MA
             matrix(i_shift,j_shift) = matrix(i_shift,j_shift) + A(i,k,1)*aux
          end do
       end if
    end do
 end do
!
 end subroutine twobody_zz
