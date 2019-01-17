! Purpose:
!    To generate the new matrix factors resulting from a logical
!    two-bit operation over a matrix product state. It picks the
!    singular values to be used according to a predetermined sequence.
!
!    INPUT: S, U, V, MA, Rtil, Rseq, NB, seq
!    OUTPUT: A, B
!
!    Define A_{ik} = U_{ij} [S_{j}]^{1/2}, j = seq(k) 
!
!    with i=1,...,2*MA and k=1,...,Rseq. Then,
!
!    [A0]_{ik} = A_{ik}, i=1,...,MA, k=1,...,Rseq
!    [A1]_{ik} = A_{i+MA,k}, i=1,...,MA, k=1,...,Rseq.
!
!    Define B_{kj} = [S_{i}]^{1/2} V_{ij}, i = seq(k)
!
!    with k=1,...,Rseq and j=1,...,2*NB. Then,
!
!    [B0]_{kj} = B_{kj}, k=1,...,Rseq, j=1,...,NB
!    [B1]_{kj} = B_{k,j+NB}, k=1,...,Rseq, j=1,...,NB.
!
!    E. Mucciolo & C. Chamon (Feb/26/2013)
!
 subroutine factorize_seq (S, U, V, MA, Rtil, Rseq, NB, seq, Atil, Btil, ifail)
!
 use double_precision_mod
!
 implicit none
!
! external variables
 integer, intent(in) :: MA, Rtil, Rseq, NB
 integer, intent(out) :: ifail
 integer, intent(in), dimension(Rseq) :: seq
 real (kind=double), intent(in), dimension(2*MA,2*MA) :: U
 real (kind=double), intent(in), dimension(2*NB,2*NB) :: V
 real (kind=double), intent(in), dimension(Rtil) :: S
 real (kind=double), intent(out), dimension(MA,Rtil,0:1) :: Atil
 real (kind=double), intent(out), dimension(Rtil,NB,0:1) :: Btil
!
! internal variables
 integer :: i, j, k, i_shift, j_shift
 real (kind=double) :: s_sqrt
!
! check for bounds
 if (maxval(seq) > Rtil) then
    ifail = 1
 else
    ifail = 0
!
! matrices A
    Atil = 0.d0
    do k = 1, Rseq
       j = seq(k)
       s_sqrt = sqrt(S(j))
       do i = 1, MA
          i_shift = i + MA
          Atil(i,k,0) = U(i,j)*s_sqrt
          Atil(i,k,1) = U(i_shift,j)*s_sqrt
       end do
    end do
!
! matrices B
    Btil = 0.d0
    do j = 1, NB
       j_shift = j + NB
       do k = 1, Rseq
          i = seq(k)
          s_sqrt = sqrt(S(i))
          Btil(k,j,0) = s_sqrt*V(i,j)
          Btil(k,j,1) = s_sqrt*V(i,j_shift)
       end do
    end do
!
 end if
!
 end subroutine factorize_seq
