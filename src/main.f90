!>
program main

  use hydrogen
  use vmatrix

  implicit none

  ! hydrogen basis parameters
  integer :: l
  double precision :: alpha
  integer :: n_b

  ! radial grid variables
  integer :: n_r
  double precision :: d_r, r_max
  double precision , allocatable :: r_grid(:)
  double precision , allocatable :: r_weights(:)

  ! momentum grid variables
  integer :: n_k
  double precision :: d_k, k_max
  double precision , allocatable :: k_grid(:)

  ! hydrogen state variables
  double precision , allocatable :: wf(:, :)
  double precision , allocatable :: energies(:)

  ! local variables
  integer :: status

end program main
