!>
program main

  use hydrogen
  use vmatrix

  implicit none

  ! hydrogen basis parameters
  integer :: l
  double precision :: alpha
  integer :: n_b

  ! scattering parameters
  double precision :: energy

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

  ! v-matrix elements
  double precision , allocatable :: V_direct(:, :)
  double precision , allocatable :: V_exchange(:, :)

  ! local variables
  integer :: status
  integer :: i_r, i_k, i_f, i_i

!> program execution
  write (*, *) "[debug] program start"

  status = 0

  ! read input
  write (*, *) "[debug] reading input"

  call read_input(l, alpha, n_b, d_r, r_max, d_k, k_max, energy, status)
  if (status /= 0) call exit(status)

  ! validate input
  write (*, *) "[debug] validating input"

  if (l < 0) status = 1
  if (alpha < 0d0) status = 1
  if (n_b < 0) status = 1
  if (d_r < 0d0) status = 1
  if (r_max < 0d0) status = 1
  if (d_k < 0d0) status = 1
  if (k_max < 0d0) status = 1
  if (energy < 0d0) status = 1
  if (status /= 0) call exit(status)

  ! scale energy from eV to atomic units
  write (*, *) "[debug] scaling energy from [eV] to [Ha]"

  energy = energy / 27.2113862459d0

  ! determine number of radial and momentum grid points
  write (*, *) "[debug] calculating size of radial and momentum grids"

  n_r = ceiling(r_max / d_r)
  if (mod(n_r, 2) == 1) n_r = n_r + 1

  n_k = ceiling(k_max / d_k)

  ! allocate memory
  write (*, *) "[debug] allocating memory"

  allocate(r_grid(n_r))
  allocate(r_weights(n_r))

  allocate(k_grid(n_k))

  allocate(wf(n_r, n_b))
  allocate(energies(n_b))

  allocate(V_direct(n_k, n_k))
  allocate(V_exchange(n_k, n_k))

  ! setup radial grid and weights
  write (*, *) "[debug] setting up radial grid and weights"

  do i_r = 1, n_r
    r_grid(i_r) = d_r * i_r
  end do

  r_weights(1:n_r:2) = 4d0 * (d_r / 3d0)
  r_weights(2:n_r:2) = 2d0 * (d_r / 3d0)

  ! setup momentum grid
  write (*, *) "[debug] setting up momentum grid"

  do i_k = 1, n_k
    k_grid(i_k) = d_k * i_k
  end do

  ! calculate hydrogen wavefunctions and energies
  write (*, *) "[debug] calculating hydrogen wavefunctions and energies"

  call hydrogen_wf(l, alpha, n_b, n_r, r_grid, wf, energies, status)
  if (status /= 0) call exit(status)

  ! iterate through final hydrogen states (1s, 2s, 3s)
  write (*, *) "[debug] iterating through transitions"

  i_i = 1
  do i_f = 1, 3
    write (*, "(a,i1,a,i1,a)") " [debug] transition ",i_i,"s -> ", i_f,"s"

    ! calculate v-matrix elements numerically
    write (*, *) "[debug] calculating direct terms"

    call direct_matrix(n_r, r_grid, r_weights, &
        (i_f == i_i), wf(:, i_i), wf(:, i_f), &
        n_k, k_grid, V_direct, status)
    if (status /= 0) call exit(status)

    write (*, *) "[debug] calculating exchange terms"

    call exchange_matrix(n_r, r_grid, r_weights, &
        energy, wf(:, i_i), wf(:, i_f), &
        n_k, k_grid, V_exchange, status)
    if (status /= 0) call exit(status)

    ! calculate analytic v-matrix elements

    ! write v-matrix elements to file
    write (*, *) "[debug] writing output"

    call write_output(n_k, k_grid, i_i, i_f, V_direct, V_exchange, status)
    if (status /= 0) call exit(status)
  end do

  write (*, *) "[debug] program finish"

end program main

! read_input
subroutine read_input (l, alpha, n_b, d_r, r_max, d_k, k_max, energy, status)
  integer , intent(out) :: l
  double precision , intent(out) :: alpha
  integer , intent(out) :: n_b
  double precision , intent(out) :: d_r
  double precision , intent(out) :: r_max
  double precision , intent(out) :: d_k
  double precision , intent(out) :: k_max
  double precision , intent(out) :: energy
  integer , intent(out) :: status
  ! local variables
  integer , parameter :: io_unit = 10
  character(len=200) :: io_file

  ! set filename
  io_file = "data/input/input.txt"

  ! open input file for reading
  open(unit=io_unit, file=trim(adjustl(io_file)), action="read", iostat=status)

  ! handle invalid open
  if (status /= 0) then
    write (*, *) "[error] ", trim(adjustl(io_file)), " could not be opened"
    return
  end if

  ! read basis parameters
  read (io_unit, *) l, alpha, n_b

  ! read radial grid parameters
  read (io_unit, *) d_r, r_max

  ! read momentum grid parameters
  read (io_unit, *) d_k, k_max

  ! read scattering parameters
  read (io_unit, *) energy

  ! close input file
  close(io_unit)

end subroutine read_input

! write_output
subroutine write_output (n_k, k_grid, i_i, i_f, V_direct, V_exchange, status)
  integer , intent(in) :: n_k
  double precision , intent(in) :: k_grid(n_k)
  integer , intent(in) :: i_i
  integer , intent(in) :: i_f
  double precision , intent(in) :: V_direct(n_k, n_k)
  double precision , intent(in) :: V_exchange(n_k, n_k)
  integer , intent(out) :: status
  ! local variables
  integer , parameter :: io_unit = 10
  character(len=200) :: io_file
  integer :: ii, jj

  ! write direct matrix elements
  write (io_file, "(a,i1,a,i1,a)") "data/output/",i_i,"s_",i_f,"s.direct.txt"

  open(unit=io_unit, file=trim(adjustl(io_file)), action="write", iostat=status)

  if (status /= 0) then
    write (*, *) "[error] ", trim(adjustl(io_file)), " could not be opened"
    return
  end if

  write (io_unit, *) "# k_grid(ii), k_grid(jj), V_direct(ii, jj)"
  do ii = 1, n_k
    do jj = 1, n_k
      write (io_unit, *) k_grid(ii), k_grid(jj), V_direct(ii, jj)
    end do
    write (io_unit, *)
  end do

  close(io_unit)

  ! write exchange matrix elements
  write (io_file, "(a,i1,a,i1,a)") "data/output/",i_i,"s_",i_f,"s.exchange.txt"

  open(unit=io_unit, file=trim(adjustl(io_file)), action="write", iostat=status)

  if (status /= 0) then
    write (*, *) "[error] ", trim(adjustl(io_file)), " could not be opened"
    return
  end if

  write (io_unit, *) "# k_grid(ii), k_grid(jj), V_exchange(ii, jj)"
  do ii = 1, n_k
    do jj = 1, n_k
      write (io_unit, *) k_grid(ii), k_grid(jj), V_exchange(ii, jj)
    end do
    write (io_unit, *)
  end do

  close(io_unit)

end subroutine write_output
