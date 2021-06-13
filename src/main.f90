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

  ! v-matrix on-shell elements
  double precision , allocatable :: on_direct(:, :)
  double precision , allocatable :: on_exchange(:, :)

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
  if (alpha < 0.0d0) status = 1
  if (n_b < 0) status = 1
  if (d_r < 0.0d0) status = 1
  if (r_max < 0.0d0) status = 1
  if (d_k < 0.0d0) status = 1
  if (k_max < 0.0d0) status = 1
  if (energy < 0.0d0) status = 1
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

  allocate(on_direct(n_k, 2))
  allocate(on_exchange(n_k, 2))

  ! setup radial grid and weights
  write (*, *) "[debug] setting up radial grid and weights"

  do i_r = 1, n_r
    r_grid(i_r) = d_r * i_r
  end do

  r_weights(1:n_r:2) = 4.0d0 * (d_r / 3.0d0)
  r_weights(2:n_r:2) = 2.0d0 * (d_r / 3.0d0)

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

    ! calculate analytic on-shell v-matrix elements
    write (*, *) "[debug] calculating on-shell elements"

    call on_elements(n_r, r_grid, r_weights, n_k, k_grid, n_b, wf, energies, &
        i_f, on_direct, on_exchange, status)
    if (status /= 0) call exit(status)

    ! write v-matrix elements to file
    write (*, *) "[debug] writing output"

    call write_output(n_k, k_grid, i_i, i_f, V_direct, V_exchange, &
        on_direct, on_exchange, status)
    if (status /= 0) call exit(status)
  end do

  write (*, *) "[debug] program finish"

end program main

! read_input
subroutine read_input (l, alpha, n_b, d_r, r_max, d_k, k_max, energy, status)
  implicit none

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
subroutine write_output (n_k, k_grid, i_i, i_f, V_direct, V_exchange, &
    on_direct, on_exchange, status)
  implicit none

  integer , intent(in) :: n_k
  double precision , intent(in) :: k_grid(n_k)
  integer , intent(in) :: i_i
  integer , intent(in) :: i_f
  double precision , intent(in) :: V_direct(n_k, n_k)
  double precision , intent(in) :: V_exchange(n_k, n_k)
  double precision , intent(in) :: on_direct(n_k, 2)
  double precision , intent(in) :: on_exchange(n_k, 2)
  integer , intent(out) :: status
  ! local variables
  integer , parameter :: io_unit = 10
  character(len=200) :: io_file
  integer :: ii, jj

  ! write direct matrix elements
  write (io_file, "(a,i1,a,i1,a)") "data/output/",i_i,"s_",i_f,"s.dir.txt"

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
  write (io_file, "(a,i1,a,i1,a)") "data/output/",i_i,"s_",i_f,"s.exc.txt"

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

  ! write on-shell direct
  write (io_file, "(a,i1,a,i1,a)") "data/output/",i_i,"s_",i_f,"s.dir.on.txt"

  open(unit=io_unit, file=trim(adjustl(io_file)), action="write", iostat=status)

  if (status /= 0) then
    write (*, *) "[error] ", trim(adjustl(io_file)), " could not be opened"
    return
  end if

  write (io_unit, *) "# k_grid(ii), on_direct(ii, 1), on_direct(ii, 2)"
  write (io_unit, *) "#           , calculated      , analytic        "
  do ii = 1, n_k
    write (io_unit, *) k_grid(ii), on_direct(ii, 1), on_direct(ii, 2)
  end do

  close(io_unit)

  ! write on-shell exchange
  write (io_file, "(a,i1,a,i1,a)") "data/output/",i_i,"s_",i_f,"s.exc.on.txt"

  open(unit=io_unit, file=trim(adjustl(io_file)), action="write", iostat=status)

  if (status /= 0) then
    write (*, *) "[error] ", trim(adjustl(io_file)), " could not be opened"
    return
  end if

  write (io_unit, *) "# k_grid(ii), on_exchange(ii, 1), on_exchange(ii, 2)"
  write (io_unit, *) "#           , calculated        , analytic          "
  do ii = 1, n_k
    write (io_unit, *) k_grid(ii), on_exchange(ii, 1), on_exchange(ii, 2)
  end do

  close(io_unit)

end subroutine write_output

! on_elements
!
! Numerically calculated and analytic on-shell elements of the direct and
! exchange terms (for 1s -> ns transitions)
subroutine on_elements (n_r, r_grid, r_weights, n_k, k_grid, &
    n_b, wf, energies, i_f, on_direct, on_exchange, status)

  use vmatrix
  implicit none

  integer , intent(in) :: n_r
  double precision , intent(in) :: r_grid(n_r)
  double precision , intent(in) :: r_weights(n_r)
  integer , intent(in) :: n_k
  double precision , intent(in) :: k_grid(n_k)
  integer , intent(in) :: n_b
  double precision , intent(in) :: wf(n_r, n_b)
  double precision , intent(in) :: energies(n_b)
  integer , intent(in) :: i_f
  double precision , intent(out) :: on_direct(n_k, 2)
  double precision , intent(out) :: on_exchange(n_k, 2)
  integer , intent(out) :: status
  ! local variables
  double precision :: k_on, ki, energy
  double precision :: k, k2, k4, k6
  integer :: ii


  ! check if arguments are valid
  if ((n_r < 1) .or. (n_k < 1) .or. (n_b < 1)) then
    status = 1
    return
  end if

  ! iterate through momentum grid
  write (*, *) "[debug] energies(i) = ", energies(1)
  write (*, *) "[debug] energies(f) = ", energies(i_f)

  do ii = 1, n_k
    ki = k_grid(ii)

    ! calculate on-shell energy
    energy = energies(1) + (0.5d0 * (ki ** 2))

    write (*, *) "[debug] ki = ", ki
    write (*, *) "[debug] energy = ", energy

    ! handle impossible transition
    if (energy < energies(i_f)) then
      on_direct(ii, :) = 0.0d0
      on_exchange(ii, :) = 0.0d0

      write (*, *) "[debug] forbidden transition "
      write (*, *)
      cycle
    end if

    ! calculate on-shell momentum
    k_on = sqrt(2.0d0*(energy - energies(i_f)))

    write (*, *) "[debug] k_on = ", k_on

    ! calculate on-shell direct term
    call direct_matrix_element(n_r, r_grid, r_weights, &
        (i_f == 1), wf(:, 1), wf(:, i_f), ki, k_on, on_direct(ii, 1), status)
    if (status /= 0) return

    ! calculate on-shell exchange term
    call exchange_matrix_element(n_r, r_grid, r_weights, &
        energy, wf(:, 1), wf(:, i_f), ki, k_on, on_exchange(ii, 1), status)
    if (status /= 0) return

    ! analytic on-shell direct term
    ! (this was painful to write out)
    k = ki
    k2 = k**2
    k4 = k**4
    k6 = k**6
    select case (i_f)
    case (1)
      on_direct(ii, 2) = - 0.25d0 * ((k2 / (k2+1.0d0)) + log(k2+1.0d0))
      on_exchange(ii, 2) = - (k2 * (k2-3.0d0)) / ((k2+1.0d0)**3)

    case (2)
      on_direct(ii, 2) = &
          (16.0d0/81.0d0) * k * (4.0d0*k2 + 3.0d0) * sqrt(8.0d0*k2 - 6.0d0) &
          / ((4.0d0*k2 + 1.0d0)**2)
      on_exchange(ii, 2) = &
          (-16.0d0/9.0d0) * k * (16.0d0*k4 - 72.0d0*k2 + 13.0d0) &
          * sqrt(8.0d0*k2 - 6.0d0) &
          / ((4.0d0*k2 + 1.0d0)**4)

    case (3)
      on_direct(ii, 2) = &
          (9.0d0*sqrt(3.0d0)/128.0d0) * k * (135.0d0*k4 + 87.0d0*k2 - 4.0d0) &
          * sqrt(9.0d0*k2 - 8.0d0) / ((9.0d0*k2 + 1.0d0)**3)
      on_exchange(ii, 2) = &
          (-9.0d0*sqrt(3.0d0)/8.0d0) * k &
          * (1701.0d0*k6 - 8208.0d0*k4 + 2325.0d0*k2 - 70.0d0)&
          * sqrt(9.0d0*k2 - 8.0d0) / ((9.0d0*k2 + 1.0d0)**5)
    end select

    write (*, *) "[debug]"
  end do

end subroutine on_elements
