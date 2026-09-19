module half_chgcar
  use half_kinds, only: dp, i32, i64
  use half_math, only: inverse3
  use half_types, only: charge_grid_t, crystal_t
#ifdef HALF_HAVE_HDF5
  use half_vaspwave,only:is_hdf5_file,read_vaspwave_h5
#endif
  implicit none
  private
  public :: read_chgcar
contains
  subroutine read_chgcar(path, crystal, charge)
    character(len=*), intent(in) :: path
    type(crystal_t), intent(out) :: crystal
    type(charge_grid_t), intent(out) :: charge
    integer :: unit, ios, i, j
    integer(i64) :: ngrid
    real(dp) :: scale
    character(len=1024) :: line
    character(len=16) :: words(64)
    integer :: nwords
    logical :: cartesian

#ifdef HALF_HAVE_HDF5
    if(is_hdf5_file(path))then;call read_vaspwave_h5(path,crystal,charge);return;end if
#endif
    open(newunit=unit, file=path, status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "HALF: cannot open CHGCAR"

    read(unit, '(A)', iostat=ios) line
    call require_read(ios, "system name")
    crystal%system_name = trim(line)
    read(unit, *, iostat=ios) scale
    call require_read(ios, "lattice scale")
    if (scale <= 0.0_dp) error stop "HALF: only positive POSCAR scale is supported"
    do i = 1, 3
      read(unit, *, iostat=ios) crystal%lattice(i,:)
      call require_read(ios, "lattice vector")
    end do
    crystal%lattice = scale*crystal%lattice
    call crystal%update_geometry()

    read(unit, '(A)', iostat=ios) line
    call require_read(ios, "species or counts")
    call split_words(line, words, nwords)
    if (nwords == 0) error stop "HALF: empty species/count line"
    if (all_integer_words(words, nwords)) then
      crystal%ntypes = nwords
      allocate(crystal%species(nwords), crystal%counts(nwords))
      do i = 1, nwords
        write(crystal%species(i), '(A,I0)') "Type", i
        read(words(i), *) crystal%counts(i)
      end do
    else
      crystal%ntypes = nwords
      allocate(crystal%species(nwords), crystal%counts(nwords))
      crystal%species = words(:nwords)
      read(unit, '(A)', iostat=ios) line
      call require_read(ios, "ion counts")
      call split_words(line, words, nwords)
      if (nwords /= crystal%ntypes .or. .not. all_integer_words(words, nwords)) &
        error stop "HALF: invalid ion counts"
      do i = 1, nwords
        read(words(i), *) crystal%counts(i)
      end do
    end if
    crystal%nions = sum(crystal%counts)

    read(unit, '(A)', iostat=ios) line
    call require_read(ios, "coordinate mode")
    if (lowercase(line(1:1)) == "s") then
      read(unit, '(A)', iostat=ios) line
      call require_read(ios, "coordinate mode after selective dynamics")
    end if
    cartesian = index("ck", lowercase(line(1:1))) > 0
    allocate(crystal%positions(crystal%nions, 3))
    do i = 1, crystal%nions
      read(unit, '(A)', iostat=ios) line
      call require_read(ios, "atomic position")
      read(line, *, iostat=ios) (crystal%positions(i,j), j=1,3)
      call require_read(ios, "atomic position values")
    end do
    if (cartesian) crystal%positions = matmul(crystal%positions, inverse3(crystal%lattice))

    do
      read(unit, '(A)', iostat=ios) line
      call require_read(ios, "FFT grid")
      if (len_trim(line) > 0) exit
    end do
    read(line, *, iostat=ios) charge%shape
    call require_read(ios, "FFT grid dimensions")
    if (any(charge%shape <= 0)) error stop "HALF: invalid FFT grid dimensions"
    ngrid = int(charge%shape(1),i64)*int(charge%shape(2),i64)*int(charge%shape(3),i64)
    allocate(charge%values(ngrid))
    read(unit, *, iostat=ios) charge%values
    call require_read(ios, "smooth charge grid")
    close(unit)
  end subroutine read_chgcar

  subroutine require_read(ios, field)
    integer, intent(in) :: ios
    character(len=*), intent(in) :: field
    if (ios /= 0) then
      write(*, '(A,A)') "HALF: failed to read ", trim(field)
      error stop
    end if
  end subroutine require_read

  subroutine split_words(line, words, count)
    character(len=*), intent(in) :: line
    character(len=16), intent(out) :: words(:)
    integer, intent(out) :: count
    integer :: first, last, n
    n = len_trim(line)
    first = 1
    count = 0
    do while (first <= n)
      do while (first <= n .and. (line(first:first) == ' ' .or. line(first:first) == achar(9)))
        first = first + 1
      end do
      if (first > n) exit
      last = first
      do while (last <= n .and. line(last:last) /= ' ' .and. line(last:last) /= achar(9))
        last = last + 1
      end do
      count = count + 1
      if (count > size(words)) error stop "HALF: too many fields"
      words(count) = line(first:last-1)
      first = last + 1
    end do
  end subroutine split_words

  logical function all_integer_words(words, count)
    character(len=*), intent(in) :: words(:)
    integer, intent(in) :: count
    integer :: i, ios, value
    all_integer_words = .true.
    do i = 1, count
      read(words(i), *, iostat=ios) value
      if (ios /= 0) then
        all_integer_words = .false.
        return
      end if
    end do
  end function all_integer_words

  pure function lowercase(character) result(lower)
    character(len=1), intent(in) :: character
    character(len=1) :: lower
    integer :: code
    code = iachar(character)
    if (code >= iachar('A') .and. code <= iachar('Z')) then
      lower = achar(code + iachar('a') - iachar('A'))
    else
      lower = character
    end if
  end function lowercase
end module half_chgcar
