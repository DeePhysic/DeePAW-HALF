program half_potcar_inspect
  use half_types, only: potcar_t
  use half_potcar, only: read_potcar
  implicit none
  type(potcar_t),allocatable::p(:)
  character(len=1024)::path
  if(command_argument_count()/=1)then
    write(*,'(A)')'usage: half-potcar-inspect POTCAR'; stop 2
  end if
  call get_command_argument(1,path); call read_potcar(trim(path),p)
  write(*,'(A,I0)')'datasets=',size(p)
  if(size(p)>0)then
    write(*,'(A,A)')'element=',trim(p(1)%element)
    write(*,'(A,ES24.16)')'zval=',p(1)%zval
    write(*,'(A,ES24.16)')'enmax=',p(1)%enmax
    write(*,'(A,I0)')'channels=',p(1)%channels
    write(*,'(A,*(I0,1X))')'lps=',p(1)%lps
    write(*,'(A,I0)')'nmax=',p(1)%nmax
    write(*,'(A,ES24.16)')'local_checksum=',sum(p(1)%psp_local)
    write(*,'(A,ES24.16)')'projector_checksum=',sum(p(1)%pspnl)
    write(*,'(A,ES24.16)')'dion_checksum=',sum(p(1)%dion)
    write(*,'(A,ES24.16)')'qpaw_checksum=',sum(p(1)%qpaw)
    write(*,'(A,ES24.16)')'radial_checksum=',sum(p(1)%rgrid)
    write(*,'(A,ES24.16)')'waves_checksum=',sum(p(1)%wae)+sum(p(1)%wps)
  end if
end program half_potcar_inspect
