#include "symbol.inc"
      MODULE profiling
      USE prec
      USE c2f_interface, ONLY: VTIME2
#ifdef NVTX
      USE c2f_interface, ONLY: NVTX_PUSH,NVTX_POP,C_NULL_CHAR
#endif
!**********************************************************************
!
!>  This module implements a small profiling utility for VASP.
!>
!>  USAGE:
!>  ======
!>
!>  To include a routine to the profiling insert at the beginning
!>
!>  #ifdef PROFILING
!>     CALL START_PROFILING(#NAME#)
!>  #endif
!>
!>  and before each RETURN and at the end of the routine
!>
!>  #ifdef PROFILING
!>     CALL STOP_PROFILING(#NAME)
!>  #endif
!>
!>  Mind that the profiling routines will mess all things up, if a
!>  call to START_PROFILING exist without a corresponding call to
!>  STOP_PROFILING and vice versa.
!>
!>  In addition you can include the calls also inside routines
!>  like it was done to profile RECURSIVE SUBROUTINES in broyden.F.
!
!**********************************************************************

      ! maximum number of timers
      INTEGER, PARAMETER :: maxtimer=100000
      INTEGER, PARAMETER :: maxlevel=50
      ! actual used timers
      INTEGER, SAVE :: used_ptimers
      ! tag for the name of the routine
      CHARACTER (LEN=40), SAVE :: name_tag(maxtimer)
      ! temporary array for storing timings, only the size of the
      ! maximum nesting level is required since two routines
      ! can not run in parallel in a nesting level
      REAL(q), SAVE :: cpu(maxlevel)
      REAL(q), SAVE :: cpu_profiling
      REAL(q), SAVE :: cpu_profiling_stop
      ! summed up timings
      REAL(q), SAVE :: sum_cpu(maxtimer)
      REAL(q), SAVE :: sum_cpu_profiling(4)
      REAL(q), SAVE :: sum_cpu_profiling_stop
      ! number of calls, nesting level and the end in the arrays of the routine
      INTEGER, SAVE :: calls(maxtimer),level(maxtimer),endsubrout(maxtimer)
      INTEGER, SAVE :: startindex(maxlevel)
      ! some indices
      INTEGER, SAVE :: actlevel,actindex,endindex
      !
!      LOGICAL,PRIVATE  :: writeout
!      INTEGER,PRIVATE  :: writeoutstart
!      INTEGER, PARAMETER, PRIVATE :: IOP=1000

      INTEGER :: PLEVEL
      REAL(q) :: PTHRESHOLD
      LOGICAL :: PFLAT
      CHARACTER (LEN=40) :: PROUTINE

      CHARACTER (LEN=40), SAVE :: opened_tag(maxlevel)
      CONTAINS

!***********************************************************************
!
!> INIT_PROFILING
!> initialises all arrays and starts the profiling for the total timing.
!
!***********************************************************************
      SUBROUTINE INIT_PROFILING
      IMPLICIT NONE

      REAL(q) :: TC
      INTEGER :: I
      !
      ! start timing of profiling
      !
      CALL VTIME2(TC)
      cpu_profiling=TC
      sum_cpu_profiling(:)=0._q
      sum_cpu_profiling_stop=0._q
      !
      ! set all arrays to zero
      !
      cpu(:)=0._q
      sum_cpu(:)=0._q
      calls(:)=0
      name_tag(:)=''
      level(:)=0
      endsubrout(:)=0
      startindex(:)=0
      !
      ! starts profiling for the whole program (total time), called in main.F
      !
      CALL VTIME2(TC)

      cpu(1)=TC
      name_tag(1)='total_time'
      level(1)=1
      calls(1)=1
      startindex(1)=1
      !
      ! update opened tag 
      !
      opened_tag(1) = 'total_time'
      !
      ! set values for indices
      !
      used_ptimers=1
      actlevel=2
      actindex=1
      endindex=1
      !
      ! stop timing of profiling
      !
      CALL VTIME2(TC)
      sum_cpu_profiling=sum_cpu_profiling+TC-cpu_profiling

#ifdef NVTX
      CALL NVTX_PUSH('total_time')
#endif

      RETURN
      END SUBROUTINE


!***********************************************************************
!
! PRINT_FULL_PROFILING
!> Prints the full profile, setting a threshold for the time is possible.
!
!***********************************************************************
      SUBROUTINE PRINT_FULL_PROFILE(IU,PTHRESHOLD,PFLAT)

      IMPLICIT NONE

      INTEGER :: IU, I, J, K, skipped
      INTEGER :: outfile
      ! threshold for the time
      REAL(q), OPTIONAL :: PTHRESHOLD
      LOGICAL, OPTIONAL :: PFLAT
      !
      ! used for flat profile
      !
      CHARACTER (LEN=40), DIMENSION(used_ptimers) :: flat_names
      REAL(q), DIMENSION(used_ptimers) :: flat_times_cpu
      INTEGER, DIMENSION(used_ptimers) :: flat_calls
      CHARACTER (LEN=40) :: flat_name
      REAL(q) :: flat_time_cpu, sum_time
      INTEGER :: used_flat_timers,flat_call
      LOGICAL :: FOUND

      outfile=iu
      skipped=0

      IF (.NOT.PFLAT) THEN
        WRITE(outfile,*)''
        WRITE(outfile,'(A,I7)')' PROFILE, used timers: ', used_ptimers
        WRITE(6,'(A,I7)')' PROFILE, used timers: ', used_ptimers
        WRITE(outfile,*)'============================='
        WRITE(outfile,*)''
        WRITE(outfile,'(A)')' index                 routine                    time                            calls level'
        WRITE(outfile,'(A)')'----------------------------------------------------------------------------------------------'
        DO I=1,used_ptimers
          IF (PRESENT(PTHRESHOLD)) THEN
            IF (sum_cpu(i)<PTHRESHOLD) THEN
              skipped=skipped+1
              CYCLE
            ENDIF
          ENDIF
          SELECT CASE (level(i))
            CASE (0)
              WRITE(outfile,1000) i,name_tag(i),sum_cpu(i),calls(i),level(i)
            CASE (1)
              WRITE(outfile,1001) i,name_tag(i),sum_cpu(i),calls(i),level(i)
            CASE (2)
              WRITE(outfile,1002) i,name_tag(i),sum_cpu(i),calls(i),level(i)
            CASE (3)
              WRITE(outfile,1003) i,name_tag(i),sum_cpu(i),calls(i),level(i)
            CASE (4)
              WRITE(outfile,1004) i,name_tag(i),sum_cpu(i),calls(i),level(i)
            CASE (5)
              WRITE(outfile,1005) i,name_tag(i),sum_cpu(i),calls(i),level(i)
            CASE (6)
              WRITE(outfile,1006) i,name_tag(i),sum_cpu(i),calls(i),level(i)
            CASE (7)
              WRITE(outfile,1007) i,name_tag(i),sum_cpu(i),calls(i),level(i)
            CASE (8)
              WRITE(outfile,1008) i,name_tag(i),sum_cpu(i),calls(i),level(i)
            CASE (9)
              WRITE(outfile,1009) i,name_tag(i),sum_cpu(i),calls(i),level(i)
            CASE (10)
              WRITE(outfile,1010) i,name_tag(i),sum_cpu(i),calls(i),level(i)
            CASE (11)
              WRITE(outfile,1011) i,name_tag(i),sum_cpu(i),calls(i),level(i)
            CASE (12)
              WRITE(outfile,1012) i,name_tag(i),sum_cpu(i),calls(i),level(i)
            CASE DEFAULT
              WRITE(outfile,1012) i,name_tag(i),sum_cpu(i),calls(i),level(i)
          END SELECT
        ENDDO
        IF (PRESENT(PTHRESHOLD)) THEN
          WRITE(outfile,*)''
          WRITE(outfile,'(A,I6,A,F10.6,A)') 'skipped ',skipped,' routines with timings under ',pthreshold,' sec'
        ENDIF
      ENDIF
      !
      ! now calculate the absolute time used in all subroutines. Therefore we have to
      ! subtract the time of all called subroutines higher by one in level within the subroutine under
      ! consideration to get the absolute time of the subroutine. Additionally we
      ! have to sum all times up.
      !
      used_flat_timers=0
      flat_names(:)=''
      flat_times_cpu(:)=0._q
      flat_calls(:)=0
      skipped=0
      DO I=1,used_ptimers
        IF (PRESENT(PTHRESHOLD)) THEN
          IF (sum_cpu(i)<PTHRESHOLD) THEN
            CYCLE
          ENDIF
        ENDIF
        flat_time_cpu=sum_cpu(i)
        !
        ! First subtract all times of subroutines with a level-1 from the time of the subroutine
        ! under consideration
        !
        DO J=I+1,endsubrout(i)
          IF (level(I)+1==level(J)) THEN
            flat_time_cpu=flat_time_cpu-sum_cpu(j)
          ENDIF
        ENDDO
        !
        ! now try to find it in the temporary arrays for output
        !
        found=.FALSE.
        DO J=1,used_flat_timers
          IF (name_tag(i)==flat_names(j)) THEN
            flat_names(j)=name_tag(i)
            flat_times_cpu(j)=flat_times_cpu(j)+flat_time_cpu
            flat_calls(j)=flat_calls(j)+calls(i)
            found=.TRUE.
            EXIT
          ENDIF
        ENDDO
        IF (.NOT.FOUND) THEN
          flat_names(used_flat_timers+1)=name_tag(i)
          flat_times_cpu(used_flat_timers+1)=flat_time_cpu
          flat_calls(used_flat_timers+1)=calls(i)
          used_flat_timers=used_flat_timers+1
        ENDIF
      ENDDO
      !
      ! sort by time
      !
      DO I=2,used_flat_timers
        DO J=1,I-1
          IF (flat_times_cpu(I)>flat_times_cpu(J)) THEN
            flat_name=flat_names(J)
            flat_call=flat_calls(J)
            flat_time_cpu=flat_times_cpu(J)
            flat_names(J)=flat_names(I)
            flat_calls(J)=flat_calls(I)
            flat_times_cpu(J)=flat_times_cpu(I)
            flat_names(I)=flat_name
            flat_calls(I)=flat_call
            flat_times_cpu(I)=flat_time_cpu
          ENDIF
        ENDDO
      ENDDO
      !
      ! print out
      !
      sum_time=0._q
      WRITE(outfile,*)''
      WRITE(outfile,*)'Flat profile'
      WRITE(outfile,*)'============'
      WRITE(outfile,*)''
      WRITE(outfile,*)'            routine name                   CPU           calls'
      WRITE(outfile,*)'---------------------------------------------------------------'
      DO I=1,used_flat_timers
        WRITE(outfile,2000)flat_names(I),flat_times_cpu(I),flat_calls(I)
        sum_time=sum_time+flat_times_cpu(I)
      ENDDO
      WRITE(outfile,*)'---------------------------------------------------------------'
      WRITE(outfile,*)' summed up times ',sum_time
      !
      ! now calculate the absolute time used in all subroutines.
      !
      used_flat_timers=0
      flat_names(:)=''
      flat_times_cpu(:)=0._q
      flat_calls(:)=0
      skipped=0
      DO I=1,used_ptimers
        IF (PRESENT(PTHRESHOLD)) THEN
          IF (sum_cpu(i)<PTHRESHOLD) THEN
            CYCLE
          ENDIF
        ENDIF
        flat_time_cpu=sum_cpu(i)
        !
        ! now try to find it in the temporary arrays for output
        !
        found=.FALSE.
        DO J=1,used_flat_timers
          IF (name_tag(i)==flat_names(j)) THEN
            flat_names(j)=name_tag(i)
            flat_times_cpu(j)=flat_times_cpu(j)+flat_time_cpu
            flat_calls(j)=flat_calls(j)+calls(i)
            found=.TRUE.
            EXIT
          ENDIF
        ENDDO
        IF (.NOT.FOUND) THEN
          flat_names(used_flat_timers+1)=name_tag(i)
          flat_times_cpu(used_flat_timers+1)=flat_time_cpu
          flat_calls(used_flat_timers+1)=calls(i)
          used_flat_timers=used_flat_timers+1
        ENDIF
      ENDDO
      !
      ! sort by time
      !
      DO I=2,used_flat_timers
        DO J=1,I-1
          IF (flat_times_cpu(I)>flat_times_cpu(J)) THEN
            flat_name=flat_names(J)
            flat_call=flat_calls(J)
            flat_time_cpu=flat_times_cpu(J)
            flat_names(J)=flat_names(I)
            flat_calls(J)=flat_calls(I)
            flat_times_cpu(J)=flat_times_cpu(I)
            flat_names(I)=flat_name
            flat_calls(I)=flat_call
            flat_times_cpu(I)=flat_time_cpu
          ENDIF
        ENDDO
      ENDDO
      !
      ! print out
      !
      WRITE(outfile,*)''
      WRITE(outfile,*)'Accumulative profile'
      WRITE(outfile,*)'===================='
      WRITE(outfile,*)''
      WRITE(outfile,*)'            routine name                   CPU           calls'
      WRITE(outfile,*)'---------------------------------------------------------------'
      DO I=1,used_flat_timers
        WRITE(outfile,2000)flat_names(I),flat_times_cpu(I),flat_calls(I)
      ENDDO
      WRITE(outfile,*)'---------------------------------------------------------------'

      IF (PRESENT(PTHRESHOLD)) THEN
        WRITE(outfile,*)''
        WRITE(outfile,'(A,I6,A,F10.6,A)') 'skipped ',skipped,' routines with timings under ',pthreshold,' sec'
        WRITE(outfile,*)'mind that the threshold is used for the total exection time of a routine, '
        WRITE(outfile,*)'not the absolute time'
      ENDIF

      WRITE(outfile,*)''
      WRITE(outfile,'(A,4F10.6,A)')'Profiling took ',sum_cpu_profiling,' seconds'
      WRITE(outfile,'(A,F10.6,A)')'Profiling took ',sum_cpu_profiling_stop,' seconds'

      RETURN

1000  FORMAT(I5,2X,A30,2X,F14.6,'                         ',I10,2X,I2)
1001  FORMAT(I5,4X,A30,2X,F14.6,'                       ',I10,2X,I2)
1002  FORMAT(I5,6X,A30,2X,F14.6,'                     ',I10,2X,I2)
1003  FORMAT(I5,8X,A30,2X,F14.6,'                   ',I10,2X,I2)
1004  FORMAT(I5,10X,A30,2X,F14.6,'                 ',I10,2X,I2)
1005  FORMAT(I5,12X,A30,2X,F14.6,'               ',I10,2X,I2)
1006  FORMAT(I5,14X,A30,2X,F14.6,'             ',I10,2X,I2)
1007  FORMAT(I5,16X,A30,2X,F14.6,'           ',I10,2X,I2)
1008  FORMAT(I5,18X,A30,2X,F14.6,'         ',I10,2X,I2)
1009  FORMAT(I5,20X,A30,2X,F14.6,'       ',I10,2X,I2)
1010  FORMAT(I5,22X,A30,2X,F14.6,'     ',I10,2X,I2)
1011  FORMAT(I5,24X,A30,2X,F14.6,'   ',I10,2X,I2)
1012  FORMAT(I5,26X,A30,2X,F14.6,' ',I10,2X,I2)

2000  FORMAT(1X,A30,4X,F14.6,4X,I8)
      END SUBROUTINE PRINT_FULL_PROFILE

#ifdef MPI
      SUBROUTINE PRINT_STAT_PROFILE(IU)
      USE mpimy
      
      IMPLICIT NONE

      INTEGER :: rank, size, ierror
      INTEGER :: IU, I, J, K, skipped
      INTEGER :: outfile
      !
      ! used for flat profile
      !
      CHARACTER (LEN=30), DIMENSION(used_ptimers) :: flat_names
      REAL(q), DIMENSION(used_ptimers) :: flat_times_cpu
      INTEGER, DIMENSION(used_ptimers) :: flat_calls
      CHARACTER (LEN=30) :: flat_name
      REAL(q) :: flat_time_cpu, sum_time
      INTEGER :: used_flat_timers,flat_call
      LOGICAL :: FOUND
      CHARACTER (LEN=30) :: name_to_process
      REAL    :: time_to_send, vmin, vmax, vavg
      REAL, ALLOCATABLE, DIMENSION(:) :: times_cpu
      INTEGER :: calls_to_send, cmin, cmax, cavg, sum_calls
      INTEGER, ALLOCATABLE, DIMENSION(:) :: num_calls
      INTEGER :: rtmin,rtmax,rcmin,rcmax ! ranks have min/max times and # of calls

      outfile=iu
      skipped=0

      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
      !
      ! now calculate the absolute time used in all subroutines. Therefore we have to
      ! subtract the time of all called subroutines higher by one in level within the subroutine under
      ! consideration to get the absolute time of the subroutine. Additionally we
      ! have to sum all times up.
      !
      used_flat_timers=0
      flat_names(:)=''
      flat_times_cpu(:)=0._q
      flat_calls(:)=0
      skipped=0
      DO I=1,used_ptimers
        flat_time_cpu=sum_cpu(i)
        !
        ! First subtract all times of subroutines with a level-1 from the time of the subroutine
        ! under consideration
        !
        DO J=I+1,endsubrout(i)
          IF (level(I)+1==level(J)) THEN
            flat_time_cpu=flat_time_cpu-sum_cpu(j)
          ENDIF
        ENDDO
        !
        ! now try to find it in the temporary arrays for output
        !
        found=.FALSE.
        DO J=1,used_flat_timers
          IF (name_tag(i)==flat_names(j)) THEN
            flat_names(j)=name_tag(i)
            flat_times_cpu(j)=flat_times_cpu(j)+flat_time_cpu
            flat_calls(j)=flat_calls(j)+calls(i)
            found=.TRUE.
            EXIT
          ENDIF
        ENDDO
        IF (.NOT.FOUND) THEN
          flat_names(used_flat_timers+1)=name_tag(i)
          flat_times_cpu(used_flat_timers+1)=flat_time_cpu
          flat_calls(used_flat_timers+1)=calls(i)
          used_flat_timers=used_flat_timers+1
        ENDIF
      ENDDO
      !
      ! sort by time
      !
      DO I=2,used_flat_timers
        DO J=1,I-1
          IF (flat_times_cpu(I)>flat_times_cpu(J)) THEN
            flat_name=flat_names(J)
            flat_call=flat_calls(J)
            flat_time_cpu=flat_times_cpu(J)
            flat_names(J)=flat_names(I)
            flat_calls(J)=flat_calls(I)
            flat_times_cpu(J)=flat_times_cpu(I)
            flat_names(I)=flat_name
            flat_calls(I)=flat_call
            flat_times_cpu(I)=flat_time_cpu
          ENDIF
        ENDDO
      ENDDO

      IF(outfile.gt.0) THEN
          WRITE(outfile,*)''
          WRITE(outfile,*)'Flat statistics profile'
          WRITE(outfile,*)'============'
          WRITE(outfile,*)''
          WRITE(outfile,*)'            routine name                                CPU time                                      # calls'
          WRITE(outfile,*)'                                      min,rank          max,rank         avg        min,rank          max,rank         avg '
          WRITE(outfile,*)'-------------------------------+-----------------+------------------+-----------+----------------+----------------+---------+'
          FLUSH(outfile)
      ENDIF

      ALLOCATE(times_cpu(size))
      ALLOCATE(num_calls(size))

      CALL  MPI_BCAST(used_flat_timers, 1, MPI_INT, 0, MPI_COMM_WORLD, ierror)
      DO I=1,used_flat_timers
         name_to_process = flat_names(I)
         time_to_send = 0.
         calls_to_send = 0
         CALL MPI_BCAST(name_to_process, 30, MPI_CHAR, 0, MPI_COMM_WORLD, ierror)
         DO J=1,used_flat_timers
            IF (name_to_process==flat_names(J)) THEN
               time_to_send = flat_times_cpu(J)
               calls_to_send = flat_calls(J)
               EXIT
            ENDIF
         ENDDO

         CALL MPI_GATHER(time_to_send, 1, MPI_REAL, times_cpu, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierror)
         CALL MPI_GATHER(calls_to_send, 1, MPI_INTEGER, num_calls, 1, MPI_INT, 0, MPI_COMM_WORLD, ierror)

         sum_time=0._q
         sum_calls=0
         vmin = times_cpu(1)
         vmax = times_cpu(1)
         cmin = num_calls(1)
         cmax = num_calls(1)
         rtmin=0; rtmax=0; rcmin=0; rcmax=0
         DO J=1, size
            sum_time = sum_time + times_cpu(J)
            sum_calls = sum_calls + num_calls(J)
            IF(vmin > times_cpu(J)) rtmin = J
            IF(vmin > times_cpu(J)) vmin = times_cpu(J)
            IF(vmax < times_cpu(J)) rtmax = J
            IF(vmax < times_cpu(J)) vmax = times_cpu(J)
            IF(cmin > num_calls(J)) rcmin = J
            IF(cmin > num_calls(J)) cmin = num_calls(J)
            IF(cmax < num_calls(J)) rcmax = J
            IF(cmax < num_calls(J)) cmax = num_calls(J)
         ENDDO
         vavg = sum_time/size
         cavg = sum_calls/size
         IF(outfile.gt.0) &
              WRITE(outfile,2001)name_to_process,vmin,rtmin,vmax,rtmax,vavg,cmin,rcmin,cmax,rcmax,cavg
      ENDDO

      IF(outfile.gt.0) THEN 
          WRITE(outfile,*)'---------------------------------+-----------+-----------+-----------+---------+---------+---------+'
          FLUSH(outfile)
      ENDIF
      DEALLOCATE(times_cpu)
      DEALLOCATE(num_calls)
2001  FORMAT(1X,A30,2X,F10.3,1X,I6,2X,F10.3,1X,I6,2X,F10.3,2X,I8,1X,I6,2X,I8,1X,I6,2X,I8,2X)
      
      END SUBROUTINE PRINT_STAT_PROFILE
#endif
!***********************************************************************
!
! PRINT_ROUTINE_PROFILING
!> Prints the full profile.
!
!***********************************************************************

      SUBROUTINE PRINT_ROUTINE_PROFILE(IU,PROUTINE,PFLAT)

      IMPLICIT NONE

      INTEGER :: IU, I, J, K, skipped
      INTEGER :: outfile
      ! name of optional output file
      LOGICAL, OPTIONAL :: PFLAT
      CHARACTER (LEN=*) :: PROUTINE
      !
      ! used for flat profile
      !
      CHARACTER (LEN=30), DIMENSION(used_ptimers) :: flat_names
      REAL(q), DIMENSION(used_ptimers) :: flat_times_cpu
      INTEGER, DIMENSION(used_ptimers) :: flat_calls
      CHARACTER (LEN=30) :: flat_name
      REAL(q) :: flat_time_cpu, sum_time
      INTEGER :: used_flat_timers,flat_call
      LOGICAL :: FOUND

      outfile=iu
      !
      ! first check if a routine name is supplied whether it exists, otherwise
      ! default to the case where no routine name is given
      !
      found=.FALSE.
      DO i =1,used_ptimers
        IF(name_tag(i)==PROUTINE) found=.TRUE.
      ENDDO

      IF (.NOT.FOUND) THEN
        WRITE(outfile,*)''
        WRITE(outfile,'(A)')'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        WRITE(outfile,'(A,A30,A)')' Routine ',PROUTINE,' not found --> defaulting back to profiling all routines'
        WRITE(outfile,'(A)')'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        WRITE(outfile,*)''
        PROUTINE=''
      ENDIF

      skipped=0

      IF (.NOT.PFLAT) THEN
        WRITE(outfile,*)''
        WRITE(outfile,'(A,A30)')' PROFILE, profiling: ', PROUTINE
        WRITE(outfile,*)'=============================================='
        WRITE(outfile,*)''
        WRITE(outfile,'(A)')' index                 routine                    time                            calls level'
        WRITE(outfile,'(A)')'----------------------------------------------------------------------------------------------'
        FOUND=.FALSE.
        DO i=1,used_ptimers
          IF ((name_tag(I)==PROUTINE).or.(PROUTINE=='')) THEN
            FOUND=.TRUE.
            DO J=i,endsubrout(i)
              SELECT CASE (level(j))
                CASE (0)
                  WRITE(outfile,1000) j,name_tag(j),sum_cpu(j),calls(j),level(j)
                CASE (1)
                  WRITE(outfile,1001) j,name_tag(j),sum_cpu(j),calls(j),level(j)
                CASE (2)
                  WRITE(outfile,1002) j,name_tag(j),sum_cpu(j),calls(j),level(j)
                CASE (3)
                  WRITE(outfile,1003) j,name_tag(j),sum_cpu(j),calls(j),level(j)
                CASE (4)
                  WRITE(outfile,1004) j,name_tag(j),sum_cpu(j),calls(j),level(j)
                CASE (5)
                  WRITE(outfile,1005) j,name_tag(j),sum_cpu(j),calls(j),level(j)
                CASE (6)
                  WRITE(outfile,1006) j,name_tag(j),sum_cpu(j),calls(j),level(j)
                CASE (7)
                  WRITE(outfile,1007) j,name_tag(j),sum_cpu(j),calls(j),level(j)
                CASE (8)
                  WRITE(outfile,1008) j,name_tag(j),sum_cpu(j),calls(j),level(j)
                CASE (9)
                  WRITE(outfile,1009) j,name_tag(j),sum_cpu(j),calls(j),level(j)
                CASE (10)
                  WRITE(outfile,1010) j,name_tag(j),sum_cpu(j),calls(j),level(j)
                CASE (11)
                  WRITE(outfile,1011) j,name_tag(j),sum_cpu(j),calls(j),level(j)
                CASE (12)
                  WRITE(outfile,1012) j,name_tag(j),sum_cpu(j),calls(j),level(j)
                CASE DEFAULT
                  WRITE(outfile,1012) j,name_tag(j),sum_cpu(j),calls(j),level(j)
              END SELECT
            ENDDO
          ENDIF
        ENDDO
      ENDIF
      !
      ! now calculate the absolute time for the routine under consideration. Therefore we have to
      ! subtract the time of all called subroutines higher by one in level within the subroutine under
      ! consideration to get the absolute time of the subroutine. Additionally we
      ! have to sum all times up.
      !
      used_flat_timers=0
      flat_names(:)=''
      flat_times_cpu(:)=0._q
      flat_calls(:)=0
      DO I=1,used_ptimers
        IF ((name_tag(i)==PROUTINE).or.(PROUTINE=='')) THEN
          DO j=i,endsubrout(i)
            flat_time_cpu=sum_cpu(j)
            !
            ! First subtract all times of subroutines with a level-1 from the time of the subroutine
            ! under consideration
            !
            DO k=j+1,endsubrout(j)
              IF (level(j)+1==level(k)) THEN
                flat_time_cpu=flat_time_cpu-sum_cpu(k)
              ENDIF
            ENDDO
            !
            ! now try to find it in the temporary arrays for output
            !
            found=.FALSE.
            DO k=1,used_flat_timers
              IF (name_tag(j)==flat_names(k)) THEN
                flat_names(k)=name_tag(j)
                flat_times_cpu(k)=flat_times_cpu(k)+flat_time_cpu
                flat_calls(k)=flat_calls(k)+calls(j)
                found=.TRUE.
                EXIT
              ENDIF
            ENDDO
            IF (.NOT.FOUND) THEN
              flat_names(used_flat_timers+1)=name_tag(j)
              flat_times_cpu(used_flat_timers+1)=flat_time_cpu
              flat_calls(used_flat_timers+1)=calls(j)
              used_flat_timers=used_flat_timers+1
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      !
      ! sort by time
      !
      DO I=2,used_flat_timers
        DO J=1,I-1
          IF (flat_times_cpu(I)>flat_times_cpu(J)) THEN
            flat_name=flat_names(J)
            flat_call=flat_calls(J)
            flat_time_cpu=flat_times_cpu(J)
            flat_names(J)=flat_names(I)
            flat_calls(J)=flat_calls(I)
            flat_times_cpu(J)=flat_times_cpu(I)
            flat_names(I)=flat_name
            flat_calls(I)=flat_call
            flat_times_cpu(I)=flat_time_cpu
          ENDIF
        ENDDO
      ENDDO
      !
      ! print out
      !
      sum_time=0._q
      WRITE(outfile,*)''
      WRITE(outfile,*)'Flat profile'
      WRITE(outfile,*)'============'
      WRITE(outfile,*)''
      WRITE(outfile,*)'            routine name                   CPU           calls'
      WRITE(outfile,*)'---------------------------------------------------------------'
      DO I=1,used_flat_timers
        WRITE(outfile,2000)flat_names(I),flat_times_cpu(I),flat_calls(I)
        sum_time=sum_time+flat_times_cpu(I)
      ENDDO
      WRITE(outfile,*)'---------------------------------------------------------------'
      WRITE(outfile,*)' summed up times ',sum_time

      WRITE(outfile,'(A,4F10.6,A)')'Profiling took ',sum_cpu_profiling,' seconds'
      WRITE(outfile,'(A,F10.6,A)')'Profiling took ',sum_cpu_profiling_stop,' seconds'

      RETURN

1000  FORMAT(I5,2X,A30,2X,F14.6,'                         ',I10,2X,I2)
1001  FORMAT(I5,4X,A30,2X,F14.6,'                       ',I10,2X,I2)
1002  FORMAT(I5,6X,A30,2X,F14.6,'                     ',I10,2X,I2)
1003  FORMAT(I5,8X,A30,2X,F14.6,'                   ',I10,2X,I2)
1004  FORMAT(I5,10X,A30,2X,F14.6,'                 ',I10,2X,I2)
1005  FORMAT(I5,12X,A30,2X,F14.6,'               ',I10,2X,I2)
1006  FORMAT(I5,14X,A30,2X,F14.6,'             ',I10,2X,I2)
1007  FORMAT(I5,16X,A30,2X,F14.6,'           ',I10,2X,I2)
1008  FORMAT(I5,18X,A30,2X,F14.6,'         ',I10,2X,I2)
1009  FORMAT(I5,20X,A30,2X,F14.6,'       ',I10,2X,I2)
1010  FORMAT(I5,22X,A30,2X,F14.6,'     ',I10,2X,I2)
1011  FORMAT(I5,24X,A30,2X,F14.6,'   ',I10,2X,I2)
1012  FORMAT(I5,26X,A30,2X,F14.6,' ',I10,2X,I2)

2000  FORMAT(1X,A30,4X,F14.6,4X,I8)
      END SUBROUTINE PRINT_ROUTINE_PROFILE


!***********************************************************************
!
! PRINT_LEVEL_PROFILING
!> Prints the profile up to maximum level depth.
!
!***********************************************************************
      SUBROUTINE PRINT_LEVEL_PROFILE(IU,PLEVEL,PFLAT)

      IMPLICIT NONE

      INTEGER :: IU, I, J, K, skipped
      INTEGER :: outfile
      ! threshold for the time
      INTEGER :: PLEVEL
      LOGICAL, OPTIONAL :: PFLAT
      !
      ! used for flat profile
      !
      CHARACTER (LEN=30), DIMENSION(used_ptimers) :: flat_names
      REAL(q), DIMENSION(used_ptimers) :: flat_times_cpu
      INTEGER, DIMENSION(used_ptimers) :: flat_calls
      CHARACTER (LEN=30) :: flat_name
      REAL(q) :: flat_time_cpu, sum_time
      INTEGER :: used_flat_timers,flat_call
      LOGICAL :: FOUND

      outfile=iu
      skipped=0

      IF (.NOT.PFLAT) THEN
        WRITE(outfile,*)''
        WRITE(outfile,'(A,I7)')' PROFILE, used timers: ', used_ptimers
        WRITE(outfile,*)'============================='
        WRITE(outfile,*)''
        WRITE(outfile,'(A)')' index                 routine                    time                            calls level'
        WRITE(outfile,'(A)')'----------------------------------------------------------------------------------------------'
        DO I=1,used_ptimers
          IF (level(i)>PLEVEL) THEN
            skipped=skipped+1
            CYCLE
          ENDIF
          SELECT CASE (level(i))
            CASE (0)
              WRITE(outfile,1000) i,name_tag(i),sum_cpu(i),calls(i),level(i)
            CASE (1)
              WRITE(outfile,1001) i,name_tag(i),sum_cpu(i),calls(i),level(i)
            CASE (2)
              WRITE(outfile,1002) i,name_tag(i),sum_cpu(i),calls(i),level(i)
            CASE (3)
              WRITE(outfile,1003) i,name_tag(i),sum_cpu(i),calls(i),level(i)
            CASE (4)
              WRITE(outfile,1004) i,name_tag(i),sum_cpu(i),calls(i),level(i)
            CASE (5)
              WRITE(outfile,1005) i,name_tag(i),sum_cpu(i),calls(i),level(i)
            CASE (6)
              WRITE(outfile,1006) i,name_tag(i),sum_cpu(i),calls(i),level(i)
            CASE (7)
              WRITE(outfile,1007) i,name_tag(i),sum_cpu(i),calls(i),level(i)
            CASE (8)
              WRITE(outfile,1008) i,name_tag(i),sum_cpu(i),calls(i),level(i)
            CASE (9)
              WRITE(outfile,1009) i,name_tag(i),sum_cpu(i),calls(i),level(i)
            CASE (10)
              WRITE(outfile,1010) i,name_tag(i),sum_cpu(i),calls(i),level(i)
            CASE (11)
              WRITE(outfile,1011) i,name_tag(i),sum_cpu(i),calls(i),level(i)
            CASE (12)
              WRITE(outfile,1012) i,name_tag(i),sum_cpu(i),calls(i),level(i)
            CASE DEFAULT
              WRITE(outfile,1012) i,name_tag(i),sum_cpu(i),calls(i),level(i)
          END SELECT
        ENDDO
        WRITE(outfile,*)''
        WRITE(outfile,'(A,I6,A,I3)') 'skipped ',skipped,' routines with a level over ',PLEVEL
      ENDIF
      !
      ! now calculate the absolute time used in all subroutines. Therefore we have to
      ! subtract the time of all called subroutines higher by one in level within the subroutine under
      ! consideration to get the absolute time of the subroutine. Additionally we
      ! have to sum all times up.
      !
      used_flat_timers=0
      flat_names(:)=''
      flat_times_cpu(:)=0._q
      flat_calls(:)=0
      skipped=0
      DO I=1,used_ptimers
          IF (level(i)>PLEVEL) THEN
            CYCLE
          ENDIF
        flat_time_cpu=sum_cpu(i)
        !
        ! First subtract all times of subroutines with a level-1 from the time of the subroutine
        ! under consideration
        !
        DO J=I+1,endsubrout(i)
          IF ((level(I)+1==level(J)).AND.(level(i)<PLEVEL)) THEN
            flat_time_cpu=flat_time_cpu-sum_cpu(j)
          ENDIF
        ENDDO
        !
        ! now try to find it in the temporary arrays for output
        !
        found=.FALSE.
        DO J=1,used_flat_timers
          IF (name_tag(i)==flat_names(j)) THEN
            flat_names(j)=name_tag(i)
            flat_times_cpu(j)=flat_times_cpu(j)+flat_time_cpu
            flat_calls(j)=flat_calls(j)+calls(i)
            found=.TRUE.
            EXIT
          ENDIF
        ENDDO
        IF (.NOT.FOUND) THEN
          flat_names(used_flat_timers+1)=name_tag(i)
          flat_times_cpu(used_flat_timers+1)=flat_time_cpu
          flat_calls(used_flat_timers+1)=calls(i)
          used_flat_timers=used_flat_timers+1
        ENDIF
      ENDDO
      !
      ! sort by time
      !
      DO I=2,used_flat_timers
        DO J=1,I-1
          IF (flat_times_cpu(I)>flat_times_cpu(J)) THEN
            flat_name=flat_names(J)
            flat_call=flat_calls(J)
            flat_time_cpu=flat_times_cpu(J)
            flat_names(J)=flat_names(I)
            flat_calls(J)=flat_calls(I)
            flat_times_cpu(J)=flat_times_cpu(I)
            flat_names(I)=flat_name
            flat_calls(I)=flat_call
            flat_times_cpu(I)=flat_time_cpu
          ENDIF
        ENDDO
      ENDDO
      !
      ! print out
      !
      sum_time=0._q
      WRITE(outfile,*)''
      WRITE(outfile,*)'Flat profile'
      WRITE(outfile,*)'============'
      WRITE(outfile,*)''
      WRITE(outfile,*)'            routine name                   CPU          calls'
      WRITE(outfile,*)'--------------------------------------------------------------'
      DO I=1,used_flat_timers
        WRITE(outfile,2000)flat_names(I),flat_times_cpu(I),flat_calls(I)
        sum_time=sum_time+flat_times_cpu(I)
      ENDDO
      WRITE(outfile,*)'---------------------------------------------------------------'
      WRITE(outfile,*)' summed up times ',sum_time

      WRITE(outfile,*)''
      WRITE(outfile,'(A,I6,A,I3)') 'skipped ',skipped,' routines with a level over ',PLEVEL
      WRITE(outfile,'(A,4F10.6,A)')'Profiling took ',sum_cpu_profiling,' seconds'
      WRITE(outfile,'(A,F10.6,A)')'Profiling took ',sum_cpu_profiling_stop,' seconds'

      RETURN

1000  FORMAT(I5,2X,A30,2X,F14.6,'                         ',I10,2X,I2)
1001  FORMAT(I5,4X,A30,2X,F14.6,'                       ',I10,2X,I2)
1002  FORMAT(I5,6X,A30,2X,F14.6,'                     ',I10,2X,I2)
1003  FORMAT(I5,8X,A30,2X,F14.6,'                   ',I10,2X,I2)
1004  FORMAT(I5,10X,A30,2X,F14.6,'                 ',I10,2X,I2)
1005  FORMAT(I5,12X,A30,2X,F14.6,'               ',I10,2X,I2)
1006  FORMAT(I5,14X,A30,2X,F14.6,'             ',I10,2X,I2)
1007  FORMAT(I5,16X,A30,2X,F14.6,'           ',I10,2X,I2)
1008  FORMAT(I5,18X,A30,2X,F14.6,'         ',I10,2X,I2)
1009  FORMAT(I5,20X,A30,2X,F14.6,'       ',I10,2X,I2)
1010  FORMAT(I5,22X,A30,2X,F14.6,'     ',I10,2X,I2)
1011  FORMAT(I5,24X,A30,2X,F14.6,'   ',I10,2X,I2)
1012  FORMAT(I5,26X,A30,2X,F14.6,' ',I10,2X,I2)

2000  FORMAT(1X,A30,4X,F14.6,4X,I8)
      END SUBROUTINE PRINT_LEVEL_PROFILE

      END MODULE


!*******************************************************************************
!
! START_PROFILING
!> Searches for already existing timings in the current level and
!> starts the timer.
!>
!> Search for the subroutine \n
!> 3 cases are possible (where 2 are actually treated in the same way here) \n
!> 1) overall first appearance of a subroutine
!>    - if we are at the end of the arrays we insert it there
!>    - due to an if clause, a call to a subroutine is triggered at a later
!>      point in the execution of the code. This subroutine has to be inserted
!>      at the appropriate place in the profiling-arrays. \n
!> 2) The routine is already listed in the profiling-arrays, update the
!>    corresponding information.
!>
!> Most important part is to keep in the indices the right values.
!> - actlevel is the actual level, so if we enter a new subroutine we
!>   increase it by one at the end, so that the value is correct for
!>   the next START_PROFILING call. (decreased when STOP_PROFILING is called.
!> - actindex is the position of the actual routine, thats it.
!> - endindex is always set to the index from which START_PROFILING should
!>   start the search in the arrays. So normally to the end of the outer
!>   subroutine (if already there), or to actindex if the routine is newly added.
!
!********************************************************************************

      SUBROUTINE START_PROFILING(TAG)
      USE profiling
      IMPLICIT NONE

      INTEGER i,j,j1
      CHARACTER (LEN=*) :: TAG
      REAL(q) TC
      LOGICAL :: FOUND
      INTEGER THREADS

!$    INTEGER id,il
!$    INTEGER, EXTERNAL :: omp_get_thread_num,omp_get_level,omp_get_ancestor_thread_num

#ifdef NVTX
      CALL NVTX_PUSH(TAG//C_NULL_CHAR)
#endif

!$    id=omp_get_thread_num()
!$    DO il=1,omp_get_level()-1
!$       id=id+omp_get_ancestor_thread_num(il)
!$    ENDDO
!$    IF (id==0) THEN
!!$OMP MASTER
      !
      ! start timing of profiling
      !
      CALL VTIME2(TC)
      cpu_profiling=TC
      !
      ! store opened_tag name
      !
      opened_tag(actlevel) = tag
      !                     
      ! inidcates whether routine already exists in profiling arrays
      !
      found=.FALSE.

      DO i=endindex,startindex(actlevel-1),-1
        IF (level(i)==actlevel) THEN
          IF (name_tag(i)==tag) THEN
            !
            ! because routines can be present at different positions in an "open"
            ! subroutine in the same nesting level, we have to find the first in the
            ! profiling arrays after an "open" subroutine
            !
            actindex=i
            found=.TRUE.
          ENDIF
        ENDIF
      ENDDO

      CALL VTIME2(TC)
      sum_cpu_profiling(1)=sum_cpu_profiling(1)+TC-cpu_profiling
      cpu_profiling=TC

      IF (found) THEN
        !
        ! now an existing entry in the profiling arrays was found, so we have to
        ! update this one and set all the other informations appropriatly
        !
        calls(actindex)=calls(actindex)+1
        startindex(actlevel)=actindex
        endindex=endsubrout(actindex)
        CALL VTIME2(TC)
        cpu(actlevel)=TC
        ! increasing the actual level is necessary if another routine is opened, for
        ! stopping the time it is irrelevant
        actlevel=actlevel+1
        !
        ! stop timing of profiling
        !
        CALL VTIME2(TC)
        sum_cpu_profiling(2)=sum_cpu_profiling(2)+TC-cpu_profiling
      ELSE
        !
        ! new subroutine is found and should be inserted at endsubrout(actindex)+1
        !
        ! now check whether we are in an existing routine or at the end of
        ! of the profiling array
        !
        IF (endindex==used_ptimers) THEN
          used_ptimers=used_ptimers+1
          endindex=used_ptimers
          actindex=endindex
          !
          ! level has to be one higher than the last "open" subroutine
          !
          level(actindex)=actlevel
          startindex(actlevel)=actindex
          name_tag(actindex)=TAG
          calls(actindex)=1
          CALL VTIME2(TC)
          cpu(actlevel)=TC
          ! increasing the actual level is necessary if another routine is opened, for
          ! stopping the time it is irrelevant
          actlevel=actlevel+1
          !
          ! stop timing of profiling
          !
          CALL VTIME2(TC)
          sum_cpu_profiling(3)=sum_cpu_profiling(3)+TC-cpu_profiling
        ELSE
          !
          ! we are adding a timer inside an existing subroutine
          !
          IF (level(actindex)>actlevel) THEN
            !
            ! update actindex: forward search for the index of the next existing entry
            ! with the same or lower level as the entry we are going to add
            !
            DO i=actindex,endindex
              IF (level(i)<=actlevel) THEN
                 actindex=i-1
                 EXIT
              ENDIF
            ENDDO
          ENDIF
          !
          ! move all entries starting from endsubrout(actindex)+1 one step up, first increase used_timers
          !
          used_ptimers=used_ptimers+1
          DO j=used_ptimers,endsubrout(actindex)+2,-1
            j1=j-1
            level(j)=level(j1)
            sum_cpu(j)=sum_cpu(j1)
            name_tag(j)=name_tag(j1)
            calls(j)=calls(j1)
            endsubrout(j)=endsubrout(j1)+1
          ENDDO
          !
          ! now insert new subroutine at endsubrout(actindex)+1
          !
          actindex=endsubrout(actindex)+1
          endindex=actindex
          endsubrout(actindex)=endindex
          !
          ! level has to be one higher than the last "open" subroutine
          !
          level(actindex)=actlevel
          startindex(actlevel)=actindex
          !
          ! set timing values to zero, since the old values are still stored there
          !
          sum_cpu(actindex)=0._q
          name_tag(actindex)=TAG
          calls(actindex)=1
          CALL VTIME2(TC)
          cpu(actlevel)=TC
          ! increasing the actual level is necessary if another routine is opened, for
          ! stopping the time it is irrelevant
          actlevel=actlevel+1
          !
          ! stop timing of profiling
          !
          CALL VTIME2(TC)
          sum_cpu_profiling(4)=sum_cpu_profiling(4)+TC-cpu_profiling
        ENDIF
      ENDIF
!!$OMP END MASTER
!$    ENDIF
      END SUBROUTINE START_PROFILING

!***********************************************************************
!
! STOP_PROFILING
!> Searches for the existing timings and stops the timer.
!
!***********************************************************************

      SUBROUTINE STOP_PROFILING(TAG)
      USE profiling
      IMPLICIT NONE

      INTEGER :: i,subindex,endsub
      CHARACTER (LEN=*) :: TAG
      REAL(q) TC
!$    INTEGER id,il
!$    INTEGER, EXTERNAL :: omp_get_thread_num,omp_get_level,omp_get_ancestor_thread_num
!$    id=omp_get_thread_num()
!$    DO il=1,omp_get_level()-1
!$       id=id+omp_get_ancestor_thread_num(il)
!$    ENDDO
!$    IF (id==0) THEN
!!$OMP MASTER
      !
      ! first stop the time before we begin to search in the arrays
      ! this time also used for timing to profiling
      !
      CALL VTIME2(TC)
      cpu_profiling_stop=TC
      !
      ! sanity check 
      !
      if( opened_tag( actlevel-1) /= tag ) then
         write(*,*) " !!! ERROR !!!, unclosed tag present: ", opened_tag( actlevel-1 )," ", tag 
      endif
      !
      ! find actual subroutine to stop timing
      !
      DO i=actindex,1,-1
        IF (name_tag(i)==tag) THEN
          subindex=i
          EXIT
        ENDIF
      ENDDO
      !
      ! set the level back to the correct value
      ! and search for endsubroutine (which
      ! gives you the array element number where the routine
      ! ends).
      !
      actlevel=level(subindex)
      endindex=used_ptimers
      endsubrout(subindex)=endindex
      DO i=actindex+1, used_ptimers
        IF(actlevel>=level(i)) then
          endsubrout(subindex)=i-1
          exit
        endif
      enddo
      !
      ! search for the end of the enclosing subroutine (is important for the
      ! next START_PROFILING call, otherwise an already existing subroutine
      ! will not be found in the corresponding level
      !
      DO i=endsubrout(subindex),used_ptimers
        IF (actlevel>level(i)) THEN
          endindex=i-1
          EXIT
        ENDIF
      ENDDO

      sum_cpu(subindex)=sum_cpu(subindex)+TC-cpu(actlevel)
      !
      ! stop timing of profiling
      !
      CALL VTIME2(TC)
      sum_cpu_profiling_stop=sum_cpu_profiling_stop+TC-cpu_profiling_stop
!!$OMP END MASTER
!$    ENDIF

100   FORMAT (A10,4X,A30,4X,'STOP ',4I9,4X,L,I9)

#ifdef NVTX
      CALL NVTX_POP()
#endif
      END SUBROUTINE STOP_PROFILING
