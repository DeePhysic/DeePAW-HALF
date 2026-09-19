#include "symbol.inc"

MODULE ml_interface_writer
   USE mpimy
   USE prec
   USE tutor, ONLY: vtutor, isError, RandomSeed, argument
   
   CONTAINS

!=======================================================================
!
! This subroutine does some post processing after the energy, forces 
! and stress was predicted in a previous step. This routine also does
! the writing to the OUTCAR file.
!
!=======================================================================


   SUBROUTINE ML_TO_VASP_WRITING(DYN, GRIDC, IO, LATT_CUR, LREMOVE_DRIFT, LWRITE_FORCE, &
                                 NSTEP, T_INFO, TOTEN, TIFOR, TSIF, LDELTADFT)
      USE base
      USE classicfields
      USE constant
      USE ini
      USE lattice
      USE poscar
      USE mgrid
      USE mpimy
      USE msymmetry
      USE poscar
      ! Input and output variables
      TYPE (dynamics)               :: DYN
      TYPE (grid_3d)                :: GRIDC             ! grid for potentials/charge
      TYPE (in_struct)              :: IO
      TYPE (latt), INTENT(IN)       :: LATT_CUR          ! Bravais lattice information
      TYPE (type_info), INTENT(IN)  :: T_INFO            ! type and position information
      LOGICAL, INTENT(IN)           :: LREMOVE_DRIFT     
      LOGICAL, INTENT(IN)           :: LWRITE_FORCE      ! If LWRITE_FORCE=.TRUE., write forces in OUTCAR file.
      INTEGER, INTENT(IN)           :: NSTEP             ! MD step
      LOGICAL, INTENT(IN)           :: LDELTADFT
      REAL(q), INTENT(INOUT)        :: TOTEN             ! total energy
      REAL(q), INTENT(INOUT)        :: TIFOR(:,:)        !(1:3,1:NIOND) ! forces on all atoms
      REAL(q), INTENT(IN)           :: TSIF(:,:)         !(1:3,1:3)      ! stress tensor
      ! Local variables
      REAL(q)   :: EV2KB
      REAL(q)   :: FAC
      INTEGER   :: I
      INTEGER   :: INIONS
      INTEGER   :: IXYZ
      INTEGER   :: JXYZ
      INTEGER   :: NOFFS
      INTEGER   :: NI
      INTEGER   :: NT
      REAL(q)   :: OFIELD_E                        ! energy contribution from order field
      REAL(q)   :: OFIELD_FOR(1:3,T_INFO%NIONS)
      REAL(q)   :: PRESS
      REAL(q)   :: pressure
      REAL(q)   :: TMP(1:6)
      REAL(q)   :: VEL(1:3)
      REAL(q)   :: VTMP(1:3)


      ! Constants.
      EV2KB = EVTOJ*1E+22_q/LATT_CUR%OMEGA
      FAC  = AMTOKG* &
             1E5_q* &
             1E5_q* &
             1E30_q* &
             1E-8_q
      ! Add the external order field, if necessary
      CALL OFIELD(GRIDC%COMM, IO, LATT_CUR, DYN, T_INFO, OFIELD_FOR, OFIELD_E)
      TIFOR=TIFOR+OFIELD_FOR
      TOTEN=TOTEN+OFIELD_E
      ! Remove drift.
      IF (DYN%IBRION/=0) THEN
         ! remove drift from the forces
         IF (LREMOVE_DRIFT) CALL SYMVEC(T_INFO%NIONS,TIFOR)
      ENDIF
      ! Only write at every DYN%NBLOCK steps
      IF (MOD(NSTEP,DYN%NBLOCK)==0) THEN
      ! Output to OUTCAR file.
         IF (IO%IU6>=0) THEN
            PRESS=(TSIF(1,1)+TSIF(2,2)+TSIF(3,3))/3._q- &
                  DYN%PSTRESS/(EVTOJ*1E22_q)*LATT_CUR%OMEGA
            IF (LDELTADFT) THEN
               WRITE(IO%IU6,29) ""
               WRITE(IO%IU6,21)      (TSIF(IXYZ,IXYZ),IXYZ=1,3),      TSIF (1,2),     TSIF(2,3),     TSIF (3,1), &
                                     (EV2KB*TSIF(IXYZ,IXYZ),IXYZ=1,3),EV2KB*TSIF(1,2),EV2KB*TSIF(2,3),EV2KB*TSIF(3,1), &
                                     PRESS*EV2KB,DYN%PSTRESS
            ELSE
               WRITE(IO%IU6,1)       (TSIF(IXYZ,IXYZ),IXYZ=1,3),      TSIF (1,2),     TSIF(2,3),     TSIF (3,1), &
                                     (EV2KB*TSIF(IXYZ,IXYZ),IXYZ=1,3),EV2KB*TSIF(1,2),EV2KB*TSIF(2,3),EV2KB*TSIF(3,1), &
                                     PRESS*EV2KB,DYN%PSTRESS
            ENDIF
            TMP = 0._q
            NOFFS = 0
            DO NT=1,T_INFO%NTYP
               DO NI=1+NOFFS,T_INFO%NITYP(NT)+NOFFS
                  VEL(1) = DYN%VEL(1,NI)/DYN%POTIM
                  VEL(2) = DYN%VEL(2,NI)/DYN%POTIM
                  VEL(3) = DYN%VEL(3,NI)/DYN%POTIM
                  CALL  DIRKAR( 1, VEL, LATT_CUR%A )
                  TMP(1) = TMP(1) + VEL(1)*VEL(2) * T_INFO%POMASS(NT)
                  TMP(2) = TMP(2) + VEL(2)*VEL(3) * T_INFO%POMASS(NT)
                  TMP(3) = TMP(3) + VEL(3)*VEL(1) * T_INFO%POMASS(NT)
                  TMP(4) = TMP(4) + VEL(1)*VEL(1) * T_INFO%POMASS(NT)
                  TMP(5) = TMP(5) + VEL(2)*VEL(2) * T_INFO%POMASS(NT)
                  TMP(6) = TMP(6) + VEL(3)*VEL(3) * T_INFO%POMASS(NT)
               ENDDO
               NOFFS = NOFFS + T_INFO%NITYP(NT)
            ENDDO
            TMP = TMP/LATT_CUR%OMEGA*FAC         ! NOW TMP IS IN KB
            WRITE(IO%IU6,'(''  kinetic pressure (ideal gas correction) = '',F9.2,'' kB'')') (TMP(4)+TMP(5)+TMP(6))/3
            pressure=(TSIF(1,1)+TSIF(2,2)+TSIF(3,3))*EV2KB
            pressure=pressure + (TMP(4)+TMP(5)+TMP(6))
            pressure=pressure/3.
            WRITE(IO%IU6,'(''  total pressure  = '',F9.2,'' kB'')') pressure
            WRITE(IO%IU6,'(''  Total+kin. '',F9.3,5F12.3)') TSIF(1,1)*EV2KB + TMP(4),  &
                                                            TSIF(2,2)*EV2KB + TMP(5),  &
                                                            TSIF(3,3)*EV2KB + TMP(6),  &
                                                            TSIF(1,2)*EV2KB + TMP(1),  &
                                                            TSIF(2,3)*EV2KB + TMP(2),  &
                                                            TSIF(3,1)*EV2KB + TMP(3)
            WRITE(IO%IU6,2) LATT_CUR%OMEGA, &
                 ((LATT_CUR%A(IXYZ,JXYZ),IXYZ=1,3),(LATT_CUR%B(IXYZ,JXYZ),IXYZ=1,3),JXYZ=1,3), &
                 (LATT_CUR%ANORM(IXYZ),IXYZ=1,3),(LATT_CUR%BNORM(IXYZ),IXYZ=1,3)
            IF(LWRITE_FORCE) THEN
               IF (LDELTADFT) THEN
                  WRITE(IO%IU6,23)
               ELSE
                  WRITE(IO%IU6,3)
               ENDIF
               DO INIONS=1, T_INFO%NIONS
                  VTMP(1)=T_INFO%POSION(1,INIONS)
                  VTMP(2)=T_INFO%POSION(2,INIONS)
                  VTMP(3)=T_INFO%POSION(3,INIONS)
                  CALL  DIRKAR(1,VTMP(1),LATT_CUR%A(1,1))
                  WRITE(IO%IU6,4) (VTMP(I),I=1,3),(TIFOR(I,INIONS),I=1,3)
               ENDDO
               VTMP(1:3)=0.0_q
               DO INIONS=1, T_INFO%NIONS
                  DO IXYZ=1, 3
                     VTMP(IXYZ)=VTMP(IXYZ)+TIFOR(IXYZ,INIONS)
                  ENDDO
               ENDDO
               WRITE(IO%IU6,5) VTMP
            ENDIF
            WRITE(IO%IU6,6)
            IF (LDELTADFT) THEN
               WRITE(IO%IU6,27) '  DELTA FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)',TOTEN,TOTEN,TOTEN
            ELSE
               WRITE(IO%IU6,7) '  ML FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)',TOTEN,TOTEN,TOTEN
            ENDIF
         ENDIF
      ENDIF
      ! Convert to enthalpy
      IF (DYN%PSTRESS/=0) THEN
         TOTEN=TOTEN+DYN%PSTRESS/(EVTOJ*1E22_q)*LATT_CUR%OMEGA
         ! Again printing only at every NBLOCK step
         IF (MOD(NSTEP,DYN%NBLOCK)==0) THEN
            IF (LDELTADFT) THEN
               IF (IO%IU6>=0) WRITE(IO%IU6,28) TOTEN,DYN%PSTRESS/(EVTOJ*1E22_q)*LATT_CUR%OMEGA
            ELSE
               IF (IO%IU6>=0) WRITE(IO%IU6,8) TOTEN,DYN%PSTRESS/(EVTOJ*1E22_q)*LATT_CUR%OMEGA
            ENDIF
         ENDIF
      ENDIF
! Format
1        FORMAT('  ML FORCE on cell =-STRESS in cart. coord. units (eV/cell)'/ &
                 '  Direction',4X,'XX', 10X,'YY', 10X,'ZZ', 10X,'XY', 10X,'YZ', 10X,'ZX'/ &
                 '  --------------------------------------------------------------------------------------'/ &
                 '  Total: ',6F12.5/ &
                 '  in kB  ',6F12.5/, &
                 '  external pressure = ',F11.2,' kB', &
                 '  Pullay stress = ',F11.2,' kB'/)
2        FORMAT( '  volume of cell :  ',F10.2/ &
                 '      direct lattice vectors',17X,'reciprocal lattice vectors'/ &
                 3(2(3X,3F13.9)/)/ &
                 '  length of vectors'/ &
                  (2(3X,3F13.9)/) /)
3        FORMAT( '  POSITION    ',35X,'TOTAL-FORCE (eV/Angst) (ML)'/ &
                 ' -----------------------------------------------------------------------------------')
4        FORMAT((3F13.5,3X,3F14.6))
5        FORMAT( ' ----------------------------------------------', &
                 '-------------------------------------',/ &
                 '    total drift:      ',20X,3F14.6)
6        FORMAT (5X, //, &
                '----------------------------------------------------', &
                '----------------------------------------------------'//)
7        FORMAT(/ &
                A/ &
                '  ---------------------------------------------------'/ &
                '  free  energy ML TOTEN  = ',F18.8,' eV'// &
                '  ML energy  without entropy=',F18.8, &
                '  ML energy(sigma->0) =',F18.8/)
8        FORMAT ('  enthalpy is ML TOTEN    = ',F18.8,' eV   P V=',F18.8/)

21        FORMAT(/'  DELTA FORCE on cell =-STRESS in cart. coord. units (eV/cell)'/ &
                 '  Direction',4X,'XX', 10X,'YY', 10X,'ZZ', 10X,'XY', 10X,'YZ', 10X,'ZX'/ &
                 '  --------------------------------------------------------------------------------------'/ &
                 '  Total: ',6F12.5/ &
                 '  in kB  ',6F12.5/, &
                 '  external pressure = ',F11.2,' kB', &
                 '  Pullay stress = ',F11.2,' kB'/)
23        FORMAT( '  POSITION    ',35X,'TOTAL-FORCE (eV/Angst) (DELTA)'/ &
                 ' -----------------------------------------------------------------------------------')
27        FORMAT(/ &
                A/ &
                '  ---------------------------------------------------'/ &
                '  free  energy DELTA TOTEN  = ',F18.8,' eV'// &
                '  DELTA energy  without entropy=',F18.8, &
                '  DELTA energy(sigma->0) =',F18.8/)
28        FORMAT ('  enthalpy is DELTA TOTEN    = ',F18.8,' eV   P V=',F18.8/)
29        FORMAT(/ &
                A/ &
                '--------------------------------------------------------------------------------------------------------')
   END SUBROUTINE ML_TO_VASP_WRITING

END MODULE ml_interface_writer
