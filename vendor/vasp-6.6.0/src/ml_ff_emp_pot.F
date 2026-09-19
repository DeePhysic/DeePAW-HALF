#if defined(MPI) || defined(MPI_CHAIN)
#include "symbol.inc"

      MODULE EMP_POT
         USE ML_FF_PREC
         USE ML_FF_CONSTANT
         USE MPI_DATA
        
      CONTAINS 


      SUBROUTINE INITIALIZE_NEIGHBOR_LIST_MODEL_POTENTIAL_FAST(EMP, PAR_SUP_HANDLE, &
                 LATTICE, INVERSE_LATTICE, POSION, NTYP, LNTYP, &
                 RCUT1, RCUT2, NEAREST_NEIGHBOR_RADIAL, NEAREST_NEIGHBOR_ANGULAR, W1, W2)
         USE ML_FF_STRUCT, ONLY : EMP_POT_TYPE, PARALLEL_SUPER, VERLET_CELL_LIST_TYPE
         IMPLICIT NONE
         TYPE (EMP_POT_TYPE), INTENT(INOUT) :: EMP
         TYPE (PARALLEL_SUPER), INTENT(IN) :: PAR_SUP_HANDLE
         REAL(q), INTENT(IN)  :: LATTICE(:, :)
         REAL(q), INTENT(IN)  :: INVERSE_LATTICE(:, :)
         REAL(q), INTENT(IN)  :: POSION(:, :)
         INTEGER, INTENT(IN)  :: NTYP
         INTEGER, INTENT(IN)  :: LNTYP(:)
         Real(q), INTENT(IN), OPTIONAL :: RCUT1
         Real(q), INTENT(IN), OPTIONAL :: RCUT2
         TYPE (VERLET_CELL_LIST_TYPE), INTENT(IN), OPTIONAL :: NEAREST_NEIGHBOR_RADIAL
         TYPE (VERLET_CELL_LIST_TYPE), INTENT(IN), OPTIONAL :: NEAREST_NEIGHBOR_ANGULAR
         REAL(q), INTENT(IN), OPTIONAL :: W1
         REAL(q), INTENT(IN), OPTIONAL :: W2

         ! Start profiling
         PROFILING_START('init_emp_pot_neighbor')
         
         CALL EMP%NEAREST_NEIGHBOR%INIT_NN_CELL_LIST(EMP%EMPPOT_RCUT, LATTICE, POSION, NTYP)
         CALL UPDATE_NEIGHBOR_LIST_MODEL_POTENTIAL_FAST(EMP, PAR_SUP_HANDLE, &
              LATTICE, INVERSE_LATTICE, POSION, LNTYP, &
              RCUT1, RCUT2, NEAREST_NEIGHBOR_RADIAL, NEAREST_NEIGHBOR_ANGULAR, W1, W2)

         ! STOP profiling
         PROFILING_STOP('init_emp_pot_neighbor')

      END SUBROUTINE INITIALIZE_NEIGHBOR_LIST_MODEL_POTENTIAL_FAST


      SUBROUTINE MAKE_NEIGHBOR_LIST_MODEL_POTENTIAL_FAST(EMP, FFM, ABN, PAR_SUP_HANDLE, &
                 RMEM_HELP1, &
                 RCUT1, RCUT2, NEAREST_NEIGHBOR_RADIAL, NEAREST_NEIGHBOR_ANGULAR, W1, W2)
         USE VERLET_CELL_LIST_MOD, ONLY : VERLET_CELL_LIST_TYPE
         USE ML_FF_STRUCT, ONLY : EMP_POT_TYPE, ABINITIONEW, PARALLEL_SUPER, FFM_PAR
         USE NEIGHBOR_HELP
         USE tutor, ONLY: vtutor
         IMPLICIT NONE
         TYPE (EMP_POT_TYPE), INTENT(INOUT) :: EMP
         TYPE (FFM_PAR), INTENT(INOUT)      :: FFM
         TYPE (ABINITIONEW), INTENT(IN) :: ABN
         TYPE (PARALLEL_SUPER), INTENT(IN) :: PAR_SUP_HANDLE
         REAL(q), INTENT(INOUT) :: RMEM_HELP1
         Real(q), INTENT(IN), OPTIONAL :: RCUT1
         Real(q), INTENT(IN), OPTIONAL :: RCUT2
         TYPE (VERLET_CELL_LIST_TYPE), INTENT(IN), OPTIONAL :: NEAREST_NEIGHBOR_RADIAL
         TYPE (VERLET_CELL_LIST_TYPE), INTENT(IN), OPTIONAL :: NEAREST_NEIGHBOR_ANGULAR
         REAL(q), INTENT(IN), OPTIONAL :: W1
         REAL(q), INTENT(IN), OPTIONAL :: W2
         ! Local variables
         INTEGER :: MROW_EST_MODEL, MCOL_EST_MODEL
         INTEGER :: INIONS_HELP, INIONS, LOCVAR_COL
 
         ! Start profiling
         PROFILING_START('emp_pot_neighbor')

         ! Determine the maximum number of neighbors
         EMP%MNEIB_EST_MODEL  =  1
         CALL UPDATE_NEIGHBOR_LIST_MODEL_POTENTIAL_FAST(EMP, PAR_SUP_HANDLE, &
              ABN%A(1:,1:,ABN%NCONF), ABN%B(1:,1:,ABN%NCONF), ABN%POSION(1:,1:,ABN%NCONF), ABN%LNTYP, &
              RCUT1, RCUT2, NEAREST_NEIGHBOR_RADIAL, NEAREST_NEIGHBOR_ANGULAR, W1, W2)
         EMP%MNEIB_EST_MODEL  =  EMP%NEAREST_NEIGHBOR%GET_MAX_NUMBER_NEIGHBORS()

         ! Set sizes for allocation
         IF (FFM%LCOUPLE) THEN
            MROW_EST_MODEL=EMP%MNEIB_EST_MODEL+1
            MCOL_EST_MODEL=CEILING(REAL(ABN%MIONS,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
         ELSE
            MROW_EST_MODEL=1
            MCOL_EST_MODEL=1
         ENDIF

         ! Allocation
         IF (ALLOCATED(EMP%LNIONS_EST_MODEL)) DEALLOCATE(EMP%LNIONS_EST_MODEL)
         ALLOCATE (EMP%LNIONS_EST_MODEL(1:MROW_EST_MODEL,1:ABN%MTYP_ALLOCATE,1:MCOL_EST_MODEL))
         RMEM_HELP1=RMEM_HELP1+REAL(MROW_EST_MODEL,q)*REAL(ABN%MTYP_ALLOCATE,q)*REAL(MCOL_EST_MODEL,q)/1.0E+06_q
 
         IF (ALLOCATED(EMP%NNEIB_EST_MODEL)) DEALLOCATE(EMP%NNEIB_EST_MODEL)
         ALLOCATE (EMP%NNEIB_EST_MODEL(1:ABN%MTYP_ALLOCATE,1:MCOL_EST_MODEL))
         RMEM_HELP1=RMEM_HELP1+REAL(ABN%MTYP_ALLOCATE,q)*REAL(MCOL_EST_MODEL,q)/1.0E+06_q

         IF (ALLOCATED(EMP%XYZ_EST_MODEL)) DEALLOCATE(EMP%XYZ_EST_MODEL)
         ALLOCATE (EMP%XYZ_EST_MODEL(1:MROW_EST_MODEL,1:3,1:ABN%MTYP_ALLOCATE,1:MCOL_EST_MODEL))
         RMEM_HELP1=RMEM_HELP1+3*REAL(MROW_EST_MODEL,q)*REAL(ABN%MTYP_ALLOCATE,q)*REAL(MCOL_EST_MODEL,q)*8/1.0E+06_q
 
         IF (ALLOCATED(EMP%XYZ_EST_MODEL_NORM)) DEALLOCATE(EMP%XYZ_EST_MODEL_NORM)
         ALLOCATE (EMP%XYZ_EST_MODEL_NORM(1:MROW_EST_MODEL,1:3,1:ABN%MTYP_ALLOCATE,1:MCOL_EST_MODEL))
         RMEM_HELP1=RMEM_HELP1+3*REAL(MROW_EST_MODEL,q)*REAL(ABN%MTYP_ALLOCATE,q)*REAL(MCOL_EST_MODEL,q)*8/1.0E+06_q
 
         IF (ALLOCATED(EMP%NORMALIZATION_XYZ_MODEL)) DEALLOCATE(EMP%NORMALIZATION_XYZ_MODEL)
         ALLOCATE (EMP%NORMALIZATION_XYZ_MODEL(1:MROW_EST_MODEL,1:ABN%MTYP_ALLOCATE,1:MCOL_EST_MODEL))
         RMEM_HELP1=RMEM_HELP1+3*REAL(MROW_EST_MODEL,q)*8/1.0E+06_q

         ! Calculate neighbor list
         DO INIONS_HELP = 1, ABN%NIONS, PAR_SUP_HANDLE%COMM_WORLD%NCPU
            INIONS = INIONS_HELP + PAR_SUP_HANDLE%COMM_WORLD%NODE_ME
            IF (INIONS.LE.ABN%NIONS) THEN
               LOCVAR_COL=CEILING(REAL(INIONS,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
               CALL UPDATE_SINGLE_NEAREST_NEIGHBOR_ARRAY(EMP%NEAREST_NEIGHBOR%NN_LIST(INIONS), &     
                                                         EMP%LNIONS_EST_MODEL(:, :, LOCVAR_COL), &
                                                         EMP%NNEIB_EST_MODEL(:, LOCVAR_COL), &
                                                         EMP%NORMALIZATION_XYZ_MODEL(:, :, LOCVAR_COL), &
                                                         EMP%XYZ_EST_MODEL(:, :, :, LOCVAR_COL), &
                                                         EMP%XYZ_EST_MODEL_NORM(:, :, :, LOCVAR_COL), &
                                                         ABN%LNTYP(INIONS), &
                                                         INIONS)
            ENDIF
         ENDDO

         ! Stop profiling
         PROFILING_STOP('emp_pot_neighbor')
      END SUBROUTINE MAKE_NEIGHBOR_LIST_MODEL_POTENTIAL_FAST


      SUBROUTINE UPDATE_NEIGHBOR_LIST_MODEL_POTENTIAL_FAST(EMP, PAR_SUP_HANDLE, &
                 LATTICE, INVERSE_LATTICE, POSION, LNTYP, &
                 RCUT1, RCUT2, NEAREST_NEIGHBOR_RADIAL, NEAREST_NEIGHBOR_ANGULAR, W1, W2)
         USE ML_FF_STRUCT, ONLY : EMP_POT_TYPE, PARALLEL_SUPER, VERLET_CELL_LIST_TYPE
         USE tutor, ONLY: vtutor
         IMPLICIT NONE
         TYPE (EMP_POT_TYPE), INTENT(INOUT) :: EMP
         TYPE (PARALLEL_SUPER), INTENT(IN) :: PAR_SUP_HANDLE
         REAL(q), INTENT(IN)  :: LATTICE(:, :)
         REAL(q), INTENT(IN)  :: INVERSE_LATTICE(:, :)
         REAL(q), INTENT(IN)  :: POSION(:, :)
         INTEGER, INTENT(IN)  :: LNTYP(:)
         Real(q), INTENT(IN), OPTIONAL :: RCUT1
         Real(q), INTENT(IN), OPTIONAL :: RCUT2
         TYPE (VERLET_CELL_LIST_TYPE), INTENT(IN), OPTIONAL :: NEAREST_NEIGHBOR_RADIAL
         TYPE (VERLET_CELL_LIST_TYPE), INTENT(IN), OPTIONAL :: NEAREST_NEIGHBOR_ANGULAR
         REAL(q), INTENT(IN), OPTIONAL :: W1
         REAL(q), INTENT(IN), OPTIONAL :: W2

         IF (PRESENT(RCUT1) .AND. PRESENT(RCUT2)) THEN
            IF ( RCUT1.GE.RCUT2 .AND. RCUT1.GE.EMP%EMPPOT_RCUT .AND. W1.GT.0.0_q ) THEN
               CALL EMP%NEAREST_NEIGHBOR%compute_nn_list(NEAREST_NEIGHBOR_RADIAL,&
                                                         PAR_SUP_HANDLE%COMM_WORLD%NODE_ME, &
                                                         PAR_SUP_HANDLE%COMM_WORLD%NCPU)
            ELSE IF ( RCUT2.GE.RCUT1 .AND. RCUT2.GE.EMP%EMPPOT_RCUT .AND. W2.GE.0.0_q ) THEN
               CALL EMP%NEAREST_NEIGHBOR%compute_nn_list(NEAREST_NEIGHBOR_ANGULAR,&
                                                         PAR_SUP_HANDLE%COMM_WORLD%NODE_ME, &
                                                         PAR_SUP_HANDLE%COMM_WORLD%NCPU)
            ELSE IF ( EMP%EMPPOT_RCUT.GE.RCUT1 .AND. EMP%EMPPOT_RCUT.GE.RCUT2 ) THEN
               CALL EMP%NEAREST_NEIGHBOR%forced_update(POSION, LNTYP, &
                                                       LATTICE, INVERSE_LATTICE,&
                                                       PAR_SUP_HANDLE%COMM_WORLD%NODE_ME, &
                                                       PAR_SUP_HANDLE%COMM_WORLD%NCPU)
            ELSE
               CALL vtutor%error("UPDATE_NEIGHBOR_LIST_MODEL_POTENTIAL: Combination of &
                            &EMPPOT_RCUT, RCUT1 and RCUT2 not valid.")
            END IF
         ELSE
            CALL EMP%NEAREST_NEIGHBOR%forced_update(POSION, LNTYP, &
                                                    LATTICE, INVERSE_LATTICE,&
                                                    PAR_SUP_HANDLE%COMM_WORLD%NODE_ME, &
                                                    PAR_SUP_HANDLE%COMM_WORLD%NCPU)
         ENDIF
      END SUBROUTINE UPDATE_NEIGHBOR_LIST_MODEL_POTENTIAL_FAST


!------------------------------------------------------------------------------------------------------------------
! Model potential calculations.
! Model potential is added to energy, forces and stress within this subroutine.
! Model potential used as an intermediate surrogate modelis calculated.
! Model potential comprised of 1)  soft repulsive potential, 2) Morse potential.
!------------------------------------------------------------------------------------------------------------------

      SUBROUTINE MODEL_POTENTIAL_FOR_TI_FAST (EMP, ABN, PAR_SUP_HANDLE, &
                 FFM, TOTEN_MB, TIFOR_MB, TSIF_MB)
        USE BASIS_SOAP
        USE ML_FF_CONSTANT
        USE ML_FF_STRUCT, ONLY: EMP_POT_TYPE, ABINITIONEW, PARALLEL_SUPER, ML_SHMEM, FF_PAR, FFM_PAR
        IMPLICIT NONE
        ! Input and output variables
        TYPE (EMP_POT_TYPE)     :: EMP
        TYPE (ABINITIONEW)      :: ABN
        TYPE (PARALLEL_SUPER)   :: PAR_SUP_HANDLE
        TYPE (FFM_PAR)          :: FFM
        REAL(q), INTENT(INOUT)  :: TOTEN_MB
        REAL(q), INTENT(INOUT)  :: TIFOR_MB(:,:) !(1:3,1:MIONS)
        REAL(q), INTENT(INOUT)  :: TSIF_MB(:,:) !(1:3,1:3)
        ! Local variables
        REAL(q)                 :: TOTEN_HELP
        REAL(q)                 :: TIFOR_HELP(1:3, 1:ABN%NIONS)
        REAL(q)                 :: TSIF_HELP(1:3, 1:3)
        INTEGER                 :: ATOM_SWITCH(1:ABN%NIONS, 1:ABN%NIONS)
        INTEGER                 :: ATOM_PAIRMAP(1:ABN%NIONS, 1:ABN%NIONS)

        ! Start profiling
        PROFILING_START('emp_pot')

        ! Make list of interaction of atoms for empirical potentials
        CALL MAKE_MAP_FOR_INTERACTIONS(EMP, FFM, ATOM_SWITCH, ATOM_PAIRMAP)

        ! Add repulsive potential
        CALL REPULSIVE_EXP_POT(EMP, PAR_SUP_HANDLE, ABN, FFM, &
                 ATOM_SWITCH, TOTEN_HELP, TIFOR_HELP, TSIF_HELP)

        ! Update of energy, forces and stress
        CALL UPDATE_ENERGY_FORCES_STRESS(ABN%NIONS, &
                        TOTEN_HELP, TIFOR_HELP, TSIF_HELP, &
                        TOTEN_MB, TIFOR_MB, TSIF_MB)

        ! Add pairwise morse potential
        CALL MORSE_POTENTIAL(EMP, PAR_SUP_HANDLE, ABN, FFM, &
                 EMP%XYZ_EST_MODEL, EMP%NNEIB_EST_MODEL, EMP%LNIONS_EST_MODEL, &
                 ATOM_SWITCH, ATOM_PAIRMAP, &
                 TOTEN_HELP, TIFOR_HELP, TSIF_HELP)

        ! Update of energy, forces and stress
        CALL UPDATE_ENERGY_FORCES_STRESS(ABN%NIONS, &
                        TOTEN_HELP, TIFOR_HELP, TSIF_HELP, &
                        TOTEN_MB, TIFOR_MB, TSIF_MB)

        ! Stop profiling
        PROFILING_STOP('emp_pot')
      END SUBROUTINE MODEL_POTENTIAL_FOR_TI_FAST


      SUBROUTINE MAKE_MAP_FOR_INTERACTIONS(EMP, FFM, ATOM_SWITCH, ATOM_PAIRMAP)
         USE ML_FF_STRUCT, ONLY: EMP_POT_TYPE, FFM_PAR
         IMPLICIT NONE
         TYPE (EMP_POT_TYPE)    :: EMP
         TYPE (FFM_PAR)         :: FFM
         INTEGER, INTENT(OUT)   :: ATOM_SWITCH(:,:)
         INTEGER, INTENT(OUT)   :: ATOM_PAIRMAP(:,:)
         ! Local variables
         INTEGER :: IPAIR_HELP

         ATOM_SWITCH = 0

         DO IPAIR_HELP = 1, EMP%MOPOT_NM
            ATOM_SWITCH(EMP%IM_COUPLE(IPAIR_HELP), EMP%JM_COUPLE(IPAIR_HELP)) = 1
            ATOM_PAIRMAP(EMP%IM_COUPLE(IPAIR_HELP), EMP%JM_COUPLE(IPAIR_HELP)) = IPAIR_HELP
         ENDDO
      END SUBROUTINE MAKE_MAP_FOR_INTERACTIONS


      SUBROUTINE PREPARE_RN0(EMPPOT_RCUT, XYZ_EST_MODEL, INEIB, RN0, &
                            FCUT1, DFCUT1)
         USE BASIS_SOAP, ONLY: RJ1
         IMPLICIT NONE 
         REAL(q), INTENT(IN)    :: EMPPOT_RCUT
         REAL(q), INTENT(IN)    :: XYZ_EST_MODEL(:,:)
         INTEGER, INTENT(IN)    :: INEIB
         REAL(q), INTENT(OUT)   :: RN0
         REAL(q), INTENT(OUT)   :: FCUT1
         REAL(q), INTENT(OUT)   :: DFCUT1
         ! Local variables
         INTEGER     :: IXYZ
         REAL(q)     :: FCUT_DUMMY(1:1)
         REAL(q)     :: DFCUT_DUMMY(1:1)
         REAL(q)     :: RN0_DUMMY(1:1)
         RN0 = 0.0_q
         DO IXYZ = 1, 3
            RN0 = RN0+XYZ_EST_MODEL(INEIB,IXYZ)**2
         ENDDO
         RN0 = DSQRT(RN0)
         RN0_DUMMY(1) = RN0
         CALL RJ1 (FCUT_DUMMY(1:), DFCUT_DUMMY(1:), 1, EMPPOT_RCUT, RN0_DUMMY(1:))
         FCUT1 = FCUT_DUMMY(1) 
         DFCUT1 = DFCUT_DUMMY(1)
      END SUBROUTINE PREPARE_RN0


      SUBROUTINE REPULSIVE_EXP_POT(EMP, PAR_SUP_HANDLE, ABN, FFM, &
                 ATOM_SWITCH, TOTEN_HELP, TIFOR_HELP, TSIF_HELP)
         USE ML_FF_STRUCT, ONLY: EMP_POT_TYPE, PARALLEL_SUPER, ABINITIONEW, FFM_PAR
         IMPLICIT NONE
         TYPE (EMP_POT_TYPE)    :: EMP
         TYPE (PARALLEL_SUPER)  :: PAR_SUP_HANDLE
         TYPE (ABINITIONEW)     :: ABN
         TYPE (FFM_PAR)         :: FFM
         INTEGER, INTENT(IN)    :: ATOM_SWITCH(:,:)
         REAL(q), INTENT(OUT)   :: TOTEN_HELP
         REAL(q), INTENT(OUT)   :: TIFOR_HELP(1:3, 1:ABN%NIONS)
         REAL(q), INTENT(OUT)   :: TSIF_HELP(1:3, 1:3)
         ! Local variables
         REAL(q)       :: RN_HELP, EXP_HELP, RN0
         INTEGER       :: INIONS, INIONS_HELP, LOCVAR_COL
         INTEGER       :: JNTYP, JNTYP0, INEIB, JNIONS
         INTEGER       :: IXYZ, JXYZ, IERR
         REAL(q)       :: FCUT1, DFCUT1, PREFACTOR, HELP_PRODUCT

         TOTEN_HELP = 0.0_q
         TIFOR_HELP = 0.0_q
         TSIF_HELP = 0.0_q
          
         DO INIONS_HELP = 1, ABN%NIONS, PAR_SUP_HANDLE%COMM_WORLD%NCPU
            INIONS = INIONS_HELP + PAR_SUP_HANDLE%COMM_WORLD%NODE_ME
            IF (INIONS.LE.ABN%NIONS) THEN
               LOCVAR_COL = CEILING(REAL(INIONS,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
               DO JNTYP = 1, ABN%NTYP
                  JNTYP0 = ABN%LADD(JNTYP)
                  DO INEIB = 2, EMP%NNEIB_EST_MODEL(JNTYP0,LOCVAR_COL)+1
                     JNIONS = EMP%LNIONS_EST_MODEL(INEIB,JNTYP0,LOCVAR_COL)
                     IF (JNIONS.LE.ABN%NIONS) THEN
                        IF ((FFM%LIONS_COUPLE(INIONS)) .AND. &
                           (.NOT.FFM%LIONS_COUPLE(JNIONS)) .AND. ATOM_SWITCH(INIONS,JNIONS).NE.1) THEN
                        ELSE
                           CYCLE
                        ENDIF
                        CALL PREPARE_RN0(EMP%EMPPOT_RCUT, EMP%XYZ_EST_MODEL(:,:,JNTYP0,LOCVAR_COL), &
                             INEIB, RN0, FCUT1, DFCUT1)
                        RN_HELP = (RN0/EMP%SRPOT_S0)**EMP%SRPOT_N0
                        EXP_HELP = DEXP(-RN_HELP)
                        ! Compute total energy.
                        TOTEN_HELP =  TOTEN_HELP + FCUT1 * EMP%SRPOT_B0 * EXP_HELP
                        ! Compute forces and stress.
                        PREFACTOR = -EMP%SRPOT_B0 * (-FCUT1 * EMP%SRPOT_N0 * &
                                    RN_HELP/RN0 * EXP_HELP + DFCUT1 * EXP_HELP) / RN0
                        DO IXYZ=1, 3
                           ! Forces
                           HELP_PRODUCT = PREFACTOR * EMP%XYZ_EST_MODEL(INEIB,IXYZ,JNTYP0,LOCVAR_COL)
                           TIFOR_HELP(IXYZ,INIONS) = TIFOR_HELP(IXYZ,INIONS) - HELP_PRODUCT
                           TIFOR_HELP(IXYZ,JNIONS) = TIFOR_HELP(IXYZ,JNIONS) + HELP_PRODUCT
                           ! Stress
                           DO JXYZ = IXYZ, 3
                              TSIF_HELP(JXYZ,IXYZ) = TSIF_HELP(JXYZ,IXYZ) + EMP%XYZ_EST_MODEL(INEIB,JXYZ,JNTYP0,LOCVAR_COL) * &
                                                     HELP_PRODUCT / ABN%OMEGA(ABN%NCONF)
                           ENDDO
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF
         ENDDO

         CALL MPI_ALLREDUCE(MPI_IN_PLACE, TOTEN_HELP, 1, MPI_REAL8, MPI_SUM, PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM, IERR)
         CALL MPI_ALLREDUCE(MPI_IN_PLACE, TIFOR_HELP, 3*ABN%NIONS, MPI_REAL8, MPI_SUM, PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM, IERR)
         CALL MPI_ALLREDUCE(MPI_IN_PLACE, TSIF_HELP, 9, MPI_REAL8, MPI_SUM, PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM, IERR)
      END SUBROUTINE REPULSIVE_EXP_POT


      SUBROUTINE MORSE_POTENTIAL(EMP, PAR_SUP_HANDLE, ABN, FFM, &
                 XYZ_EST_MODEL, NNEIB_EST_MODEL, LNIONS_EST_MODEL, &
                 ATOM_SWITCH, ATOM_PAIRMAP, &
                 TOTEN_HELP, TIFOR_HELP, TSIF_HELP)
         USE ML_FF_STRUCT, ONLY: EMP_POT_TYPE, PARALLEL_SUPER, ABINITIONEW, FFM_PAR
         IMPLICIT NONE
         TYPE (EMP_POT_TYPE)    :: EMP
         TYPE (PARALLEL_SUPER)  :: PAR_SUP_HANDLE
         TYPE (ABINITIONEW)     :: ABN
         TYPE (FFM_PAR)         :: FFM
         REAL(q), INTENT(IN)    :: XYZ_EST_MODEL(:,:,:,:)
         INTEGER, INTENT(IN)    :: NNEIB_EST_MODEL(:,:)
         INTEGER, INTENT(IN)    :: LNIONS_EST_MODEL(:,:,:)
         INTEGER, INTENT(IN)    :: ATOM_SWITCH(:,:)
         INTEGER, INTENT(IN)    :: ATOM_PAIRMAP(:,:)
         REAL(q), INTENT(OUT)   :: TOTEN_HELP
         REAL(q), INTENT(OUT)   :: TIFOR_HELP(1:3, 1:ABN%NIONS)
         REAL(q), INTENT(OUT)   :: TSIF_HELP(1:3, 1:3)
         ! Local variables
         REAL(q)       :: RN_HELP, EXP_HELP, RN0
         INTEGER       :: INIONS, INIONS_HELP, LOCVAR_COL
         INTEGER       :: JNTYP, JNTYP0, INEIB, JNIONS
         INTEGER       :: IXYZ, JXYZ, IERR, IPAIR
         REAL(q)       :: FCUT1, DFCUT1, PREFACTOR, HELP_PRODUCT

         TOTEN_HELP = 0.0_q
         TIFOR_HELP = 0.0_q
         TSIF_HELP = 0.0_q
          
         DO INIONS_HELP = 1, ABN%NIONS, PAR_SUP_HANDLE%COMM_WORLD%NCPU
            INIONS = INIONS_HELP+PAR_SUP_HANDLE%COMM_WORLD%NODE_ME
            IF (INIONS.LE.ABN%NIONS) THEN
               LOCVAR_COL = CEILING(REAL(INIONS,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
               DO JNTYP = 1, ABN%NTYP
                  JNTYP0 = ABN%LADD(JNTYP)
                  DO INEIB = 2, EMP%NNEIB_EST_MODEL(JNTYP0,LOCVAR_COL)+1
                     JNIONS = EMP%LNIONS_EST_MODEL(INEIB,JNTYP0,LOCVAR_COL)
                     IF (JNIONS.LE.ABN%NIONS) THEN
                        IF (ATOM_SWITCH(INIONS,JNIONS) .EQ. 1) THEN
                            IPAIR = ATOM_PAIRMAP(INIONS, JNIONS)
                            CALL PREPARE_RN0(EMP%EMPPOT_RCUT, EMP%XYZ_EST_MODEL(:,:,JNTYP0,LOCVAR_COL), &
                                 INEIB, RN0, FCUT1, DFCUT1)
                            EXP_HELP = DEXP(-EMP%MOPOT_RKM(IPAIR) * (RN0-EMP%MOPOT_RM(IPAIR)))
                            ! Compute total energy.
                            TOTEN_HELP = TOTEN_HELP + FCUT1 * EMP%MOPOT_DM(IPAIR) * ((1.0_q-EXP_HELP)**2-1.0_q)
                            ! Compute forces and stress.
                            PREFACTOR = EMP%MOPOT_DM(IPAIR) / RN0 * (-DFCUT1*((1.0_q-EXP_HELP)**2-1.0_q) - &
                                           2.0_q * EMP%MOPOT_RKM(IPAIR) * FCUT1 * (1.0_q-EXP_HELP)*EXP_HELP)
                            DO IXYZ=1, 3
                               HELP_PRODUCT = PREFACTOR * EMP%XYZ_EST_MODEL(INEIB,IXYZ,JNTYP0,LOCVAR_COL)
                               TIFOR_HELP(IXYZ, INIONS) = TIFOR_HELP(IXYZ, INIONS) - HELP_PRODUCT
                               TIFOR_HELP(IXYZ, JNIONS) = TIFOR_HELP(IXYZ, JNIONS) + HELP_PRODUCT
                               DO JXYZ=IXYZ, 3
                                  TSIF_HELP(JXYZ, IXYZ) = TSIF_HELP(JXYZ, IXYZ) + HELP_PRODUCT * &
                                                          EMP%XYZ_EST_MODEL(INEIB,JXYZ,JNTYP0,LOCVAR_COL) / &
                                                          ABN%OMEGA(ABN%NCONF)
                               ENDDO
                            ENDDO
                        ENDIF
                     ENDIF 
                  ENDDO
               ENDDO
            ENDIF
         ENDDO

         CALL MPI_ALLREDUCE(MPI_IN_PLACE, TOTEN_HELP, 1, MPI_REAL8, MPI_SUM, PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM, IERR)
         CALL MPI_ALLREDUCE(MPI_IN_PLACE, TIFOR_HELP, 3*ABN%NIONS, MPI_REAL8, MPI_SUM, PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM, IERR)
         CALL MPI_ALLREDUCE(MPI_IN_PLACE, TSIF_HELP, 9, MPI_REAL8, MPI_SUM, PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM, IERR)
      END SUBROUTINE MORSE_POTENTIAL


      SUBROUTINE UPDATE_ENERGY_FORCES_STRESS(NIONS, &
                        TOTEN_IN, TIFOR_IN, TSIF_IN, &
                        TOTEN_UPDATE, TIFOR_UPDATE, TSIF_UPDATE)
        IMPLICIT NONE
        INTEGER, INTENT(IN)   :: NIONS
        REAL(q), INTENT(IN)   :: TOTEN_IN
        REAL(q), INTENT(IN)   :: TIFOR_IN(:,:)
        REAL(q), INTENT(IN)   :: TSIF_IN(:,:)
        REAL(q), INTENT(INOUT) :: TOTEN_UPDATE
        REAL(q), INTENT(INOUT) :: TIFOR_UPDATE(:,:)
        REAL(q), INTENT(INOUT) :: TSIF_UPDATE(:,:)
        ! Local variables
        INTEGER :: IXYZ, JXYZ, INIONS

        TOTEN_UPDATE = TOTEN_UPDATE + TOTEN_IN
        DO INIONS = 1, NIONS
           DO IXYZ=1, 3
              TIFOR_UPDATE(IXYZ, INIONS) = TIFOR_UPDATE(IXYZ, INIONS) + TIFOR_IN(IXYZ, INIONS)
           ENDDO
        ENDDO
        DO IXYZ=1, 3
           DO JXYZ=IXYZ, 3
              TSIF_UPDATE(JXYZ, IXYZ) = TSIF_UPDATE(JXYZ, IXYZ) + TSIF_IN(JXYZ, IXYZ)
              ! Symmetrization of stress tensor
              TSIF_UPDATE(IXYZ, JXYZ) = TSIF_UPDATE(IXYZ, JXYZ) + TSIF_IN(JXYZ, IXYZ)
           ENDDO
        ENDDO
      END SUBROUTINE UPDATE_ENERGY_FORCES_STRESS

      SUBROUTINE DEALLOCATE_EMP_POT(EMP)
         USE ML_FF_STRUCT, ONLY : PARALLEL_SUPER, EMP_POT_TYPE
         IMPLICIT NONE
         TYPE (EMP_POT_TYPE) :: EMP

         IF (ALLOCATED(EMP%LNIONS_EST_MODEL)) DEALLOCATE(EMP%LNIONS_EST_MODEL)
         IF (ALLOCATED(EMP%NNEIB_EST_MODEL)) DEALLOCATE(EMP%NNEIB_EST_MODEL)
         IF (ALLOCATED(EMP%XYZ_EST_MODEL)) DEALLOCATE(EMP%XYZ_EST_MODEL)
         IF (ALLOCATED(EMP%XYZ_EST_MODEL_NORM)) DEALLOCATE(EMP%XYZ_EST_MODEL_NORM)
         IF (ALLOCATED(EMP%NORMALIZATION_XYZ_MODEL)) DEALLOCATE(EMP%NORMALIZATION_XYZ_MODEL)
         
         IF (ALLOCATED(EMP%IM_COUPLE)) DEALLOCATE(EMP%IM_COUPLE)
         IF (ALLOCATED(EMP%JM_COUPLE)) DEALLOCATE(EMP%JM_COUPLE)
         IF (ALLOCATED(EMP%MOPOT_DM)) DEALLOCATE(EMP%MOPOT_DM)
         IF (ALLOCATED(EMP%MOPOT_RM)) DEALLOCATE(EMP%MOPOT_RM)
         IF (ALLOCATED(EMP%MOPOT_RKM)) DEALLOCATE(EMP%MOPOT_RKM)
         IF (ALLOCATED(EMP%MOPOT_IJM)) DEALLOCATE(EMP%MOPOT_IJM)
      END SUBROUTINE DEALLOCATE_EMP_POT

      END MODULE EMP_POT
#endif
