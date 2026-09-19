#include "symbol.inc"

!********************************************************************************
!> This module contains routines for a mixed stochastic-deterministic compression
!> of virtual bands as described in PRL 132, 086401 (2024)
!********************************************************************************

MODULE msdgw
   USE prec
   USE base
   USE ini
   USE mpimy
   USE main_mpi
   USE fileio
   USE wave_struct_def
   USE chi_glb 
   USE constant 
   USE tutor
   USE string, ONLY: str 
   USE msdgw_struct
   IMPLICIT NONE

   PRIVATE
   PUBLIC :: MSD_COMPRESSION
   PUBLIC :: MSD_COMPRESSED
   REAL(q), PARAMETER :: EPS = 1.E-6_q
   LOGICAL :: COMPRESSION_CALLED = .FALSE. 

   TYPE compression_map
      INTEGER :: NB_START = 0  !< last frozen band index
      INTEGER :: NB_FINAL = 0  !< last relevant band index 
      INTEGER :: NB_ORIG  = 0  !< original number of bands that are compressed
      INTEGER :: NB_DEGEN = 0  !< number of degenerate bands 
      INTEGER :: NB_TOT = 0    !< total number of bands after compression
      REAL(q) :: F = 0         !< average energy ratio F 
      TYPE( band_ebin ), ALLOCATABLE :: EBIN(:)
   END TYPE compression_map

   INTEGER :: IOUT=16
   CHARACTER(LEN=6) :: FM1='(32I8)'
   CHARACTER(LEN=11) :: FM2='(I3,32F8.3)'

   CONTAINS

   !> returns current status of compression
   FUNCTION MSD_COMPRESSED()
      LOGICAL :: MSD_COMPRESSED
      MSD_COMPRESSED = COMPRESSION_CALLED
   END FUNCTION MSD_COMPRESSED
 
   !*****************************************************************************
   !> main calling routine of module. compresses virtual (un-occupied) bands by
   !> sub-dividing virtual states lying away from the Fermi energy into a 
   !> smaller number of degenerate bands with equal energy and stochastic eigen
   !> vectors. This routine is called prior to a GW calculation to reduce the 
   !> computational cost of setting up the Green's function and related 
   !> quantities
   !>  
   !> @param[in]  WDES (waveded) wavefunction descriptor
   !> @param[in]  W (wavespin) original wavefunction with eigenvalues and states
   !> @param[in]  LATT_CUR (latt) Bravais lattice description
   !> @param[in]  INFO (info_struct) information about settings of job
   !> @param[in]  T_INFO (type_info) information about atoms in unit cell
   !> @param[in]  IO (in_struct) in- out handle 
   !*****************************************************************************
   SUBROUTINE MSD_COMPRESSION( WDES, W, LATT_CUR, INFO, T_INFO, IO )
      USE wave, ONLY: GET_NBANDS_SERIAL
      USE random_seeded, ONLY: pcg_generator, NEW_PCG_GENERATOR, RANDOM_NUMBER
      TYPE (wavedes), INTENT(INOUT) :: WDES
      TYPE (wavespin), INTENT(INOUT) :: W
      TYPE (latt), INTENT(IN) :: LATT_CUR 
      TYPE (info_struct), INTENT(IN) :: INFO
      TYPE (type_info), INTENT(IN) ::  T_INFO
      TYPE (in_struct), INTENT(IN) :: IO
      !> local variables
      TYPE( compression_map ), ALLOCATABLE :: CMAP(:,:)
      TYPE( compression_map ), ALLOCATABLE :: CMAP_OCC(:,:)
      TYPE(pcg_generator) ::  R !< random number generator 
      INTEGER :: NBANDS_DEFAULT
      INTEGER :: NODE_ME = 1
#ifdef MPI
      NODE_ME = W%WDES%COMM%NODE_ME 
#endif
      !> number of degenerate states per k-point must be set 
      !> or constant energy ratio for each degenerate energy 
      IF ( .NOT. ( MSDGW_NS > 0 .OR. MSDGW_F > EPS ) )  THEN
         RETURN 
      ENDIF

      !> obtain default number bands, if NBANDS has not been set in INCAR
      NBANDS_DEFAULT=GET_NBANDS_SERIAL(WDES%LNONCOLLINEAR,INFO%ISPIN,INFO%NELECT,T_INFO)
      NBANDS_DEFAULT=((NBANDS_DEFAULT+WDES%NB_PAR-1)/WDES%NB_PAR)*WDES%NB_PAR

      !> set number of protected bands to two times the default number of bands 
      IF( MSDGW_NP < 0 ) THEN
         MSDGW_NP = 2 * NBANDS_DEFAULT 
      ENDIF

      !> some sanity checks for selected parameter
      IF( MSDGW_NP > WDES%NB_TOT ) THEN
          CALL VTUTOR%ERROR( "Failed to set default number of protected bands for "//&
             "pseudo-deterministic compression of virtual states, because: "//&
             "MSDGW_NP, NBANDS_DEFAULT, NB_TOT = "//str([ MSDGW_NP, NBANDS_DEFAULT, WDES%NB_TOT])//&
             ". NBANDS is probably too small for this job." )
      ENDIF

      !> allocate compression map handles
      IF( ALLOCATED( CMAP ) )THEN
         CALL VTUTOR%BUG( "compression map handle already allocated", __FILE__, __LINE__)
      ENDIF
      ALLOCATE( CMAP( W%WDES%NKPTS, W%WDES%ISPIN ))


      !> optionally, subdivide energy interval such that energy ratio is constant among bins
      IF( MSDGW_F_OCC > 0) THEN
         ALLOCATE( CMAP_OCC( W%WDES%NKPTS, W%WDES%ISPIN ))
         IF( NODE_ME == 1) THEN
            CALL SET_CONST_F_CMAP_OCC( W, CMAP_OCC, IO )
         ENDIF
      ENDIF
      !> set sensible defaults for MSDGW_NP and MSDGW_NNXI 
      !> based on requested MSDGW_NS or MSDGW_F. 
      !> Latter overwrites MDDGW_NS
      IF( MSDGW_F > EPS ) THEN
         !> subdivide energy interval such that energy ratio is constant among bins
         IF( NODE_ME == 1) THEN
            IF( MSDGW_F_OCC > 0 ) THEN
               CALL SET_CONST_F_CMAP( W, CMAP, IO, CMAP_OCC )
            ELSE
               CALL SET_CONST_F_CMAP( W, CMAP, IO )
            ENDIF
         ENDIF 
      ELSEIF( MSDGW_NS> 0) THEN
         !> set sensible parameter based on selected number of stochastic subspaces 
         CALL PARAM_FROM_NS(WDES, INFO, T_INFO, IO )
         !> check parameter for consistency 
         CALL CHECK_PARAM( WDES, IO ) 
         !> determine compression map using stochastic subdivision of bands 
         !> based on distribution of energies 
         IF( NODE_ME == 1) THEN
            CALL SET_CUMULANT_CMAP( W, CMAP, IO )
         ENDIF 
      ELSE
         CALL VTUTOR%BUG( "I should never go through this line", __FILE__, __LINE__)
      ENDIF

      CALL SYNC_CMAP( CMAP )
      IF( MSDGW_F_OCC > 0 ) CALL SYNC_CMAP( CMAP_OCC )
       
      !> initialize random number generator 
      R = NEW_PCG_GENERATOR( INT( MSDGW_SEED, KIND=qi8) , 1_qi8 )

      !> perform actual compression of WAVECAR here
      !> map bands to bins, also do the same for the plane-wave coefficients 
      IF( MSDGW_F_OCC > 0 ) THEN
         CALL COMPRESS_BANDS( WDES, W, CMAP, LATT_CUR, R, IO, CMAP_OCC )
      ELSE
         CALL COMPRESS_BANDS( WDES, W, CMAP, LATT_CUR, R, IO )
      ENDIF

      !> free compression map 
      CALL FREE_COMPRESSION_MAP( CMAP )
      DEALLOCATE( CMAP )
      IF( MSDGW_F_OCC > 0 ) THEN
         CALL FREE_COMPRESSION_MAP( CMAP_OCC )
         DEALLOCATE( CMAP_OCC )
      ENDIF


      !> set compression 
      COMPRESSION_CALLED = .TRUE. 

      CONTAINS 

      !> deallocate compression handle 
      SUBROUTINE FREE_COMPRESSION_MAP( CMAP )
         TYPE( compression_map ), INTENT(INOUT) :: CMAP(:,:)
         INTEGER :: N, I, ISP, NK 
#ifdef debug_msdgw
         DO ISP = 1, SIZE( CMAP, 2)
            DO NK = 1, SIZE(CMAP, 1)
               DO N = 1, CMAP(NK,ISP)%NB_DEGEN
                  WRITE(500+NK,'(F20.10, 1000I5)') CMAP(NK,ISP)%EBIN(N)%EBAR,(CMAP(NK,ISP)%EBIN(N)%IBAND(I), I=1, CMAP(NK,ISP)%EBIN(N)%NBANDS)
               ENDDO
            ENDDO
         ENDDO
#endif 
         DO ISP = 1, SIZE( CMAP, 2)
            DO NK = 1, SIZE( CMAP, 1)
               DO N = 1, CMAP( NK, ISP)%NB_DEGEN
                  DEALLOCATE( CMAP( NK, ISP )%EBIN( N ) %IBAND )
               ENDDO
               DEALLOCATE( CMAP( NK, ISP )%EBIN )
            ENDDO
         ENDDO
      END  SUBROUTINE FREE_COMPRESSION_MAP

      !> deallocate compression handle 
      SUBROUTINE SYNC_CMAP( CMAP )
         TYPE( compression_map ), INTENT(INOUT) :: CMAP(:,:)
         INTEGER :: N, I, ISP, NK 
         DO ISP = 1, SIZE( CMAP, 2)
            DO NK = 1, SIZE(CMAP, 1)
               !> broadcast number of energy bins 
               CALLMPI( M_bcast_i( W%WDES%COMM, CMAP( NK, ISP)%NB_START, 1))
               CALLMPI( M_bcast_i( W%WDES%COMM, CMAP( NK, ISP)%NB_FINAL, 1))
               CALLMPI( M_bcast_i( W%WDES%COMM, CMAP( NK, ISP)%NB_ORIG, 1))
               CALLMPI( M_bcast_i( W%WDES%COMM, CMAP( NK, ISP)%NB_DEGEN, 1))
               CALLMPI( M_bcast_i( W%WDES%COMM, CMAP( NK, ISP)%NB_TOT, 1))
               CALLMPI( M_bcast_d( W%WDES%COMM, CMAP( NK, ISP)%F, 1))

               !> allocate bin arrays, if they haven't been allocated yet
               IF( .NOT. ALLOCATED( CMAP( NK, ISP)%EBIN ) ) THEN
                  ALLOCATE( CMAP( NK, ISP )%EBIN(  CMAP( NK, ISP)%NB_DEGEN ) )
               ENDIF

               !> now broadcast the information about each bin from root rank 
               DO I = 1, CMAP( NK, ISP)%NB_DEGEN
                  !> like the number of bands that fall into this bin
                  CALLMPI( M_bcast_i( W%WDES%COMM, CMAP( NK, ISP)%EBIN(I)%NBANDS, 1))
                  !> energy of the bin 
                  CALLMPI( M_bcast_d( W%WDES%COMM, CMAP( NK, ISP)%EBIN(I)%EBAR, 1))
                  !> energy window of the bin 
                  CALLMPI( M_bcast_d( W%WDES%COMM, CMAP( NK, ISP)%EBIN(I)%F, 1))
                  !> number of pseudo bands for this bin (can be smaller than requested one)
                  CALLMPI( M_bcast_i( W%WDES%COMM, CMAP( NK, ISP)%EBIN(I)%NXI, 1))

                  !> also allocate array that stores which band goes into this bin 
                  IF( .NOT. ALLOCATED(CMAP( NK, ISP)%EBIN(I)%IBAND)) THEN
                      ALLOCATE(CMAP( NK, ISP)%EBIN(I)%IBAND(CMAP( NK, ISP)%EBIN(I)%NBANDS ) )
                  ENDIF
                  !> broadcast these band indices now 
                  CALLMPI( M_bcast_i( W%WDES%COMM, CMAP( NK, ISP)%EBIN(I)%IBAND,CMAP( NK, ISP)%EBIN(I)%NBANDS))
               ENDDO
            ENDDO
         ENDDO
      END SUBROUTINE SYNC_CMAP


   END SUBROUTINE MSD_COMPRESSION
   
   !> helper routine: sets sensible values for global variables 
   !> based on selected number of stochastic subspaces 
   !> @param[in]  WDES (wavedes) wavefunction descriptor
   !> @param[in]  W (wavespin) original wavefunction with eigenvalues and states
   !> @param[in]  LATT_CUR (latt) Bravais lattice description
   !> @param[in]  INFO (info_struct) information about settings of job
   !> @param[in]  T_INFO (type_info) information about atoms in unit cell
   !> @param[in]  IO (in_struct) in- out handle 
   !> @param[out] MSDGW_NP number of protected bands 
   !> @param[out] MSDGW_NS number of stochastic subspaces
   !> @param[out] MSDGW_NXI number of degenerate states for each subspace 
   SUBROUTINE PARAM_FROM_NS(WDES, INFO, T_INFO, IO)
      TYPE (wavedes), INTENT(IN) :: WDES
      TYPE (info_struct), INTENT(IN) :: INFO
      TYPE (type_info), INTENT(IN) ::  T_INFO
      TYPE (in_struct), INTENT(IN) :: IO
      !local variables
      INTEGER :: ISP
      INTEGER :: NK
      INTEGER :: NB_MAX 
      REAL(q) :: RATIO, AVG

      !> check number of subspaces and compute average compression for each k-point 
      IF( IO%IU0>=0 ) WRITE(*,997)
      IF( IO%IU6>=0 ) WRITE(IO%IU6,997)
      AVG=0._q
      DO ISP = 1, WDES%ISPIN
         DO NK = 1, WDES%NKPTS
            IF( MSDGW_NS*MSDGW_NXI + MSDGW_NP > WDES%NB_TOTK(NK, ISP) ) THEN
                CALL VTUTOR%ERROR( "MSDGW_NS*MSDGW_NXI + MSDGW_NP ("//str( [MSDGW_NS, MSDGW_NXI, MSDGW_NP])//&
                   ") is larger than total number of bands NBANDS ("//str(WDES%NB_TOTK(NK,ISP))//") &
                   &at k-point and spin: "//str([NK, ISP] )//". MSDGW_NS and/or MSDGW_NXI is probably &
                   &too large for this job." )
            ENDIF
            !> compression ratio for this k-point in per cent is 
            RATIO = (1-( MSDGW_NP+MSDGW_NS*MSDGW_NXI )/REAL( WDES%NB_TOTK(NK,ISP), q ) )*100
            !> add it to the average 
            AVG = AVG + RATIO
            IF( IO%IU0>=0 .AND. IO%NWRITE>2) WRITE(  *   ,998) RATIO 
            IF( IO%IU6>=0 .AND. IO%NWRITE>2) WRITE(IO%IU6,998) RATIO 
         ENDDO
      ENDDO
      AVG = AVG / (WDES%ISPIN * WDES%NKPTS )

      !> print selected settings 
      NB_MAX = MAXVAL( WDES%NB_TOTK(1:WDES%NKPTS,:) )
      IF( IO%IU0>=0 ) WRITE(*,999) MSDGW_NP, MSDGW_NS, MSDGW_NXI , MSDGW_NP+MSDGW_NS*MSDGW_NXI, AVG
      IF( IO%IU6>=0 ) WRITE(IO%IU6,999) MSDGW_NP, MSDGW_NS, MSDGW_NXI , MSDGW_NP+MSDGW_NS*MSDGW_NXI, AVG

997   FORMAT( ' Mixed stochastic-deterministic compression of virtual states:')
998   FORMAT( '    compression ratio in %:    ', F6.2 )
999   FORMAT( '    number of protected bands: ', I6 ,/,&
              '    number of degenerate bands:', I6, /,&
              '    degeneracy of compr. bands:', I6, /,&
              '    total number of bands used:', I6, /,&
              '    avg. compression ratio (%):', F6.2 )
   END SUBROUTINE PARAM_FROM_NS

   !> helper routine that determines subdivisions of energy intervall 
   !> in virtual space 
   SUBROUTINE SET_CUMULANT_CMAP_OLD( W, CMAP, IO) 
      USE m_unirnk, ONLY: UNIRNK
      TYPE(wavespin)               :: W
      TYPE (in_struct), INTENT(IN) :: IO
      TYPE (compression_map), INTENT(OUT) :: CMAP(:,:)
      !> local
      INTEGER, PARAMETER :: NBINS_MIN = 1001 
      REAL(q), ALLOCATABLE :: RCOUNT(: ) 
      REAL(q), ALLOCATABLE :: RCOUNT_F(: ) 
      REAL(q), ALLOCATABLE :: F( : ) 
      INTEGER, ALLOCATABLE :: IDX_F(:,:,:)
      INTEGER, ALLOCATABLE :: IDX(:,:,:)
      INTEGER :: NBINS 
      INTEGER :: NB_START, NB_FINAL, NB
      INTEGER :: NK, ISP
      INTEGER :: I,J, N
      REAL(q), ALLOCATABLE :: R(:)
      REAL(q) :: RMIN, RMAX 
      REAL(q) :: RI, RI_
      INTEGER :: NODE_ME = 1
      REAL(q) :: E_BAR( MSDGW_NS )
      INTEGER :: ICOUNT( MSDGW_NS )
      INTEGER :: NS

      !> stores the indices of energies that are in one bin
      INTEGER, ALLOCATABLE :: IDX_IN_BIN(:,:)
     
#ifdef MPI 
      NODE_ME = W%WDES%COMM%NODE_ME 
#endif
      PROFILING_START( 'set_cumulant_cmap_old')

      !> number of subdivisions
      NBINS = MAX( NBINS_MIN, MSDGW_NS+1 )
      ALLOCATE( RCOUNT( NBINS) )
      ALLOCATE( RCOUNT_F( MSDGW_NS ) )
      ALLOCATE( F( NBINS ) )
      !> stores bin index in which energy falls into 
      ALLOCATE( IDX( W%WDES%NB_TOT-MSDGW_NP+1, W%WDES%NKPTS, W%WDES%ISPIN ) )
      ALLOCATE( IDX_F( NBINS , W%WDES%NKPTS, W%WDES%ISPIN ) )
      IDX = 0 
      IDX_F = 0 
      !> setup cummulant of energy interval
      !> do this for each k-point and spin channel seperately 
      DO ISP = 1, W%WDES%ISPIN
      DO NK = 1, W%WDES%NKPTS

         !> count the "relevant" number of elements in array 
         !> if more bands are used than provided by ENCUT,
         !> the redundant bands have a large energy eigenvalue beyond 1e4
         !> do not consider these bands here: 
         NB_FINAL = 0 
         DO I = 1, W%WDES%NB_TOTK( NK, ISP )
            IF( ABS( W%CELTOT( I, NK, ISP) ) > 9.e3 ) CYCLE
            NB_FINAL= NB_FINAL + 1
         ENDDO
         NB_START = MSDGW_NP

         !> this is the total number of energies considered 
         NB = NB_FINAL - MSDGW_NP
      
         !> store index and number of compressed bands for this k-point 
         !> last frozed band index 
         CMAP(NK,ISP)%NB_START= NB_START
         !> last band considered for compression
         CMAP(NK,ISP)%NB_FINAL= NB_FINAL
         !> number of bands considered for compression 
         CMAP(NK,ISP)%NB_ORIG = NB
         !> store these energies in R
         ALLOCATE( R(  NB ))
         DO I = 1, NB
            R( I ) = W%CELTOT( I + NB_START, NK, ISP )
         ENDDO

         !> obtain their minimal and maximal value, easy since array is already ordered
         RMIN = R(1)-EPS
         RMAX = R(NB)+EPS
         !IF( IO%IU0>=0 ) WRITE(*,*) RMIN, RMAX, NB, MSDGW_NP, NB_FINAL 

         !> usually one has NB >> MSDGW_NS  
         !> equidistant subdivision is probably not as good as an non-uniform subdivision
         !> so we subdivide the energy interval between RMIN, RMAX into MSDGW_NS sub intervals
         !> the points are selected based on the cummulant of the distribution of energies in 
         !> [RMIN, RMAX]. For this purpose we have to build this distribution by subdividing 
         !> the intervals into NBINS equi-distant intervals of width (RMAX-RMIN)/NBINS and 
         !> count the number of energies occuring in these bins. This yields the probability 
         !> distribution RCOUNT(I) for each bin. 
         !> 
         RCOUNT = 0 
         F=0
         !> loop over all bins 
         DO I = 1, NBINS-1
            RI = RMIN + (I-1) * ( RMAX - RMIN ) / REAL(NBINS-1,q)
            RI_= RMIN + I * ( RMAX - RMIN ) / REAL(NBINS-1,q)
            !> find energy that falls into this bin 
            DO N = 1, NB
               !> increase counter if energy falls into this bin
               IF ( RI < R( N ) .AND. R( N ) < RI_ ) THEN
                  RCOUNT( I ) = RCOUNT( I ) + (1._q/NB)
                  !> also store index into which this CELTOT falls into
                  IF ( IDX( N, NK, ISP)/=0 .AND. I < NBINS ) THEN
                     CALL VTUTOR%BUG( "double assignment of entry: "//str( [ I, N, IDX(N,NK,ISP) ])//&
                        ","//str([ RI, RI_, R(N) ]), __FILE__, __LINE__ )
                  ENDIF
                  IDX( N, NK, ISP ) = I 
               ENDIF
            ENDDO
            !> set cummulant 
            F(I+1) = F(I) + RCOUNT(I)
         ENDDO
         !> consistency check, RCOUNT_F should add up to number of datapoints
         IF( ABS( SUM( RCOUNT ) - 1._q ) >1.E-6_q )THEN
            CALL VTUTOR%BUG( "inconsistent bin search "//str( [SUM(RCOUNT), 1._q ])//&
               " at k-point: "//str([ NK, ISP]), __FILE__, __LINE__ )
         ENDIF 
         !> also cummulant must add to one
         IF( ABS( F(NBINS) -1 ) > 1.E-6_q )THEN
            CALL VTUTOR%BUG( "cummulant not setup correctly "//str( F(NBINS)), __FILE__, __LINE__ )
         ENDIF 
#ifdef debug_msdgw
         CALL DUMP_RARRAY( 'eps',NK, F  )
         IF ( NODE_ME == 1 ) THEN 
            DO N = 1, NB
               WRITE(100+NK,'(I5, F20.10)') IDX(N,NK,ISP), R(N)
            ENDDO
         ENDIF
#endif

         RCOUNT_F = 0 
         !> now loop over requested intersections
         DO I = 1, MSDGW_NS
            !> sub divide cummulant range into MSDGW_NS bins of same size
            RI = (I-1)/ REAL(MSDGW_NS-1,q)
            RI_= I/ REAL(MSDGW_NS-1,q)
            !> loop over all cummulant bins 
            DO N = 1, NBINS
               !> find cummulant bin that belongs to this sub-interval 
               IF( RI <= F( N ) .AND. F(N) <= RI_ ) THEN
                  !> increase counter and store index 
                  RCOUNT_F( I ) = RCOUNT_F( I ) + 1 
                  !> store cummulant bin index 
                  IDX_F( N, NK, ISP ) = I 
               ENDIF 
            ENDDO
         ENDDO 

#ifdef debug_msdgw
         CALL DUMP_RARRAY( 'p_eps',NK, RCOUNT_F, 0._q, 1._q  )
         IF ( NODE_ME == 1 ) THEN 
            DO N = 1, NBINS
               WRITE(200+NK,'(I5, F20.10)') IDX_F(N,NK,ISP)
            ENDDO
         ENDIF
#endif 
         !> consistency check, RCOUNT_F should add up to number of datapoints
         IF( ABS( SUM( RCOUNT_F ) - NBINS ) >1.E-6_q )THEN

            CALL VTUTOR%BUG( "inconsistent bin search "//str( [SUM(RCOUNT_F), REAL(NBINS,q) ])// &
               " at k-point: "//str([ NK, ISP]), __FILE__, __LINE__ )
         ENDIF 
         !> now build energy averages for each bin
         E_BAR = 0 
         ICOUNT = 0 
         DO I = 1, NBINS  
            DO J = 1, NB  
               !> pick correct energy index for this bin 
               IF ( IDX( J, NK, ISP) == I  ) THEN
                  !> obtain index of requested sub division
                  N = IDX_F( I, NK, ISP)
                  !> build average of all energies in this bin
                  E_BAR(N) = E_BAR( N )  + R( J ) 
                  !> store the number of energies in this bin
                  ICOUNT(N) = ICOUNT( N ) + 1 
               ENDIF 
            ENDDO
         ENDDO

         IF( MAXVAL( ICOUNT)==0 ) THEN
            CALL VTUTOR%BUG( "all bins are empty", __FILE__, __LINE__ )
         ENDIF

         !> now figure out which energies are in each bin 
         !> for this we require an array that can store 
         !> this information. Thus we need to know the bin
         !> with the maximum number of energies 
         ALLOCATE( IDX_IN_BIN( MSDGW_NS,  MAXVAL( ICOUNT ) ) )
         IDX_IN_BIN = 0 
         ICOUNT = 0 
         DO I = 1, NBINS  
            DO J = 1, NB  
               !> pick correct energy index for this bin 
               IF ( IDX( J, NK, ISP) == I  ) THEN
                  !> obtain index of requested sub division
                  N = IDX_F( I, NK, ISP)
                  !> store energy index for this bin 
                  ICOUNT(N) = ICOUNT( N ) + 1 
                  IDX_IN_BIN( N , ICOUNT(N) ) = J 
               ENDIF 
            ENDDO
         ENDDO
         !> it is still possible that some bins are empty 
         !> thus determine the actual number of NS
         NS =  0 
         DO I = 1, MSDGW_NS
            IF( ICOUNT( I) > 0 ) THEN
               NS = NS + 1
               E_BAR(I) = E_BAR(I)/ICOUNT(I)
            ENDIF
         ENDDO

         !> allocate compression handle for this k-point 
         CMAP( NK, ISP )%NB_DEGEN = NS 
         !> number of bands compressed 
         CMAP(NK,ISP)%NB_TOT = CMAP(NK,ISP)%NB_START + NS*MSDGW_NXI
         !> round to proper number of bands when running in parallel 
         CMAP( NK, ISP)%NB_TOT = ((CMAP( NK, ISP)%NB_TOT+ W%WDES%NB_PAR-1)/W%WDES%NB_PAR)* W%WDES%NB_PAR
         ALLOCATE( CMAP( NK, ISP )%EBIN( NS ) )
         !> loop over all bins 
         N=0
         DO I = 1, MSDGW_NS
            !> empty bins are ignored 
            IF( ICOUNT( I) ==0 ) CYCLE 
            !> loop over all energy bins
            N=N+1
            !> allocate energy bin handle 
            CMAP( NK, ISP)%EBIN(N)%NBANDS = ICOUNT(I)
            !> set averaged energy 
            CMAP( NK, ISP)%EBIN(N)%EBAR = E_BAR(I)
            !> set all band indices of original CELTOT array that fall into this bin 
            ALLOCATE( CMAP(NK,ISP)%EBIN(N)%IBAND( ICOUNT(I)))
            CMAP( NK, ISP)%EBIN(N)%IBAND(1:ICOUNT(I))=IDX_IN_BIN(I,1:ICOUNT(I))
         ENDDO

            
         !> determine energy window for this bin (delta E)
         IF ( IO%IU0>=0 .AND. IO%NWRITE>3 ) WRITE(   *  ,'(A,2I4)' )" Energy ratios for k-point: ", NK, ISP
         IF ( IO%IU6>=0 .AND. IO%NWRITE>3 ) WRITE(IO%IU6,'(A,2I4)' )" Energy ratios for k-point: ", NK, ISP
         CMAP( NK, ISP )%F = 0 
         DO I = 1, CMAP( NK, ISP )%NB_DEGEN
            CMAP( NK, ISP)%EBIN(I)%F = & 
               R( CMAP( NK, ISP)%EBIN(I)%IBAND( CMAP( NK, ISP)%EBIN(I)%NBANDS )  ) - &
                 R( CMAP( NK, ISP)%EBIN(I)%IBAND(1) ) 
            CMAP( NK, ISP)%EBIN(I)%F = CMAP( NK, ISP)%EBIN(I)%F / CMAP( NK, ISP)%EBIN(I)%EBAR 
            !> determine energy ratio 
            IF ( IO%IU0>=0 .AND. IO%NWRITE>3 ) WRITE(   *  , 889) I, CMAP( NK, ISP)%EBIN(I)%F   
            IF ( IO%IU6>=0 .AND. IO%NWRITE>3 ) WRITE(IO%IU6, 889) I, CMAP( NK, ISP)%EBIN(I)%F   
            CMAP( NK, ISP )%F = CMAP( NK, ISP )%F + CMAP( NK, ISP)%EBIN(I)%F 
         ENDDO
         CMAP( NK, ISP )%F = CMAP( NK, ISP )%F /  CMAP(NK,ISP)%NB_DEGEN
           
         DEALLOCATE( IDX_IN_BIN )
         DEALLOCATE( R ) 

         !> run unit test 
         CALL UNIT_TEST_CMAP_SETUP( CMAP( NK, ISP), NB )

      ENDDO  ! NK
      ENDDO  ! ISP

      DEALLOCATE( IDX_F )
      DEALLOCATE( IDX )
      DEALLOCATE( F )
      DEALLOCATE( RCOUNT_F )
      DEALLOCATE( RCOUNT )
889   FORMAT( '    energy ratio for bin :    ', I5, F6.2 )


      PROFILING_STOP( 'set_cumulant_cmap_old')

   END SUBROUTINE SET_CUMULANT_CMAP_OLD

   !> helper routine that determines subdivisions of energy intervall 
   !> in virtual space 
   SUBROUTINE SET_CUMULANT_CMAP( W, CMAP, IO) 
      USE m_unirnk, ONLY: UNIRNK
      TYPE(wavespin)               :: W
      TYPE (in_struct), INTENT(IN) :: IO
      TYPE (compression_map), INTENT(OUT) :: CMAP(:,:)
      !> local
      INTEGER, PARAMETER :: NBINS_GRAIN = 1
      INTEGER :: NB_FINAL
      INTEGER :: NK, ISP
      INTEGER :: I, J, K, N,  NB
      INTEGER :: IDX1
      INTEGER :: NS
      INTEGER :: NODE_ME = 1
      REAL(q) :: RI, RI_
      REAL(q), ALLOCATABLE :: R(:)
      INTEGER, ALLOCATABLE :: ICOUNT(:)
      INTEGER, ALLOCATABLE :: IDX(:)
      INTEGER, ALLOCATABLE :: JDX(:)
      REAL(q), ALLOCATABLE :: F(:)
      INTEGER :: NBINS 
#ifdef MPI 
      NODE_ME = W%WDES%COMM%NODE_ME 
#endif
      PROFILING_START( 'set_cumulant_cmap')


      !> setup cummulant of energy interval
      !> do this for each k-point and spin channel seperately 
      DO ISP = 1, W%WDES%ISPIN
      DO NK = 1, W%WDES%NKPTS

         !> count the "relevant" number of elements in array 
         !> if more bands are used than provided by ENCUT,
         !> the redundant bands have a large energy eigenvalue, usually beyond 1E4
         !> do not consider these bands here: 
         NB_FINAL = 0 
         DO I = 1, W%WDES%NB_TOTK( NK, ISP )
            IF( ABS( W%CELTOT( I, NK, ISP) ) > 9.E3_q) CYCLE
            NB_FINAL= NB_FINAL + 1
         ENDDO
         !> this is the total number of energies considered 
         NB = NB_FINAL - MSDGW_NP
      
         !> store index and number of compressed bands for this k-point 
         !> last frozed band index 
         CMAP(NK,ISP)%NB_START= MSDGW_NP
         !> last band considered for compression
         CMAP(NK,ISP)%NB_FINAL= NB_FINAL
         !> number of bands considered for compression 
         CMAP(NK,ISP)%NB_ORIG = NB
         !> store these energies in R
         ALLOCATE( R(  NB ))
         DO I = 1, NB
            R( I ) = W%CELTOT( I + MSDGW_NP, NK, ISP )
         ENDDO

         !> set up cummulant 
         !> sampling bins 
         NBINS =NBINS_GRAIN*MSDGW_NS
         ALLOCATE( ICOUNT( NBINS) )
         ALLOCATE( IDX( NB ) )
         ALLOCATE( F( NBINS ) )
         !
         ICOUNT = 0 
         IDX = 0 
         IDX1 = 1
         !> find out which band falls into which energy bin
         DO N = 1, NB
            !> find energy bin 
            bin: DO I = IDX1, NBINS 
               RI = R(1) + I*(R(NB) - R(1) )/ REAL( NBINS , q)  + EPS
               IF( R(N) < RI ) THEN    
                  !> set bin index for this band 
                  IDX( N ) = I 
                  !> count number of energies in this bin 
                  ICOUNT( I) = ICOUNT( I ) + 1 

                  !> increase bin index for next loop 
                  IDX1 = I
                  EXIT bin
               ENDIF
            ENDDO bin
         ENDDO
         !> sanity check, ICOUNT must sum to number of bands
         IF( SUM( ICOUNT ) /= NB) THEN
            CALL VTUTOR%BUG( "sum of bands in bins inconsistent with number of bands "//&
               str( [SUM(ICOUNT), NB]) , __FILE__, __LINE__ ) 
         ENDIF

         !> set cumulant 
         F=0
         F(1) = ICOUNT(1)
         DO I = 2, NBINS
            F(I) = F(I-1) + ICOUNT(I)
         ENDDO
         !> normalize cummulant 
         F = F / F( NBINS )

         !> at this state, cumulant has been set up, 
         !> also IDX stores which band goes into which cumulant bin
         !> now determine number of non-empty bins
         ALLOCATE( JDX( NB ) )
         CALL UNIRNK( IDX, JDX, NS )

#ifdef debug_msdgw
         WRITE(*, *) "distinct bins:", NS , MSDGW_NP,  NB,  R(1), R(NB)
         CALL DUMP_RARRAY( 'cum', NK, F, 0._q, 1._q  )
         CALL DUMP_RARRAY( 'eps', NK, R )
         CALL DUMP_IARRAY( 'idx', NK, IDX )
         CALL DUMP_IIARRAY( 'ibin', NK, IDX,JDX )
         CALL DUMP_IARRAY( 'bin', NK, ICOUNT )
#endif         
      
         !> allocate bins 
         ALLOCATE( CMAP( NK, ISP)%EBIN( NS ))
         !> allocate compression handle for this k-point 
         CMAP( NK, ISP )%NB_DEGEN = NS 
         !> number of bands compressed 
         CMAP(NK,ISP)%NB_TOT = CMAP(NK,ISP)%NB_START + NS*MSDGW_NXI
         !> round to proper number of bands when running in parallel 
         CMAP( NK, ISP)%NB_TOT = ((CMAP( NK, ISP)%NB_TOT+ W%WDES%NB_PAR-1)/W%WDES%NB_PAR)* W%WDES%NB_PAR
         J=0
         IDX1 = 1
         DO I = 1, NBINS
            !> skip empty bins 
            IF( ICOUNT(I) == 0  ) CYCLE
            J=J+1
            !> set number of bands for this bin 
            CMAP( NK, ISP )%EBIN( J )%NBANDS  = ICOUNT( I )
            CMAP( NK, ISP)%EBIN( J )%EBAR = 0
            CMAP( NK, ISP)%EBIN( J )%F = 0
            ALLOCATE( CMAP( NK, ISP )%EBIN( J )%IBAND( ICOUNT( I ) ) ) 
            !> set band indices that belong into this bin
            K=0
            DO N = IDX1, NB
               IF( IDX( N ) == I ) THEN
                  K=K+1
                  CMAP( NK, ISP )%EBIN( J )%IBAND( K ) = N
                  CMAP( NK, ISP)%EBIN( J )%EBAR = CMAP( NK, ISP)%EBIN( J )%EBAR + R( N )
               ENDIF
            ENDDO
            !> set start index for next loop 
            IDX1 = CMAP( NK, ISP )%EBIN( J )%IBAND( K ) + 1 
         ENDDO

         IF ( J /= NS ) THEN
            CALL VTUTOR%BUG( "inconsistent bins:"//str( [NS, J]), __FILE__, __LINE__ )
         ENDIF

         !> normalization
         CMAP( NK, ISP )%F = 0 
         IF ( IO%IU0>=0 .AND. IO%NWRITE>3 ) WRITE(   *  ,'(A,2I4)' )" Energy ratios for k-point: ", NK, ISP
         IF ( IO%IU6>=0 .AND. IO%NWRITE>3 ) WRITE(IO%IU6,'(A,2I4)' )" Energy ratios for k-point: ", NK, ISP
         DO J = 1, CMAP( NK, ISP )%NB_DEGEN
            CMAP( NK, ISP)%EBIN( J )%EBAR = CMAP( NK, ISP)%EBIN( J )%EBAR/ &
                CMAP( NK, ISP)%EBIN( J )%NBANDS
            !> energy ratio 
            RI = R( CMAP( NK, ISP)%EBIN( J )%IBAND(1) )
            RI_ = R( CMAP( NK, ISP)%EBIN( J )%IBAND(CMAP( NK, ISP)%EBIN( J )%NBANDS ) )
            CMAP( NK, ISP)%EBIN( J )%F = (RI_-RI)/CMAP( NK, ISP)%EBIN( J )%EBAR
            IF ( IO%IU0>=0 .AND. IO%NWRITE>3 ) WRITE(   *  , 889) J, CMAP( NK, ISP)%EBIN(J)%F   
            IF ( IO%IU6>=0 .AND. IO%NWRITE>3 ) WRITE(IO%IU6, 889) J, CMAP( NK, ISP)%EBIN(J)%F   
            CMAP( NK, ISP )%F = CMAP( NK, ISP )%F + CMAP( NK, ISP)%EBIN(J)%F 
         ENDDO
         CMAP( NK, ISP )%F = CMAP( NK, ISP )%F /  CMAP(NK,ISP)%NB_DEGEN

#ifdef debug_msdgw
         CALL DUMP_RARRAY( 'ebar', NK, CMAP(NK,ISP)%EBIN(:)%EBAR )
         CALL DUMP_IARRAY( 'ibands', NK, CMAP(NK,ISP)%EBIN(:)%NBANDS )
#endif
         IF( IO%NWRITE > 2) THEN
            CALL DUMP_IRARRAY( 'hist', NK, CMAP(NK,ISP)%EBIN(:)%NBANDS, CMAP(NK,ISP)%EBIN(:)%EBAR  )
         ENDIF

         !> run unit test 
         CALL UNIT_TEST_CMAP_SETUP( CMAP( NK, ISP), NB )

         DEALLOCATE( JDX )
         DEALLOCATE( F )
         DEALLOCATE( IDX )
         DEALLOCATE( ICOUNT )
         DEALLOCATE( R )

      ENDDO  ! NK
      ENDDO  ! ISP
889   FORMAT( '    energy ratio for bin :    ', I5, F6.2 )

      PROFILING_STOP( 'set_cumulant_cmap')

   END SUBROUTINE SET_CUMULANT_CMAP

   !> helper routine that determines subdivisions of energy intervall 
   !> in virtual space based on constant energy ratio, follows 
   !> BerkeleyGW implementation
   SUBROUTINE SET_CONST_F_CMAP( W, CMAP, IO, CMAP_OCC) 
      TYPE(wavespin)               :: W
      TYPE (compression_map), INTENT(OUT) :: CMAP(:,:)
      TYPE (in_struct), INTENT(IN) :: IO
      TYPE (compression_map), INTENT(IN), OPTIONAL :: CMAP_OCC(:,:)
      !> local
      INTEGER :: NB_FINAL
      INTEGER :: NK, ISP
      INTEGER :: I,J, NB
      INTEGER :: NODE_ME = 1
      INTEGER :: NS
      INTEGER :: NB_OUT
      INTEGER :: NB_TOT
      REAL(q), ALLOCATABLE :: R(:)
#ifdef MPI 
      NODE_ME = W%WDES%COMM%NODE_ME 
#endif
      PROFILING_START( 'set_const_f_cmap')

      !>  loop over each k-point separately 
      DO ISP = 1, W%WDES%ISPIN
      DO NK = 1, W%WDES%NKPTS

         !> count the "relevant" number of elements in array 
         !> if more bands are used than provided by ENCUT,
         !> the redundant bands have a large energy eigenvalue, usually beyond 1E4
         !> do not consider these bands here: 
         NB_FINAL = 0 
         DO I = 1, W%WDES%NB_TOTK( NK, ISP )
            IF( ABS( W%CELTOT( I, NK, ISP) ) > 9.E3_q) CYCLE
            NB_FINAL= NB_FINAL + 1
         ENDDO
         !> this is the total number of energies considered 
         NB = NB_FINAL - MSDGW_NP
      
         !> store index and number of compressed bands for this k-point 
         !> last frozed band index 
         CMAP(NK,ISP)%NB_START= MSDGW_NP
         !> last band considered for compression
         CMAP(NK,ISP)%NB_FINAL= NB_FINAL
         !> number of bands considered for compression 
         CMAP(NK,ISP)%NB_ORIG = NB
         !> store these energies in R
         ALLOCATE( R(  NB ))
         DO I = 1, NB
            R( I ) = W%CELTOT( I + MSDGW_NP, NK, ISP )
         ENDDO

         CMAP(NK,ISP)%F = MSDGW_F 
         !> determine number of bins first
         CALL SET_E_BINS( NB, NS, R, MSDGW_F )
         !WRITE(*,*)' bands in total for this k-point: ', NK, ISP, NB, NB + MSDGW_NP
         
         !> allocate bins 
         ALLOCATE( CMAP( NK, ISP)%EBIN( NS ))

         !> set bins now 
         CALL SET_E_BINS( NB, NS, R, MSDGW_F, CMAP(NK,ISP), MSDGW_NXI )
         !> define number of degenerate bands / MSDGW_NXI 
         CMAP( NK, ISP )%NB_DEGEN = SIZE( CMAP(NK,ISP)%EBIN, 1)
         !> total number of bands, after compression
         NB_TOT = CMAP(NK,ISP)%NB_START !+ CMAP(NK,ISP)%NB_DEGEN* MSDGW_NXI 
         DO I = 1,  CMAP( NK, ISP )%NB_DEGEN
            NB_TOT = NB_TOT + CMAP(NK,ISP)%EBIN(I)%NXI 
         ENDDO

         !> optionally, take also compressed, occupied states into account
         IF( PRESENT( CMAP_OCC )) THEN
            !> this is the number of states, after compression including compressed, occupied states 
            !NB_TOT = MSDGW_NP - CMAP_OCC(NK,ISP)%NB_FINAL + CMAP_OCC( NK, ISP)%NB_TOT & 
            !       + CMAP(NK,ISP)%NB_DEGEN* MSDGW_NXI 
            NB_TOT = MSDGW_NP - CMAP_OCC(NK,ISP)%NB_FINAL + CMAP_OCC( NK, ISP)%NB_TOT 
            DO I = 1, CMAP(NK,ISP)%NB_DEGEN
               NB_TOT = NB_TOT + CMAP(NK,ISP)%EBIN(I)%NXI 
            ENDDO
#ifdef debug_msdgw
            WRITE(*,*) '-------------------------------------------', NK
            WRITE(*,*) 'frozen:        ', MSDGW_NP - CMAP_OCC( NK, ISP)%NB_TOT + 1 
            WRITE(*,*) 'compressed occ:', CMAP_OCC(NK,ISP)%NB_DEGEN* MSDGW_NXI_OCC, CMAP_OCC(NK,ISP)%NB_DEGEN
            WRITE(*,*) 'compressed    :', CMAP(NK,ISP)%NB_DEGEN* MSDGW_NXI, CMAP(NK,ISP)%NB_DEGEN
            WRITE(*,*) 'occ inc. NP:   ', MSDGW_NP - CMAP_OCC(NK,ISP)%NB_FINAL + CMAP_OCC( NK, ISP)%NB_TOT
            WRITE(*,*) 'total:        ', NB_OUT
#endif
         ENDIF
         !> round to proper number of bands when running in parallel 
         NB_OUT = ((NB_TOT+ W%WDES%NB_PAR-1)/W%WDES%NB_PAR)* W%WDES%NB_PAR
         CMAP( NK, ISP)%NB_TOT = NB_OUT 
         DEALLOCATE( R ) 
           
         !> run unit test 
         CALL UNIT_TEST_CMAP_SETUP( CMAP( NK, ISP), NB )

         IF( IO%NWRITE > 2) THEN
            CALL DUMP_IRARRAY( 'hist', NK, CMAP(NK,ISP)%EBIN(:)%NBANDS, CMAP(NK,ISP)%EBIN(:)%EBAR  )
         ENDIF
      ENDDO  ! NK
      ENDDO  ! ISP
      PROFILING_STOP( 'set_const_f_cmap')
   END SUBROUTINE SET_CONST_F_CMAP

   !> helper routine that determines subdivisions of energy intervall 
   !> in virtual space based on constant energy ratio, follows 
   !> BerkeleyGW implementation
   SUBROUTINE SET_CONST_F_CMAP_OCC( W, CMAP, IO) 
      TYPE(wavespin)               :: W
      TYPE (in_struct), INTENT(IN) :: IO
      TYPE (compression_map), INTENT(OUT) :: CMAP(:,:)
      !> local
      INTEGER :: NB_FINAL
      INTEGER :: NK, ISP
      INTEGER :: I,J, NB
      INTEGER :: NODE_ME = 1
      INTEGER :: NS
      INTEGER :: NB_OUT
      INTEGER :: NB_TOT
      REAL(q), ALLOCATABLE :: R(:)
#ifdef MPI 
      NODE_ME = W%WDES%COMM%NODE_ME 
#endif
      PROFILING_START( 'set_const_f_cmap_occ')

      !>  loop over each k-point separately 
      DO ISP = 1, W%WDES%ISPIN
      DO NK = 1, W%WDES%NKPTS

         !> count the "relevant" number of elements in array 
         NB_FINAL = 0 
         DO I = 1, W%WDES%NB_TOTK( NK, ISP )
            !> exit on first unoccupied band 
            IF( ABS( W%FERTOT( I, NK, ISP) ) < 1.E-8_q) EXIT
            NB_FINAL= NB_FINAL + 1
         ENDDO
         !> this is the total number of energies considered 
         NB = NB_FINAL - MSDGW_NP_OCC
         IF( NB < 0 ) THEN
            CALL VTUTOR%ERROR("MSDGW_NP_OCC larger than number of occupied states: "//&
               str([ MSDGW_NP_OCC, NB_FINAL])//". Reduce MSDGW_NP_OCC will fix this problem.")
         ENDIF
      
         !> store index and number of compressed bands for this k-point 
         !> last frozed band index 
         CMAP(NK,ISP)%NB_START= 1
         !> last band considered for compression
         CMAP(NK,ISP)%NB_FINAL= NB_FINAL
         !> number of bands considered for compression 
         CMAP(NK,ISP)%NB_ORIG = NB
         !> store these energies in R
         ALLOCATE( R(  NB ))
         DO I = 1, NB
            R( I ) = W%CELTOT( I, NK, ISP )
         ENDDO

         CMAP(NK,ISP)%F = MSDGW_F_OCC
         !> determine number of bins first
         CALL SET_E_BINS( NB, NS, R, MSDGW_F_OCC )
         !WRITE(*,*)' bins for k-point: ', NK, ISP, NS, NB, NS + MSDGW_NP_OCC, NB_FINAL

         !> allocate bins 
         ALLOCATE( CMAP( NK, ISP)%EBIN( NS ))

         !> set bins now 
         CALL SET_E_BINS( NB, NS, R, MSDGW_F_OCC, CMAP(NK,ISP), MSDGW_NXI_OCC )
         !> define number of degenerate, occupied  bands / MSDGW_NXI_OCC
         CMAP( NK, ISP )%NB_DEGEN = SIZE( CMAP(NK,ISP)%EBIN, 1)
         !> total number of occupied bands, after compression
         NB_TOT = MSDGW_NP_OCC !+ CMAP(NK,ISP)%NB_DEGEN* MSDGW_NXI_OCC
         DO I = 1,  CMAP( NK, ISP )%NB_DEGEN
            NB_TOT = NB_TOT + CMAP(NK,ISP)%EBIN(I)%NXI 
         ENDDO

         !> this number should not exceed the number of occupied states 
         IF( NB_TOT > CMAP(NK,ISP)%NB_FINAL) THEN
            CALL VTUTOR%ERROR("number of compressed, occupied states &
               & MSDGW_NP_OCC + NB_DGEN*MSDGW_NXI_OCC="//&
               str([ MSDGW_NP_OCC, CMAP(NK,ISP)%NB_DEGEN, MSDGW_NXI_OCC])//&
               " exceeds the number of occupied states ("//str( NB_FINAL)//")&
               &. Reduce MSDGW_NXI_OCC, MSDGW_NP_OCC or increase MSDGW_F_OCC." )
         ENDIF

         CMAP( NK, ISP)%NB_TOT = NB_TOT 
         DEALLOCATE( R ) 

         !> run unit test 
         CALL UNIT_TEST_CMAP_SETUP( CMAP( NK, ISP), NB )

         IF( IO%NWRITE > 2) THEN
            CALL DUMP_IRARRAY( 'hist_occ', NK, CMAP(NK,ISP)%EBIN(:)%NBANDS, CMAP(NK,ISP)%EBIN(:)%EBAR  )
         ENDIF
      ENDDO  ! NK
      ENDDO  ! ISP
      PROFILING_STOP( 'set_const_f_cmap_occ')
      
   END SUBROUTINE SET_CONST_F_CMAP_OCC

   !> helper routine that returns number of energy bins, 
   !> optionally, CMAP is set as well 
   SUBROUTINE SET_E_BINS( NB, NS, R, F, CMAP, NXI)
      INTEGER, INTENT(IN) :: NB
      INTEGER, INTENT(INOUT) :: NS
      REAL(q), INTENT(IN) :: R(:) !< energies 
      REAL(q), INTENT(IN) :: F !< requested energy ratio 
      TYPE(compression_map), INTENT(INOUT), OPTIONAL ::  CMAP
      INTEGER, INTENT(IN), OPTIONAL :: NXI
      !> local variables
      INTEGER, PARAMETER :: ITER_MAX=100000
      INTEGER :: IDX1, IDX2
      INTEGER :: I, J
      REAL(q) :: E1, E2, DE

      !> number of pseudo bands must be passed too 
      IF( PRESENT( CMAP ) .AND. .NOT. PRESENT( NXI ) ) THEN
          CALL VTUTOR%BUG("Wrong usuage of this routine, NXI not passed", __FILE__, __LINE__ )
      ENDIF

      NS = 0
      IDX1 = 1
      !> set cmap as well 
      !> first get an estimate for the number of subdivisions 
      DO J = 1, ITER_MAX
         E1 = R( IDX1 ) 
         DE = ABS( E1* F )
         E2 = E1 + DE 
         !> find last energy smaller than E2 
         IDX2 = GET_IDX2(IDX1+1, NB, E2, R)
         !> last index reached 
         IF( IDX2 < 1 ) THEN
            WRITE(*,*)' last bin?', IDX1, IDX2, R( IDX1 ), J
            EXIT
         ENDIF
         !> increase number of bins
         NS = NS +1 

         !> set compression map, optionally 
         IF( PRESENT( CMAP ) ) THEN
            !> set number of indices for this bin 
            CMAP%EBIN( NS )%NBANDS =  IDX2 - IDX1 + 1 
            !> number of pseudo bands should be 
            CMAP%EBIN( NS )%NXI = MIN( CMAP%EBIN( NS )%NBANDS, NXI )

            !> set average energy and ratio for this bin
            CMAP%EBIN(NS)%EBAR = 0 
            CMAP%EBIN(NS)%F = 0 
            !> set band indices for this bin 
            ALLOCATE( CMAP%EBIN(NS)%IBAND(CMAP%EBIN( NS )%NBANDS))
            DO I = 1, CMAP%EBIN(NS)%NBANDS 
               CMAP%EBIN(NS)%IBAND( I ) = I + IDX1-1 
               CMAP%EBIN(NS)%EBAR = CMAP%EBIN(NS)%EBAR + R( I+IDX1-1 )
            ENDDO
            CMAP%EBIN(NS)%EBAR = CMAP%EBIN(NS)%EBAR/CMAP%EBIN(NS)%NBANDS
            CMAP%EBIN(NS)%F = (E2 - E1)/CMAP%EBIN(NS)%EBAR
            !WRITE(*,'(A,I4,A,3I5,A,2F10.4, 100I4)')' bin', NS,":", IDX1, IDX2, IDX2-IDX1+1, " :",CMAP%EBIN(NS)%EBAR, CMAP%EBIN(NS)%F, CMAP%EBIN(NS)%IBAND(1: CMAP%EBIN(NS)%NBANDS)
            !WRITE(*,'(A,I4,A,3I5,A,2F10.4, 100I4)')' bin', NS,":", IDX1, IDX2, IDX2-IDX1+1, " :",CMAP%EBIN(NS)%EBAR, CMAP%EBIN(NS)%F
         ENDIF

         !> increase starting index for next bin
         IDX1 = IDX2 + 1 

         !> in case last bin has been set, no need to continue iteration
         IF( IDX1 > NB ) EXIT 
      ENDDO
      IF( J > ITER_MAX ) THEN
         CALL VTUTOR%BUG( "ITER_MAX="//str( ITER_MAX)//" reached before all bins have been set", __FILE__, __LINE__)
      ENDIF

      CONTAINS 

      !> returns last index with energy below E2 
      !> If E2 is larger than last energy, last index is returned 
      FUNCTION GET_IDX2(IDX1, NB, E2, R)
         INTEGER :: GET_IDX2 
         INTEGER, INTENT(IN) :: IDX1, NB
         REAL(q), INTENT(IN) :: E2
         REAL(q), INTENT(IN) :: R(:)
         INTEGER :: I

         GET_IDX2 = -1

         !> maybe last bin has smaller energy window
         IF( R(NB) < E2  ) THEN
            GET_IDX2 = NB 
            RETURN
         ENDIF

         !> look for last index with energy below E2 
         DO I = IDX1, NB
           IF (R(I)>E2) THEN
             GET_IDX2 = I - 1
             !>optionally, add last index to prevent R(NB) forming its own bin
             IF( I == NB) THEN
                GET_IDX2 = NB
             ENDIF
             RETURN
           ENDIF
         ENDDO
      END FUNCTION GET_IDX2

   END SUBROUTINE SET_E_BINS

   !> simple unit test subroutine 
   SUBROUTINE UNIT_TEST_CMAP_SETUP( CMAP, NB )
      USE m_unirnk, ONLY: UNIRNK
      TYPE( compression_map ), INTENT(IN) :: CMAP
      INTEGER, INTENT(IN) :: NB 
      !> local variables
      INTEGER :: N, I, J
      INTEGER, ALLOCATABLE :: IDX(:)
      INTEGER, ALLOCATABLE :: JDX(:)

      !> 
      !> unit test: all energies must fall into 
      !> 
      !> one bin. This is equivalent that IDX_IN_BIN contains
      !> all indices, at least once
      N = 0 
      DO I = 1, CMAP%NB_DEGEN
         DO J = 1, CMAP%EBIN(I)%NBANDS
            IF (  CMAP%EBIN(I)%IBAND(J) <1 ) CYCLE 
            !> also, check if we are out of bounds
            IF (  CMAP%EBIN(I)%IBAND(J) > NB ) THEN
               CALL VTUTOR%BUG( "index in bin out of bounds: "//str( [ I,J,CMAP%EBIN(I)%IBAND(J), NB]), __FILE__, __LINE__ )
            ENDIF
            !> count number of indices 
            N = N + 1 
         ENDDO
      ENDDO
      IF ( N /= NB ) THEN
         CALL VTUTOR%BUG( "inconsistent number of indices in bins found: "//str( [ N, NB]), __FILE__, __LINE__ )
      ENDIF
      !> second test, no double entries are allowed
      ALLOCATE( IDX( N ))
      ALLOCATE( JDX( N ))
      N=0
      DO I = 1, CMAP%NB_DEGEN
         DO J = 1, CMAP%EBIN(I)%NBANDS
            IF (  CMAP%EBIN(I)%IBAND(J) <1 ) CYCLE 
            !> count number of indices 
            N = N + 1 
            !> populate test array with distributed indices 
            IDX( N ) = CMAP%EBIN(I)%IBAND(J)
         ENDDO
      ENDDO
      CALL UNIRNK( IDX, JDX, N )
      IF( N < NB ) THEN
         CALL VTUTOR%BUG( "double entry in bins found: "//str( [ N, NB]), __FILE__, __LINE__ )
      ENDIF

      DEALLOCATE( JDX )
      DEALLOCATE( IDX )
      
   END SUBROUTINE UNIT_TEST_CMAP_SETUP
      

   !> helper routine, checks parameter for consistency 
   !> @param[in]  WDES (wavedes) wavefunction decriptor
   !> @param[in]  IO (in_struct) in- out handle 
   SUBROUTINE CHECK_PARAM( WDES, IO)
      TYPE(wavedes), INTENT(IN) :: WDES
      TYPE(in_struct), INTENT(IN) :: IO
      !> local  
      INTEGER :: NK, ISP
      
      !> check if number of compressed bands is not too large
      DO ISP = 1, WDES%ISPIN
         DO NK = 1, WDES%NKPTS
            !> number of compressed bands is 
            IF( WDES%NB_TOTK( NK, ISP )  < MSDGW_NP+ MSDGW_NXI*MSDGW_NS ) THEN
               CALL VTUTOR%ERROR( "Number of compressed bands ("//str(MSDGW_NP+ MSDGW_NXI*MSDGW_NS)//&
                  ") exceeds actual number of bands ("//str(WDES%NB_TOTK(NK,ISP))//" at k-point:"//&
                  str([NK,ISP])//". You can reduce MSDGW_NS or MSDGW_NXI or both and try again."//&
                  str([ MSDGW_NP, MSDGW_NS, MSDGW_NXI]))
            ENDIF
         ENDDO
      ENDDO
   END SUBROUTINE CHECK_PARAM

   !> helper routine: performs the actual compression of WAVECAR
   !> based on the information stored in CMAP 
   !> @param[inout] W (wavespin)
   !> changes following attributes 
   !> W%WDES%NB_TOT
   !> W%CELTOT
   !> W%FERTOT
   !> W%CW
   SUBROUTINE COMPRESS_BANDS( WDES, W, CMAP, LATT_CUR, R, IO, CMAP_OCC )
      USE dfast
      USE wave_high
      USE scala, ONLY: INIT_SCALAAWARE , LscaAWARE, INIT_scala
      USE fock, ONLY : RESETUP_FOCK_WDES, RESETUP_FOCK
      USE random_seeded, ONLY: pcg_generator
      TYPE(wavedes), INTENT(INOUT) :: WDES
      TYPE(wavespin), INTENT(INOUT) :: W 
      TYPE(compression_map), INTENT(IN) :: CMAP(:,:)
      TYPE (latt), INTENT(IN) :: LATT_CUR 
      TYPE(pcg_generator), INTENT(INOUT)::  R
      TYPE(in_struct), INTENT(IN) :: IO 
      TYPE(compression_map), OPTIONAL, INTENT(IN) :: CMAP_OCC(:,:)
      !> local variables
      TYPE (wavespin)    WC !< compressed wavefunction
      TYPE(wavedes1) WDES1
      TYPE(wavefuna) WA
      GDEF, ALLOCATABLE :: U(:,:)
      INTEGER :: NK, ISP , I
      INTEGER :: NODE_ME = 1 
      INTEGER :: NB_COMP
      INTEGER :: NB_OCC
      INTEGER :: NB_ORIG
      REAL(q) :: F_I 
      REAL(q) :: F_AVG_K, F_AVG
      INTEGER :: NXI_ADD 
#ifdef MPI 
      NODE_ME = W%WDES%COMM%NODE_ME 
#endif
      PROFILING_START( 'compress_bands')

      !> determine new number of bands 
      NB_COMP = -1 
      DO ISP = 1, W%WDES%ISPIN
         DO NK = 1, W%WDES%NKPTS
            NB_COMP=MAX(NB_COMP, CMAP( NK, ISP)%NB_TOT) 
         ENDDO
      ENDDO

      !> todo allocate with reduced number of states 
      CALL ALLOCW( W%WDES, WC )
      !> copy wavefunction 
      WC%WDES  =W%WDES
      WC%CW    =W%CW
      WC%CPROJ =W%CPROJ 
      WC%FERWE =W%FERWE
      WC%AUX   =W%AUX
      WC%OVER_BAND=W%OVER_BAND
      WC%CELEN =W%CELEN
      WC%FERTOT=W%FERTOT
      WC%AUXTOT=W%AUXTOT
      WC%CELTOT=W%CELTOT

      !>
      !> set up rotation matrix for LINCOM
      !>
      F_AVG = 0 
      DO ISP = 1, W%WDES%ISPIN
         DO NK = 1, W%WDES%NKPTS
            IF (LscaAWARE) CALL INIT_scala(WDES%COMM_KIN,WDES%NB_TOTK(NK,ISP))
            ALLOCATE( U(W%WDES%NB_TOTK( NK, ISP),W%WDES%NB_TOTK( NK, ISP)) )
            IF( PRESENT( CMAP_OCC)) THEN
               CALL SET_PROJECTION( WC, NK, ISP, CMAP(NK,ISP), R, U, IO%IU0, CMAP_OCC(NK,ISP))
            ELSE
               CALL SET_PROJECTION( WC, NK, ISP, CMAP(NK,ISP), R, U, IO%IU0)
            ENDIF
            !rotate the Wavefunction 
            CALL SETWDES(W%WDES,WDES1,NK)
            WA=ELEMENTS(W,WDES1,ISP)
            IF (W%WDES%DO_REDIS) THEN
               CALL REDIS_PROJ(WDES1, W%WDES%NBANDS, WC%CPROJ(1,1,NK,ISP))
               CALL REDIS_PW  (WDES1, W%WDES%NBANDS, WC%CW   (1,1,NK,ISP))
            ENDIF
            ! build linear combinations of the wave functions in
            ! accordance with the effective rotation matrix U
            IF( LscaAWARE ) THEN
               CALL LINCOM_DISTRI('F',WA%CW_RED,WA%CPROJ_RED,U(1,1), &
                &   W%WDES%NB_TOTK(NK,ISP), &
                &   WDES1%NPL_RED,WDES1%NPRO_RED,W%WDES%NRPLWV_RED,W%WDES%NPROD_RED,W%WDES%NB_TOT, &
                &   W%WDES%COMM_KIN,NBLK)
            ELSE
               CALL LINCOM('F',WA%CW_RED,WA%CPROJ_RED,U(1,1), &
                &   W%WDES%NB_TOTK(NK,ISP),W%WDES%NB_TOTK(NK,ISP), &
                &   WDES1%NPL_RED,WDES1%NPRO_RED,W%WDES%NRPLWV_RED,W%WDES%NPROD_RED,W%WDES%NB_TOT, &
                &   WA%CW_RED,WA%CPROJ_RED)
            ENDIF

            ! redistribute back over bands
            IF (W%WDES%DO_REDIS) THEN
               CALL REDIS_PROJ(WDES1, WC%WDES%NBANDS, WC%CPROJ(1,1,NK,ISP))
               CALL REDIS_PW  (WDES1, WC%WDES%NBANDS, WC%CW   (1,1,NK,ISP))
            ENDIF

            DEALLOCATE( U ) 
            IF ( IO%IU0>=0 .AND. IO%NWRITE>2 ) WRITE(*, 888) NK, ISP, W%WDES%NB_TOTK(NK,ISP), &
                 CMAP(NK,ISP)%NB_TOT, REAL(CMAP(NK,ISP)%NB_TOT)/REAL(W%WDES%NB_TOTK(NK,ISP))*100
            IF ( IO%IU6>=6 .AND. IO%NWRITE>2 ) WRITE(IO%IU6, 888) NK, ISP, W%WDES%NB_TOTK(NK,ISP), &
                 CMAP(NK,ISP)%NB_TOT, REAL(CMAP(NK,ISP)%NB_TOT)/REAL(W%WDES%NB_TOTK(NK,ISP))*100

            !> dump avg. energy ratio , optionally for each bin 
            F_AVG_K = 0
            DO I = 1, CMAP(NK,ISP)%NB_DEGEN
               F_I = CMAP(NK,ISP)%EBIN(I)%F
               F_AVG_K = F_AVG_K + F_I 
               IF ( IO%IU0>=0 .AND. IO%NWRITE>3 ) WRITE(   *  , 889) I, F_I
               IF ( IO%IU6>=0 .AND. IO%NWRITE>3 ) WRITE(IO%IU6, 889) I, F_I
            ENDDO
            F_AVG_K = F_AVG_K / CMAP(NK,ISP)%NB_DEGEN
            IF ( IO%IU0>=0 .AND. IO%NWRITE>2 ) WRITE(   *  , 890) F_AVG_K
            IF ( IO%IU6>=0 .AND. IO%NWRITE>2 ) WRITE(IO%IU6, 890) F_AVG_K

            F_AVG = F_AVG + F_AVG_K * W%WDES%WTKPT( NK )

#ifdef debug_msdgw
            WRITE(*,*)'Re[ U x C_nG(1:'//str(IOUT)//')]'
            WRITE(*,FM1)( I ,I=1,IOUT)
            DO I = 1, IOUT
               WRITE(*,FM2)I, REAL(WC%CW(1:IOUT,I,1,1))
            ENDDO
            WRITE(101,'(2F20.10)')WC%CW(:,:,1,1)
!            WRITE(*,*)'Im( C_nG(1:'//str(IOUT)//')'
!            WRITE(*,FM1)( I ,I=1,IOUT)
!            DO I = 1, IOUT
!               WRITE(*,FM2)I, AIMAG(WC%CW(1:IOUT,I,1,1))
!            ENDDO
            STOP
#endif

         ENDDO
      ENDDO
      IF ( IO%IU0>=0 ) WRITE(   *  , 891) F_AVG*100 
      IF ( IO%IU6>=0 ) WRITE(IO%IU6, 891) F_AVG*100 

      !> reduce number of bands 
      WC%WDES%NB_TOT = NB_COMP

      IF( MSDGW_LWAVE ) THEN
         !> dump WAVECAR to WAVECAR.psd
         CALL OUTWAV(IO, WC%WDES, WC, LATT_CUR, W%EFERMI(1), EXT='psd' )
         IF ( IO%IU0>=0 ) WRITE(*, *) WC%WDES%NB_TOT, "bands written to WAVECAR.psd"
         IF ( IO%IU6>=0 ) WRITE(IO%IU6, *) WC%WDES%NB_TOT, "bands written to WAVECAR.psd"
      ENDIF
      IF ( IO%IU0>=0 ) WRITE(*, 892) WC%WDES%NB_TOT
      IF ( IO%IU6>=0 ) WRITE(IO%IU6, 892) WC%WDES%NB_TOT

      !> delete original wavefunction 
      CALL DEALLOCW(W)

      !> overwrite original wavedes
      WDES=WC%WDES
      WDES%NB_TOT=NB_COMP
      WDES%NBANDS=NB_COMP/WC%WDES%NB_PAR
      CALL INIT_SCALAAWARE( WDES%NB_TOT, WDES%NRPLWV, WDES%COMM_KIN )
      
      NULLIFY(WDES%NB_TOTK)
      ALLOCATE(WDES%NB_TOTK(WC%WDES%NKDIM,2))
! set the maximum number of bands k-point dependent
      DO NK=1,WDES%NKPTS
         IF (WDES%LGAMMA) THEN
            WDES%NB_TOTK(NK,:)=MIN(WDES%NB_TOT,WDES%NPLWKP_TOT(NK)*2-1)
         ELSE
            WDES%NB_TOTK(NK,:)=MIN(WDES%NB_TOT,WDES%NPLWKP_TOT(NK))
         ENDIF
      ENDDO
      CALL RESETUP_FOCK_WDES(WDES, LATT_CUR, LATT_CUR, -1)

      CALL ALLOCW(WDES,W)

      CALL DUMP_ALLOCATE(IO%IU6)

      W%FERTOT=0
      W%CELTOT=0
      ! copy data back to work array
      W%CW(:,1:WDES%NBANDS,:,:)    =WC%CW(:,1:WDES%NBANDS,:,:)
      W%CPROJ(:,1:WDES%NBANDS,:,:) =WC%CPROJ(:,1:WDES%NBANDS,:,:)
      W%CELTOT(1:WDES%NB_TOT,:,:)=WC%CELTOT(1:WDES%NB_TOT,:,:)
      W%FERTOT(1:WDES%NB_TOT,:,:)=WC%FERTOT(1:WDES%NB_TOT,:,:)

      !> remove temporary wavefunction
      CALL DEALLOCW( WC )

888   FORMAT( '    compression details for k-point, spin:', 2I4,/,&
              '      original number of bands:  ', I6, /,&
              '      bands after compression:   ', I6, /,&
              '      compression ratio (%):     ', F6.2 )
889   FORMAT( '       energy ratio F for bin:   ', I6, F6.3 )
890   FORMAT( '    energy ratio for k-point: ', F6.3 )
891   FORMAT( '    => avg. energy ratio F (%):', F6.2 )
892   FORMAT( '    bands after compression:   ', I6 )
      PROFILING_STOP( 'compress_bands')
      CONTAINS


      !> helper routine: sets the projection matrix that compresses bands 
      !> @param[in]  WDES (wavedes) wavefunction descriptor
      !> @param[in]  NK (integer) k-point
      !> @param[in]  ISP (integer) spin
      !> @param[in]  CMAP (compression_map) compression info
      !> @param[out] U (GDEF) projection matrix
      !> @param[in]  IU0 (integer) stdout unit 
      SUBROUTINE SET_PROJECTION( WC, NK, ISP, CMAP, R, U, IU0, CMAP_OCC)
         USE random_seeded, ONLY: pcg_generator, RANDOM_NUMBER
         TYPE(wavespin), INTENT(INOUT) :: WC
         INTEGER, INTENT(IN) :: NK
         INTEGER, INTENT(IN) :: ISP
         TYPE(compression_map), INTENT(IN) :: CMAP 
         TYPE(pcg_generator), INTENT(INOUT) :: R !< random number generator 
         GDEF, INTENT(INOUT) :: U(:,:)
         INTEGER, INTENT(IN) :: IU0 
         TYPE(compression_map), INTENT(IN), OPTIONAL :: CMAP_OCC 
         !> local variables 
         INTEGER :: I, I_, J, N, NXI, M
         INTEGER :: NB_TOT
         INTEGER :: NB_OUT
         COMPLEX(q) :: CPHASE 
         REAL(q):: RNORM
         REAL(q):: RAND
         PROFILING_START( 'set_projection|compress_bands')

         !> consistency check 
         IF( WC%WDES%NB_TOTK( NK, ISP )  < CMAP%NB_TOT ) THEN
            CALL VTUTOR%ERROR( "Number of compressed bands ("//str(CMAP%NB_TOT)//&
               ") exceeds actual number of bands ("//str(WC%WDES%NB_TOTK(NK,ISP))//" at k-point:"//&
               str([NK,ISP])//". You can reduce MSDGW_NS or MSDGW_NXI or both and try again."//&
               str([CMAP%NB_START,CMAP%NB_DEGEN,MSDGW_NXI]))
         ENDIF

         !> initialize rotation matrix 
         U =  0
         J=0
         I_=0
         IF ( PRESENT( CMAP_OCC )) THEN
            I_ = CMAP_OCC%NB_DEGEN*MSDGW_NXI_OCC
            J  = CMAP_OCC%NB_ORIG
         ENDIF
         !> keep frozen 
         DO I = 1 + J , CMAP%NB_START 
            I_ = I_+1 
            U( I, I_ ) = 1 
         ENDDO

         !> round to number of cores now 
         !> it assumes that every degenerate energy has the same number of pseudo bands 
         NB_TOT = CMAP%NB_START 
         DO N = 1, CMAP%NB_DEGEN !< loop over all degenerate energies
             NB_TOT = NB_TOT +  CMAP%EBIN(N)%NXI 
         ENDDO

         I_=0
         !> overwrite in case of compressed, occupied states 
         IF ( PRESENT( CMAP_OCC )) THEN
            !> now determine compressed bands 
            !RNORM = SQRT(1._q/MSDGW_NXI_OCC )
            DO N = 1, CMAP_OCC%NB_DEGEN !< loop over all degenerate energies
               RNORM = SQRT(1._q/CMAP_OCC%EBIN(N)%NXI )
               !WRITE(*,'(I4,":", 50I3)')N, CMAP_OCC%EBIN(N)%IBAND
   
               DO NXI = 1, MSDGW_NXI_OCC
                  !M =  CMAP_OCC%NB_START + MSDGW_NXI_OCC * (N-1) + NXI-1
                  M =  CMAP_OCC%NB_START + CMAP_OCC%EBIN(N)%NXI * (N-1) + NXI-1
                  !> set degenerate energies 
                  WC%CELTOT( M , NK, ISP) = CMAP_OCC%EBIN(N)%EBAR
                  !> normalization factor 
                  !RNORM = SQRT( REAL(CMAP_OCC%EBIN(N)%NBANDS ,q)/CMAP_OCC%EBIN(N)%NXI )
                  DO I = 1, CMAP_OCC%EBIN( N )%NBANDS !< loop over all original bands that have this energy
                     J = CMAP_OCC%EBIN( N)%IBAND( I )
                     !> determine random phase for each degenerate band 
                     CALL RANDOM_NUMBER(R, RAND )
                     CPHASE = EXP( CMPLX( 0._q, 2*PI * RAND, q) ) * RNORM
                     !> multiply with number of states in this bin 
                     U( J, M) = CPHASE 
                     !WRITE(*,*) J, M, ABS(CPHASE), RAND 
                  ENDDO
               ENDDO
            ENDDO
            !> set shift for occupied states 
            !NB_TOT = CMAP%NB_START -  CMAP_OCC%NB_FINAL + CMAP_OCC%NB_TOT & 
            !       + CMAP%NB_DEGEN* MSDGW_NXI 

            NB_TOT = CMAP%NB_START -  CMAP_OCC%NB_FINAL + CMAP_OCC%NB_TOT 
            DO N = 1, CMAP%NB_DEGEN
               NB_TOT = NB_TOT + CMAP%EBIN(N)%NXI
            ENDDO
            I_ = -CMAP_OCC%NB_FINAL + CMAP_OCC%NB_TOT
         ENDIF

         !> number of bands after compression: should be divisible by the number of ranks
         NB_OUT = ((NB_TOT+ WC%WDES%NB_PAR-1)/WC%WDES%NB_PAR)* WC%WDES%NB_PAR
         !IF(IU0>=0 ) WRITE(*,*) NB_TOT, NB_OUT, CMAP%NB_TOT, NB_OUT - NB_TOT 
         
         !> now determine compressed bands 
         DO N = 1, CMAP%NB_DEGEN !< loop over all degenerate energies
            RNORM = SQRT(1._q/CMAP%EBIN(N)%NXI )

            !> add additional pseudo bands to round to number of bands chosen 
            !> when running in parallel
            NXI_ADD =0
            IF( N==CMAP%NB_DEGEN ) NXI_ADD= NB_OUT - NB_TOT

            !DO NXI = 1, MSDGW_NXI + NXI_ADD
            !   M =  CMAP%NB_START + MSDGW_NXI * (N-1) + NXI  + I_
            DO NXI = 1, CMAP%EBIN(N)%NXI + NXI_ADD
               M =  CMAP%NB_START + CMAP%EBIN(N)%NXI * (N-1) + NXI  + I_
               !> set degenerate energies 
               WC%CELTOT( M , NK, ISP) = CMAP%EBIN(N)%EBAR
               !> normalization factor 
               !RNORM = SQRT( REAL(CMAP%EBIN(N)%NBANDS ,q)/CMAP%EBIN(N)%NXI )
               DO I = 1, CMAP%EBIN( N )%NBANDS !< loop over all original bands that have this energy
                  J = CMAP%EBIN( N)%IBAND( I ) + CMAP%NB_START
                  !> determine random phase for each degenerate band 
                  CALL RANDOM_NUMBER(R, RAND )
                  CPHASE = EXP( CMPLX( 0._q, 2*PI * RAND, q) ) * RNORM
                  !> multiply with number of states in this bin 
                  U( J, M) = CPHASE 
                  !WRITE(*,*) J, M, ABS(CPHASE), CPHASE, RAND
               ENDDO
            ENDDO
         ENDDO

#ifdef debug_msdgw
         !> dump U here 
         IF( IU0>=0 )THEN
            WRITE(*,*)'Re(U)'
               WRITE(*,FM1)( I ,I=1,IOUT)
            DO I = 1, IOUT
               WRITE(*,FM2)I, REAL(U(1:IOUT,I))
            ENDDO
!            WRITE(*,*)'Im(U)'
!            DO I = 1, IOUT
!               WRITE(*,FM2)I,AIMAG(U(I,1:IOUT))
!            ENDDO
!
            WRITE(*,*)'Re( C_nG(1:'//str(IOUT)//')'
            WRITE(*,FM1)( I ,I=1,IOUT)
            DO I = 1, IOUT
               WRITE(*,FM2)I, REAL(WC%CW(1:IOUT,I,1,1))
            ENDDO
            WRITE(100,'(2F20.10)')WC%CW(:,:,1,1)
!            WRITE(*,*)'Im( C_nG(1:'//str(IOUT)//')'
!            WRITE(*,FM1)( I ,I=1,IOUT)
!            DO I = 1, IOUT
!               WRITE(*,FM2)I, AIMAG(WC%CW(1:IOUT,I,1,1))
!            ENDDO
         ENDIF
#endif
         PROFILING_STOP( 'set_projection|compress_bands')
      END SUBROUTINE SET_PROJECTION
   END SUBROUTINE COMPRESS_BANDS

   SUBROUTINE DUMP_RARRAY( C, NK, RCOUNT, RMIN, RMAX )
      USE main_mpi, ONLY: COMM_WORLD
      CHARACTER(LEN=*)    :: C
      INTEGER, INTENT(IN) :: NK
      REAL(q), INTENT(IN) :: RCOUNT(:)
      REAL(q), INTENT(IN), OPTIONAL :: RMIN, RMAX 
      !> local variables
      CHARACTER(LEN=16)   :: EXT
      INTEGER :: IUNIT = 101 
      INTEGER :: I
      INTEGER :: NBINS
      REAL(q) :: RI, RI_
      INTEGER :: NODE_ME = 1
#ifdef MPI
      NODE_ME = COMM_WORLD%NODE_ME 
#endif

      NBINS = SIZE( RCOUNT )
      WRITE( EXT, * ) NK
   
      IF ( NODE_ME == 1 ) THEN 
         OPEN( NEWUNIT=IUNIT, FILE = C//'-'//TRIM(ADJUSTL(EXT))//'.dat', STATUS='REPLACE' )
   
         IF( PRESENT( RMIN )) THEN
            DO I = 1, NBINS
               RI = RMIN + (I-1) * ( RMAX - RMIN ) / REAL(NBINS-1)
               WRITE(IUNIT,'(2F20.10)') RI , RCOUNT( I )
            ENDDO
         ELSE
            DO I = 1, NBINS
               WRITE(IUNIT,'(2F20.10)') RCOUNT( I )
            ENDDO
         ENDIF
         CLOSE( IUNIT ) 
      ENDIF
   END SUBROUTINE DUMP_RARRAY

   SUBROUTINE DUMP_IARRAY( C, NK, RCOUNT, RMIN, RMAX )
      USE main_mpi, ONLY: COMM_WORLD
      CHARACTER(LEN=*)    :: C
      INTEGER, INTENT(IN) :: NK
      INTEGER, INTENT(IN) :: RCOUNT(:)
      REAL(q), INTENT(IN), OPTIONAL :: RMIN, RMAX 
      !> local variables
      CHARACTER(LEN=16)   :: EXT
      INTEGER :: IUNIT = 101 
      INTEGER :: I
      INTEGER :: NBINS
      REAL(q) :: RI, RI_
      INTEGER :: NODE_ME = 1
#ifdef MPI
      NODE_ME = COMM_WORLD%NODE_ME 
#endif

      NBINS = SIZE( RCOUNT )
      WRITE( EXT, * ) NK
   
      IF ( NODE_ME == 1 ) THEN 
         OPEN( NEWUNIT=IUNIT, FILE = C//'-'//TRIM(ADJUSTL(EXT))//'.dat', STATUS='REPLACE' )
   
         IF( PRESENT( RMIN )) THEN
            DO I = 1, NBINS
               RI = RMIN + (I-1) * ( RMAX - RMIN ) / REAL(NBINS-1)
               RI_= RMIN + I * ( RMAX - RMIN ) / REAL(NBINS-1)
               WRITE(IUNIT,'(F20.10, I20)') RI + (RI_-RI)/2, RCOUNT( I )
            ENDDO
         ELSE
            DO I = 1, NBINS
               WRITE(IUNIT,'(I20)') RCOUNT( I )
            ENDDO
         ENDIF
         CLOSE( IUNIT ) 
      ENDIF
   END SUBROUTINE DUMP_IARRAY

   SUBROUTINE DUMP_IRARRAY( C, NK, IX, Y)
      USE main_mpi, ONLY: COMM_WORLD
      CHARACTER(LEN=*)    :: C
      INTEGER, INTENT(IN) :: NK
      INTEGER, INTENT(IN) :: IX(:)
      REAL(q), INTENT(IN) :: Y(:)
      !> local variables
      CHARACTER(LEN=16)   :: EXT
      INTEGER :: IUNIT = 101 
      INTEGER :: I
      INTEGER :: NBINS
      REAL(q) :: RI, RI_
      INTEGER :: NODE_ME = 1
#ifdef MPI
      NODE_ME = COMM_WORLD%NODE_ME 
#endif

      NBINS = SIZE( Y )
      WRITE( EXT, * ) NK
   
      IF ( NODE_ME == 1 ) THEN 
         OPEN( NEWUNIT=IUNIT, FILE = C//'-'//TRIM(ADJUSTL(EXT))//'.dat', STATUS='REPLACE' )
         DO I = 1, NBINS
            WRITE(IUNIT,'(I20, F20.10)') IX(I), Y( I )
         ENDDO
         CLOSE( IUNIT ) 
      ENDIF
   END SUBROUTINE DUMP_IRARRAY

   SUBROUTINE DUMP_RRARRAY( C, NK, X, Y)
      USE main_mpi, ONLY: COMM_WORLD
      CHARACTER(LEN=*)    :: C
      INTEGER, INTENT(IN) :: NK
      REAL(q), INTENT(IN) :: X(:)
      REAL(q), INTENT(IN) :: Y(:)
      !> local variables
      CHARACTER(LEN=16)   :: EXT
      INTEGER :: IUNIT = 101 
      INTEGER :: I
      INTEGER :: NBINS
      REAL(q) :: RI, RI_
      INTEGER :: NODE_ME = 1
#ifdef MPI
      NODE_ME = COMM_WORLD%NODE_ME 
#endif

      NBINS = SIZE( X )
      WRITE( EXT, * ) NK
   
      IF ( NODE_ME == 1 ) THEN 
         OPEN( NEWUNIT=IUNIT, FILE = C//'-'//TRIM(ADJUSTL(EXT))//'.dat', STATUS='REPLACE' )
         DO I = 1, NBINS
            WRITE(IUNIT,'(2F20.10)')X(I), Y( I )
         ENDDO
         CLOSE( IUNIT ) 
      ENDIF
   END SUBROUTINE DUMP_RRARRAY

   SUBROUTINE DUMP_IIARRAY( C, NK, X, Y)
      USE main_mpi, ONLY: COMM_WORLD
      CHARACTER(LEN=*)    :: C
      INTEGER, INTENT(IN) :: NK
      INTEGER, INTENT(IN) :: X(:)
      INTEGER, INTENT(IN) :: Y(:)
      !> local variables
      CHARACTER(LEN=16)   :: EXT
      INTEGER :: IUNIT = 101 
      INTEGER :: I
      INTEGER :: NBINS
      REAL(q) :: RI, RI_
      INTEGER :: NODE_ME = 1
#ifdef MPI
      NODE_ME = COMM_WORLD%NODE_ME 
#endif

      NBINS = SIZE( X )
      WRITE( EXT, * ) NK
   
      IF ( NODE_ME == 1 ) THEN 
         OPEN( NEWUNIT=IUNIT, FILE = C//'-'//TRIM(ADJUSTL(EXT))//'.dat', STATUS='REPLACE' )
         DO I = 1, NBINS
            WRITE(IUNIT,'(2I20)')X(I), Y( I )
         ENDDO
         CLOSE( IUNIT ) 
      ENDIF
   END SUBROUTINE DUMP_IIARRAY

END MODULE msdgw
