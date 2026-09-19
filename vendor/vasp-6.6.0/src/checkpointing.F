#include "symbol.inc"

MODULE checkpointing
   USE prec
   USE tutor, ONLY: vtutor
   USE mpimy, only: communic
   USE string, ONLY: str, lowercase
#ifdef VASP_HDF5
   USE hdf5, only: HID_T
   USE vhdf5_base
   USE vhdf5, only: vh5_write_version
   IMPLICIT NONE
#else
   IMPLICIT NONE
   ! Define HID_T as regular integer when HDF5 is not available
   INTEGER, PARAMETER :: HID_T = KIND(0)
#endif

#ifdef realmode
   LOGICAL, PARAMETER :: IS_REALMODE = .TRUE.
#else
   LOGICAL, PARAMETER :: IS_REALMODE = .FALSE.
#endif
#ifdef gammareal
   LOGICAL, PARAMETER :: IS_GAMMAREAL = .TRUE.
#else
   LOGICAL, PARAMETER :: IS_GAMMAREAL = .FALSE.
#endif

   TYPE STRING_CONTAINER
      CHARACTER(:), ALLOCATABLE :: S
   END TYPE

   ! Checkpoint modes for general checkpoint settings
   ENUM, BIND(C)
      ENUMERATOR :: CHECKPOINT_MODE = 0
      !> Continue from existing checkpoint or start new
      ENUMERATOR :: CHECKPOINT_MODE_CONTINUE
      !> Reset checkpoint and start anew
      ENUMERATOR :: CHECKPOINT_MODE_RESET
      !> Prepare a split checkpointing scheme
      ENUMERATOR :: CHECKPOINT_MODE_PREPARE
      !> Run a single indenpendent calculation
      !> as part of a split checkpointing scheme
      ENUMERATOR :: CHECKPOINT_MODE_SINGLE
      !> Mode for compatibility
      !> No checkpoint logic will be executed
      !> No checkpoint IO will be performed
      ENUMERATOR :: CHECKPOINT_MODE_NONE
  END ENUM

   ! General checkpoint settings struct
   TYPE checkpoint_settings_struct
      INTEGER(KIND(CHECKPOINT_MODE)) :: MODE = CHECKPOINT_MODE_RESET
   CONTAINS
      PROCEDURE :: DO_CONTINUE => CHECKPOINT_SETTINGS_DO_CONTINUE
      PROCEDURE :: DO_PREPARE => CHECKPOINT_SETTINGS_DO_PREPARE
      PROCEDURE :: DO_SINGLE => CHECKPOINT_SETTINGS_DO_SINGLE
      PROCEDURE :: DO_OVERWRITE_CHECKPOINT => CHECKPOINT_SETTINGS_DO_OVERWRITE_CHECKPOINT
      PROCEDURE :: DO_WRITE_CHECKPOINT => CHECKPOINT_SETTINGS_DO_WRITE_CHECKPOINT
      PROCEDURE :: DOES_FINISH => CHECKPOINT_SETTINGS_DOES_FINISH
   END TYPE

   ! Checkpoint file handler struct
   TYPE checkpoint_struct
      !> Path to checkpoint file
      CHARACTER(:), ALLOCATABLE, PRIVATE :: FILE_PATH
      !> Whether checkpoint has been initialized
      LOGICAL, PRIVATE :: IS_INITIALIZED = .FALSE.
      !> Whether checkpointing is enabled (can be false even after initialization)
      LOGICAL, PRIVATE :: ENABLED = .FALSE.
      !> Whether a data group is currently open for writing
      LOGICAL, PRIVATE :: DATA_IS_OPEN = .FALSE.
      !> Whether to use subdirectories for checkpoint files
      !> This is used for splitting the calculations into multiple separate runs
      LOGICAL, PRIVATE :: USE_SUBDIRS = .FALSE.
      !> Whether checkpoint is in reading mode
      LOGICAL, PRIVATE :: IS_READING = .FALSE.
      !> Number of available data groups that hold actual completed data
      INTEGER, PRIVATE :: DATA_COUNT = 0
      !> Subdirectory file paths for split checkpointing
      TYPE(STRING_CONTAINER), ALLOCATABLE, PRIVATE :: SUB_FILE_PATHS(:)
   CONTAINS
      PROCEDURE :: INIT => CHECKPOINT_INIT
      PROCEDURE :: INIT_DISABLED => CHECKPOINT_INIT_DISABLED
      PROCEDURE :: OPEN_NEW_DATA_GROUP => CHECKPOINT_OPEN_NEW_DATA_GROUP
      PROCEDURE :: CLOSE_DATA_GROUP => CHECKPOINT_CLOSE_DATA_GROUP
      PROCEDURE :: WRITE_DATA => CHECKPOINT_WRITE_DATA
      PROCEDURE :: READ_DATA => CHECKPOINT_READ_DATA
      PROCEDURE :: WRITE_METADATA => CHECKPOINT_WRITE_METADATA
      PROCEDURE :: WRITE_SUBDIR_METADATA => CHECKPOINT_WRITE_SUBDIR_METADATA
      PROCEDURE :: READ_METADATA => CHECKPOINT_READ_METADATA
      PROCEDURE :: SHOULD_READ => CHECKPOINT_SHOULD_READ
      PROCEDURE :: USES_SUBDIRS => CHECKPOINT_USES_SUBDIRS
   END TYPE

   ! Generic interfaces for checkpoint operations
   INTERFACE
      SUBROUTINE CHECKPOINT_HANDLER(GROUP_ID)
         IMPORT :: HID_T
         INTEGER(HID_T), INTENT(IN) :: GROUP_ID
      END SUBROUTINE
   END INTERFACE

   INTERFACE
      SUBROUTINE CHECKPOINT_HANDLER_IDX(GROUP_ID, IDX)
         IMPORT :: HID_T
         INTEGER(HID_T), INTENT(IN) :: GROUP_ID
         INTEGER, INTENT(IN) :: IDX
      END SUBROUTINE
   END INTERFACE

CONTAINS

   !*********************************************************************
   !
   ! Methods for checkpoint_settings_struct
   !
   !*********************************************************************

   !> Check if we should continue checkpointing.
   !> This is true for both CONTINUE and RESET modes, as both involve writing checkpoint data
   !> (RESET will overwrite existing checkpoint, while CONTINUE will append to it).
   PURE LOGICAL FUNCTION CHECKPOINT_SETTINGS_DO_CONTINUE(THIS) RESULT(DOES)
      CLASS(checkpoint_settings_struct), INTENT(IN) :: THIS
      DOES = THIS%MODE == CHECKPOINT_MODE_CONTINUE .OR. THIS%MODE == CHECKPOINT_MODE_RESET
   END FUNCTION

   !> Check if we should prepare a split checkpointing scheme.
   !> This is true only for PREPARE mode.
   PURE LOGICAL FUNCTION CHECKPOINT_SETTINGS_DO_PREPARE(THIS) RESULT(DOES)
      CLASS(checkpoint_settings_struct), INTENT(IN) :: THIS
      DOES = THIS%MODE == CHECKPOINT_MODE_PREPARE
   END FUNCTION

   !> Check if we should run a single independent calculation.
   !> This is true only for SINGLE mode.
   PURE LOGICAL FUNCTION CHECKPOINT_SETTINGS_DO_SINGLE(THIS) RESULT(DOES)
      CLASS(checkpoint_settings_struct), INTENT(IN) :: THIS
      DOES = THIS%MODE == CHECKPOINT_MODE_SINGLE
   END FUNCTION

   !> Check if we should overwrite existing checkpoint data.
   !> This is true for modes that involve writing checkpoint data and
   !> where we want to start fresh (RESET, PREPARE, SINGLE).
   PURE LOGICAL FUNCTION CHECKPOINT_SETTINGS_DO_OVERWRITE_CHECKPOINT(THIS) RESULT(DOES)
      CLASS(checkpoint_settings_struct), INTENT(IN) :: THIS
      DOES = &
         THIS%MODE == CHECKPOINT_MODE_RESET .OR. &
         THIS%MODE == CHECKPOINT_MODE_PREPARE .OR. &
         THIS%MODE == CHECKPOINT_MODE_SINGLE
   END FUNCTION

   !> Check if we should write any kind of checkpoint data.
   PURE LOGICAL FUNCTION CHECKPOINT_SETTINGS_DO_WRITE_CHECKPOINT(THIS) RESULT(DOES)
      CLASS(checkpoint_settings_struct), INTENT(IN) :: THIS
      DOES = THIS%MODE == CHECKPOINT_MODE_CONTINUE .OR. THIS%DO_OVERWRITE_CHECKPOINT()
   END FUNCTION

   !> Check if the current checkpoint mode implies that the calculation will finish in a single run.
   !> What "finish" means here is that all the checkpointed information will be available by the end.
   !> So SINGLE and PREPARE modes will not finish as they do not compute all data in one go.
   PURE LOGICAL FUNCTION CHECKPOINT_SETTINGS_DOES_FINISH(THIS) RESULT(DOES)
      CLASS(checkpoint_settings_struct), INTENT(IN) :: THIS
      DOES = .NOT. (THIS%MODE == CHECKPOINT_MODE_SINGLE .OR. THIS%MODE == CHECKPOINT_MODE_PREPARE)
   END FUNCTION

   !*********************************************************************
   !
   ! Methods for checkpoint_struct
   !
   !*********************************************************************

   !> @brief Initializes the checkpoint structure
   !>
   !> This subroutine initializes the checkpoint structure by opening or
   !> creating the checkpoint file based on the specified mode. It handles
   !> both reading and writing modes, checks for existing checkpoint data,
   !> and sets up the internal state of the checkpoint structure.
   !> The behavior depends on the mode:
   !> - In 'w' (write) mode, it creates a new checkpoint file, overwriting
   !>   any existing one.
   !> - In 'a' (append) mode, it checks if the checkpoint file exists.
   !>   If it does, it opens it for reading; if not, it creates a new file
   !>   for writing. When opening an existing file, it verifies the
   !>   integrity of the checkpoint data and counts the number of
   !>   completed data groups.
   !> The communicator is only used to check for file existence on the
   !> I/O node and broadcast the result to all nodes.
   SUBROUTINE CHECKPOINT_INIT(THIS, FILE_PATH, MODE, COMM)
      USE mpimy, ONLY: communic

      !> Checkpoint structure to initialize
      CLASS(checkpoint_struct), INTENT(OUT) :: THIS
      !> Path to checkpoint file
      CHARACTER(*), INTENT(IN) :: FILE_PATH
      !> Mode for opening the checkpoint file ('w' or 'a')
      CHARACTER(1), INTENT(IN) :: MODE
      !> Communicator for checking file existence
      TYPE(communic), INTENT(IN) :: COMM

      INTEGER(HID_T) :: FILE_ID
      INTEGER :: IDX
      CHARACTER(:), ALLOCATABLE :: FILE_NAME
      CHARACTER(:), ALLOCATABLE :: SUBDIR_PREFIX

#ifdef VASP_HDF5
      IF (THIS%DATA_IS_OPEN) CALL vtutor%bug( &
         'While trying to initialize checkpoint: Data group already open', __FILE__, __LINE__)

      THIS%FILE_PATH = FILE_PATH
      THIS%DATA_COUNT = 0
      THIS%USE_SUBDIRS = .FALSE.
      THIS%IS_READING = .FALSE.
      IF (ALLOCATED(THIS%SUB_FILE_PATHS)) DEALLOCATE(THIS%SUB_FILE_PATHS)

      SELECT CASE (lowercase(MODE))
      CASE ('w')
         VH5_CHECK(vh5_file_create_or_overwrite(FILE_PATH, FILE_ID))
      CASE ('a')
         IF (COMM%NODE_ME == COMM%IONODE) INQUIRE(FILE=FILE_PATH, EXIST=THIS%IS_READING)
         CALLMPI(M_bcast_l(COMM, THIS%IS_READING, 1))
         IF (THIS%IS_READING) THEN
            VH5_CHECK(vh5_file_open_read(FILE_PATH, FILE_ID))
         ELSE
            VH5_CHECK(vh5_file_create(FILE_PATH, FILE_ID))
         ENDIF
      CASE DEFAULT
         CALL vtutor%bug('Unknown mode ' // MODE // ' in checkpoint init', __FILE__, __LINE__)
      END SELECT

      ! Here, we perform a few checks if we are in reading mode
      ! to ensure that the checkpoint data is valid
      ! In particular, we need the number of completed steps
      IF (THIS%IS_READING) THEN
         THIS%USE_SUBDIRS = vh5_exists(FILE_ID, 'subdir_prefix')

         IF (THIS%USE_SUBDIRS) THEN
            ! In the subdirectory scheme, we require additional metadata
            VH5_CHECK(vh5_read(FILE_ID, 'subdir_prefix', SUBDIR_PREFIX))
            VH5_CHECK(vh5_read(FILE_ID, 'total_count', THIS%DATA_COUNT))

            IDX = SCAN(FILE_PATH, '/', BACK=.TRUE.)
            FILE_NAME = FILE_PATH(IDX+1:)

            ALLOCATE(THIS%SUB_FILE_PATHS(THIS%DATA_COUNT))
            DO IDX = 1, THIS%DATA_COUNT
               THIS%SUB_FILE_PATHS(IDX)%S = SUBDIR_PREFIX // str(IDX) // '/' // FILE_NAME
            ENDDO
         ENDIF

         CALL VALIDATE()
         VH5_CHECK(vh5_file_close(FILE_ID))
      ELSE
         ! If we are writing, write the version info
         CALL vh5_write_version(FILE_ID)
         VH5_CHECK(vh5_file_close_writing(FILE_ID))
      ENDIF

      THIS%IS_INITIALIZED = .TRUE.
      THIS%ENABLED = .TRUE.

   CONTAINS

      !> @brief Validates an existing checkpoint file
      !>
      !> This subroutine validates an existing checkpoint file by checking the
      !> integrity of the checkpoint data and ensuring compatibility with the
      !> current executable. It also determines the number of completed data
      !> groups by discarding any incomplete ones.
      !>
      !> For the validation of VASP versions, the following table summarizes
      !> the expected behavior based on the executable type that created the
      !> checkpoint and the current executable type:
      !>
      !> | File \ Exe | gam   | std   | ncl   |
      !> |------------|-------|-------|-------|
      !> | gam        | OK    | alert | alert |
      !> | std        | error | OK    | alert |
      !> | ncl        | error | error | OK    |
      !>
      !> - OK: checkpoint can be read without issues.
      !> - alert: checkpoint can be read, but the executable types differ.
      !>          An alert can be emitted but execution continues.
      !> - error: checkpoint cannot be read due to incompatibility.
      !>          The program will stop with an error message.
      SUBROUTINE VALIDATE()
         INTEGER(HID_T) :: SUB_FILE_ID
         LOGICAL :: IS_OK, EXISTS
         CHARACTER(:), ALLOCATABLE :: PATH
         LOGICAL :: FILE_REALMODE, FILE_GAMMAREAL
         INTEGER :: COMPLEXITY, FILE_COMPLEXITY
         CHARACTER(:), ALLOCATABLE :: ALERT_MESSAGE

         COMPLEXITY = GET_EXE_COMPLEXITY(IS_REALMODE, IS_GAMMAREAL)
         ALERT_MESSAGE = ''

         IF (THIS%USE_SUBDIRS) THEN
            DO IDX = 1, THIS%DATA_COUNT
               PATH = THIS%SUB_FILE_PATHS(IDX)%S

               INQUIRE(FILE=PATH, EXIST=EXISTS)

               IF (.NOT. EXISTS) CALL vtutor%error( &
                  'While trying to read checkpoint: Checkpoint file "' // PATH // &
                  '" does not exist, but required for restart. ' // &
                  'Please (re)run the affected single calculation.')

               VH5_CHECK(vh5_file_open_read(PATH, SUB_FILE_ID))

               IF (vh5_read(SUB_FILE_ID, 'data-1/completed', IS_OK) /= 0) &
                  CALL vtutor%error( &
                  'While trying to read checkpoint: Checkpoint file "' // PATH // &
                  '" is corrupted: missing or incomplete data group. ' // &
                  'Please (re)run the affected single calculation.')

               IF (.NOT. IS_OK) CALL vtutor%error( &
                  'While trying to read checkpoint: Checkpoint file "' // PATH // &
                  '" is corrupted: data group not completed. ' // &
                  'Please (re)run the affected single calculation.')

               VH5_CHECK(vh5_read(SUB_FILE_ID, 'data-1/realmode', FILE_REALMODE))
               VH5_CHECK(vh5_read(SUB_FILE_ID, 'data-1/gammareal', FILE_GAMMAREAL))
               FILE_COMPLEXITY = GET_EXE_COMPLEXITY(FILE_REALMODE, FILE_GAMMAREAL)

               IF (FILE_COMPLEXITY > COMPLEXITY) THEN
                  CALL COMPLEXITY_ERROR(COMPLEXITY, FILE_COMPLEXITY, PATH)
               ELSE IF (FILE_COMPLEXITY < COMPLEXITY) THEN
                  ALERT_MESSAGE = ALERT_MESSAGE // PATH // ' '
               ENDIF

               VH5_CHECK(vh5_file_close(SUB_FILE_ID))
            ENDDO

            CALL COMPLEXITY_ALERT(COMPLEXITY, ALERT_MESSAGE)
         ELSE
            DO
               PATH = 'data-'//str(THIS%DATA_COUNT + 1)
               IF (vh5_read(FILE_ID, PATH//'/completed', IS_OK) /= 0) EXIT
               IF (.NOT. IS_OK) EXIT
               THIS%DATA_COUNT = THIS%DATA_COUNT + 1

               VH5_CHECK(vh5_read(FILE_ID, PATH//'/realmode', FILE_REALMODE))
               VH5_CHECK(vh5_read(FILE_ID, PATH//'/gammareal', FILE_GAMMAREAL))
               FILE_COMPLEXITY = GET_EXE_COMPLEXITY(FILE_REALMODE, FILE_GAMMAREAL)

               IF (FILE_COMPLEXITY > COMPLEXITY) THEN
                  CALL COMPLEXITY_ERROR(COMPLEXITY, FILE_COMPLEXITY, THIS%FILE_PATH)
               ELSE IF (FILE_COMPLEXITY < COMPLEXITY) THEN
                  ALERT_MESSAGE = ALERT_MESSAGE // str(THIS%DATA_COUNT) // ' '
               ENDIF
            ENDDO

            CALL COMPLEXITY_ALERT(COMPLEXITY, ALERT_MESSAGE, THIS%FILE_PATH)
         ENDIF
      END SUBROUTINE

      SUBROUTINE COMPLEXITY_ERROR(COMPLEXITY, FILE_COMPLEXITY, PATH)
         INTEGER, INTENT(IN) :: COMPLEXITY
         INTEGER, INTENT(IN) :: FILE_COMPLEXITY
         CHARACTER(*), INTENT(IN) :: PATH
         CALL vtutor%error( &
            'While trying to read checkpoint: Checkpoint file "' // PATH // &
            '" was created by the ' // GET_EXE_NAME(FILE_COMPLEXITY) // &
            ', which cannot be read by the ' // GET_EXE_NAME(COMPLEXITY) // &
            '. Please restart using a compatible executable or regenerate the checkpoint.')
      END SUBROUTINE

      SUBROUTINE COMPLEXITY_ALERT(COMPLEXITY, ALERT_MESSAGE, PATH)
         INTEGER, INTENT(IN) :: COMPLEXITY
         CHARACTER(*), INTENT(IN) :: ALERT_MESSAGE
         CHARACTER(*), OPTIONAL, INTENT(IN) :: PATH

         CHARACTER(:), ALLOCATABLE :: CREATED_BY, STR1, STR2

         IF (ALERT_MESSAGE == '') RETURN

         CREATED_BY = GET_EXE_NAME(1)
         IF (COMPLEXITY == 3) CREATED_BY = CREATED_BY // ' or ' // GET_EXE_NAME(2)

         IF (PRESENT(PATH)) THEN
            STR1 = 'Some steps in checkpoint file "' // PATH // '"'
            STR2 = 'steps'
         ELSE
            STR1 = 'Some checkpoint files'
            STR2 = 'files'
         ENDIF

         CALL vtutor%alert( &
            'While trying to read checkpoint: ' // STR1 // ' were created by the ' // &
            CREATED_BY // ' but you are using the ' // GET_EXE_NAME(COMPLEXITY) // &
            '. The checkpoint data can still be read, but take extra care to ensure ' // &
            'that the input parameters and settings give the same results.' // &
            NEW_LINE('A') // NEW_LINE('A') // 'Affected ' // STR2 // ': ' // ALERT_MESSAGE)
      END SUBROUTINE
#endif
   END SUBROUTINE

   !> @ brief Initializes the checkpoint structure in a disabled state
   !>
   !> After calling this subroutine, the checkpoint structure will
   !> be initialized but permanently disabled.
   !> Any subsequent calls to checkpoint methods will be no-ops,
   !> and no file I/O will occur.
   SUBROUTINE CHECKPOINT_INIT_DISABLED(THIS)
      CLASS(checkpoint_struct), INTENT(OUT) :: THIS

      ! Stub initialization - checkpoint is permanently disabled
      ! No file I/O, minimal setup for a no-op checkpoint
      THIS%ENABLED = .FALSE.
      THIS%IS_INITIALIZED = .TRUE.
      THIS%DATA_IS_OPEN = .FALSE.
      THIS%USE_SUBDIRS = .FALSE.
      THIS%IS_READING = .FALSE.
      THIS%DATA_COUNT = 0
      THIS%FILE_PATH = ''
   END SUBROUTINE

   !> @brief Opens a new data group for writing checkpoint data
   !>
   !> This subroutine opens a new data group in the checkpoint file for writing.
   !> It handles both the case where we are using subdirectories
   !> (for split checkpointing) and the case where we are not.
   !>
   !> A data group in this context refers to a group in the HDF5 file that
   !> represents a single completed step. When using subdirectories, each step
   !> is stored in a separate file, and the data group is created in that file.
   !>
   !> The subroutine also automatically switches to writing mode if we were
   !> previously in reading mode.
   SUBROUTINE CHECKPOINT_OPEN_NEW_DATA_GROUP(THIS)
      CLASS(checkpoint_struct), INTENT(INOUT) :: THIS

      INTEGER(HID_T) :: FILE_ID, DATA_ID
      INTEGER :: IGNORE

#ifdef VASP_HDF5
      IF (.NOT. THIS%IS_INITIALIZED) CALL vtutor%bug( &
         'While trying to open data group: Checkpoint not initialized', __FILE__, __LINE__)
      IF (.NOT. THIS%ENABLED) RETURN
      IF (THIS%DATA_IS_OPEN) CALL vtutor%bug( &
         'While trying to open data group: Data group already open', __FILE__, __LINE__)

      THIS%IS_READING = .FALSE.
      THIS%DATA_COUNT = THIS%DATA_COUNT + 1

      VH5_CHECK(vh5_file_open_readwrite(THIS%FILE_PATH, FILE_ID))
      ! Delete any existing data group with the same name
      IGNORE = vh5_delete(FILE_ID, 'data-'//str(THIS%DATA_COUNT))
      ! Create new data group
      VH5_CHECK(vh5_group_create(FILE_ID, 'data-'//str(THIS%DATA_COUNT), DATA_ID))
      VH5_CHECK(vh5_write(DATA_ID, 'realmode', IS_REALMODE))
      VH5_CHECK(vh5_write(DATA_ID, 'gammareal', IS_GAMMAREAL))
      ! Close so we minimize risk of corruption
      VH5_CHECK(vh5_group_close_writing(DATA_ID))
      VH5_CHECK(vh5_file_close_writing(FILE_ID))

      THIS%DATA_IS_OPEN = .TRUE.
#endif
   END SUBROUTINE

   !> @brief Closes the currently open data group
   !>
   !> This subroutine closes the currently open data group in the checkpoint file.
   !> It also marks the data group as completed by writing a flag, and flushes
   !> the file to ensure data integrity.
   SUBROUTINE CHECKPOINT_CLOSE_DATA_GROUP(THIS)
      CLASS(checkpoint_struct), INTENT(INOUT) :: THIS

      INTEGER(HID_T) :: FILE_ID

#ifdef VASP_HDF5
      IF (.NOT. THIS%IS_INITIALIZED) CALL vtutor%bug( &
         'While trying to close data group: Checkpoint not initialized', __FILE__, __LINE__)
      IF (.NOT. THIS%ENABLED) RETURN
      IF (.NOT. THIS%DATA_IS_OPEN) CALL vtutor%bug( &
         'While trying to close data group: Data group not open', __FILE__, __LINE__)
      IF (THIS%IS_READING) CALL vtutor%bug( &
         'While trying to close data group: Checkpoint is in read mode', __FILE__, __LINE__)

      VH5_CHECK(vh5_file_open_readwrite(THIS%FILE_PATH, FILE_ID))
      VH5_CHECK(vh5_write(FILE_ID, 'data-'//str(THIS%DATA_COUNT)//'/completed', .TRUE.))
      VH5_CHECK(vh5_file_close_writing(FILE_ID))

      THIS%DATA_IS_OPEN = .FALSE.
#endif
   END SUBROUTINE

   !> @brief Writes data to the checkpoint file using a provided handler
   !>
   !> This subroutine writes data to the checkpoint file by calling a
   !> user-provided handler subroutine. The handler is responsible for writing
   !> the actual data to the provided group ID. This method handles the logic
   !> of opening a new group for writing and ensures that the checkpoint is in
   !> the correct state for writing.
   SUBROUTINE CHECKPOINT_WRITE_DATA(THIS, NAME, WRITER)
      CLASS(checkpoint_struct), INTENT(INOUT) :: THIS !< Checkpoint structure
      CHARACTER(*), INTENT(IN) :: NAME !< Name of the data group to write
      PROCEDURE(CHECKPOINT_HANDLER) :: WRITER !< Writer handler

      INTEGER(HID_T) :: FILE_ID, GROUP_ID

#ifdef VASP_HDF5
      IF (.NOT. THIS%IS_INITIALIZED) CALL vtutor%bug( &
         'While trying to write data group "'//NAME//'": Checkpoint not initialized', __FILE__, __LINE__)
      IF (.NOT. THIS%ENABLED) RETURN
      IF (.NOT. THIS%DATA_IS_OPEN) CALL vtutor%bug( &
         'While trying to write data group "'//NAME//'": Data group not open (not in write mode)', __FILE__, __LINE__)
      IF (THIS%IS_READING) CALL vtutor%bug( &
         'While trying to write data group "'//NAME//'": Checkpoint is in read mode', __FILE__, __LINE__)

      VH5_CHECK(vh5_file_open_readwrite(THIS%FILE_PATH, FILE_ID))
      VH5_CHECK(vh5_group_create(FILE_ID, 'data-'//str(THIS%DATA_COUNT)//'/'//NAME, GROUP_ID))
      CALL WRITER(GROUP_ID)
      VH5_CHECK(vh5_group_close_writing(GROUP_ID))
      VH5_CHECK(vh5_file_close_writing(FILE_ID))
#endif
   END SUBROUTINE

   !> @brief Reads data from the checkpoint file using a provided handler
   !>
   !> This subroutine reads data from the checkpoint file by calling a
   !> user-provided handler subroutine. The handler is responsible for reading
   !> the actual data from the provided group ID. This method handles the logic
   !> of opening the group for reading and ensures that the checkpoint is in
   !> the correct state for reading. The reader handler also receives the index
   !> of the data group being read, as the code loops over all available data
   !> groups when reading.
   SUBROUTINE CHECKPOINT_READ_DATA(THIS, NAME, READER)
      CLASS(checkpoint_struct), INTENT(IN) :: THIS !< Checkpoint structure
      CHARACTER(*), INTENT(IN) :: NAME !< Name of the data group to read
      PROCEDURE(CHECKPOINT_HANDLER_IDX) :: READER !< Reader handler

      INTEGER(HID_T) :: FILE_ID, GROUP_ID
      INTEGER :: IDX
      CHARACTER(:), ALLOCATABLE :: PATH
      LOGICAL :: IS_OK

#ifdef VASP_HDF5
      IF (.NOT. THIS%IS_INITIALIZED) CALL vtutor%bug( &
         'While trying to read data group "'//NAME//'": Checkpoint not initialized', __FILE__, __LINE__)
      IF (.NOT. THIS%ENABLED) RETURN
      IF (THIS%DATA_IS_OPEN) CALL vtutor%bug( &
         'While trying to read data group "'//NAME//'": Data group is open (not in read mode)', __FILE__, __LINE__)
      IF (.NOT. THIS%IS_READING) CALL vtutor%bug( &
         'While trying to read data group "'//NAME//'": Checkpoint is not in read mode', __FILE__, __LINE__)

      IF (THIS%USE_SUBDIRS) THEN
         DO IDX = 1, THIS%DATA_COUNT
            ! Open file in subdirectory
            PATH = THIS%SUB_FILE_PATHS(IDX)%S
            VH5_CHECK(vh5_file_open_read(PATH, FILE_ID))
            ! Open group in file
            IS_OK = vh5_group_open(FILE_ID, 'data-1/'//NAME, GROUP_ID) == 0
            IF (.NOT. IS_OK) CALL MISSING_ERROR()
            ! Perform reading using handler
            CALL READER(GROUP_ID, IDX)
            VH5_CHECK(vh5_group_close(GROUP_ID))
            VH5_CHECK(vh5_file_close(FILE_ID))
         END DO
      ELSE
         ! Open main checkpoint file
         PATH = THIS%FILE_PATH
         VH5_CHECK(vh5_file_open_read(PATH, FILE_ID))
         DO IDX = 1, THIS%DATA_COUNT
            ! Open group for reading
            IS_OK = vh5_group_open(FILE_ID, 'data-'//str(IDX)//'/'//NAME, GROUP_ID) == 0
            IF (.NOT. IS_OK) CALL MISSING_ERROR()
            ! Perform reading using handler
            CALL READER(GROUP_ID, IDX)
            VH5_CHECK(vh5_group_close(GROUP_ID))
         ENDDO
         VH5_CHECK(vh5_file_close(FILE_ID))
      END IF

   CONTAINS
      SUBROUTINE MISSING_ERROR()
         CALL vtutor%error( &
            'While trying to read checkpoint: Data group "' // NAME // &
            '" is missing in checkpoint file "' // PATH // '". Most likely, ' // &
            'you are requesting a feature that was not checkpointed in the ' // &
            'previous run. Please check that your input settings are ' // &
            'consistent with the previous run or regenerate the checkpoint ' // &
            'file with your current settings.')
      END SUBROUTINE
#endif
   END SUBROUTINE

   !> @brief Writes metadata to the checkpoint file using a provided handler
   !>
   !> This subroutine writes metadata to the checkpoint file by calling a
   !> user-provided handler subroutine. The handler is responsible for writing
   !> the actual metadata to the provided group ID. This method handles the
   !> logic of creating a "metadata" group for writing and ensures that the
   !> checkpoint is in the correct state for writing metadata.
   SUBROUTINE CHECKPOINT_WRITE_METADATA(THIS, WRITER)
      CLASS(checkpoint_struct), INTENT(INOUT) :: THIS !< Checkpoint structure
      PROCEDURE(CHECKPOINT_HANDLER) :: WRITER !< Writer handler

      INTEGER(HID_T) :: FILE_ID, GROUP_ID

#ifdef VASP_HDF5
      IF (.NOT. THIS%IS_INITIALIZED) CALL vtutor%bug( &
         'While trying to write metadata: Checkpoint not initialized', __FILE__, __LINE__)
      IF (.NOT. THIS%ENABLED) RETURN
      IF (THIS%DATA_IS_OPEN) CALL vtutor%bug( &
         'While trying to write metadata: Data group is still open', __FILE__, __LINE__)
      IF (THIS%IS_READING) CALL vtutor%bug( &
         'While trying to write metadata: Checkpoint is in read mode', __FILE__, __LINE__)

      VH5_CHECK(vh5_file_open_readwrite(THIS%FILE_PATH, FILE_ID))
      VH5_CHECK(vh5_group_create(FILE_ID, 'metadata', GROUP_ID))
      CALL WRITER(GROUP_ID)
      VH5_CHECK(vh5_group_close_writing(GROUP_ID))
      VH5_CHECK(vh5_file_close_writing(FILE_ID))
#endif
   END SUBROUTINE

   !> @brief Reads metadata from the checkpoint file using a provided handler
   !>
   !> This subroutine reads metadata from the checkpoint file by calling a
   !> user-provided handler subroutine. The handler is responsible for reading
   !> the actual metadata from the provided group ID. This method handles the
   !> logic of opening the "metadata" group for reading and ensures that the
   !> checkpoint is in the correct state for reading metadata.
   !>
   !> The metadata group is expected to contain general information about the
   !> checkpointing setup, such as the total number of steps, subdirectory
   !> prefix, or other relevant information that is not specific to individual
   !> data groups.
   !>
   !> It is a good idea to check the integrity of the metadata upon reading,
   !> to ensure that, e.g., input parameters match between the current run
   !> and the checkpoint data.
   SUBROUTINE CHECKPOINT_READ_METADATA(THIS, READER)
      CLASS(checkpoint_struct), INTENT(INOUT) :: THIS !< Checkpoint structure
      PROCEDURE(CHECKPOINT_HANDLER) :: READER !< Reader handler

      INTEGER(HID_T) :: FILE_ID, GROUP_ID

#ifdef VASP_HDF5
      IF (.NOT. THIS%IS_INITIALIZED) CALL vtutor%bug( &
         'While trying to read metadata: Checkpoint not initialized', __FILE__, __LINE__)
      IF (.NOT. THIS%ENABLED) RETURN
      IF (THIS%DATA_IS_OPEN) CALL vtutor%bug( &
         'While trying to read metadata: Data group is still open', __FILE__, __LINE__)
      IF (.NOT. THIS%IS_READING) CALL vtutor%bug( &
         'While trying to read metadata: Checkpoint is not in read mode', __FILE__, __LINE__)

      VH5_CHECK(vh5_file_open_read(THIS%FILE_PATH, FILE_ID))
      VH5_CHECK(vh5_group_open(FILE_ID, 'metadata', GROUP_ID))
      CALL READER(GROUP_ID)
      VH5_CHECK(vh5_group_close(GROUP_ID))
      VH5_CHECK(vh5_file_close(FILE_ID))
#endif
   END SUBROUTINE

   !> @brief Writes subdirectory metadata to the checkpoint file
   !>
   !> This subroutine writes metadata related to subdirectory checkpointing
   !> to the checkpoint file. It takes the subdirectory prefix and the total
   !> count of subdirectories as input and writes them to the checkpoint file.
   !> This information is crucial for correctly reading the checkpoint data
   !> when combining results from multiple subdirectory files.
   SUBROUTINE CHECKPOINT_WRITE_SUBDIR_METADATA(THIS, PREFIX, TOTAL_COUNT)
      CLASS(checkpoint_struct), INTENT(INOUT) :: THIS !< Checkpoint structure
      CHARACTER(*), INTENT(IN) :: PREFIX !< Prefix for subdirectories
      INTEGER, INTENT(IN) :: TOTAL_COUNT !< Total number of subdirectories

      INTEGER(HID_T) :: FILE_ID

#ifdef VASP_HDF5
      IF (.NOT. THIS%IS_INITIALIZED) CALL vtutor%bug( &
         'While writing subdir metadata: Checkpoint not initialized', __FILE__, __LINE__)
      IF (.NOT. THIS%ENABLED) RETURN
      IF (THIS%DATA_IS_OPEN) CALL vtutor%bug( &
         'While writing subdir metadata: Data group is open', __FILE__, __LINE__)
      IF (THIS%IS_READING) CALL vtutor%bug( &
         'While writing subdir metadata: Checkpoint is in read mode', __FILE__, __LINE__)

      VH5_CHECK(vh5_file_open_readwrite(THIS%FILE_PATH, FILE_ID))
      VH5_CHECK(vh5_write(FILE_ID, 'subdir_prefix', PREFIX))
      VH5_CHECK(vh5_write(FILE_ID, 'total_count', TOTAL_COUNT))
      VH5_CHECK(vh5_file_close_writing(FILE_ID))
#endif
   END SUBROUTINE

   !> @brief Checks if the checkpoint is in reading mode
   !>
   !> You should use this to check if you should read metadata or actual
   !> data from the checkpoint.
   PURE LOGICAL FUNCTION CHECKPOINT_SHOULD_READ(THIS) RESULT(SHOULD)
      CLASS(checkpoint_struct), INTENT(IN) :: THIS
      SHOULD = THIS%IS_READING
   END FUNCTION

   !> @brief Checks if the checkpoint is using subdirectories for data groups
   !>
   !> You can use this to, e.g., provide output messages that are consistent
   !> with the checkpointing scheme being used. It is not necessary to check
   !> this when reading or writing checkpoint data, as logic is handled
   !> internally in the checkpoint methods.
   PURE LOGICAL FUNCTION CHECKPOINT_USES_SUBDIRS(THIS) RESULT(USES)
      CLASS(checkpoint_struct), INTENT(IN) :: THIS
      USES = THIS%USE_SUBDIRS
   END FUNCTION

   !> @brief Determines the complexity level of the executable based on input parameters
   !>
   !> The complexity level is determined based on whether we are in real mode and
   !> whether we are in gamma-only real mode. The levels are defined as follows:
   !> - 1: gamma-only executable (vasp_gam)
   !> - 2: standard executable (vasp_std)
   !> - 3: noncollinear executable (vasp_ncl)
   INTEGER FUNCTION GET_EXE_COMPLEXITY(REAL_MODE, GAMMA_REAL) RESULT(COMPLEXITY)
      LOGICAL, INTENT(IN) :: REAL_MODE, GAMMA_REAL

      IF (.NOT. REAL_MODE) THEN
         COMPLEXITY = 3
      ELSE IF (GAMMA_REAL) THEN
         COMPLEXITY = 1
      ELSE
         COMPLEXITY = 2
      ENDIF
   END FUNCTION

   !> @brief Returns the name of the executable corresponding to a given complexity level
   !>
   !> The names are defined as follows:
   !> - 1: gamma-only executable (vasp_gam)
   !> - 2: standard executable (vasp_std)
   !> - 3: noncollinear executable (vasp_ncl)
   FUNCTION GET_EXE_NAME(COMPLEXITY) RESULT(NAME)
      INTEGER, INTENT(IN) :: COMPLEXITY
      CHARACTER(:), ALLOCATABLE :: NAME

      SELECT CASE (COMPLEXITY)
      CASE (1)
         NAME = 'gamma-only executable (vasp_gam)'
      CASE (2)
         NAME = 'standard executable (vasp_std)'
      CASE (3)
         NAME = 'noncollinear executable (vasp_ncl)'
      CASE DEFAULT
         CALL vtutor%bug('Invalid complexity level: ' // str(COMPLEXITY), __FILE__, __LINE__)
      END SELECT
   END FUNCTION

END MODULE checkpointing
