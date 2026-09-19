! Fortran "template" workaround from incar_reader.inc.
#if RANK == 0
#   define RANK_SUFFIX
#elif RANK == 1
#   define RANK_SUFFIX (:)
#elif RANK == 2
#   define RANK_SUFFIX (:,:)
#elif RANK == 3
#   define RANK_SUFFIX (:,:,:)
#endif

#if defined(TYPE_INTEGER)
#   define TYPE_VAL INTEGER
#elif defined(TYPE_LOGICAL)
#   define TYPE_VAL LOGICAL
#elif defined(TYPE_REAL)
#   define TYPE_VAL REAL(q)
#elif defined(TYPE_CHARACTER)
#   define TYPE_VAL CHARACTER (len=:), ALLOCATABLE
#endif

#ifdef INCLUDE_READ_TAG_BODY
!==============================================================================
!> Search and assign given ML tag (+ alternatives) in contents of INCAR file.
!>
!> Example usage to get a READ_TAG function for 1-dimensional integer arrays:
!> ~~~
!>     #define INCLUDE_READ_TAG_BODY
!>     #define TYPE_REAL
!>     #define RANK 1
!>        INTEGER FUNCTION READ_TAG_REAL1 &
!>     #include "ml_interface.inc"
!>        END FUNCTION READ_TAG_REAL1
!>     #undef RANK
!>     #undef TYPE_REAL
!>     #undef INCLUDE_READ_TAG_BODY
!> ~~~
!> @attention Don't forget the '&' at the end of the function definition!
!>
!> @returns NARGS ... Number of elements found for the given tag.
!==============================================================================
!  INTEGER FUNCTION READ_TAG_...       &
                            (INCAR,    &
                             TAG_LIST, &
                             PREFIX,   &
                             TAG,      &
                             VARIABLE, &
                             DEFAULT_VALUE) RESULT(NARGS)
      USE incar_reader, only: INCAR_FILE, PROCESS_INCAR, COUNT_ELEMENTS
      USE ml_ff_taglist
#if defined(TYPE_CHARACTER)
      USE ml_ff_string, only: LOWERCASE
#endif
      IMPLICIT NONE

      !> String representation of INCAR file.
      TYPE(INCAR_FILE), INTENT(INOUT) :: INCAR
      !> List with all tag states.
      TYPE(TAGLIST),    INTENT(INOUT) :: TAG_LIST
      !> INCAR tag prefix (in case of multiple ML instances).
      CHARACTER(len=*), INTENT(IN)    :: PREFIX
      !> INCAR tag to be processed.
      CHARACTER(len=*), INTENT(IN)    :: TAG
      !> Variable where to store default or INCAR tag argument.
      TYPE_VAL,          INTENT(INOUT) :: VARIABLE RANK_SUFFIX
      !> Default value for variable.
      TYPE_VAL,          INTENT(IN)    :: DEFAULT_VALUE

      ! Alternative tags
      CHARACTER(len=:), ALLOCATABLE   :: ALT_TAGS(:)
      INTEGER                         :: I

      ! Set default value.
      VARIABLE = DEFAULT_VALUE
      CALL TAG_LIST%SET_STATE(TAG, TL_DEFAULT)

      ! Collect alternative tag list.
      ALT_TAGS = TAGLIST_GET_ALT_TAGS(TAG)

      NARGS = 0
      ! Loop over alternative tags (includes current "preferred" tag as
      ! first entry).
      DO I = 1, SIZE(ALT_TAGS)
         ! Check for per-instance INCAR tag first.
         NARGS = COUNT_ELEMENTS(INCAR, ADJUSTL(TRIM(PREFIX // ALT_TAGS(I))))
         IF (NARGS > 0) THEN
            CALL PROCESS_INCAR(INCAR,                                &
                               ADJUSTL(TRIM(PREFIX // ALT_TAGS(I))), &
                               VARIABLE)
            EXIT ! Skip alternative tags if one is already found.
         END IF

         ! Now comes the same check without prefix.
         NARGS = COUNT_ELEMENTS(INCAR, ADJUSTL(TRIM(ALT_TAGS(I))))
         IF (NARGS > 0) THEN
            CALL PROCESS_INCAR(INCAR, ADJUSTL(TRIM(ALT_TAGS(I))), VARIABLE)
            EXIT ! Skip alternative tags if one is already found.
         END IF
      END DO

#if defined(TYPE_CHARACTER)
      ! Automatically lower the case of all string-type values.
      VARIABLE = LOWERCASE(VARIABLE)
#endif

      ! If found anywhere, set tag list state accordingly.
      IF (NARGS > 0) THEN
         IF (I > 1) THEN
            CALL TAG_LIST%SET_STATE(TAG, TL_INCAR_ALT, I)
         ELSE
            CALL TAG_LIST%SET_STATE(TAG, TL_INCAR, I)
         END IF
      END IF

!  END FUNCTION READ_TAG_...
#endif

#undef TYPE_VAL
