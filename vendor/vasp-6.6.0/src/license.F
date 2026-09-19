#include "symbol.inc"
!> This module provides functionality to check the validity of a license key.
!>
!> The license key is read from a file and parsed into a structured format. The
!> license key contains information about the user, the license id, the valid VASP
!> version, and the expiry dates. The module provides functions to check the
!> validity of the license key based on the current date, the VASP version, and
!> the expected license id and user id.
!>
#ifdef LICENSE
module license_check

    use, intrinsic :: iso_c_binding, only: c_char, c_loc, c_int, c_long, c_null_char, c_ptr
    use version, only: date_t

    implicit none

    integer, parameter :: not_set = -1
    type(date_t), parameter :: not_set_date = date_t(year = not_set, month = not_set, day = not_set)
    character(len=*), parameter :: invalid_user_name = "invalid user name"
    character(len=*), parameter :: invalid_user_id = "invalid user id"
    character(len=*), parameter :: hash_salt = "sfKDSLtv5Q#p"  ! salt used when hashing user name
    character(len=*), parameter :: license_id = LICENSE  ! this variable is replace by the preprocessor
    character(len=*), parameter :: revoked_keys_path = REVOKED_KEYS_PATH  ! this variable is replaced by the preprocessor

    type version_number
        integer :: major = not_set
        integer :: minor = not_set
        integer :: patch = not_set
    end type version_number

    type license_data
        character(len=:), allocatable :: user_name
        !< The name of the user
        character(len=:), allocatable :: user_hash
        !< The hash of the user data
        character(len=:), allocatable :: license_id
        !< The id of the license e.g. of the HPC center
        character(len=:), allocatable :: license_text
        !< The license text without the signature
        character(len=:), allocatable :: signature
        !< The signature of the license text (base64 encoded)
        type(version_number) valid_vasp_version
        !< The version restriction for VASP (major, minor, patch)
        type(date_t) :: key_expiry = not_set_date
        !< The expiry date of the license key
        type(date_t) :: license_expiry = not_set_date
        !< The expiry date of the license
    contains
        procedure :: check_license_signature
        procedure :: check_license_date
        procedure :: check_license_version
        procedure :: check_license_id
        procedure :: check_user_hash
        procedure :: check_revoked_keys_signature
        procedure :: check_revoked_keys
    end type license_data

    type key_data
        character(len=:), allocatable :: license_id
        !< The id of the license e.g. of the HPC center
        character(len=:), allocatable :: list_hashes
        !< A string containing all hashes that are concerned
        type(date_t) :: created_at = not_set_date
        !< The date when the information about the keys was generated
        character(len=:), allocatable :: content
        !< The content of the key data file without the signature
        character(len=:), allocatable :: signature
        !< The signature of the license text (base64 encoded) to check the content
    end type key_data

    !> A struct to hold a C-style string with its length
    type, bind(c) :: string_with_length
        type(c_ptr) :: text
        integer(c_int) :: length
    end type string_with_length

    !> A struct to hold a Fortran string and its corresponding C-style string
    type string_with_buffer
        type(string_with_length) :: content
        character(len=:, kind=c_char), allocatable :: buffer
    end type string_with_buffer

    enum, bind(c)
        enumerator :: LICENSE_OK = 0
        enumerator :: LICENSE_WARNING_KEY_ALMOST_EXPIRED
        enumerator :: LICENSE_ERROR_KEY_NOT_FOUND
        enumerator :: LICENSE_ERROR_KEY_READING
        enumerator :: LICENSE_ERROR_KEY_PARSING
        enumerator :: LICENSE_ERROR_KEY_INVALID
        enumerator :: LICENSE_ERROR_KEY_EXPIRED
        enumerator :: LICENSE_ERROR_KEY_REVOKED
        enumerator :: LICENSE_ERROR_KEY_ID
        enumerator :: LICENSE_ERROR_REVOKED_READING
        enumerator :: LICENSE_ERROR_REVOKED_PARSING
        enumerator :: LICENSE_ERROR_REVOKED_INVALID
        enumerator :: LICENSE_ERROR_REVOKED_EXPIRED
        enumerator :: LICENSE_ERROR_REVOKED_ID
        enumerator :: LICENSE_ERROR_LICENSE_EXPIRED
        enumerator :: LICENSE_ERROR_VERSION
        enumerator :: LICENSE_ERROR_USER_NAME
        enumerator :: LICENSE_ERROR_USER_ID
        enumerator :: LICENSE_ERROR_USER_HASH
        enumerator :: LICENSE_ERROR_DATE
    end enum

    interface operator(<)
        module procedure earlier_date, earlier_version
    end interface

    interface
        function verify_signature(message, signature) bind(c, name="verify_signature")
            import :: c_int, string_with_length
            integer(c_int) :: verify_signature
            type(string_with_length), intent(in) :: message, signature
        end function verify_signature

        function get_user_name_c(user_name) bind(c, name="get_user_name")
            import :: c_int, string_with_length
            integer(c_int) :: get_user_name_c
            type(string_with_length), intent(out) :: user_name
        end function get_user_name_c

        function get_user_id_c() bind(c, name="get_user_id")
            import :: c_long
            integer(c_long) :: get_user_id_c
        end function get_user_id_c

        function get_current_epoch() bind(c, name="get_current_epoch")
            import :: c_long
            integer(c_long) :: get_current_epoch
        end function get_current_epoch

        subroutine hash_user_data_c(user_data, user_hash) bind(c, name="hash_user_data")
            import :: string_with_length
            type(string_with_length), intent(in) :: user_data
            type(string_with_length), intent(out) :: user_hash
        end subroutine hash_user_data_c
    end interface

contains

    !> Check the license and stop the program if the license is invalid.
    !>
    !> This subroutine checks the license on the I/O node and broadcasts the result
    !> to all other nodes. If the license is invalid, an appropriate error message
    !> is printed and the program is terminated.
    subroutine check_license_and_stop_if_invalid(io, comm, release_date)
        use base, only: in_struct
        use mpimy
        use tutor, only: vtutor
        type(in_struct), intent(in) :: io
        !< I/O structure for printing messages
        type(communic), intent(in) :: comm
        !< Communicator structure for MPI communication
        type(date_t) :: release_date
        !< The release date of this version of VASP.
        character(len=:), allocatable :: additional_info
        integer error
        !
        error = LICENSE_OK
        if (comm%ionode == comm%node_me) then
            write(io%iu0,'(1x,a)', advance="no") "Checking license..."
            error = conduct_license_check(release_date, additional_info)
            if (error == LICENSE_OK .OR. error == LICENSE_WARNING_KEY_ALMOST_EXPIRED) then
                write(io%iu0,'(a,/)') " OK"
            else
                write(io%iu0,'(a,/)') " FAILED"
            end if
        else
            additional_info = ""
        end if
        CALLMPI( M_bcast_i(comm, error, 1) )
        select case (error)
        case (LICENSE_OK)
            return
        case (LICENSE_WARNING_KEY_ALMOST_EXPIRED)
            call vtutor%alert("VASP license keys are time limited to ensure that the data in the &
                &VASP portal is accurately reflected. Your license key will expire on " // &
                additional_info // ". To avoid calculations stopping after that day, please &
                &use the VASP portal to update your license key as soon as possible.")
        case (LICENSE_ERROR_KEY_NOT_FOUND)
            call vtutor%error("No license file found. Please provide a valid license file. &
                &VASP will look in the following places:\n&
                & * The path set by the environment variable VASP_LICENSE_FILE.\n&
                & * The file `vasp_license` in the current working directory.\n&
                & * The file `vasp_license` in $HOME/.vasp/.")
        case (LICENSE_ERROR_KEY_READING)
            call vtutor%error("Could not read the license file. Please check the permissions of file " &
                // additional_info // ".")
        case (LICENSE_ERROR_KEY_PARSING)
            call vtutor%error("Could not parse the license file. VASP read the following content:\n" &
                // additional_info // "\nPlease check that this contains user_name, user_hash, &
                &license_id, valid_vasp_version, key_expiry, license_expiry, and signature.")
        case (LICENSE_ERROR_KEY_INVALID)
            call vtutor%error("The signature of the license is invalid. This indicates that the &
                &license file has been modified after its creation. Please create a new license " &
                // "and replace the one you are currently using.")
        case (LICENSE_ERROR_KEY_EXPIRED)
            call vtutor%error("VASP license keys are time limited to ensure that the data in the &
                &VASP portal is accurately reflected. Since your license key expired, VASP will &
                &stop the execution. If you still have an active VASP license all you need to do &
                &is to create a new license key with an updated expiry date.")
        case (LICENSE_ERROR_KEY_REVOKED)
            call vtutor%error("VASP license key revoked! The license key currently in use (" // &
                additional_info // ") is present on this site's revocation list. To resume &
                &calculations, please provide a valid, active license key. If you require &
                &assistance or did not authorize this revocation, please contact your site &
                &administrator or reach out to us at https://vasp.at/forum.")
        case (LICENSE_ERROR_KEY_ID)
            call vtutor%error("The license is not valid for this site. " // additional_info &
                // " Please generate a specific license key for this site.")
        case (LICENSE_ERROR_REVOKED_READING)
            call vtutor%error("Could not read the file of revoked keys. VASP tried to read the &
                &following file:\n" // revoked_keys_path // "\nPlease make sure this file exists &
                &and is readable.")
        case (LICENSE_ERROR_REVOKED_PARSING)
            call vtutor%error("Could not parse the file of revoked keys. VASP tried to read the &
                &following file:\n" // revoked_keys_path // "\nPlease make sure this file contains &
                &license_id, revoked_keys, created_at, and signature.")
        case (LICENSE_ERROR_REVOKED_INVALID)
            call vtutor%error("The signature of the revoked keys file is invalid. This indicates &
                &that the file has been modified after its creation or the public key is outdated. &
                &Please get in touch with the person that installed VASP and report the issue in &
                &the forum (https://vasp.at/forum).")
        case (LICENSE_ERROR_REVOKED_EXPIRED)
            call vtutor%error("The file of revoked keys is outdated. It must not be older than &
                &one month. Please ask your admin to keep this file up to date with a cronjob.")
        case (LICENSE_ERROR_REVOKED_ID)
            call vtutor%error("VASP was setup inconsistently. " // additional_info &
                // " Please ask your system administrator to fix this.")
        case (LICENSE_ERROR_LICENSE_EXPIRED)
            call vtutor%error("Your VASP license expired before the release date of this version &
                &of VASP. You are may not use this version of VASP unless you update your VASP &
                &license. Please get in touch with licensing@vasp.at. If you have already updated &
                &your VASP license, your license key does not reflect this yet. Please use the &
                &VASP portal to create a new license key.")
        case (LICENSE_ERROR_VERSION)
            call vtutor%error("The license is not valid for this version of VASP. Please update &
                &the license key or renew your VASP license.")
        case (LICENSE_ERROR_USER_NAME)
            call vtutor%error("Could not determine the user name on this system.")
        case (LICENSE_ERROR_USER_ID)
            call vtutor%error("Could not determine the user id on this system.")
        case (LICENSE_ERROR_USER_HASH)
            call vtutor%error("The license is not valid for this user. This error can occur &
                &if the user name or id does not agree with the settings provided during &
                &generation of the key.\nVASP determined " // additional_info // ". Please &
                &generate a new key with this data.")
        case (LICENSE_ERROR_DATE)
            call vtutor%error("Could not determine the date on this system.")
        case default
            call vtutor%error("An unknown error occurred during license check.")
        end select
        !
    end subroutine check_license_and_stop_if_invalid

    !> Conduct the license check.
    !>
    !> This function reads the license key from a file, parses it, and checks its
    !> validity. If the license is valid, LICENSE_OK is returned. Otherwise, an
    !> error code is returned.
    integer function conduct_license_check(release_date, additional_info) result (error)
        use string, only: str, string_from_file
        type(date_t), intent(in) :: release_date
        !< The release date of this version of VASP.
        character(len=:), allocatable, intent(out) :: additional_info
        !< Additional information useful for constructing error messages
        character(len=:), allocatable :: license_file, license_content, revoked_keys_content
        character(len=:), allocatable :: user_name, user_id, user_hash
        type(license_data) license
        type(key_data) revoked_keys
        !
        error = find_license_key(license_file)
        if (error /= LICENSE_OK) return
        !
        license_content = string_from_file(license_file, error)
        if (error /= 0) then
            error = LICENSE_ERROR_KEY_READING
            additional_info = license_file
            return
        end if
        !
        additional_info = license_content
        error = parse_license_key(license_content, license)
        if (error /= LICENSE_OK) return
        !
        revoked_keys_content = string_from_file(revoked_keys_path, error)
        if (error /= 0) then
            error = LICENSE_ERROR_REVOKED_READING
            additional_info = revoked_keys_path
            return
        end if
        !
        error = parse_revoked_keys(revoked_keys_content, revoked_keys)
        if (error /= LICENSE_OK) return
        !
        error = get_user_name(user_name)
        if (error /= LICENSE_OK) return
        !
        error = get_user_id(user_id)
        if (error /= LICENSE_OK) return
        !
        user_hash = hash_user_data(user_name, user_id, license%license_id)
        !
        error = check_license(license, user_hash, release_date, revoked_keys)
        if (error == LICENSE_ERROR_KEY_ID) then
            additional_info = "This license is valid for `" // license%license_id &
                // "` and not for `" // license_id // "`."
        else if (error == LICENSE_ERROR_KEY_REVOKED) then
            additional_info = license_file
        else if (error == LICENSE_ERROR_REVOKED_ID) then
            additional_info = "This license is compiled with `" // license_id &
                // "` but the revoked keys are for `" // revoked_keys%license_id // "`."
        else if (error == LICENSE_ERROR_USER_HASH) then
            additional_info = "user_name: " // user_name // ", user_id: " // user_id
        else if (error == LICENSE_WARNING_KEY_ALMOST_EXPIRED) then
            additional_info = str(license%key_expiry%year) // "-" // str(license%key_expiry%month, "(i2.2)") &
                // "-" // str(license%key_expiry%day, "(i2.2)")
        end if
        !
    end function conduct_license_check

    !> Find the license key file.
    !>
    !> The function first checks the environment variable VASP_LICENSE_FILE.
    !> If it is set and points to an existing file, that file is used.
    !> Otherwise, the function checks for a file named "vasp_license" in the
    !> current working directory. If that file exists, it is used.
    !> Finally, the function checks for a file named "vasp_license" in the
    !> directory "$HOME/.vasp/". If that file exists, it is used.
    !> If none of these files exist, the function returns an empty string.
    integer function find_license_key(filepath) result (error)
        character(len=*), parameter :: default_filename = "vasp_license"
        character(len=:), allocatable, intent(out) :: filepath
        character(len=512) :: path
        logical file_exists
        !
        error = LICENSE_OK
        call get_environment_variable("VASP_LICENSE_FILE", path)
        if (len_trim(path) > 0) then
            filepath = trim(path)
            return
        end if
        !
        inquire(file=default_filename, exist=file_exists)
        if (file_exists) then
            filepath = default_filename
            return
        end if
        !
        call get_environment_variable("HOME", path)
        if (len_trim(path) > 0) then
            filepath = trim(path) // "/.vasp/" // default_filename
            inquire(file=filepath, exist=file_exists)
            if (file_exists) return
        end if
        !
        filepath = ""
        error = LICENSE_ERROR_KEY_NOT_FOUND
        !
    end function find_license_key

    !> Parse a license key from a string.
    !>
    !> The license key is expected to be in a format similar to the INCAR file.
    !> The following keys are expected:
    !> - user_name: Name of the user
    !> - user_hash: Hash of the user data
    !> - license_id: License of the HPC center
    !> - signature: Signature string (base64 encoded)
    !> - valid_vasp_version: Version restriction for VASP (major, minor, patch)
    !> - key_expiry: Expiry date of the license key (year, month, day)
    !> - license_expiry: Expiry date of the license (year, month, day)
    !>
    integer function parse_license_key(license_content, license) result (error)
        use incar_reader, only: incar_file, incar_from_string, process_incar
        character(len=*), intent(in) :: license_content
        !< The content of the license key as a string.
        type(license_data), intent(out) :: license
        !< The parsed license key data.
        !
        integer begin_signature
        type(incar_file) parsed_license
        !
        error = LICENSE_OK
        !
        ! parse signature
        begin_signature = index(license_content, "signature")
        if (begin_signature == 0) then
            error = LICENSE_ERROR_KEY_PARSING
            return
        end if
        parsed_license = incar_from_string(license_content(begin_signature:))
        call process_incar(parsed_license, "signature", license%signature)
        if (.not.allocated(license%signature)) error = LICENSE_ERROR_KEY_PARSING
        !
        ! parse remaining license content
        license%license_text = license_content(:begin_signature-1)
        parsed_license = incar_from_string(license%license_text)
        !
        call process_incar(parsed_license, "user_name", license%user_name)
        if (.not.allocated(license%user_name)) error = LICENSE_ERROR_KEY_PARSING
        call process_incar(parsed_license, "user_hash", license%user_hash)
        if (.not.allocated(license%user_hash)) error = LICENSE_ERROR_KEY_PARSING
        call process_incar(parsed_license, "license_id", license%license_id)
        if (.not.allocated(license%license_id)) error = LICENSE_ERROR_KEY_PARSING
        !
        call process_incar(parsed_license, "valid_vasp_version/major", license%valid_vasp_version%major)
        if (license%valid_vasp_version%major == not_set) error = LICENSE_ERROR_KEY_PARSING
        call process_incar(parsed_license, "valid_vasp_version/minor", license%valid_vasp_version%minor)
        if (license%valid_vasp_version%minor == not_set) error = LICENSE_ERROR_KEY_PARSING
        call process_incar(parsed_license, "valid_vasp_version/patch", license%valid_vasp_version%patch)
        if (license%valid_vasp_version%patch == not_set) error = LICENSE_ERROR_KEY_PARSING
        !
        call process_incar(parsed_license, "key_expiry/year", license%key_expiry%year)
        if (license%key_expiry%year == not_set) error = LICENSE_ERROR_KEY_PARSING
        call process_incar(parsed_license, "key_expiry/month", license%key_expiry%month)
        if (license%key_expiry%month == not_set) error = LICENSE_ERROR_KEY_PARSING
        call process_incar(parsed_license, "key_expiry/day", license%key_expiry%day)
        if (license%key_expiry%day == not_set) error = LICENSE_ERROR_KEY_PARSING
        !
        call process_incar(parsed_license, "license_expiry/year", license%license_expiry%year)
        if (license%license_expiry%year == not_set) error = LICENSE_ERROR_KEY_PARSING
        call process_incar(parsed_license, "license_expiry/month", license%license_expiry%month)
        if (license%license_expiry%month == not_set) error = LICENSE_ERROR_KEY_PARSING
        call process_incar(parsed_license, "license_expiry/day", license%license_expiry%day)
        if (license%license_expiry%day == not_set) error = LICENSE_ERROR_KEY_PARSING
        !
    end function parse_license_key

    !> Parse a revoked keys data from a string.
    !>
    !> The revoked keys data is expected to be in a format similar to the INCAR file.
    !> The following keys are expected:
    !> - license_id: License of the HPC center
    !> - revoked_keys: A string containing all revoked keys
    !> - signature: Signature string (base64 encoded)
    !> - create_at: Date when the key data was generated (year, month, day)
    !>
    integer function parse_revoked_keys(revoked_keys_content, revoked_keys) result (error)
        use incar_reader, only: incar_file, incar_from_string, process_incar
        character(len=*), intent(in) :: revoked_keys_content
        !< The content of the revoked keys file as a string.
        type(key_data), intent(out) :: revoked_keys
        !< The parsed revoked keys data.
        !
        integer begin_signature
        type(incar_file) parsed_key_data
        !
        error = LICENSE_OK
        !
        ! parse signature
        begin_signature = index(revoked_keys_content, "signature")
        if (begin_signature == 0) then
            error = LICENSE_ERROR_REVOKED_PARSING
            return
        end if
        parsed_key_data = incar_from_string(revoked_keys_content(begin_signature:))
        call process_incar(parsed_key_data, "signature", revoked_keys%signature)
        if (.not.allocated(revoked_keys%signature)) error = LICENSE_ERROR_REVOKED_PARSING
        !
        ! parse remaining license content
        revoked_keys%content = revoked_keys_content(:begin_signature-1)
        parsed_key_data = incar_from_string(revoked_keys%content)
        !
        call process_incar(parsed_key_data, "license_id", revoked_keys%license_id)
        if (.not.allocated(revoked_keys%license_id)) error = LICENSE_ERROR_REVOKED_PARSING
        call process_incar(parsed_key_data, "revoked_keys", revoked_keys%list_hashes)
        if (.not.allocated(revoked_keys%list_hashes)) error = LICENSE_ERROR_REVOKED_PARSING
        !
        call process_incar(parsed_key_data, "created_at/year", revoked_keys%created_at%year)
        if (revoked_keys%created_at%year == not_set) error = LICENSE_ERROR_REVOKED_PARSING
        call process_incar(parsed_key_data, "created_at/month", revoked_keys%created_at%month)
        if (revoked_keys%created_at%month == not_set) error = LICENSE_ERROR_REVOKED_PARSING
        call process_incar(parsed_key_data, "created_at/day", revoked_keys%created_at%day)
        if (revoked_keys%created_at%day == not_set) error = LICENSE_ERROR_REVOKED_PARSING
        !
    end function parse_revoked_keys

    !> Check whether the license is valid.
    !>
    !> This function checks the license signature, the expiry dates, the version
    !> restriction, the license id, and the user id. If any of these checks fail, an
    !> appropriate error code is returned. If all checks pass, LICENSE_OK is returned.
    !> The user_hash is passed as an argument to check against the license.
    integer function check_license(license, user_hash, release_date, revoked_keys) result (error)
        use version, only: major, minor, patch
        class(license_data), intent(in) :: license
        !< The license to check.
        character(len=*), intent(in) :: user_hash
        !< The user hash to check against.
        type(date_t), intent(in) :: release_date
        !< The release date of this version of VASP.
        type(key_data), intent(in) :: revoked_keys
        !< stores information which keys are revoked and should not be allowed to use VASP
        type(date_t) :: current_date
        character(len=8) :: date_string
        !
        error = license%check_license_signature()
        if (error /= LICENSE_OK) return
        !
        error = license%check_license_version(version_number(major, minor, patch))
        if (error /= LICENSE_OK) return
        !
        error = license%check_license_id(license_id)
        if (error /= LICENSE_OK) return
        !
        error = license%check_user_hash(user_hash)
        if (error /= LICENSE_OK) return
        !
        error = license%check_revoked_keys_signature(revoked_keys)
        if (error /= LICENSE_OK) return
        !
        error = license%check_revoked_keys(revoked_keys)
        if (error /= LICENSE_OK) return
        !
        ! we check the date last so that we can raise a warning for almost expired licenses
        error = get_current_date(current_date)
        if (error /= LICENSE_OK) return
        error = license%check_license_date(current_date, release_date, revoked_keys%created_at)
        !
    end function check_license

    !> Check whether the license signature is valid.
    integer function check_license_signature(license) result (error)
        class(license_data), intent(in) :: license
        !< check the signature of this license
        integer(c_int) :: openssl_result
        type(string_with_buffer) message, signature
        !
        error = LICENSE_OK
        call get_c_string(license%license_text, message)
        call get_c_string(license%signature, signature)
        openssl_result = verify_signature(message%content, signature%content)
        if (openssl_result /= 0) error = LICENSE_ERROR_KEY_INVALID
        !
    end function check_license_signature

    !> Check whether the signature of the revoked keys is valid
    integer function check_revoked_keys_signature(license, revoked_keys) result (error)
        class(license_data), intent(in) :: license
        !< license information of the user
        type(key_data), intent(in) :: revoked_keys
        !< check the signature of this type
        integer(c_int) :: openssl_result
        type(string_with_buffer) message, signature
        !
        error = LICENSE_OK
        call get_c_string(revoked_keys%content, message)
        call get_c_string(revoked_keys%signature, signature)
        openssl_result = verify_signature(message%content, signature%content)
        if (openssl_result /= 0) error = LICENSE_ERROR_REVOKED_INVALID
        !
    end function check_revoked_keys_signature

    !> Check whether the license date is valid
    !>
    !> For a license date to be valid, the following conditions must be fulfilled
    !> * the key must not have been expired yet
    !> * the license must be valid for longer than the release date of VASP
    !> * the revoked keys file must be recent (not older than a month)
    !> To counteract spoofing of the current date, the current date is increased
    !> to the day the revoked keys file was generated if that is larger
    integer function check_license_date(license, current_date_in, release_date, revoked_keys_date) result (error)
        class(license_data), intent(in) :: license
        !< check the date of this license
        type(date_t), intent(in) :: current_date_in
        !< the current date when running the code
        type(date_t), intent(in) :: release_date
        !< the release date when this version was published
        type(date_t), intent(in) :: revoked_keys_date
        !< the date when the revoked keys file was created
        type(date_t) current_date, in_a_month
        !
        error = LICENSE_OK
        if (revoked_keys_date < current_date_in) then
            current_date = current_date_in
        else
            current_date = revoked_keys_date
        end if
        if (license%key_expiry < next_month(current_date)) &
            error = LICENSE_WARNING_KEY_ALMOST_EXPIRED
        if (next_month(revoked_keys_date) < current_date) &
            error = LICENSE_ERROR_REVOKED_EXPIRED
        if (license%key_expiry < current_date) &
            error = LICENSE_ERROR_KEY_EXPIRED
        if (license%license_expiry < release_date) &
            error = LICENSE_ERROR_LICENSE_EXPIRED
        !
    end function check_license_date

    logical function earlier_date(left_date, right_date)
        type(date_t), intent(in) :: left_date, right_date
        earlier_date = left_date%year < right_date%year
        if (left_date%year == right_date%year) &
            earlier_date = left_date%month < right_date%month
        if (left_date%year == right_date%year .and. left_date%month == right_date%month) &
            earlier_date = left_date%day < right_date%day
    end function earlier_date

    type(date_t) function next_month(current_date)
        type(date_t), intent(in) :: current_date
        next_month = current_date
        if (next_month%month < 12) then
            next_month%month = next_month%month + 1
        else
            next_month%month = 1
            next_month%year = next_month%year + 1
        end if
    end function next_month

    !> Check whether the licensed version is valid, i.e., the version number of VASP is
    !> below the last version in the key.
    integer function check_license_version(license, current_version) result (error)
        class(license_data), intent(in) :: license
        !< check the version of this license
        type(version_number), intent(in) :: current_version
        !< the current version number of VASP
        !
        error = LICENSE_OK
        if (license%valid_vasp_version < current_version) &
            error = LICENSE_ERROR_VERSION
    end function check_license_version

    logical function earlier_version(left_version, right_version)
        type(version_number), intent(in) :: left_version, right_version
        earlier_version = left_version%major < right_version%major
        if (left_version%major == right_version%major) &
            earlier_version = left_version%minor < right_version%minor
        if (left_version%major == right_version%major .and. left_version%minor == right_version%minor) &
            earlier_version = left_version%patch < right_version%patch
    end function earlier_version

    !> Check whether the license has the same id like the code it was compiled with
    integer function check_license_id(license, expected_license_id) result (error)
        class(license_data), intent(in) :: license
        !< check the license id of this license
        character(len=*), intent(in) :: expected_license_id
        !< the license id the code was compiled with
        !
        error = LICENSE_OK
        if (license%license_id /= expected_license_id) &
            error = LICENSE_ERROR_KEY_ID
    end function check_license_id

    !> Check whether the license has the same user id like the expected one
    integer function check_user_hash(license, expected_user_hash) result (error)
        class(license_data), intent(in) :: license
        !< check the user hash of this license
        character(len=*), intent(in) :: expected_user_hash
        !< the expected user hash
        !
        error = LICENSE_OK
        if (license%user_hash /= expected_user_hash) &
            error = LICENSE_ERROR_USER_HASH
    end function check_user_hash

    !> Check whether the license id of the license and the revoked keys agree
    !> up to a suffix and the user_hash is not on the list of revoked keys
    integer function check_revoked_keys(license, revoked_keys) result (error)
        class(license_data), intent(in) :: license
        !< check the license id and user hash of this license
        type(key_data), intent(in) :: revoked_keys
        !< the data of the revoked user keys
        character(len=:), allocatable :: license_id
        integer index_
        !
        error = LICENSE_OK
        index_ = index(license%license_id, "@")
        if (index_ > 0) then
            ! trim suffix
            license_id = license%license_id(:index_-1)
        else
            license_id = license%license_id
        end if
        if (license_id /= revoked_keys%license_id) &
            error = LICENSE_ERROR_REVOKED_ID
        index_ = index(revoked_keys%list_hashes, license%user_hash)
        if (index_ > 0) &
            error = LICENSE_ERROR_KEY_REVOKED
    end function check_revoked_keys

    !> Get the user name from the system using a C function.
    !> If the user name cannot be obtained, return "invalid user name".
    integer function get_user_name(user_name) result (error)
        character(len=:), allocatable :: user_name
        type(string_with_buffer) :: user_name_c
        !
        call new_string_with_buffer(256, user_name_c)
        error = get_user_name_c(user_name_c%content)
        if (error == 0) then
            user_name = get_f_string(user_name_c)
        else
            error = LICENSE_ERROR_USER_NAME
            user_name = invalid_user_name
        end if
        !
    end function get_user_name

    !> Get the user id from the system using a C function.
    !> If the user id cannot be obtained, return "invalid user id".
    integer function get_user_id(user_id) result (error)
        use string, only: str
        character(len=:), allocatable :: user_id
        integer(c_long) user_id_c
        !
        user_id_c = get_user_id_c()
        if (user_id_c >= 0) then
            error = LICENSE_OK
            user_id = str(user_id_c)
        else
            error = LICENSE_ERROR_USER_ID
            user_id = invalid_user_id
        end if
        !
    end function get_user_id

    integer function get_current_date(current_date) result (error)
        type(date_t), intent(out) :: current_date
        integer(c_long) current_epoch
        integer days, days_in_year, days_in_month
        logical is_leap_year
        !
        current_epoch = get_current_epoch()
        if (current_epoch < 0) then
            error = LICENSE_ERROR_DATE
            return
        end if
        !
        error = LICENSE_OK
        days = current_epoch / 86400
        current_date = date_t(year=1970, month=1, day=1)
        is_leap_year = .false.
        do while (days > 0)
            if (is_leap_year) then
                days_in_year = 366
            else
                days_in_year = 365
            end if
            if (days >= days_in_year) then
                current_date%year = current_date%year + 1
                is_leap_year = (mod(current_date%year, 4) == 0) .and. &
                    (.not.mod(current_date%year, 100) == 0 .or. mod(current_date%year, 400) == 0)
                days = days - days_in_year
                cycle
            end if
            select case (current_date%month)
            case (1, 3, 5, 7, 8, 10, 12)
                days_in_month = 31
            case (4, 6, 9, 11)
                days_in_month = 30
            case (2)
                if (is_leap_year) then
                    days_in_month = 29
                else
                    days_in_month = 28
                end if
            end select
            if (days >= days_in_month) then
                current_date%month = current_date%month + 1
                days = days - days_in_month
            else
                current_date%day = current_date%day + days
                days = 0
            end if
        end do
        !
    end function get_current_date

    !> Hash the user data (user name, user id, and license id) using SHA-256.
    function hash_user_data(user_name, user_id, license_id) result (hash)
        character(len=*), intent(in) :: user_name, user_id, license_id
        character(len=:), allocatable :: user_string, hash
        type(string_with_buffer) user_data, hash_c
        integer(c_int) length
        !
        user_string = user_name // "-" // hash_salt // "-" // user_id // "-" // license_id
        call get_c_string(user_string, user_data)
        call new_string_with_buffer(256, hash_c)
        call hash_user_data_c(user_data%content, hash_c%content)
        hash = get_f_string(hash_c)
        !
    end function hash_user_data

    !> Create a string_with_length object with a given length.
    subroutine new_string_with_buffer(length, output)
        integer, intent(in) :: length
        type(string_with_buffer), intent(out), target :: output
        output%buffer = repeat("_", length) // c_null_char
        output%content%text = c_loc(output%buffer)
        output%content%length = length
    end subroutine new_string_with_buffer

    !> Convert a string_with_length object to a Fortran string.
    function get_f_string(c_string) result (f_string)
        type(string_with_buffer), intent(in) :: c_string
        character(len=:), allocatable :: f_string
        f_string = c_string%buffer(:c_string%content%length)
    end function get_f_string

    !> Convert a Fortran string to a C string (null-terminated).
    subroutine get_c_string(f_string, c_string)
        character(len=*), intent(in) :: f_string
        type(string_with_buffer), intent(out) :: c_string
        call new_string_with_buffer(len(f_string), c_string)
        c_string%buffer(:len(f_string)) = f_string
    end subroutine get_c_string

end module license_check

#else

module license_check

    implicit none

contains

    !> If the code is compiled without -DLICENSE, no license check is conducted
    subroutine check_license_and_stop_if_invalid(io, comm, release_date)
        use base, only: in_struct
        use mpimy
        use version, only: date_t
        type(in_struct), intent(in) :: io
        type(communic), intent(in) :: comm
        type(date_t), intent(in) :: release_date
    end subroutine check_license_and_stop_if_invalid

end module license_check
#endif
