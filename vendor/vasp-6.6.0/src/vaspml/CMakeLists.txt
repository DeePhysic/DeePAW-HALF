
cmake_minimum_required( VERSION 3.0 )

project( vasp_ml )


## define language and standards
set(CMAKE_CXX_STANDARD 17)
set( CMAKE_Fortran_COMPILER gfortran )
set( CMAKE_CPP_COMPILER mpic++)
find_package(MPI REQUIRED)
enable_language(Fortran)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)





## define flags for release and debugg... and so on
## define those things for bothh c++ and fortran
#if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
#    set(CMAKE_Fortran_FLAGS         "-O3")
#    set(CMAKE_Fortran_FLAGS_DEBUG   "-g -traceback")
#endif()
#
#if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
#    set( CMAKE_Fortran_FLAGS         "-O3 -fopenmp -Wunused" )
#    set( CMAKE_Fortran_FLAGS_DEBUG   "-O3 -fopenmp -g3 -fbounds-check -Wmaybe-uninitialized -Wuninitialized" )
#    set( CMAKE_Fortran_FLAGS_RELEASE "-O3 -fopenmp" )
#endif()








# switch on compiler optimization if no build type is given
if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Release )
endif()



## set the release flags
set( CMAKE_CXX_FLAGS_RELEASE "-O3 -lgfortran" )
set( CMAKE_Fortran_FLAGS_RELEASE "-O3 -ffree-form" )

# build special function and interpolation library
add_subdirectory( src/libvaspml_f )
include_directories( src/libvaspml_f )

#adding preprocessor tags
add_compile_definitions( VASPML_USE_CBLAS VASPML_DEBUG_LEVEL=3 )



## include eigen directory
include_directories( Eigen "/fsc/home/jona/source/eigen-3.4.0" )
include_directories( libvasp_ml "src/libvaspml/" )
include_directories( boost "/fsc/home/jona/source/boost_1_81_0" )
include_directories( cblas $ENV{OPENBLAS_ROOT}/include )
link_directories( openblas $ENV{OPENBLAS_ROOT}/lib )


file(GLOB SOURCES
    src/libvaspml/*.cpp
)

set(INCS fortran_math fortran_asa Eigen libvasp_ml clas_inc)
set(LLIBS fortran_math fortran_asa openblas)

#compile a file
function(compile testCase path)
   add_executable( ${testCase} ${path}${testCase}.cpp ${SOURCES})
   target_include_directories( ${testCase} PRIVATE ${INCS})
   target_link_libraries( ${testCase} ${LLIBS} MPI::MPI_CXX)
endfunction()

# Define a function to extract filename without path and extension
function(EXTRACT_FILENAME INPUT_STRING OUTPUT_VARIABLE)
    # Remove the path
    #    string(REPLACE "./src/applications/" "" TEMP_STRING "${INPUT_STRING}")
    # Find the last occurrence of the substring "/"
    string(FIND "${INPUT_STRING}" "/" LAST_SLASH_REVERSE_OFFSET REVERSE)
    # Extract the substring after the last occurrence of "/"
    math(EXPR LAST_SLASH_REVERSE_OFFSET "${LAST_SLASH_REVERSE_OFFSET} + 1")
    string(SUBSTRING "${INPUT_STRING}" ${LAST_SLASH_REVERSE_OFFSET} -1 TEMP_STRING)
    # Remove the extension
    string(REPLACE ".cpp" "" TEMP_STRING "${TEMP_STRING}")
    # Set the output variable
    set(${OUTPUT_VARIABLE} "${TEMP_STRING}" PARENT_SCOPE)
endfunction()


function(compile_list path)
  file(GLOB file_list "${path}*.cpp")
  foreach( file ${file_list} )
     extract_filename("${file}" RESULT)
     # Output the result
     message("Compiling ${RESULT}")
     compile( "${RESULT}" "${path}")
  endforeach()
endfunction()


## add libraries to link against here
## this creates the main executable
add_executable( vasp_ml ./src/applications/test_energyPrediction.cpp ${SOURCES})
target_include_directories( vasp_ml PRIVATE ${INCS})
target_link_libraries( vasp_ml ${LLIBS} MPI::MPI_CXX)



set( path "./src/applications/" )
compile_list( "${path}" )



# fill proper vasp colors here
add_custom_command( TARGET vasp_ml POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --blue
                    " _    _____   _____ ____  __  _____  " 
                    COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --blue
                    "| |  / /   | / ___// __ \\/  |/  / /  "
                    COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan
                    "| | / / /| | \\__ \\/ /_/ / /|_/ / /   "
                    COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --green
                    "| |/ / ___ |___/ / ____/ /  / / /___ "
                    COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --red
                    "|___/_/  |_/____/_/   /_/  /_/_____/"
                    VERBATIM )
