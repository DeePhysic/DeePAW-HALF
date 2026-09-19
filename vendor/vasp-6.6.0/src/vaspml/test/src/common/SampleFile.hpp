#ifndef SAMPLEFILE_HPP
#define SAMPLEFILE_HPP

#include "SemanticVersion.hpp"
#include "types.hpp"

namespace vaspml
{

struct SampleFile
{
    SampleFile(String group, String id);

    /**********************************************************************************************
     * Relative path from test execution directory to "data" directory.
     **********************************************************************************************/
    String          pathPrefix;
    /**********************************************************************************************
     * Group to which the file belongs, specified as part of the full path below the "data"
     * directory before the version number.
     *
     * Examples:
     * * Force field files: "ff"
     * * Training data files: "sets"
     **********************************************************************************************/
    String          group;
    /**********************************************************************************************
     * The remainder of the file path after the group specifier.
     *
     * The first subdirectory specifies the VASP version which created the file.
     *
     * Examples:
     * * Force field: "6.4.1/refit/ML_FF.MAPbI3"
     * * Training data: "6.5.1/train/ML_AB.MAPbI3"
     *
     * Everything below the version number depends on the group, i.e., is handled by the subclasses.
     **********************************************************************************************/
    String          id;
    /**********************************************************************************************
     * Full path of the file relative to the test execution directory.
     **********************************************************************************************/
    String          path;
    String          system;
    SemanticVersion version;
};

} //namespace vaspml

#endif
