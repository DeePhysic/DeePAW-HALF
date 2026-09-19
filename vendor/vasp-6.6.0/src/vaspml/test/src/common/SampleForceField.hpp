#ifndef SAMPLEFORCEFIELD_HPP
#define SAMPLEFORCEFIELD_HPP

#include "SampleFile.hpp"
#include "types.hpp"

namespace vaspml
{

struct Record;

struct SampleForceField : public SampleFile
{
    SampleForceField(String id);
    Record& load();

    String mode;
    bool   fast;
    /**********************************************************************************************
     * The actual force field record (enclosed with shared_ptr).
     **********************************************************************************************/
    ShRec ff;
};

} //namespace vaspml

#endif
