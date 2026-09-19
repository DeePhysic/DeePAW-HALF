#ifndef SAMPLETRAININGDATA_HPP
#define SAMPLETRAININGDATA_HPP

#include "SampleFile.hpp"

#include "types.hpp"

namespace vaspml
{

struct Record;

struct SampleTrainingData : public SampleFile
{
    SampleTrainingData(String id);
    Record& load();

    String mode;
    /**********************************************************************************************
     * The actual data set (enclosed with shared_ptr).
     **********************************************************************************************/
    ShRec set;
};

} //namespace vaspml

#endif
