#include "Frame.hpp"

#include "DescriptorSHS2.hpp"
#include "DescriptorSHS3.hpp"
#include "DescriptorSHS3ReducedLinElem.hpp"
#include "ParallelEnvironment.hpp"
#include "Record.hpp"
#include "Timer.hpp"
#include "TypeMap.hpp"
#include "setup.hpp"

using namespace vaspml;

MapString vaspml::dataMapFrame()
{
    return MapString{
        {"structure",             "ShRec"},
        {"2-body-neighbor-list",  "ShRec"},
        {"3-body-neighbor-list",  "ShRec"},
        {"2-body-SHS2",           "ShRec"},
        {"3-body-SHS2",           "ShRec"},
        {"3-body-SHS3",           "ShRec"},
        {"3-body-SHS3LinElement", "ShRec"}
    };
}

Frame::Frame(ShRec record) : data(assignOrMakeRecord(record, dataMapFrame()))
{}

void Frame::init(const Record& inputParameters, const BasisFunctionMap& basisFunctions)
{
    neighborLists = setup::createNeighborLists(
        {
            {"2-body", inputParameters.cget<Real>("ML_RCUT1")},
            {"3-body", inputParameters.cget<Real>("ML_RCUT2")}
    },
        data);

    structure = std::make_shared<Structure>(data->get<ShRec>("structure"));
    forceFieldTypes = std::make_shared<Vec1String>(inputParameters.cget<Vec1String>("types"));

    switch (inputParameters.cget<Int>("ML_DESC_TYPE"))
    {
    case 0:
        this->descriptor3BodyType = Descriptor3BodyType::StandardDescriptor;
        break;
    case 1:
        this->descriptor3BodyType = Descriptor3BodyType::LinearScalingDescriptor;
        break;
    }
    descriptors = setup::createDescriptors(inputParameters,
                                           *basisFunctions.at("3-body"),
                                           descriptor3BodyType,
                                           data);

    set_basisFunctionsSHS2(basisFunctions);

    return;
}

void Frame::initRefit(const Record&           incar,
                      const Vec1String&       types,
                      const BasisFunctionMap& basisFunctions)
{
    neighborLists = setup::createNeighborLists(
        {
            {"2-body", incar.cget<Real>("ML_RCUT1")},
            {"3-body", incar.cget<Real>("ML_RCUT2")}
    },
        data);

    forceFieldTypes = std::make_shared<Vec1String>(types);
    switch (incar.cget<Int>("ML_DESC_TYPE"))
    {
    case 0:
        this->descriptor3BodyType = Descriptor3BodyType::StandardDescriptor;
        break;
    case 1:
        this->descriptor3BodyType = Descriptor3BodyType::LinearScalingDescriptor;
        break;
    }
    descriptors = setup::createDescriptors(incar,
                                           types.size(),
                                           *basisFunctions.at("3-body"),
                                           descriptor3BodyType,
                                           data);

    set_basisFunctionsSHS2(basisFunctions);

    return;
}

void Frame::set_basisFunctionsSHS2(const BasisFunctionMap& basisFunctions)
{
    for (const auto& x : neighborLists)
    {
        String          descKey = "SHS2-" + x.first;
        DescriptorSHS2* desc = static_cast<DescriptorSHS2*>(descriptors[descKey].get());
        desc->set_basisFunctions(basisFunctions.at(x.first));
    }

    return;
}

void Frame::prepare_DescriptorSHS3()
{
    std::shared_ptr<DescriptorSHS2> desc2 =
        std::static_pointer_cast<DescriptorSHS2>(descriptors["SHS2-3-body"]);
    DescriptorSHS3* desc3 = static_cast<DescriptorSHS3*>(descriptors["SHS3-3-body"].get());
    desc3->set_parameters(desc2, *typeMap);
    if (!global::parallel.off()) desc3->make_lookUpTables(*typeMap);

    return;
}

void Frame::update()
{
    VASPML_PROFILING_START("Frame::DescriptorSHS2");
    for (auto& [key, item] : neighborLists)
    {
        String descKey = "SHS2-" + key;
        if (descriptors[descKey]->get_weight() <= 0) continue;
        item->compute_centralAtomIndexPerType();
        DescriptorSHS2* desc = static_cast<DescriptorSHS2*>(descriptors[descKey].get());
        desc->updatePairCoefficients(item);
        desc->computeVaspCoefficientsFromPairCoefficients();
        if (!global::parallel.off()) desc->resize_LookUpTablesForceCompute(*typeMap);
        desc->allocate_typeIndexArray(typeMap->countStructureTypes());
    }
    VASPML_PROFILING_STOP("Frame::DescriptorSHS2");
    VASPML_PROFILING_START("Frame::DescriptorSHS3");
    // compute SHS3 and SHS3LinElem descriptor
    std::shared_ptr<DescriptorSHS2> desc2 =
        std::static_pointer_cast<DescriptorSHS2>(descriptors["SHS2-3-body"]);
    if (descriptor3BodyType == Descriptor3BodyType::StandardDescriptor)
    {
        DescriptorSHS3* desc3 = static_cast<DescriptorSHS3*>(descriptors["SHS3-3-body"].get());
        desc3->computeSHS3(desc2, *typeMap);
    }
    else if (descriptor3BodyType == Descriptor3BodyType::LinearScalingDescriptor)
    {
        DescriptorSHS3ReducedLinElem* desc3 =
            static_cast<DescriptorSHS3ReducedLinElem*>(descriptors["SHS3-3-body"].get());
        desc3->computeSHS3(desc2, *typeMap);
    }
    VASPML_PROFILING_STOP("Frame::DescriptorSHS3");

    return;
}

void Frame::updateNeighborLists(const String& fileName, const Real unitFactor)
{
    structure->readPoscar(fileName);
    structure->convertUnits(unitFactor);
    structure->cartesianToDirect();
    structureTypes = structure->get_typeNames();

    typeMap = std::make_shared<TypeMap>(*forceFieldTypes, structureTypes);
    // update neighbor lists
    for (auto& list : neighborLists)
    {
        const String descKey = "SHS2-" + list.first;
        if (descriptors[descKey]->get_weight() <= 0) continue;
        list.second->computeNearestNeighborsDirectCoordinates(*structure);
        descriptors[descKey]->set_neighborList(list.second);
    }

    return;
}

void Frame::updateNeighborLists(Structure& structure, const Vec1Int& distributed)
{
    structure.cartesianToDirect();
    structureTypes = structure.get_typeNames();
    typeMap = std::make_shared<TypeMap>(*forceFieldTypes, structureTypes);
    // update neighbor lists
    for (auto& list : neighborLists)
    {
        const String descKey = "SHS2-" + list.first;
        if (descriptors[descKey]->get_weight() <= 0) continue;
        if (distributed.empty()) list.second->computeNearestNeighborsDirectCoordinates(structure);
        else list.second->computeNearestNeighborsDirectCoordinates(structure, distributed);
        descriptors[descKey]->set_neighborList(list.second);
    }

    return;
}

void Frame::update(Structure& structure, const Vec1Int& distributed)
{
    updateNeighborLists(structure, distributed);
    update();

    return;
}

void Frame::update(const String& fileName, const Real unitFactor)
{
    updateNeighborLists(fileName, unitFactor);
    update();

    return;
}

void Frame::update(const std::shared_ptr<TypeMap>& typeMap)
{
    this->typeMap = typeMap;
    update();

    return;
}

const DescriptorMap& Frame::get_descriptors() const
{
    return descriptors;
}

std::shared_ptr<const TypeMap> Frame::get_typeMap() const
{
    return typeMap;
}

const NeighborListMap& Frame::get_neighborLists() const
{
    return neighborLists;
}

std::shared_ptr<Structure> Frame::get_structure() const
{
    return structure;
}
