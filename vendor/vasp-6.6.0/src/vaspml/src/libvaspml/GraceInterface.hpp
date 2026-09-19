#ifndef GRACEINTERFACE_HPP
#define GRACEINTERFACE_HPP

#include "NeighborListGrace.hpp"
#include "Record.hpp"
#include "Structure.hpp"
#include "nearest_neighbor.hpp"
#include "rec.hpp"
#include "types.hpp"

#ifdef VASPML_ENABLE_GRACE
#include <cppflow/cppflow.h>
#else
#include "cppflow_stub.hpp"
#endif

#include <algorithm>  // for max
#include <cstdlib>    // for getenv
#include <filesystem> // for path, operator/, exists, is_directory
#include <fstream>    // for basic_ifstream, ifstream
#include <iterator>   // for istreambuf_iterator, operator!=
#include <stdexcept>  // for runtime_error
#include <tuple>      // for tuple

namespace vaspml
{

class GraceInterface
{
    struct GraceMetaData
    {
        std::filesystem::path model_path;
        CutoffMap             type_cutoff;
        IndexMap              types2indices;
        GraceMetaData(std::string path) : model_path(path), type_cutoff(), types2indices()
        {
            // normalize model path
            if (!std::filesystem::exists(model_path))
            {
#ifndef _WIN32
                auto home = std::getenv("HOME");
#else
                auto home = std::getenv("USERPROFILE");
#endif
                model_path = std::filesystem::path(home) / ".cache" / "grace" / model_path;
                if (!std::filesystem::exists(model_path))
                {
                    throw std::runtime_error("ERROR: Path to GRACE model does not exist! "
                                             "Current path is \""
                                             + model_path.native() + "\".");
                }
            }
            if (!std::filesystem::is_directory(model_path))
            {
                throw std::runtime_error(
                    "ERROR: Path to GRACE model must point to a directory! "
                    "Check that ML_GRACE_MODEL points to the correct location! "
                    "Current path is \""
                    + model_path.native() + "\".");
            }

            // read meta json
            std::ifstream jsonFile((model_path / "metadata.json").native());
            if (!jsonFile.is_open())
            {
                throw std::runtime_error("ERROR: metadata.json file does not exist for this model! "
                                         "Convert from yaml file and try again.");
            }
            std::string jsonText((std::istreambuf_iterator<char>(jsonFile)),
                                 std::istreambuf_iterator<char>());

            Record meta;
            rec::fromJson(jsonText, meta);

            auto    elements = meta.get<Vec1String>("chemical_symbols");
            int32_t i = 0;
            for (auto e : elements)
            {
                types2indices[e] = i;
                i += 1;
            }
            if (meta.typeOf("cutoff") == "Int")
            {
                double cutoff = static_cast<double>(meta.get<int32_t>("cutoff"));
                for (auto e : types2indices) { type_cutoff[e.first] = cutoff / 2; }
            }
            else if (meta.typeOf("cutoff") == "Real")
            {
                double cutoff = meta.get<double>("cutoff");
                for (auto e : types2indices) { type_cutoff[e.first] = cutoff / 2; }
            }
            else
            {
                throw std::runtime_error("NOTIMPLEMENTED: Per element cutoffs not parsed yet!");
            }
        };
    };

    struct GraceResult
    {
        double              energy = 0.0;
        std::vector<double> energy_per_atom;
        std::vector<double> forces;
        std::vector<double> stress;
        GraceResult() {};
        GraceResult([[maybe_unused]] std::vector<cppflow::tensor> output)
        {
#ifdef VASPML_ENABLE_GRACE
            energy = *static_cast<double*>(TF_TensorData(output[1].get_tensor().get()));
            energy_per_atom = output[0].get_data<double>();
            // TODO: drop padded atom properly
            energy_per_atom.resize(energy_per_atom.size() - 1);
            forces = output[2].get_data<double>();
            forces.resize(forces.size() - 3);
            stress = output[3].get_data<double>();
#endif
        };
    } result;

    cppflow::model model;

  public:
    NeighborListGrace neighbors;
    CutoffMap         type_cutoff;
    IndexMap          types2indices;

    GraceInterface(GraceMetaData meta, double padding) :
        model(meta.model_path.native()),
        neighbors(padding, meta.types2indices, meta.type_cutoff),
        type_cutoff(meta.type_cutoff),
        types2indices(meta.types2indices) {};

    GraceInterface(std::string path, double padding) :
        GraceInterface(GraceMetaData(path), padding) {};

    void updateNeighborList(Vec1String typeNames, NearestNeighborNSquare vaspNeighbors)
    {
        neighbors.copyFromVaspml(typeNames, vaspNeighbors);
    }

    void updateNeighborList(Structure structure, NearestNeighborNSquare vaspNeighbors)
    {
        neighbors.copyFromVaspml(structure, vaspNeighbors);
    }

    void update()
    {
        auto                     inputs = neighbors.exportToList("serving_default_");
        std::vector<std::string> output_names({
            "StatefulPartitionedCall:0", // atomic_energy [nat,1]
            "StatefulPartitionedCall:1", // total_energy [-1, 1]
            "StatefulPartitionedCall:2", // total_f [n_at, 3]
            "StatefulPartitionedCall:3", // virial [6]
        });
        result = GraceResult(model(inputs, output_names));
    }

    double getEnergy() const { return result.energy; }

    std::vector<double> getEnergyPerAtom() const { return result.energy_per_atom; }

    std::vector<double> getForces() const { return result.forces; }

    std::vector<double> getStressTensor() const { return result.stress; }

    double getMaximumCutoff() const
    {
        double mrc = 0.0;
        for (auto p1 : type_cutoff)
        {
            for (auto p2 : type_cutoff) { mrc = std::max(mrc, p1.second + p2.second); }
        }
        return mrc;
    }
};

} //namespace vaspml

#endif
