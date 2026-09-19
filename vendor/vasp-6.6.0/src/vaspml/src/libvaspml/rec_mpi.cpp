#include "MlMPI.hpp"

#include "rec_mpi.hpp"

#include "Record.hpp"
#include "Serializer.hpp"
#include "Tutor.hpp"
#include "debug.hpp"
#include "math.hpp"
#include "rec.hpp"

void vaspml::rec::mpi::bcast(Record& data, MlMPI& mlmpi, const Int root)
{
    Buffer buffer;
    Int    bufferSize = 0;
    if (mlmpi.get_rank() == root)
    {
        rec::toBuffer(data, buffer);
        bufferSize = buffer.size() * sizeof(Buffer::value_type);
    }
    mlmpi.bcast(bufferSize, 0);
    if (mlmpi.get_rank() != root) buffer.resize(bufferSize);
    mlmpi.bcast(buffer.data(), bufferSize, root);
    if (mlmpi.get_rank() != root) rec::fromBuffer(buffer, data);

    return;
}

vaspml::Vec1ShRec vaspml::rec::mpi::allGatherVVector(Record& data, MlMPI& mlmpi)
{
    Buffer buffer;
    rec::toBuffer(data, buffer);

    Vec1Int recCounts(mlmpi.get_numberRanks());
    recCounts[mlmpi.get_rank()] = buffer.size();

    mlmpi.allReduceSum(recCounts.data(), recCounts.data(), mlmpi.get_numberRanks());
    const Vec1Int offsets = math::partialSum0(recCounts);
    const Int     totalElements = math::sumVector(recCounts);

    Buffer totalBuffer(totalElements);

    [[maybe_unused]] Int info = MPI_Allgatherv(buffer.data(),
                                               buffer.size(),
                                               MPI_BYTE,
                                               totalBuffer.data(),
                                               recCounts.data(),
                                               offsets.data(),
                                               MPI_BYTE,
                                               mlmpi.get_communicator());
    VASPML_DEBUG_L1(
        if (info != 0)
        {
            global::tutor.error("ERROR: " + flf(VASPML_FLF)
                                + " MPI_Allgatherv exited with error code " + std::to_string(info)
                                + ".\n");
        }
    );
    Vec1ShRec totalRecord(mlmpi.get_numberRanks());
    for (Int i = 0; i < mlmpi.get_numberRanks(); i++)
    {
        totalRecord[i] = std::make_shared<Record>();
        vaspml::rec::detail::deserialize(totalBuffer, *(totalRecord[i]), offsets[i]);
    }

    return totalRecord;
}

vaspml::ShRec vaspml::rec::mpi::allGatherV(Record& data, MlMPI& mlmpi)
{
    Buffer buffer;
    rec::toBuffer(data, buffer);
    Vec1Int recCounts(mlmpi.get_numberRanks());
    recCounts[mlmpi.get_rank()] = buffer.size();
    mlmpi.allReduceSum(recCounts.data(), recCounts.data(), mlmpi.get_numberRanks());
    const Vec1Int offsets = math::partialSum0(recCounts);
    const Int     totalElements = math::sumVector(recCounts);
    Buffer        totalBuffer(totalElements);

    [[maybe_unused]] Int info = MPI_Allgatherv(buffer.data(),
                                               buffer.size(),
                                               MPI_BYTE,
                                               totalBuffer.data(),
                                               recCounts.data(),
                                               offsets.data(),
                                               MPI_BYTE,
                                               mlmpi.get_communicator());
    VASPML_DEBUG_L1(
        if (info != 0)
        {
            global::tutor.error("ERROR: " + flf(VASPML_FLF)
                                + " MPI_Allgatherv exited with error code " + std::to_string(info)
                                + ".\n");
        }
    );
    ShRec totalRecord;
    totalRecord = std::make_shared<Record>();
    vaspml::rec::detail::deserialize(totalBuffer, *totalRecord, offsets[0]);
    for (Int i = 1; i < mlmpi.get_numberRanks(); i++)
    {
        Record temp;
        vaspml::rec::detail::deserialize(totalBuffer, temp, offsets[i]);
        rec::merge(*totalRecord, temp);
    }

    return totalRecord;
}

vaspml::ShRec vaspml::rec::mpi::gatherV(Record& data, MlMPI& mlmpi, const Int root)
{
    Buffer buffer;
    rec::toBuffer(data, buffer);

    // Prepare to gather buffer sizes at root
    Vec1Int recCounts;
    if (mlmpi.get_rank() == root) recCounts.resize(mlmpi.get_numberRanks(), 0);

    // Send local buffer size to root
    Int sendCount = buffer.size();
    mlmpi.gather(&sendCount, recCounts.data(), 1, 1, root

    );
    Vec1Int offsets;
    Int     totalElements = 0;
    if (mlmpi.get_rank() == root)
    {
        offsets = math::partialSum0(recCounts);
        totalElements = math::sumVector(recCounts);
    }
    Buffer totalBuffer;
    if (mlmpi.get_rank() == root) totalBuffer.resize(totalElements);
    MPI_Gatherv(buffer.data(),
                sendCount,
                MPI_BYTE,
                totalBuffer.data(),
                recCounts.data(),
                offsets.data(),
                MPI_BYTE,
                root,
                mlmpi.get_communicator());
    if (mlmpi.get_rank() == root)
    {
        ShRec totalRecord = std::make_shared<Record>();
        vaspml::rec::detail::deserialize(totalBuffer, *totalRecord, offsets[0]);
        for (Int i = 1; i < mlmpi.get_numberRanks(); i++)
        {
            Record temp;
            vaspml::rec::detail::deserialize(totalBuffer, temp, offsets[i]);
            rec::merge(*totalRecord, temp);
        }
        return totalRecord;
    }
    else { return nullptr; }
}
