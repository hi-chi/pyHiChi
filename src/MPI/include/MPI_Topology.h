#pragma once

#include <vector>

#include "Vectors.h"
#include "MPI_Utilities.h"

namespace mpi
{

// A class that describes abstracted mpi topology
class Topology final
{
public:

    // MPI topology looping type
    enum class LoopType
    {
        None,
        LoopX, LoopY, LoopZ,
        LoopXY, LoopYZ, LoopXZ,
        LoopXYZ
    };

public:
    // Constructing an abstracted MPI topology based on the number of nodes (ranks) and the desired topology type
    // `sections` - Number of sections for each axis
    Topology(const pfc::Int3& sections, LoopType loop_type, int node_count);
    Topology(const pfc::Int3& sections, bool loop_mask[3], int node_count);

    // Frees mpi communicator
    ~Topology();

private:

    // Number of sections for each axis
    pfc::Int3 sections_;
    // Topology looping type
    bool loop[3] = {false, false, false};
    // Number of nodes in this topology
    int size;

    // MPI Communicator for processes that participate in this topology
    MPI_Comm communicator = MPI_COMM_NULL;

    // Separates a topology communicator from MPI_COMM_WORLD
    void initCommunicator();

    // Converts a 3D index into a 1D index
    // Returns MPI_INVALID_RANK if at least one of coordinates is MPI_INVALID_RANK
    int convertIndexToRank(int x, int y, int z) const;
    // Converts a 3D index into a 1D index
    // Returns MPI_INVALID_RANK if at least one of coordinates is MPI_INVALID_RANK
    int convertIndexToRank(pfc::Int3 index) const { return convertIndexToRank(index.x, index.y, index.z); }

    // Converts a rank index into a 3D index
    // Return MPI_INVALID_SECTION if rank is invalid
    pfc::Int3 convertRankToIndex(int rank) const;

public:

    static bool doesLoopOverX(const LoopType& loop_type);
    static bool doesLoopOverY(const LoopType& loop_type);
    static bool doesLoopOverZ(const LoopType& loop_type);

    // Returns sections
    const pfc::Int3& getSections() const { return sections_; }

    // Returns the type of looping used for this topology
    void getLoopMask(bool out_loop_mask[3]) const;

    // Returns the number of nodes in this topology
    int getSize() const { return size; }

    // Returns rank of the node, that is responsible for processing the specified section
    // Returns MPI_INVALID_RANK if at least one of coordinates is MPI_INVALID_RANK
    int getResponsibleNode(const pfc::Int3& section) const { return convertIndexToRank(section); }

    // Returns a 3D index of the section, assigned to a node with the specified rank
    // Return MPI_INVALID_SECTION if rank is invalid
    pfc::Int3 getNodeSection(int rank) const { return convertRankToIndex(rank); }

    // Returns rank of the node's neighbor in the specified direction
    // MPI_INVALID_RANK if no neighbor exists
    int getNeighbor(int node_rank, Direction direction) const;

    // Returns rank of the node's neighbor with the specified offset
    // MPI_INVALID_RANK if no neighbor exists
    int getNeighbor(int node_rank, const pfc::Int3& offset) const;


    // MPI Communicator for processes that participate in this topology
    // NOTE: Depends on the rank this method is executed at
    MPI_Comm getTopologyCommunicator() { return communicator; }

    // Whether this rank participates in topology or not
    // NOTE: Depends on the rank this method is executed at
    bool isValidOnThisRank() { return communicator != MPI_COMM_NULL; }
};

}