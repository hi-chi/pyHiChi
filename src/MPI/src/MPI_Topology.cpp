#include "MPI_Topology.h"

namespace mpi
{
    
// Constructing an abstracted MPI topology based on the number of nodes (ranks) and the desired topology type
Topology::Topology(const pfc::Int3& sections, LoopType loop_type, int node_count) : sections_(sections)
{
    // Checking if topology is possible
    if (node_count < sections.x * sections.y * sections.z)
        throw std::runtime_error("MPI TOPOLOGY : Impossible topology for the specified number of nodes");

    size = sections.x * sections.y * sections.z;

    loop[0] = doesLoopOverX(loop_type);
    loop[1] = doesLoopOverY(loop_type);
    loop[2] = doesLoopOverZ(loop_type);


    initCommunicator();
}
Topology::Topology(const pfc::Int3& sections, bool loop_mask[3], int node_count) : sections_(sections)
{
    // Checking if topology is possible
    if (node_count < sections.x * sections.y * sections.z)
        throw std::runtime_error("MPI TOPOLOGY : Impossible topology for the specified number of nodes");

    size = sections.x * sections.y * sections.z;

    for (int i = 0; i < 3; i++)
        loop[i] = loop_mask[i];

    initCommunicator();
}

// Frees mpi communicator
Topology::~Topology()
{
    if (communicator != MPI_COMM_NULL)
        MPI_Comm_free(&communicator);
}


// Separates a topology communicator from MPI_COMM_WORLD
void Topology::initCommunicator()
{
    int mpi;
    MPI_Initialized(&mpi);
    if (!mpi) return;

    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    // Whether this rank belongs to the topology
    int color = (mpi_rank < sections_.x * sections_.y * sections_.z) ? 1 : MPI_UNDEFINED;
    int key = mpi_rank;

    MPI_Comm_split(MPI_COMM_WORLD, color, key, &communicator);
}

// Converts a 3D index into a 1D index
// Returns MPI_INVALID_RANK if at least one of coordinates is MPI_INVALID_RANK
int Topology::convertIndexToRank(int x, int y, int z) const
{
    if (x == MPI_INVALID_RANK || y == MPI_INVALID_RANK || z == MPI_INVALID_RANK)
        return MPI_INVALID_RANK;

    return x * sections_.y * sections_.z + y * sections_.z + z;
}

 // Converts a rank index into a 3D index
pfc::Int3 Topology::convertRankToIndex(int rank) const
{
    if (rank < 0 || rank >= size)
        return MPI_INVALID_SECTION;
    
    pfc::Int3 index;

    index.x = rank / (sections_.y * sections_.z);
    index.y = (rank / sections_.z) % sections_.y;
    index.z = rank % sections_.z;
    
    return index;
}

// Returns the type of looping used for this topology
void Topology::getLoopMask(bool out_loop_mask[3]) const
{
    for (int i = 0; i < 3; i++)
        out_loop_mask[i] = loop[i];
}

bool Topology::doesLoopOverX(const LoopType& loop_type)
{
    return loop_type == LoopType::LoopX 
        || loop_type == LoopType::LoopXY 
        || loop_type == LoopType::LoopXZ 
        || loop_type == LoopType::LoopXYZ;
}
bool Topology::doesLoopOverY(const LoopType& loop_type)
{
    return loop_type == LoopType::LoopY
        || loop_type == LoopType::LoopXY 
        || loop_type == LoopType::LoopYZ 
        || loop_type == LoopType::LoopXYZ;
}
bool Topology::doesLoopOverZ(const LoopType& loop_type)
{
    return loop_type == LoopType::LoopZ
        || loop_type == LoopType::LoopXZ 
        || loop_type == LoopType::LoopYZ 
        || loop_type == LoopType::LoopXYZ;
}


// Returns rank of the node's neighbor in the specified direction
// MPI_INVALID_RANK if no neighbor exists
int Topology::getNeighbor(int node_rank, Direction direction) const
{
    return getNeighbor(node_rank, mpi::DirectionToOffset[static_cast<int>(direction)]);
}

// Returns rank of the node's neighbor with the specified offset
// MPI_INVALID_RANK if no neighbor exists
int Topology::getNeighbor(int node_rank, const pfc::Int3& offset) const
{
    if (node_rank == MPI_INVALID_RANK)
        return MPI_INVALID_RANK;

    pfc::Int3 index = convertRankToIndex(node_rank);
    if (index == MPI_INVALID_SECTION)
        return MPI_INVALID_RANK;

    pfc::Int3 neighbor_index = index + offset;
    
    // Validating and looping neighbor_index
    neighbor_index.x = (neighbor_index.x >= 0) ? (neighbor_index.x) : ( loop[0] ? (sections_.x + neighbor_index.x) : MPI_INVALID_RANK );
    neighbor_index.x = (neighbor_index.x < sections_.x) ? (neighbor_index.x) : ( loop[0] ? (neighbor_index.x - sections_.x) : MPI_INVALID_RANK );

    neighbor_index.y = (neighbor_index.y >= 0) ? (neighbor_index.y) : ( loop[1] ? (sections_.y + neighbor_index.y) : MPI_INVALID_RANK );
    neighbor_index.y = (neighbor_index.y < sections_.y) ? (neighbor_index.y) : ( loop[1] ? (neighbor_index.y - sections_.y) : MPI_INVALID_RANK );

    neighbor_index.z = (neighbor_index.z >= 0) ? (neighbor_index.z) : ( loop[2] ? (sections_.z + neighbor_index.z) : MPI_INVALID_RANK );
    neighbor_index.z = (neighbor_index.z < sections_.z) ? (neighbor_index.z) : ( loop[2] ? (neighbor_index.z - sections_.z) : MPI_INVALID_RANK );

    // Converting neighbor_index to rank
    return convertIndexToRank(neighbor_index);
}

}