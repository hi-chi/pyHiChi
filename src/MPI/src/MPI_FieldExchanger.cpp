#include "MPI_FieldExchanger.h"
#include "MPI_Topology.h"
#include "MPI_FieldUtils.h"
#include "MPI_Utilities.h"

namespace mpi
{

// Constructor that will create types for sending and recieving 
FieldExchanger::FieldExchanger(const pfc::Int3& grid_num_cells, int num_external_cells)
{
    // Lower limit for num_external_cells (for cases when num_external_cells is 0, for example in spectral solvers)
    if (num_external_cells < 1)
        num_external_cells = 1;

    // Creating types
    for (int d = 0; d < 6; d++)
    {
        Direction direction = static_cast<Direction>(d);

        // Send
        MPI_Datatype send_type;
        FieldUtils::defineTransmission_Send(send_type, grid_num_cells, num_external_cells, direction);
        sendTypes.push_back(send_type);

        // Recv
        MPI_Datatype recv_type;
        FieldUtils::defineTransmission_Recv(recv_type, grid_num_cells, num_external_cells, direction);
        recvTypes.push_back(recv_type);
    }
}

// Frees created mpi types
FieldExchanger::~FieldExchanger()
{
    int mpiInit;
    MPI_Initialized(&mpiInit);

    int mpiFinalized;
    MPI_Finalized(&mpiFinalized);

    if (mpiInit && !mpiFinalized)
        for (int i = 0; i < 6; i++)
        {
            MPI_Type_free(&sendTypes[i]);
            MPI_Type_free(&recvTypes[i]);
        }
}


// Performs echange sequence iteration for node "rank" in the specified topology
// BLOCKING
// MUST BE CALLED BY EVERY PROCESS IN THE TOPOLOGY!
void FieldExchanger::performExchangeOverAxis(pfc::FP* data, pfc::CoordinateEnum axis, const class Topology& topology, int rank, MPI_Comm communicator)
{
    int posDirection = static_cast<int>(axis) * 2;
    int negDirection = posDirection + 1;
    for (int dir = posDirection; dir <= negDirection; dir++)
    {
        Direction direction = static_cast<Direction>(dir);

        // Neighbors
        int send_neighbor = topology.getNeighbor(rank, direction);
        int recv_neighbor = topology.getNeighbor(rank, invertDirection(direction));

        // Sync nodes
        MPI_Barrier(communicator);

        MPI_Request request;
        MPI_Status status;

        // Sending
        if (send_neighbor != MPI_INVALID_RANK)
            MPI_Isend(data, 1, sendTypes[dir], send_neighbor, 0, communicator, &request);
        
        // Receiving
        if (recv_neighbor != MPI_INVALID_RANK)
            MPI_Recv(data, 1, recvTypes[dir], recv_neighbor, 0, communicator, &status);

        // Waiting for send operation to complete
        if (send_neighbor != MPI_INVALID_RANK)
            MPI_Wait(&request, &status);
    }

    // Sync nodes
    MPI_Barrier(communicator);
}

// Performs echange sequence iteration for node "rank" in the specified topology
// MUST BE CALLED BY EVERY PROCESS IN THE TOPOLOGY AT THE SAME TIME!
void FieldExchanger::performExchangeSequence(pfc::FP* data, const Topology& topology, int rank, MPI_Comm communicator)
{
    performExchangeOverAxis(data, pfc::CoordinateEnum::x, topology, rank, communicator);
    performExchangeOverAxis(data, pfc::CoordinateEnum::y, topology, rank, communicator);
    performExchangeOverAxis(data, pfc::CoordinateEnum::z, topology, rank, communicator);
}

}