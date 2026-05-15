#include "gtest/gtest.h"
#include <map>
#include <random>

#include "Constants.h"
#include "FieldSolver.h"
#include "Grid.h"
#include "PmlFdtd.h"
#include "FieldBoundaryConditionFdtd.h"
#include "FieldGeneratorFdtd.h"
#include "Vectors.h"

#include "HiChi_MPI.h"

#define INVALID mpi::MPI_INVALID_RANK

using GridType = pfc::YeeGrid;
using Int3 = pfc::Int3;
using FP3 = pfc::FP3;

TEST(MPI_FieldExchanger, CanCreateFieldExchanger)
{
    Int3 num_cells = Int3(3, 3, 3);
    int num_external_cells = 1;

    ASSERT_NO_THROW( mpi::FieldExchanger exchanger(num_cells, num_external_cells) );
}

void initializeGrid(std::unique_ptr<GridType>& grid, int rank) 
{
    srand(rank);

    Int3 begin = Int3(0, 0, 0);
    Int3 end = grid->numCells;

    for (int i = begin.x; i < end.x; i++)
        for (int j = begin.y; j < end.y; j++)
            for (int k = begin.z; k < end.z; k++)
                grid->Ex(i, j, k) = static_cast<pfc::FP>(rand());
}

void zeroGrid(std::unique_ptr<GridType>& grid) 
{
    Int3 begin = Int3(0, 0, 0);
    Int3 end = grid->numCells;

    for (int i = begin.x; i < end.x; i++)
        for (int j = begin.y; j < end.y; j++)
            for (int k = begin.z; k < end.z; k++)
                grid->Ex(i, j, k) = 0.0;
}

bool isExternalCellAlongDirection(std::unique_ptr<GridType>& grid, Int3 cell, mpi::Direction dir, int num_external)
{
    const Int3& size = grid->numCells;

    if (dir == mpi::Direction::pX && cell.x >= size.x - num_external)
        return true;

    if (dir == mpi::Direction::nX && cell.x < num_external)
        return true;

    if (dir == mpi::Direction::pY && cell.y >= size.y - num_external)
        return true;

    if (dir == mpi::Direction::nY && cell.y < num_external)
        return true;

    if (dir == mpi::Direction::pZ && cell.z >= size.z - num_external)
        return true;

    if (dir == mpi::Direction::nZ && cell.z < num_external)
        return true;

    return false;
}

// Int3 mapRecvExternalCellToSendInternal(std::unique_ptr<GridType>& grid, Int3 cell)
// {
//     const Int3& size = grid->numCells;
//     const Int3& num_internal = grid->numInternalCells;
//     int num_external = grid->numExternalCells;

//     Int3 result;

//     result.x = (cell.x < num_external) ? (size.x - 2 * num_external + cell.x) : ((cell.x >= size.x - num_external) ? (2 * num_external - (size.x - cell.x)) : (cell.x) );
//     result.y = (cell.y < num_external) ? (size.y - 2 * num_external + cell.y) : ((cell.y >= size.y - num_external) ? (2 * num_external - (size.y - cell.y)) : (cell.y) );
//     result.z = (cell.z < num_external) ? (size.z - 2 * num_external + cell.z) : ((cell.z >= size.z - num_external) ? (2 * num_external - (size.z - cell.z)) : (cell.z) );

//     return result;
// }

Int3 mapRecvExternalCellToSendInternal(std::unique_ptr<GridType>& grid, Int3 cell, mpi::Direction recv_direction, int num_external)
{
    const Int3& size = grid->numCells;
    Int3 result = cell;

    switch(recv_direction) 
    {
        case mpi::Direction::nX:
            if (cell.x < num_external)
                result.x = size.x - 2 * num_external + cell.x;
            break;

        case mpi::Direction::pX:
            if (cell.x >= size.x - num_external)
                result.x = 2 * num_external - (size.x - cell.x);
            break;

        case mpi::Direction::nY:
            if (cell.y < num_external)
                result.y = size.y - 2 * num_external + cell.y;
            break;

        case mpi::Direction::pY:
            if (cell.y >= size.y - num_external)
                result.y = 2 * num_external - (size.y - cell.y);
            break;

        case mpi::Direction::nZ:
            if (cell.z < num_external)
                result.z = size.z - 2 * num_external + cell.z;
            break;

        case mpi::Direction::pZ:
            if (cell.z >= size.z - num_external)
                result.z = 2 * num_external - (size.z - cell.z);
            break;
    }
    
    return result;
}

void debugPrintGrid_xyplane(std::unique_ptr<GridType>& grid, int z)
{
    Int3 size = grid->numCells;

    for (int j = 0; j < size.y; j++)
    {
        std::string line = "";
        for (int i = 0; i < size.x; i++) 
            line += std::to_string(grid->Ex(i, j, z)) + ", ";
            //line += (isExternalCellAlongDirection(grid, Int3(i, j, z), mpi::Direction::nZ)) ? std::to_string(mapRecvExternalCellToSendInternal(grid, Int3(i, j, z)).z) : " 0 ";
        std::cout << line << std::endl;
    }

    std::cout << std::endl << std::endl;
}

void testExchange(
    Int3 topologySize, 
    mpi::Topology::LoopType loopType,
    Int3 gridSize,
    FP3 minCoords,
    FP3 maxCoords,
    FP3 gridStep,
    int numExternalCells)
{
    MPI_Barrier(MPI_COMM_WORLD);

    // Checking mpi size
    int real_size;
    MPI_Comm_size(MPI_COMM_WORLD, &real_size);

    int size = topologySize.x * topologySize.y * topologySize.z;
    if (real_size < size)
        GTEST_SKIP() << "Not enough processes " << real_size << "/" << size;

    // Running
    mpi::Topology topology(topologySize, loopType, size);

    if (!topology.isValidOnThisRank())
    {
        EXPECT_TRUE(true);
        return;
    }

    int rank;
    MPI_Comm_rank(topology.getTopologyCommunicator(), &rank);

    std::unique_ptr<GridType> grid;
    grid.reset(new GridType(gridSize, minCoords, gridStep, gridSize));
    initializeGrid(grid, rank);

    mpi::FieldExchanger exchanger(grid->numCells, numExternalCells);

    MPI_Barrier(topology.getTopologyCommunicator());
    exchanger.performExchangeSequence(grid->Ex.getData(), topology, rank, topology.getTopologyCommunicator());
    MPI_Barrier(topology.getTopologyCommunicator());

    // Check
    for (int r = 0; r < real_size; r++)
    {
        for (int dir = 0; dir < 6; dir++)
        {
            mpi::Direction direction = static_cast<mpi::Direction>(dir);
            int neighbor = topology.getNeighbor(rank, direction);
            int inverse_neighbor = topology.getNeighbor(rank, mpi::invertDirection(direction));

            // Exchange full grid data

            std::unique_ptr<GridType> neighbor_grid;
            neighbor_grid.reset(new GridType(gridSize, minCoords, gridStep, gridSize));
            zeroGrid(neighbor_grid);

            MPI_Request request;
            MPI_Status status;

            MPI_Barrier(topology.getTopologyCommunicator());
            if (neighbor != mpi::MPI_INVALID_RANK)
                MPI_Isend(grid->Ex.getData(), grid->Ex.toVector().size(), MPI_DOUBLE, neighbor, 0, topology.getTopologyCommunicator(), &request);

            if (inverse_neighbor != mpi::MPI_INVALID_RANK)
                MPI_Recv(neighbor_grid->Ex.getData(), neighbor_grid->Ex.toVector().size(), MPI_DOUBLE, inverse_neighbor, 0, topology.getTopologyCommunicator(), &status);

            if (neighbor != mpi::MPI_INVALID_RANK)
                MPI_Wait(&request, &status);

            MPI_Barrier(topology.getTopologyCommunicator());

            if (inverse_neighbor != mpi::MPI_INVALID_RANK)
            {
                // Assertion
                for (int x = 0; x < grid->numCells.x; x++)
                    for (int y = 0; y < grid->numCells.y; y++)
                        for (int z = 0; z < grid->numCells.z; z++)
                        {
                            Int3 cell(x, y, z);
                            if(isExternalCellAlongDirection(grid, cell, mpi::invertDirection(direction), numExternalCells))
                            {
                                Int3 neighbor_cell = mapRecvExternalCellToSendInternal(grid, cell, mpi::invertDirection(direction), numExternalCells);
                                bool check = abs(grid->Ex(cell) - neighbor_grid->Ex(neighbor_cell)) < 1e-5;
                                EXPECT_TRUE(check);
                            }
                        }
            }
        }
    }

    MPI_Barrier(topology.getTopologyCommunicator());
}





TEST(MPI_FieldExchanger, Exhange_Uniform_OneExternalCell_2_2_1)
{
    Int3 topologySize = pfc::Int3(2, 2, 1);
    mpi::Topology::LoopType loopType = mpi::Topology::LoopType::None;

    Int3 gridSize = Int3(6, 6, 6);
    FP3 minCoords = FP3(0, 0, 0);
    FP3 maxCoords = FP3(1, 1, 1);
    FP3 gridStep = (maxCoords - minCoords) / (FP3)gridSize;
    int numExternalCells = 1;

    testExchange(topologySize, loopType, gridSize, minCoords, maxCoords, gridStep, numExternalCells);
}

TEST(MPI_FieldExchanger, Exhange_Uniform_TwoExternalCells_2_2_1)
{
    MPI_Barrier(MPI_COMM_WORLD);

    Int3 topologySize = pfc::Int3(2, 2, 1);
    mpi::Topology::LoopType loopType = mpi::Topology::LoopType::None;

    Int3 gridSize = Int3(6, 6, 6);
    FP3 minCoords = FP3(0, 0, 0);
    FP3 maxCoords = FP3(1, 1, 1);
    FP3 gridStep = (maxCoords - minCoords) / (FP3)gridSize;
    int numExternalCells = 2;

    testExchange(topologySize, loopType, gridSize, minCoords, maxCoords, gridStep, numExternalCells);
}


TEST(MPI_FieldExchanger, Exhange_Uniform_OneExternalCell_3_3_3)
{
    Int3 topologySize = pfc::Int3(3, 3, 3);
    mpi::Topology::LoopType loopType = mpi::Topology::LoopType::None;

    Int3 gridSize = Int3(6, 6, 6);
    FP3 minCoords = FP3(0, 0, 0);
    FP3 maxCoords = FP3(1, 1, 1);
    FP3 gridStep = (maxCoords - minCoords) / (FP3)gridSize;
    int numExternalCells = 1;

    testExchange(topologySize, loopType, gridSize, minCoords, maxCoords, gridStep, numExternalCells);
}

TEST(MPI_FieldExchanger, Exhange_Uniform_TwoExternalCells_3_3_3)
{
    MPI_Barrier(MPI_COMM_WORLD);

    Int3 topologySize = pfc::Int3(3, 3, 3);
    mpi::Topology::LoopType loopType = mpi::Topology::LoopType::None;

    Int3 gridSize = Int3(6, 6, 6);
    FP3 minCoords = FP3(0, 0, 0);
    FP3 maxCoords = FP3(1, 1, 1);
    FP3 gridStep = (maxCoords - minCoords) / (FP3)gridSize;
    int numExternalCells = 2;

    testExchange(topologySize, loopType, gridSize, minCoords, maxCoords, gridStep, numExternalCells);
}


TEST(MPI_FieldExchanger, Exhange_Uniform_OneExternalCell_3_4_5)
{
    MPI_Barrier(MPI_COMM_WORLD);

    Int3 topologySize = pfc::Int3(3, 4, 5);
    mpi::Topology::LoopType loopType = mpi::Topology::LoopType::None;

    Int3 gridSize = Int3(6, 6, 6);
    FP3 minCoords = FP3(0, 0, 0);
    FP3 maxCoords = FP3(1, 1, 1);
    FP3 gridStep = (maxCoords - minCoords) / (FP3)gridSize;
    int numExternalCells = 1;

    testExchange(topologySize, loopType, gridSize, minCoords, maxCoords, gridStep, numExternalCells);
}

TEST(MPI_FieldExchanger, Exhange_Uniform_TwoExternalCells_3_4_5)
{
    MPI_Barrier(MPI_COMM_WORLD);

    Int3 topologySize = pfc::Int3(3, 4, 5);
    mpi::Topology::LoopType loopType = mpi::Topology::LoopType::None;

    Int3 gridSize = Int3(6, 6, 6);
    FP3 minCoords = FP3(0, 0, 0);
    FP3 maxCoords = FP3(1, 1, 1);
    FP3 gridStep = (maxCoords - minCoords) / (FP3)gridSize;
    int numExternalCells = 2;

    testExchange(topologySize, loopType, gridSize, minCoords, maxCoords, gridStep, numExternalCells);
}