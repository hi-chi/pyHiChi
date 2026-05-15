#include "gtest/gtest.h"
#include <map>

#include "Grid.h"
#include "Constants.h"

#include "HiChi_MPI.h"

using FP = pfc::FP;
using FP3 = pfc::FP3;
using Int3 = pfc::Int3;
using GridType = pfc::YeeGrid;

#define NUM_EXTERNAL_CELLS 1

FP testFunc(const FP3& coords)
{
    return sin(coords.x + 10) * 3 + cos(coords.y * 3 - 3.14) / 5 - pow(coords.z, 2);
}

int getResponsibleSubGrid(  Int3& outIndexOrigin,
                            const Int3& index, const FP3& minCoords, const FP3& gridStep, 
                            std::vector<int> divisions[3], const mpi::Topology& topology)
{
    Int3 section = Int3(0, 0, 0);
    for (int div_x = 0; div_x < divisions[0].size(); div_x++)
        if (index.x >= divisions[0][div_x])
        {
            section.x++;
        }

    for (int div_y = 0; div_y < divisions[1].size(); div_y++)
        if (index.y >= divisions[1][div_y])
        {
            section.y++;
        }

    for (int div_z = 0; div_z < divisions[2].size(); div_z++)
        if (index.z >= divisions[2][div_z])
        {
            section.z++;
        }

    outIndexOrigin.x = (section.x == 0) ? 0 : divisions[0][section.x - 1];
    outIndexOrigin.y = (section.y == 0) ? 0 : divisions[1][section.y - 1];
    outIndexOrigin.z = (section.z == 0) ? 0 : divisions[2][section.z - 1];

    //outIndexOrigin

    return topology.getResponsibleNode(section);
}


void debugPrintGrid_xzplane(std::unique_ptr<GridType>& grid, int y)
{
    Int3 size = grid->numCells;

    for (int j = 0; j < size.x; j++)
    {
        std::string line = "";
        for (int i = 0; i < size.z; i++)
            line += std::to_string(grid->Ex(j, y, i)) + ", ";
        std::cout << line << std::endl;
    }

    std::cout << std::endl << std::endl;
}


TEST(MPI_GridSlicer, GridSlicerIsValid_1)
{
    // Parameters
    FP3 minCoords = FP3(0, 0, 0);
    Int3 gridSize = Int3(5, 10, 15);
    FP3 gridStep = FP3(1, 1, 1);
    FP3 maxCoords = gridStep * gridSize;

    Int3 sections = Int3(2, 2, 3);
    std::vector<int> divisions[3] = {
        {2},
        {5},
        {5, 10}
    };
    int nodeCount = sections.x * sections.y * sections.z;
    mpi::Topology topology(sections, mpi::Topology::LoopType::None, nodeCount);


    // Initializing main grid
    std::unique_ptr<GridType> mainGrid;
    mainGrid.reset(new GridType(gridSize, minCoords, gridStep, gridSize));

    // Filling with values
    Int3 begin = Int3(0, 0, 0);
    Int3 end = mainGrid->numCells;

    for (int i = begin.x; i < end.x; i++)
        for (int j = begin.y; j < end.y; j++)
            for (int k = begin.z; k < end.z; k++) 
            {
                FP3 coords = mainGrid->ExPosition(i, j, k);
                mainGrid->Ex(i, j, k) = testFunc(coords);
            }

    //debugPrintGrid_xzplane(mainGrid, 0);

    // Initializing sub grids
    std::vector<std::unique_ptr<GridType>> subGrids;
    subGrids.resize(nodeCount);
    for (int r = 0; r < nodeCount; r++)
    {
        FP3 localMinCoords;
        Int3 localGridSize;

        // Slicing
        mpi::GridSlicer::getSubGridParameters(  localMinCoords, localGridSize, 
                                                minCoords, gridSize, gridStep, 
                                                divisions, r, topology);

        subGrids[r].reset(new GridType(localGridSize, localMinCoords, gridStep, localGridSize));

        std::unique_ptr<GridType>& grid = subGrids[r];

        // Filling with values
        Int3 begin = Int3(0, 0, 0);
        Int3 end = grid->numCells;

        for (int i = begin.x; i < end.x; i++)
            for (int j = begin.y; j < end.y; j++)
                for (int k = begin.z; k < end.z; k++) 
                {
                    FP3 coords = grid->ExPosition(i, j, k);
                    grid->Ex(i, j, k) = testFunc(coords);
                }

        //std::cout << r << std::endl;
        //debugPrintGrid_xzplane(grid, 0);
    }

    // Validation
    for (int x = 0; x < gridSize.x; x++)
        for (int y = 0; y < gridSize.y; y++)
            for (int z = 0; z < gridSize.z; z++)
            {
                Int3 index(x, y, z);
                FP3 coords = mainGrid->ExPosition(index.x, index.y, index.z);
                FP expected = testFunc(coords);

                Int3 local_index_origin;
                int r = getResponsibleSubGrid(local_index_origin, index, minCoords, gridStep, divisions, topology);

                //std::cout << index.x << ", " << index.z << ", " << r << std::endl;

                Int3 localIndex = index - local_index_origin;
                FP actual = subGrids[r]->Ex(localIndex);

                EXPECT_NEAR(expected, actual, 0.0001);
            }
}



TEST(MPI_GridSlicer, GridSlicerIsValid_2)
{
    // Parameters
    Int3 gridSize = Int3(32, 8, 8);
    FP3 minCoords = FP3(0, 0, 0);
    FP3 maxCoords = pfc::constants::c * (FP3)gridSize;
    FP3 gridStep = (maxCoords - minCoords) / (FP3)gridSize;

    Int3 sections = Int3(2, 2, 1);
    std::vector<int> divisions[3] = {
        {gridSize.x / 2},
        {gridSize.y / 2},
        {}
    };
    int nodeCount = sections.x * sections.y * sections.z;
    mpi::Topology topology(sections, mpi::Topology::LoopType::None, nodeCount);


    // Initializing main grid
    std::unique_ptr<GridType> mainGrid;
    mainGrid.reset(new GridType(gridSize, minCoords, gridStep, gridSize));

    // Filling with values
    Int3 begin = Int3(0, 0, 0);
    Int3 end = mainGrid->numCells;

    for (int i = begin.x; i < end.x; i++)
        for (int j = begin.y; j < end.y; j++)
            for (int k = begin.z; k < end.z; k++) 
            {
                FP3 coords = mainGrid->ExPosition(i, j, k);
                mainGrid->Ex(i, j, k) = testFunc(coords);
            }


    // Initializing sub grids
    std::vector<std::unique_ptr<GridType>> subGrids;
    subGrids.resize(nodeCount);
    for (int r = 0; r < nodeCount; r++)
    {
        FP3 localMinCoords;
        Int3 localGridSize;

        // Slicing
        mpi::GridSlicer::getSubGridParameters(  localMinCoords, localGridSize, 
                                                minCoords, gridSize, gridStep, 
                                                divisions, r, topology);

        subGrids[r].reset(new GridType(localGridSize, localMinCoords, gridStep, localGridSize));

        std::unique_ptr<GridType>& grid = subGrids[r];

        // Filling with values
        Int3 begin = Int3(0, 0, 0);
        Int3 end = grid->numCells;

        for (int i = begin.x; i < end.x; i++)
            for (int j = begin.y; j < end.y; j++)
                for (int k = begin.z; k < end.z; k++) 
                {
                    FP3 coords = grid->ExPosition(i, j, k);
                    grid->Ex(i, j, k) = testFunc(coords);
                }
    }

    // Validation
    for (int x = 0; x < gridSize.x; x++)
        for (int y = 0; y < gridSize.y; y++)
            for (int z = 0; z < gridSize.z; z++)
            {
                Int3 index(x, y, z);
                FP3 coords = mainGrid->ExPosition(index.x, index.y, index.z);
                FP expected = testFunc(coords);

                Int3 local_index_origin;
                int r = getResponsibleSubGrid(local_index_origin, index, minCoords, gridStep, divisions, topology);

                //std::cout << index.x << ", " << index.z << ", " << r << std::endl;

                Int3 localIndex = index - local_index_origin;
                FP actual = subGrids[r]->Ex(localIndex);

                EXPECT_NEAR(expected, actual, 0.0001);
            }
}