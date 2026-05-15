#include "gtest/gtest.h"
#include <map>

#include "HiChi_MPI.h"

#define INVALID mpi::MPI_INVALID_RANK

using Int3 = pfc::Int3;

// TOPOLOGY
TEST(MPI_Topology, CanCreateTopology)
{
    Int3 sections(3, 4, 5);
    mpi::Topology::LoopType loopType = mpi::Topology::LoopType::None;

    ASSERT_NO_THROW( auto topology = mpi::Topology(sections, loopType, sections.x * sections.y * sections.z) );
}

TEST(MPI_Topology, ThrowsWhenNotEnoughNodes)
{
    Int3 sections(3, 4, 5);
    mpi::Topology::LoopType loopType = mpi::Topology::LoopType::None;

    ASSERT_ANY_THROW( auto topology = mpi::Topology(sections, loopType, 1) );
}

TEST(MPI_Topology, CorrectXLineTopologyNoLoop)
{
    Int3 sections(3, 1, 1);
    mpi::Topology::LoopType loopType = mpi::Topology::LoopType::None;
    int size = sections.x * sections.y * sections.z;

    auto topology = mpi::Topology(sections, loopType, size);

    std::vector<std::map<pfc::Int3, int>> neighbor_matrix = {
        // Rank 0
        {
            {Int3( 1,  0,  0), 1},
            {Int3(-1,  0,  0), INVALID},
            {Int3( 0,  1,  0), INVALID},
            {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), INVALID},
            {Int3( 0,  0, -1), INVALID}
        },

        // Rank 1
        {
            {Int3( 1,  0,  0), 2},
            {Int3(-1,  0,  0), 0},
            {Int3( 0,  1,  0), INVALID},
            {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), INVALID},
            {Int3( 0,  0, -1), INVALID}
        },

        // Rank 2
        {
            {Int3( 1,  0,  0), INVALID},
            {Int3(-1,  0,  0), 1},
            {Int3( 0,  1,  0), INVALID},
            {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), INVALID},
            {Int3( 0,  0, -1), INVALID}
        }
    };

    for (int r = 0; r < size; r++)
        for (auto& offset_neighbor : neighbor_matrix[r])
        {
            Int3 offset = offset_neighbor.first;
            int neighbor = offset_neighbor.second;

            EXPECT_EQ(topology.getNeighbor(r, offset), neighbor);
        }
}

TEST(MPI_Topology, CorrectYLineTopologyNoLoop)
{
    Int3 sections(1, 3, 1);
    mpi::Topology::LoopType loopType = mpi::Topology::LoopType::None;
    int size = sections.x * sections.y * sections.z;

    auto topology = mpi::Topology(sections, loopType, size);

    std::vector<std::map<pfc::Int3, int>> neighbor_matrix = {
        // Rank 0
        {
            {Int3( 1,  0,  0), INVALID},
            {Int3(-1,  0,  0), INVALID},
            {Int3( 0,  1,  0), 1},
            {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), INVALID},
            {Int3( 0,  0, -1), INVALID}
        },

        // Rank 1
        {
            {Int3( 1,  0,  0), INVALID},
            {Int3(-1,  0,  0), INVALID},
            {Int3( 0,  1,  0), 2},
            {Int3( 0, -1,  0), 0},
            {Int3( 0,  0,  1), INVALID},
            {Int3( 0,  0, -1), INVALID}
        },

        // Rank 2
        {
            {Int3( 1,  0,  0), INVALID},
            {Int3(-1,  0,  0), INVALID},
            {Int3( 0,  1,  0), INVALID},
            {Int3( 0, -1,  0), 1},
            {Int3( 0,  0,  1), INVALID},
            {Int3( 0,  0, -1), INVALID}
        }
    };

    for (int r = 0; r < size; r++)
        for (auto& offset_neighbor : neighbor_matrix[r])
        {
            Int3 offset = offset_neighbor.first;
            int neighbor = offset_neighbor.second;

            EXPECT_EQ(topology.getNeighbor(r, offset), neighbor);
        }
}

TEST(MPI_Topology, CorrectZLineTopologyNoLoop)
{
    Int3 sections(1, 1, 3);
    mpi::Topology::LoopType loopType = mpi::Topology::LoopType::None;
    int size = sections.x * sections.y * sections.z;

    auto topology = mpi::Topology(sections, loopType, size);

    std::vector<std::map<pfc::Int3, int>> neighbor_matrix = {
        // Rank 0
        {
            {Int3( 1,  0,  0), INVALID},
            {Int3(-1,  0,  0), INVALID},
            {Int3( 0,  1,  0), INVALID},
            {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), 1},
            {Int3( 0,  0, -1), INVALID}
        },

        // Rank 1
        {
            {Int3( 1,  0,  0), INVALID},
            {Int3(-1,  0,  0), INVALID},
            {Int3( 0,  1,  0), INVALID},
            {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), 2},
            {Int3( 0,  0, -1), 0}
        },

        // Rank 2
        {
            {Int3( 1,  0,  0), INVALID},
            {Int3(-1,  0,  0), INVALID},
            {Int3( 0,  1,  0), INVALID},
            {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), INVALID},
            {Int3( 0,  0, -1), 1}
        }
    };

    for (int r = 0; r < size; r++)
        for (auto& offset_neighbor : neighbor_matrix[r])
        {
            Int3 offset = offset_neighbor.first;
            int neighbor = offset_neighbor.second;

            EXPECT_EQ(topology.getNeighbor(r, offset), neighbor);
        }
}

TEST(MPI_Topology, CorrectVolumeTopologyNoLoop)
{
    Int3 sections(4, 3, 2);
    mpi::Topology::LoopType loopType = mpi::Topology::LoopType::None;
    int size = sections.x * sections.y * sections.z;

    auto topology = mpi::Topology(sections, loopType, size);

    std::vector<std::map<pfc::Int3, int>> neighbor_matrix = {
        // Rank 0
        {
            {Int3( 1,  0,  0), 6}, {Int3(-1,  0,  0), INVALID},
            {Int3( 0,  1,  0), 2}, {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), 1}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 1
        {
            {Int3( 1,  0,  0), 7}, {Int3(-1,  0,  0), INVALID},
            {Int3( 0,  1,  0), 3}, {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 0}
        },

        // Rank 2
        {
            {Int3( 1,  0,  0), 8}, {Int3(-1,  0,  0), INVALID},
            {Int3( 0,  1,  0), 4}, {Int3( 0, -1,  0), 0},
            {Int3( 0,  0,  1), 3}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 3
        {
            {Int3( 1,  0,  0), 9}, {Int3(-1,  0,  0), INVALID},
            {Int3( 0,  1,  0), 5}, {Int3( 0, -1,  0), 3},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 2}
        },

        // Rank 4
        {
            {Int3( 1,  0,  0), 10}, {Int3(-1,  0,  0), INVALID},
            {Int3( 0,  1,  0), INVALID}, {Int3( 0, -1,  0), 2},
            {Int3( 0,  0,  1), 5}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 5
        {
            {Int3( 1,  0,  0), 11}, {Int3(-1,  0,  0), INVALID},
            {Int3( 0,  1,  0), INVALID}, {Int3( 0, -1,  0), 3},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 4}
        },

        // Rank 6
        {
            {Int3( 1,  0,  0), 12}, {Int3(-1,  0,  0), 0},
            {Int3( 0,  1,  0), 8}, {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), 7}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 7
        {
            {Int3( 1,  0,  0), 13}, {Int3(-1,  0,  0), 1},
            {Int3( 0,  1,  0), 9}, {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 6}
        },

        // Rank 8
        {
            {Int3( 1,  0,  0), 14}, {Int3(-1,  0,  0), 2},
            {Int3( 0,  1,  0), 10}, {Int3( 0, -1,  0), 6},
            {Int3( 0,  0,  1), 9}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 9
        {
            {Int3( 1,  0,  0), 15}, {Int3(-1,  0,  0), 3},
            {Int3( 0,  1,  0), 11}, {Int3( 0, -1,  0), 7},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 8}
        },

        // Rank 10
        {
            {Int3( 1,  0,  0), 16}, {Int3(-1,  0,  0), 4},
            {Int3( 0,  1,  0), INVALID}, {Int3( 0, -1,  0), 8},
            {Int3( 0,  0,  1), 11}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 11
        {
            {Int3( 1,  0,  0), 17}, {Int3(-1,  0,  0), 5},
            {Int3( 0,  1,  0), INVALID}, {Int3( 0, -1,  0), 9},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 10}
        },


        // Rank 12
        {
            {Int3( 1,  0,  0), 18}, {Int3(-1,  0,  0), 6},
            {Int3( 0,  1,  0), 14}, {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), 13}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 13
        {
            {Int3( 1,  0,  0), 19}, {Int3(-1,  0,  0), 7},
            {Int3( 0,  1,  0), 15}, {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 12}
        },

        // Rank 14
        {
            {Int3( 1,  0,  0), 20}, {Int3(-1,  0,  0), 8},
            {Int3( 0,  1,  0), 16}, {Int3( 0, -1,  0), 12},
            {Int3( 0,  0,  1), 15}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 15
        {
            {Int3( 1,  0,  0), 21}, {Int3(-1,  0,  0), 9},
            {Int3( 0,  1,  0), 17}, {Int3( 0, -1,  0), 13},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 14}
        },

        // Rank 16
        {
            {Int3( 1,  0,  0), 22}, {Int3(-1,  0,  0), 10},
            {Int3( 0,  1,  0), INVALID}, {Int3( 0, -1,  0), 14},
            {Int3( 0,  0,  1), 17}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 17
        {
            {Int3( 1,  0,  0), 23}, {Int3(-1,  0,  0), 11},
            {Int3( 0,  1,  0), INVALID}, {Int3( 0, -1,  0), 15},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 16}
        },

        // Rank 18
        {
            {Int3( 1,  0,  0), INVALID}, {Int3(-1,  0,  0), 12},
            {Int3( 0,  1,  0), 20}, {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), 19}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 19
        {
            {Int3( 1,  0,  0), INVALID}, {Int3(-1,  0,  0), 13},
            {Int3( 0,  1,  0), 21}, {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 18}
        },

        // Rank 20
        {
            {Int3( 1,  0,  0), INVALID}, {Int3(-1,  0,  0), 14},
            {Int3( 0,  1,  0), 22}, {Int3( 0, -1,  0), 18},
            {Int3( 0,  0,  1), 21}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 21
        {
            {Int3( 1,  0,  0), INVALID}, {Int3(-1,  0,  0), 15},
            {Int3( 0,  1,  0), 23}, {Int3( 0, -1,  0), 19},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 20}
        },

        // Rank 22
        {
            {Int3( 1,  0,  0), INVALID}, {Int3(-1,  0,  0), 16},
            {Int3( 0,  1,  0), INVALID}, {Int3( 0, -1,  0), 20},
            {Int3( 0,  0,  1), 23}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 23
        {
            {Int3( 1,  0,  0), INVALID}, {Int3(-1,  0,  0), 17},
            {Int3( 0,  1,  0), INVALID}, {Int3( 0, -1,  0), 21},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 22}
        },
    };

    for (int r = 0; r < size; r++)
        for (auto& offset_neighbor : neighbor_matrix[r])
        {
            Int3 offset = offset_neighbor.first;
            int neighbor = offset_neighbor.second;

            EXPECT_EQ(topology.getNeighbor(r, offset), neighbor);
        }
}


TEST(MPI_Topology, LoopX)
{
    Int3 sections(4, 3, 2);
    mpi::Topology::LoopType loopType = mpi::Topology::LoopType::LoopX;
    int size = sections.x * sections.y * sections.z;

    auto topology = mpi::Topology(sections, loopType, size);

    std::vector<std::map<pfc::Int3, int>> neighbor_matrix = {
        // Rank 0
        {
            {Int3( 1,  0,  0), 6}, {Int3(-1,  0,  0), 18},
            {Int3( 0,  1,  0), 2}, {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), 1}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 1
        {
            {Int3( 1,  0,  0), 7}, {Int3(-1,  0,  0), 19},
            {Int3( 0,  1,  0), 3}, {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 0}
        },

        // Rank 2
        {
            {Int3( 1,  0,  0), 8}, {Int3(-1,  0,  0), 20},
            {Int3( 0,  1,  0), 4}, {Int3( 0, -1,  0), 0},
            {Int3( 0,  0,  1), 3}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 3
        {
            {Int3( 1,  0,  0), 9}, {Int3(-1,  0,  0), 21},
            {Int3( 0,  1,  0), 5}, {Int3( 0, -1,  0), 3},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 2}
        },

        // Rank 4
        {
            {Int3( 1,  0,  0), 10}, {Int3(-1,  0,  0), 22},
            {Int3( 0,  1,  0), INVALID}, {Int3( 0, -1,  0), 2},
            {Int3( 0,  0,  1), 5}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 5
        {
            {Int3( 1,  0,  0), 11}, {Int3(-1,  0,  0), 23},
            {Int3( 0,  1,  0), INVALID}, {Int3( 0, -1,  0), 3},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 4}
        },

        // Rank 6
        {
            {Int3( 1,  0,  0), 12}, {Int3(-1,  0,  0), 0},
            {Int3( 0,  1,  0), 8}, {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), 7}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 7
        {
            {Int3( 1,  0,  0), 13}, {Int3(-1,  0,  0), 1},
            {Int3( 0,  1,  0), 9}, {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 6}
        },

        // Rank 8
        {
            {Int3( 1,  0,  0), 14}, {Int3(-1,  0,  0), 2},
            {Int3( 0,  1,  0), 10}, {Int3( 0, -1,  0), 6},
            {Int3( 0,  0,  1), 9}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 9
        {
            {Int3( 1,  0,  0), 15}, {Int3(-1,  0,  0), 3},
            {Int3( 0,  1,  0), 11}, {Int3( 0, -1,  0), 7},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 8}
        },

        // Rank 10
        {
            {Int3( 1,  0,  0), 16}, {Int3(-1,  0,  0), 4},
            {Int3( 0,  1,  0), INVALID}, {Int3( 0, -1,  0), 8},
            {Int3( 0,  0,  1), 11}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 11
        {
            {Int3( 1,  0,  0), 17}, {Int3(-1,  0,  0), 5},
            {Int3( 0,  1,  0), INVALID}, {Int3( 0, -1,  0), 9},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 10}
        },


        // Rank 12
        {
            {Int3( 1,  0,  0), 18}, {Int3(-1,  0,  0), 6},
            {Int3( 0,  1,  0), 14}, {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), 13}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 13
        {
            {Int3( 1,  0,  0), 19}, {Int3(-1,  0,  0), 7},
            {Int3( 0,  1,  0), 15}, {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 12}
        },

        // Rank 14
        {
            {Int3( 1,  0,  0), 20}, {Int3(-1,  0,  0), 8},
            {Int3( 0,  1,  0), 16}, {Int3( 0, -1,  0), 12},
            {Int3( 0,  0,  1), 15}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 15
        {
            {Int3( 1,  0,  0), 21}, {Int3(-1,  0,  0), 9},
            {Int3( 0,  1,  0), 17}, {Int3( 0, -1,  0), 13},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 14}
        },

        // Rank 16
        {
            {Int3( 1,  0,  0), 22}, {Int3(-1,  0,  0), 10},
            {Int3( 0,  1,  0), INVALID}, {Int3( 0, -1,  0), 14},
            {Int3( 0,  0,  1), 17}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 17
        {
            {Int3( 1,  0,  0), 23}, {Int3(-1,  0,  0), 11},
            {Int3( 0,  1,  0), INVALID}, {Int3( 0, -1,  0), 15},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 16}
        },

        // Rank 18
        {
            {Int3( 1,  0,  0), 0}, {Int3(-1,  0,  0), 12},
            {Int3( 0,  1,  0), 20}, {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), 19}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 19
        {
            {Int3( 1,  0,  0), 1}, {Int3(-1,  0,  0), 13},
            {Int3( 0,  1,  0), 21}, {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 18}
        },

        // Rank 20
        {
            {Int3( 1,  0,  0), 2}, {Int3(-1,  0,  0), 14},
            {Int3( 0,  1,  0), 22}, {Int3( 0, -1,  0), 18},
            {Int3( 0,  0,  1), 21}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 21
        {
            {Int3( 1,  0,  0), 3}, {Int3(-1,  0,  0), 15},
            {Int3( 0,  1,  0), 23}, {Int3( 0, -1,  0), 19},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 20}
        },

        // Rank 22
        {
            {Int3( 1,  0,  0), 4}, {Int3(-1,  0,  0), 16},
            {Int3( 0,  1,  0), INVALID}, {Int3( 0, -1,  0), 20},
            {Int3( 0,  0,  1), 23}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 23
        {
            {Int3( 1,  0,  0), 5}, {Int3(-1,  0,  0), 17},
            {Int3( 0,  1,  0), INVALID}, {Int3( 0, -1,  0), 21},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 22}
        },
    };

    for (int r = 0; r < size; r++)
        for (auto& offset_neighbor : neighbor_matrix[r])
        {
            Int3 offset = offset_neighbor.first;
            int neighbor = offset_neighbor.second;

            EXPECT_EQ(topology.getNeighbor(r, offset), neighbor);
        }
}

TEST(MPI_Topology, LoopY)
{
    Int3 sections(4, 3, 2);
    mpi::Topology::LoopType loopType = mpi::Topology::LoopType::LoopY;
    int size = sections.x * sections.y * sections.z;

    auto topology = mpi::Topology(sections, loopType, size);

    std::vector<std::map<pfc::Int3, int>> neighbor_matrix = {
        // Rank 0
        {
            {Int3( 1,  0,  0), 6}, {Int3(-1,  0,  0), INVALID},
            {Int3( 0,  1,  0), 2}, {Int3( 0, -1,  0), 4},
            {Int3( 0,  0,  1), 1}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 1
        {
            {Int3( 1,  0,  0), 7}, {Int3(-1,  0,  0), INVALID},
            {Int3( 0,  1,  0), 3}, {Int3( 0, -1,  0), 5},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 0}
        },

        // Rank 2
        {
            {Int3( 1,  0,  0), 8}, {Int3(-1,  0,  0), INVALID},
            {Int3( 0,  1,  0), 4}, {Int3( 0, -1,  0), 0},
            {Int3( 0,  0,  1), 3}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 3
        {
            {Int3( 1,  0,  0), 9}, {Int3(-1,  0,  0), INVALID},
            {Int3( 0,  1,  0), 5}, {Int3( 0, -1,  0), 3},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 2}
        },

        // Rank 4
        {
            {Int3( 1,  0,  0), 10}, {Int3(-1,  0,  0), INVALID},
            {Int3( 0,  1,  0), 0}, {Int3( 0, -1,  0), 2},
            {Int3( 0,  0,  1), 5}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 5
        {
            {Int3( 1,  0,  0), 11}, {Int3(-1,  0,  0), INVALID},
            {Int3( 0,  1,  0), 1}, {Int3( 0, -1,  0), 3},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 4}
        },

        // Rank 6
        {
            {Int3( 1,  0,  0), 12}, {Int3(-1,  0,  0), 0},
            {Int3( 0,  1,  0), 8}, {Int3( 0, -1,  0), 10},
            {Int3( 0,  0,  1), 7}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 7
        {
            {Int3( 1,  0,  0), 13}, {Int3(-1,  0,  0), 1},
            {Int3( 0,  1,  0), 9}, {Int3( 0, -1,  0), 11},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 6}
        },

        // Rank 8
        {
            {Int3( 1,  0,  0), 14}, {Int3(-1,  0,  0), 2},
            {Int3( 0,  1,  0), 10}, {Int3( 0, -1,  0), 6},
            {Int3( 0,  0,  1), 9}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 9
        {
            {Int3( 1,  0,  0), 15}, {Int3(-1,  0,  0), 3},
            {Int3( 0,  1,  0), 11}, {Int3( 0, -1,  0), 7},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 8}
        },

        // Rank 10
        {
            {Int3( 1,  0,  0), 16}, {Int3(-1,  0,  0), 4},
            {Int3( 0,  1,  0), 6}, {Int3( 0, -1,  0), 8},
            {Int3( 0,  0,  1), 11}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 11
        {
            {Int3( 1,  0,  0), 17}, {Int3(-1,  0,  0), 5},
            {Int3( 0,  1,  0), 7}, {Int3( 0, -1,  0), 9},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 10}
        },


        // Rank 12
        {
            {Int3( 1,  0,  0), 18}, {Int3(-1,  0,  0), 6},
            {Int3( 0,  1,  0), 14}, {Int3( 0, -1,  0), 16},
            {Int3( 0,  0,  1), 13}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 13
        {
            {Int3( 1,  0,  0), 19}, {Int3(-1,  0,  0), 7},
            {Int3( 0,  1,  0), 15}, {Int3( 0, -1,  0), 17},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 12}
        },

        // Rank 14
        {
            {Int3( 1,  0,  0), 20}, {Int3(-1,  0,  0), 8},
            {Int3( 0,  1,  0), 16}, {Int3( 0, -1,  0), 12},
            {Int3( 0,  0,  1), 15}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 15
        {
            {Int3( 1,  0,  0), 21}, {Int3(-1,  0,  0), 9},
            {Int3( 0,  1,  0), 17}, {Int3( 0, -1,  0), 13},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 14}
        },

        // Rank 16
        {
            {Int3( 1,  0,  0), 22}, {Int3(-1,  0,  0), 10},
            {Int3( 0,  1,  0), 12}, {Int3( 0, -1,  0), 14},
            {Int3( 0,  0,  1), 17}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 17
        {
            {Int3( 1,  0,  0), 23}, {Int3(-1,  0,  0), 11},
            {Int3( 0,  1,  0), 13}, {Int3( 0, -1,  0), 15},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 16}
        },

        // Rank 18
        {
            {Int3( 1,  0,  0), INVALID}, {Int3(-1,  0,  0), 12},
            {Int3( 0,  1,  0), 20}, {Int3( 0, -1,  0), 22},
            {Int3( 0,  0,  1), 19}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 19
        {
            {Int3( 1,  0,  0), INVALID}, {Int3(-1,  0,  0), 13},
            {Int3( 0,  1,  0), 21}, {Int3( 0, -1,  0), 23},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 18}
        },

        // Rank 20
        {
            {Int3( 1,  0,  0), INVALID}, {Int3(-1,  0,  0), 14},
            {Int3( 0,  1,  0), 22}, {Int3( 0, -1,  0), 18},
            {Int3( 0,  0,  1), 21}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 21
        {
            {Int3( 1,  0,  0), INVALID}, {Int3(-1,  0,  0), 15},
            {Int3( 0,  1,  0), 23}, {Int3( 0, -1,  0), 19},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 20}
        },

        // Rank 22
        {
            {Int3( 1,  0,  0), INVALID}, {Int3(-1,  0,  0), 16},
            {Int3( 0,  1,  0), 18}, {Int3( 0, -1,  0), 20},
            {Int3( 0,  0,  1), 23}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 23
        {
            {Int3( 1,  0,  0), INVALID}, {Int3(-1,  0,  0), 17},
            {Int3( 0,  1,  0), 19}, {Int3( 0, -1,  0), 21},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 22}
        },
    };

    for (int r = 0; r < size; r++)
        for (auto& offset_neighbor : neighbor_matrix[r])
        {
            Int3 offset = offset_neighbor.first;
            int neighbor = offset_neighbor.second;

            EXPECT_EQ(topology.getNeighbor(r, offset), neighbor);
        }
}

TEST(MPI_Topology, LoopZ)
{
    Int3 sections(4, 3, 2);
    mpi::Topology::LoopType loopType = mpi::Topology::LoopType::LoopZ;
    int size = sections.x * sections.y * sections.z;

    auto topology = mpi::Topology(sections, loopType, size);

    std::vector<std::map<pfc::Int3, int>> neighbor_matrix = {
        // Rank 0
        {
            {Int3( 1,  0,  0), 6}, {Int3(-1,  0,  0), INVALID},
            {Int3( 0,  1,  0), 2}, {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), 1}, {Int3( 0,  0, -1), 1}
        },

        // Rank 1
        {
            {Int3( 1,  0,  0), 7}, {Int3(-1,  0,  0), INVALID},
            {Int3( 0,  1,  0), 3}, {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), 0}, {Int3( 0,  0, -1), 0}
        },

        // Rank 2
        {
            {Int3( 1,  0,  0), 8}, {Int3(-1,  0,  0), INVALID},
            {Int3( 0,  1,  0), 4}, {Int3( 0, -1,  0), 0},
            {Int3( 0,  0,  1), 3}, {Int3( 0,  0, -1), 3}
        },

        // Rank 3
        {
            {Int3( 1,  0,  0), 9}, {Int3(-1,  0,  0), INVALID},
            {Int3( 0,  1,  0), 5}, {Int3( 0, -1,  0), 3},
            {Int3( 0,  0,  1), 2}, {Int3( 0,  0, -1), 2}
        },

        // Rank 4
        {
            {Int3( 1,  0,  0), 10}, {Int3(-1,  0,  0), INVALID},
            {Int3( 0,  1,  0), INVALID}, {Int3( 0, -1,  0), 2},
            {Int3( 0,  0,  1), 5}, {Int3( 0,  0, -1), 5}
        },

        // Rank 5
        {
            {Int3( 1,  0,  0), 11}, {Int3(-1,  0,  0), INVALID},
            {Int3( 0,  1,  0), INVALID}, {Int3( 0, -1,  0), 3},
            {Int3( 0,  0,  1), 4}, {Int3( 0,  0, -1), 4}
        },

        // Rank 6
        {
            {Int3( 1,  0,  0), 12}, {Int3(-1,  0,  0), 0},
            {Int3( 0,  1,  0), 8}, {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), 7}, {Int3( 0,  0, -1), 7}
        },

        // Rank 7
        {
            {Int3( 1,  0,  0), 13}, {Int3(-1,  0,  0), 1},
            {Int3( 0,  1,  0), 9}, {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), 6}, {Int3( 0,  0, -1), 6}
        },

        // Rank 8
        {
            {Int3( 1,  0,  0), 14}, {Int3(-1,  0,  0), 2},
            {Int3( 0,  1,  0), 10}, {Int3( 0, -1,  0), 6},
            {Int3( 0,  0,  1), 9}, {Int3( 0,  0, -1), 9}
        },

        // Rank 9
        {
            {Int3( 1,  0,  0), 15}, {Int3(-1,  0,  0), 3},
            {Int3( 0,  1,  0), 11}, {Int3( 0, -1,  0), 7},
            {Int3( 0,  0,  1), 8}, {Int3( 0,  0, -1), 8}
        },

        // Rank 10
        {
            {Int3( 1,  0,  0), 16}, {Int3(-1,  0,  0), 4},
            {Int3( 0,  1,  0), INVALID}, {Int3( 0, -1,  0), 8},
            {Int3( 0,  0,  1), 11}, {Int3( 0,  0, -1), 11}
        },

        // Rank 11
        {
            {Int3( 1,  0,  0), 17}, {Int3(-1,  0,  0), 5},
            {Int3( 0,  1,  0), INVALID}, {Int3( 0, -1,  0), 9},
            {Int3( 0,  0,  1), 10}, {Int3( 0,  0, -1), 10}
        },


        // Rank 12
        {
            {Int3( 1,  0,  0), 18}, {Int3(-1,  0,  0), 6},
            {Int3( 0,  1,  0), 14}, {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), 13}, {Int3( 0,  0, -1), 13}
        },

        // Rank 13
        {
            {Int3( 1,  0,  0), 19}, {Int3(-1,  0,  0), 7},
            {Int3( 0,  1,  0), 15}, {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), 12}, {Int3( 0,  0, -1), 12}
        },

        // Rank 14
        {
            {Int3( 1,  0,  0), 20}, {Int3(-1,  0,  0), 8},
            {Int3( 0,  1,  0), 16}, {Int3( 0, -1,  0), 12},
            {Int3( 0,  0,  1), 15}, {Int3( 0,  0, -1), 15}
        },

        // Rank 15
        {
            {Int3( 1,  0,  0), 21}, {Int3(-1,  0,  0), 9},
            {Int3( 0,  1,  0), 17}, {Int3( 0, -1,  0), 13},
            {Int3( 0,  0,  1), 14}, {Int3( 0,  0, -1), 14}
        },

        // Rank 16
        {
            {Int3( 1,  0,  0), 22}, {Int3(-1,  0,  0), 10},
            {Int3( 0,  1,  0), INVALID}, {Int3( 0, -1,  0), 14},
            {Int3( 0,  0,  1), 17}, {Int3( 0,  0, -1), 17}
        },

        // Rank 17
        {
            {Int3( 1,  0,  0), 23}, {Int3(-1,  0,  0), 11},
            {Int3( 0,  1,  0), INVALID}, {Int3( 0, -1,  0), 15},
            {Int3( 0,  0,  1), 16}, {Int3( 0,  0, -1), 16}
        },

        // Rank 18
        {
            {Int3( 1,  0,  0), INVALID}, {Int3(-1,  0,  0), 12},
            {Int3( 0,  1,  0), 20}, {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), 19}, {Int3( 0,  0, -1), 19}
        },

        // Rank 19
        {
            {Int3( 1,  0,  0), INVALID}, {Int3(-1,  0,  0), 13},
            {Int3( 0,  1,  0), 21}, {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), 18}, {Int3( 0,  0, -1), 18}
        },

        // Rank 20
        {
            {Int3( 1,  0,  0), INVALID}, {Int3(-1,  0,  0), 14},
            {Int3( 0,  1,  0), 22}, {Int3( 0, -1,  0), 18},
            {Int3( 0,  0,  1), 21}, {Int3( 0,  0, -1), 21}
        },

        // Rank 21
        {
            {Int3( 1,  0,  0), INVALID}, {Int3(-1,  0,  0), 15},
            {Int3( 0,  1,  0), 23}, {Int3( 0, -1,  0), 19},
            {Int3( 0,  0,  1), 20}, {Int3( 0,  0, -1), 20}
        },

        // Rank 22
        {
            {Int3( 1,  0,  0), INVALID}, {Int3(-1,  0,  0), 16},
            {Int3( 0,  1,  0), INVALID}, {Int3( 0, -1,  0), 20},
            {Int3( 0,  0,  1), 23}, {Int3( 0,  0, -1), 23}
        },

        // Rank 23
        {
            {Int3( 1,  0,  0), INVALID}, {Int3(-1,  0,  0), 17},
            {Int3( 0,  1,  0), INVALID}, {Int3( 0, -1,  0), 21},
            {Int3( 0,  0,  1), 22}, {Int3( 0,  0, -1), 22}
        },
    };

    for (int r = 0; r < size; r++)
        for (auto& offset_neighbor : neighbor_matrix[r])
        {
            Int3 offset = offset_neighbor.first;
            int neighbor = offset_neighbor.second;

            EXPECT_EQ(topology.getNeighbor(r, offset), neighbor);
        }
}

TEST(MPI_Topology, LoopXY)
{
    Int3 sections(4, 3, 2);
    mpi::Topology::LoopType loopType = mpi::Topology::LoopType::LoopXY;
    int size = sections.x * sections.y * sections.z;

    auto topology = mpi::Topology(sections, loopType, size);

    std::vector<std::map<pfc::Int3, int>> neighbor_matrix = {
        // Rank 0
        {
            {Int3( 1,  0,  0), 6}, {Int3(-1,  0,  0), 18},
            {Int3( 0,  1,  0), 2}, {Int3( 0, -1,  0), 4},
            {Int3( 0,  0,  1), 1}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 1
        {
            {Int3( 1,  0,  0), 7}, {Int3(-1,  0,  0), 19},
            {Int3( 0,  1,  0), 3}, {Int3( 0, -1,  0), 5},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 0}
        },

        // Rank 2
        {
            {Int3( 1,  0,  0), 8}, {Int3(-1,  0,  0), 20},
            {Int3( 0,  1,  0), 4}, {Int3( 0, -1,  0), 0},
            {Int3( 0,  0,  1), 3}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 3
        {
            {Int3( 1,  0,  0), 9}, {Int3(-1,  0,  0), 21},
            {Int3( 0,  1,  0), 5}, {Int3( 0, -1,  0), 3},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 2}
        },

        // Rank 4
        {
            {Int3( 1,  0,  0), 10}, {Int3(-1,  0,  0), 22},
            {Int3( 0,  1,  0), 0}, {Int3( 0, -1,  0), 2},
            {Int3( 0,  0,  1), 5}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 5
        {
            {Int3( 1,  0,  0), 11}, {Int3(-1,  0,  0), 23},
            {Int3( 0,  1,  0), 1}, {Int3( 0, -1,  0), 3},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 4}
        },

        // Rank 6
        {
            {Int3( 1,  0,  0), 12}, {Int3(-1,  0,  0), 0},
            {Int3( 0,  1,  0), 8}, {Int3( 0, -1,  0), 10},
            {Int3( 0,  0,  1), 7}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 7
        {
            {Int3( 1,  0,  0), 13}, {Int3(-1,  0,  0), 1},
            {Int3( 0,  1,  0), 9}, {Int3( 0, -1,  0), 11},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 6}
        },

        // Rank 8
        {
            {Int3( 1,  0,  0), 14}, {Int3(-1,  0,  0), 2},
            {Int3( 0,  1,  0), 10}, {Int3( 0, -1,  0), 6},
            {Int3( 0,  0,  1), 9}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 9
        {
            {Int3( 1,  0,  0), 15}, {Int3(-1,  0,  0), 3},
            {Int3( 0,  1,  0), 11}, {Int3( 0, -1,  0), 7},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 8}
        },

        // Rank 10
        {
            {Int3( 1,  0,  0), 16}, {Int3(-1,  0,  0), 4},
            {Int3( 0,  1,  0), 6}, {Int3( 0, -1,  0), 8},
            {Int3( 0,  0,  1), 11}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 11
        {
            {Int3( 1,  0,  0), 17}, {Int3(-1,  0,  0), 5},
            {Int3( 0,  1,  0), 7}, {Int3( 0, -1,  0), 9},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 10}
        },


        // Rank 12
        {
            {Int3( 1,  0,  0), 18}, {Int3(-1,  0,  0), 6},
            {Int3( 0,  1,  0), 14}, {Int3( 0, -1,  0), 16},
            {Int3( 0,  0,  1), 13}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 13
        {
            {Int3( 1,  0,  0), 19}, {Int3(-1,  0,  0), 7},
            {Int3( 0,  1,  0), 15}, {Int3( 0, -1,  0), 17},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 12}
        },

        // Rank 14
        {
            {Int3( 1,  0,  0), 20}, {Int3(-1,  0,  0), 8},
            {Int3( 0,  1,  0), 16}, {Int3( 0, -1,  0), 12},
            {Int3( 0,  0,  1), 15}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 15
        {
            {Int3( 1,  0,  0), 21}, {Int3(-1,  0,  0), 9},
            {Int3( 0,  1,  0), 17}, {Int3( 0, -1,  0), 13},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 14}
        },

        // Rank 16
        {
            {Int3( 1,  0,  0), 22}, {Int3(-1,  0,  0), 10},
            {Int3( 0,  1,  0), 12}, {Int3( 0, -1,  0), 14},
            {Int3( 0,  0,  1), 17}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 17
        {
            {Int3( 1,  0,  0), 23}, {Int3(-1,  0,  0), 11},
            {Int3( 0,  1,  0), 13}, {Int3( 0, -1,  0), 15},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 16}
        },

        // Rank 18
        {
            {Int3( 1,  0,  0), 0}, {Int3(-1,  0,  0), 12},
            {Int3( 0,  1,  0), 20}, {Int3( 0, -1,  0), 22},
            {Int3( 0,  0,  1), 19}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 19
        {
            {Int3( 1,  0,  0), 1}, {Int3(-1,  0,  0), 13},
            {Int3( 0,  1,  0), 21}, {Int3( 0, -1,  0), 23},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 18}
        },

        // Rank 20
        {
            {Int3( 1,  0,  0), 2}, {Int3(-1,  0,  0), 14},
            {Int3( 0,  1,  0), 22}, {Int3( 0, -1,  0), 18},
            {Int3( 0,  0,  1), 21}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 21
        {
            {Int3( 1,  0,  0), 3}, {Int3(-1,  0,  0), 15},
            {Int3( 0,  1,  0), 23}, {Int3( 0, -1,  0), 19},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 20}
        },

        // Rank 22
        {
            {Int3( 1,  0,  0), 4}, {Int3(-1,  0,  0), 16},
            {Int3( 0,  1,  0), 18}, {Int3( 0, -1,  0), 20},
            {Int3( 0,  0,  1), 23}, {Int3( 0,  0, -1), INVALID}
        },

        // Rank 23
        {
            {Int3( 1,  0,  0), 5}, {Int3(-1,  0,  0), 17},
            {Int3( 0,  1,  0), 19}, {Int3( 0, -1,  0), 21},
            {Int3( 0,  0,  1), INVALID}, {Int3( 0,  0, -1), 22}
        },
    };

    for (int r = 0; r < size; r++)
        for (auto& offset_neighbor : neighbor_matrix[r])
        {
            Int3 offset = offset_neighbor.first;
            int neighbor = offset_neighbor.second;

            EXPECT_EQ(topology.getNeighbor(r, offset), neighbor);
        }
}

TEST(MPI_Topology, LoopYZ)
{
    Int3 sections(4, 3, 2);
    mpi::Topology::LoopType loopType = mpi::Topology::LoopType::LoopYZ;
    int size = sections.x * sections.y * sections.z;

    auto topology = mpi::Topology(sections, loopType, size);

    std::vector<std::map<pfc::Int3, int>> neighbor_matrix = {
        // Rank 0
        {
            {Int3( 1,  0,  0), 6}, {Int3(-1,  0,  0), INVALID},
            {Int3( 0,  1,  0), 2}, {Int3( 0, -1,  0), 4},
            {Int3( 0,  0,  1), 1}, {Int3( 0,  0, -1), 1}
        },

        // Rank 1
        {
            {Int3( 1,  0,  0), 7}, {Int3(-1,  0,  0), INVALID},
            {Int3( 0,  1,  0), 3}, {Int3( 0, -1,  0), 5},
            {Int3( 0,  0,  1), 0}, {Int3( 0,  0, -1), 0}
        },

        // Rank 2
        {
            {Int3( 1,  0,  0), 8}, {Int3(-1,  0,  0), INVALID},
            {Int3( 0,  1,  0), 4}, {Int3( 0, -1,  0), 0},
            {Int3( 0,  0,  1), 3}, {Int3( 0,  0, -1), 3}
        },

        // Rank 3
        {
            {Int3( 1,  0,  0), 9}, {Int3(-1,  0,  0), INVALID},
            {Int3( 0,  1,  0), 5}, {Int3( 0, -1,  0), 3},
            {Int3( 0,  0,  1), 2}, {Int3( 0,  0, -1), 2}
        },

        // Rank 4
        {
            {Int3( 1,  0,  0), 10}, {Int3(-1,  0,  0), INVALID},
            {Int3( 0,  1,  0), 0}, {Int3( 0, -1,  0), 2},
            {Int3( 0,  0,  1), 5}, {Int3( 0,  0, -1), 5}
        },

        // Rank 5
        {
            {Int3( 1,  0,  0), 11}, {Int3(-1,  0,  0), INVALID},
            {Int3( 0,  1,  0), 1}, {Int3( 0, -1,  0), 3},
            {Int3( 0,  0,  1), 4}, {Int3( 0,  0, -1), 4}
        },

        // Rank 6
        {
            {Int3( 1,  0,  0), 12}, {Int3(-1,  0,  0), 0},
            {Int3( 0,  1,  0), 8}, {Int3( 0, -1,  0), 10},
            {Int3( 0,  0,  1), 7}, {Int3( 0,  0, -1), 7}
        },

        // Rank 7
        {
            {Int3( 1,  0,  0), 13}, {Int3(-1,  0,  0), 1},
            {Int3( 0,  1,  0), 9}, {Int3( 0, -1,  0), 11},
            {Int3( 0,  0,  1), 6}, {Int3( 0,  0, -1), 6}
        },

        // Rank 8
        {
            {Int3( 1,  0,  0), 14}, {Int3(-1,  0,  0), 2},
            {Int3( 0,  1,  0), 10}, {Int3( 0, -1,  0), 6},
            {Int3( 0,  0,  1), 9}, {Int3( 0,  0, -1), 9}
        },

        // Rank 9
        {
            {Int3( 1,  0,  0), 15}, {Int3(-1,  0,  0), 3},
            {Int3( 0,  1,  0), 11}, {Int3( 0, -1,  0), 7},
            {Int3( 0,  0,  1), 8}, {Int3( 0,  0, -1), 8}
        },

        // Rank 10
        {
            {Int3( 1,  0,  0), 16}, {Int3(-1,  0,  0), 4},
            {Int3( 0,  1,  0), 6}, {Int3( 0, -1,  0), 8},
            {Int3( 0,  0,  1), 11}, {Int3( 0,  0, -1), 11}
        },

        // Rank 11
        {
            {Int3( 1,  0,  0), 17}, {Int3(-1,  0,  0), 5},
            {Int3( 0,  1,  0), 7}, {Int3( 0, -1,  0), 9},
            {Int3( 0,  0,  1), 10}, {Int3( 0,  0, -1), 10}
        },


        // Rank 12
        {
            {Int3( 1,  0,  0), 18}, {Int3(-1,  0,  0), 6},
            {Int3( 0,  1,  0), 14}, {Int3( 0, -1,  0), 16},
            {Int3( 0,  0,  1), 13}, {Int3( 0,  0, -1), 13}
        },

        // Rank 13
        {
            {Int3( 1,  0,  0), 19}, {Int3(-1,  0,  0), 7},
            {Int3( 0,  1,  0), 15}, {Int3( 0, -1,  0), 17},
            {Int3( 0,  0,  1), 12}, {Int3( 0,  0, -1), 12}
        },

        // Rank 14
        {
            {Int3( 1,  0,  0), 20}, {Int3(-1,  0,  0), 8},
            {Int3( 0,  1,  0), 16}, {Int3( 0, -1,  0), 12},
            {Int3( 0,  0,  1), 15}, {Int3( 0,  0, -1), 15}
        },

        // Rank 15
        {
            {Int3( 1,  0,  0), 21}, {Int3(-1,  0,  0), 9},
            {Int3( 0,  1,  0), 17}, {Int3( 0, -1,  0), 13},
            {Int3( 0,  0,  1), 14}, {Int3( 0,  0, -1), 14}
        },

        // Rank 16
        {
            {Int3( 1,  0,  0), 22}, {Int3(-1,  0,  0), 10},
            {Int3( 0,  1,  0), 12}, {Int3( 0, -1,  0), 14},
            {Int3( 0,  0,  1), 17}, {Int3( 0,  0, -1), 17}
        },

        // Rank 17
        {
            {Int3( 1,  0,  0), 23}, {Int3(-1,  0,  0), 11},
            {Int3( 0,  1,  0), 13}, {Int3( 0, -1,  0), 15},
            {Int3( 0,  0,  1), 16}, {Int3( 0,  0, -1), 16}
        },

        // Rank 18
        {
            {Int3( 1,  0,  0), INVALID}, {Int3(-1,  0,  0), 12},
            {Int3( 0,  1,  0), 20}, {Int3( 0, -1,  0), 22},
            {Int3( 0,  0,  1), 19}, {Int3( 0,  0, -1), 19}
        },

        // Rank 19
        {
            {Int3( 1,  0,  0), INVALID}, {Int3(-1,  0,  0), 13},
            {Int3( 0,  1,  0), 21}, {Int3( 0, -1,  0), 23},
            {Int3( 0,  0,  1), 18}, {Int3( 0,  0, -1), 18}
        },

        // Rank 20
        {
            {Int3( 1,  0,  0), INVALID}, {Int3(-1,  0,  0), 14},
            {Int3( 0,  1,  0), 22}, {Int3( 0, -1,  0), 18},
            {Int3( 0,  0,  1), 21}, {Int3( 0,  0, -1), 21}
        },

        // Rank 21
        {
            {Int3( 1,  0,  0), INVALID}, {Int3(-1,  0,  0), 15},
            {Int3( 0,  1,  0), 23}, {Int3( 0, -1,  0), 19},
            {Int3( 0,  0,  1), 20}, {Int3( 0,  0, -1), 20}
        },

        // Rank 22
        {
            {Int3( 1,  0,  0), INVALID}, {Int3(-1,  0,  0), 16},
            {Int3( 0,  1,  0), 18}, {Int3( 0, -1,  0), 20},
            {Int3( 0,  0,  1), 23}, {Int3( 0,  0, -1), 23}
        },

        // Rank 23
        {
            {Int3( 1,  0,  0), INVALID}, {Int3(-1,  0,  0), 17},
            {Int3( 0,  1,  0), 19}, {Int3( 0, -1,  0), 21},
            {Int3( 0,  0,  1), 22}, {Int3( 0,  0, -1), 22}
        },
    };

    for (int r = 0; r < size; r++)
        for (auto& offset_neighbor : neighbor_matrix[r])
        {
            Int3 offset = offset_neighbor.first;
            int neighbor = offset_neighbor.second;

            EXPECT_EQ(topology.getNeighbor(r, offset), neighbor);
        }
}

TEST(MPI_Topology, LoopXZ)
{
    Int3 sections(4, 3, 2);
    mpi::Topology::LoopType loopType = mpi::Topology::LoopType::LoopXZ;
    int size = sections.x * sections.y * sections.z;

    auto topology = mpi::Topology(sections, loopType, size);

    std::vector<std::map<pfc::Int3, int>> neighbor_matrix = {
        // Rank 0
        {
            {Int3( 1,  0,  0), 6}, {Int3(-1,  0,  0), 18},
            {Int3( 0,  1,  0), 2}, {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), 1}, {Int3( 0,  0, -1), 1}
        },

        // Rank 1
        {
            {Int3( 1,  0,  0), 7}, {Int3(-1,  0,  0), 19},
            {Int3( 0,  1,  0), 3}, {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), 0}, {Int3( 0,  0, -1), 0}
        },

        // Rank 2
        {
            {Int3( 1,  0,  0), 8}, {Int3(-1,  0,  0), 20},
            {Int3( 0,  1,  0), 4}, {Int3( 0, -1,  0), 0},
            {Int3( 0,  0,  1), 3}, {Int3( 0,  0, -1), 3}
        },

        // Rank 3
        {
            {Int3( 1,  0,  0), 9}, {Int3(-1,  0,  0), 21},
            {Int3( 0,  1,  0), 5}, {Int3( 0, -1,  0), 3},
            {Int3( 0,  0,  1), 2}, {Int3( 0,  0, -1), 2}
        },

        // Rank 4
        {
            {Int3( 1,  0,  0), 10}, {Int3(-1,  0,  0), 22},
            {Int3( 0,  1,  0), INVALID}, {Int3( 0, -1,  0), 2},
            {Int3( 0,  0,  1), 5}, {Int3( 0,  0, -1), 5}
        },

        // Rank 5
        {
            {Int3( 1,  0,  0), 11}, {Int3(-1,  0,  0), 23},
            {Int3( 0,  1,  0), INVALID}, {Int3( 0, -1,  0), 3},
            {Int3( 0,  0,  1), 4}, {Int3( 0,  0, -1), 4}
        },

        // Rank 6
        {
            {Int3( 1,  0,  0), 12}, {Int3(-1,  0,  0), 0},
            {Int3( 0,  1,  0), 8}, {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), 7}, {Int3( 0,  0, -1), 7}
        },

        // Rank 7
        {
            {Int3( 1,  0,  0), 13}, {Int3(-1,  0,  0), 1},
            {Int3( 0,  1,  0), 9}, {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), 6}, {Int3( 0,  0, -1), 6}
        },

        // Rank 8
        {
            {Int3( 1,  0,  0), 14}, {Int3(-1,  0,  0), 2},
            {Int3( 0,  1,  0), 10}, {Int3( 0, -1,  0), 6},
            {Int3( 0,  0,  1), 9}, {Int3( 0,  0, -1), 9}
        },

        // Rank 9
        {
            {Int3( 1,  0,  0), 15}, {Int3(-1,  0,  0), 3},
            {Int3( 0,  1,  0), 11}, {Int3( 0, -1,  0), 7},
            {Int3( 0,  0,  1), 8}, {Int3( 0,  0, -1), 8}
        },

        // Rank 10
        {
            {Int3( 1,  0,  0), 16}, {Int3(-1,  0,  0), 4},
            {Int3( 0,  1,  0), INVALID}, {Int3( 0, -1,  0), 8},
            {Int3( 0,  0,  1), 11}, {Int3( 0,  0, -1), 11}
        },

        // Rank 11
        {
            {Int3( 1,  0,  0), 17}, {Int3(-1,  0,  0), 5},
            {Int3( 0,  1,  0), INVALID}, {Int3( 0, -1,  0), 9},
            {Int3( 0,  0,  1), 10}, {Int3( 0,  0, -1), 10}
        },


        // Rank 12
        {
            {Int3( 1,  0,  0), 18}, {Int3(-1,  0,  0), 6},
            {Int3( 0,  1,  0), 14}, {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), 13}, {Int3( 0,  0, -1), 13}
        },

        // Rank 13
        {
            {Int3( 1,  0,  0), 19}, {Int3(-1,  0,  0), 7},
            {Int3( 0,  1,  0), 15}, {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), 12}, {Int3( 0,  0, -1), 12}
        },

        // Rank 14
        {
            {Int3( 1,  0,  0), 20}, {Int3(-1,  0,  0), 8},
            {Int3( 0,  1,  0), 16}, {Int3( 0, -1,  0), 12},
            {Int3( 0,  0,  1), 15}, {Int3( 0,  0, -1), 15}
        },

        // Rank 15
        {
            {Int3( 1,  0,  0), 21}, {Int3(-1,  0,  0), 9},
            {Int3( 0,  1,  0), 17}, {Int3( 0, -1,  0), 13},
            {Int3( 0,  0,  1), 14}, {Int3( 0,  0, -1), 14}
        },

        // Rank 16
        {
            {Int3( 1,  0,  0), 22}, {Int3(-1,  0,  0), 10},
            {Int3( 0,  1,  0), INVALID}, {Int3( 0, -1,  0), 14},
            {Int3( 0,  0,  1), 17}, {Int3( 0,  0, -1), 17}
        },

        // Rank 17
        {
            {Int3( 1,  0,  0), 23}, {Int3(-1,  0,  0), 11},
            {Int3( 0,  1,  0), INVALID}, {Int3( 0, -1,  0), 15},
            {Int3( 0,  0,  1), 16}, {Int3( 0,  0, -1), 16}
        },

        // Rank 18
        {
            {Int3( 1,  0,  0), 0}, {Int3(-1,  0,  0), 12},
            {Int3( 0,  1,  0), 20}, {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), 19}, {Int3( 0,  0, -1), 19}
        },

        // Rank 19
        {
            {Int3( 1,  0,  0), 1}, {Int3(-1,  0,  0), 13},
            {Int3( 0,  1,  0), 21}, {Int3( 0, -1,  0), INVALID},
            {Int3( 0,  0,  1), 18}, {Int3( 0,  0, -1), 18}
        },

        // Rank 20
        {
            {Int3( 1,  0,  0), 2}, {Int3(-1,  0,  0), 14},
            {Int3( 0,  1,  0), 22}, {Int3( 0, -1,  0), 18},
            {Int3( 0,  0,  1), 21}, {Int3( 0,  0, -1), 21}
        },

        // Rank 21
        {
            {Int3( 1,  0,  0), 3}, {Int3(-1,  0,  0), 15},
            {Int3( 0,  1,  0), 23}, {Int3( 0, -1,  0), 19},
            {Int3( 0,  0,  1), 20}, {Int3( 0,  0, -1), 20}
        },

        // Rank 22
        {
            {Int3( 1,  0,  0), 4}, {Int3(-1,  0,  0), 16},
            {Int3( 0,  1,  0), INVALID}, {Int3( 0, -1,  0), 20},
            {Int3( 0,  0,  1), 23}, {Int3( 0,  0, -1), 23}
        },

        // Rank 23
        {
            {Int3( 1,  0,  0), 5}, {Int3(-1,  0,  0), 17},
            {Int3( 0,  1,  0), INVALID}, {Int3( 0, -1,  0), 21},
            {Int3( 0,  0,  1), 22}, {Int3( 0,  0, -1), 22}
        },
    };

    for (int r = 0; r < size; r++)
        for (auto& offset_neighbor : neighbor_matrix[r])
        {
            Int3 offset = offset_neighbor.first;
            int neighbor = offset_neighbor.second;

            EXPECT_EQ(topology.getNeighbor(r, offset), neighbor);
        }
}

TEST(MPI_Topology, LoopXYZ)
{
    Int3 sections(4, 3, 2);
    mpi::Topology::LoopType loopType = mpi::Topology::LoopType::LoopXYZ;
    int size = sections.x * sections.y * sections.z;

    auto topology = mpi::Topology(sections, loopType, size);

    std::vector<std::map<pfc::Int3, int>> neighbor_matrix = {
        // Rank 0
        {
            {Int3( 1,  0,  0), 6}, {Int3(-1,  0,  0), 18},
            {Int3( 0,  1,  0), 2}, {Int3( 0, -1,  0), 4},
            {Int3( 0,  0,  1), 1}, {Int3( 0,  0, -1), 1}
        },

        // Rank 1
        {
            {Int3( 1,  0,  0), 7}, {Int3(-1,  0,  0), 19},
            {Int3( 0,  1,  0), 3}, {Int3( 0, -1,  0), 5},
            {Int3( 0,  0,  1), 0}, {Int3( 0,  0, -1), 0}
        },

        // Rank 2
        {
            {Int3( 1,  0,  0), 8}, {Int3(-1,  0,  0), 20},
            {Int3( 0,  1,  0), 4}, {Int3( 0, -1,  0), 0},
            {Int3( 0,  0,  1), 3}, {Int3( 0,  0, -1), 3}
        },

        // Rank 3
        {
            {Int3( 1,  0,  0), 9}, {Int3(-1,  0,  0), 21},
            {Int3( 0,  1,  0), 5}, {Int3( 0, -1,  0), 3},
            {Int3( 0,  0,  1), 2}, {Int3( 0,  0, -1), 2}
        },

        // Rank 4
        {
            {Int3( 1,  0,  0), 10}, {Int3(-1,  0,  0), 22},
            {Int3( 0,  1,  0), 0}, {Int3( 0, -1,  0), 2},
            {Int3( 0,  0,  1), 5}, {Int3( 0,  0, -1), 5}
        },

        // Rank 5
        {
            {Int3( 1,  0,  0), 11}, {Int3(-1,  0,  0), 23},
            {Int3( 0,  1,  0), 1}, {Int3( 0, -1,  0), 3},
            {Int3( 0,  0,  1), 4}, {Int3( 0,  0, -1), 4}
        },

        // Rank 6
        {
            {Int3( 1,  0,  0), 12}, {Int3(-1,  0,  0), 0},
            {Int3( 0,  1,  0), 8}, {Int3( 0, -1,  0), 10},
            {Int3( 0,  0,  1), 7}, {Int3( 0,  0, -1), 7}
        },

        // Rank 7
        {
            {Int3( 1,  0,  0), 13}, {Int3(-1,  0,  0), 1},
            {Int3( 0,  1,  0), 9}, {Int3( 0, -1,  0), 11},
            {Int3( 0,  0,  1), 6}, {Int3( 0,  0, -1), 6}
        },

        // Rank 8
        {
            {Int3( 1,  0,  0), 14}, {Int3(-1,  0,  0), 2},
            {Int3( 0,  1,  0), 10}, {Int3( 0, -1,  0), 6},
            {Int3( 0,  0,  1), 9}, {Int3( 0,  0, -1), 9}
        },

        // Rank 9
        {
            {Int3( 1,  0,  0), 15}, {Int3(-1,  0,  0), 3},
            {Int3( 0,  1,  0), 11}, {Int3( 0, -1,  0), 7},
            {Int3( 0,  0,  1), 8}, {Int3( 0,  0, -1), 8}
        },

        // Rank 10
        {
            {Int3( 1,  0,  0), 16}, {Int3(-1,  0,  0), 4},
            {Int3( 0,  1,  0), 6}, {Int3( 0, -1,  0), 8},
            {Int3( 0,  0,  1), 11}, {Int3( 0,  0, -1), 11}
        },

        // Rank 11
        {
            {Int3( 1,  0,  0), 17}, {Int3(-1,  0,  0), 5},
            {Int3( 0,  1,  0), 7}, {Int3( 0, -1,  0), 9},
            {Int3( 0,  0,  1), 10}, {Int3( 0,  0, -1), 10}
        },


        // Rank 12
        {
            {Int3( 1,  0,  0), 18}, {Int3(-1,  0,  0), 6},
            {Int3( 0,  1,  0), 14}, {Int3( 0, -1,  0), 16},
            {Int3( 0,  0,  1), 13}, {Int3( 0,  0, -1), 13}
        },

        // Rank 13
        {
            {Int3( 1,  0,  0), 19}, {Int3(-1,  0,  0), 7},
            {Int3( 0,  1,  0), 15}, {Int3( 0, -1,  0), 17},
            {Int3( 0,  0,  1), 12}, {Int3( 0,  0, -1), 12}
        },

        // Rank 14
        {
            {Int3( 1,  0,  0), 20}, {Int3(-1,  0,  0), 8},
            {Int3( 0,  1,  0), 16}, {Int3( 0, -1,  0), 12},
            {Int3( 0,  0,  1), 15}, {Int3( 0,  0, -1), 15}
        },

        // Rank 15
        {
            {Int3( 1,  0,  0), 21}, {Int3(-1,  0,  0), 9},
            {Int3( 0,  1,  0), 17}, {Int3( 0, -1,  0), 13},
            {Int3( 0,  0,  1), 14}, {Int3( 0,  0, -1), 14}
        },

        // Rank 16
        {
            {Int3( 1,  0,  0), 22}, {Int3(-1,  0,  0), 10},
            {Int3( 0,  1,  0), 12}, {Int3( 0, -1,  0), 14},
            {Int3( 0,  0,  1), 17}, {Int3( 0,  0, -1), 17}
        },

        // Rank 17
        {
            {Int3( 1,  0,  0), 23}, {Int3(-1,  0,  0), 11},
            {Int3( 0,  1,  0), 13}, {Int3( 0, -1,  0), 15},
            {Int3( 0,  0,  1), 16}, {Int3( 0,  0, -1), 16}
        },

        // Rank 18
        {
            {Int3( 1,  0,  0), 0}, {Int3(-1,  0,  0), 12},
            {Int3( 0,  1,  0), 20}, {Int3( 0, -1,  0), 22},
            {Int3( 0,  0,  1), 19}, {Int3( 0,  0, -1), 19}
        },

        // Rank 19
        {
            {Int3( 1,  0,  0), 1}, {Int3(-1,  0,  0), 13},
            {Int3( 0,  1,  0), 21}, {Int3( 0, -1,  0), 23},
            {Int3( 0,  0,  1), 18}, {Int3( 0,  0, -1), 18}
        },

        // Rank 20
        {
            {Int3( 1,  0,  0), 2}, {Int3(-1,  0,  0), 14},
            {Int3( 0,  1,  0), 22}, {Int3( 0, -1,  0), 18},
            {Int3( 0,  0,  1), 21}, {Int3( 0,  0, -1), 21}
        },

        // Rank 21
        {
            {Int3( 1,  0,  0), 3}, {Int3(-1,  0,  0), 15},
            {Int3( 0,  1,  0), 23}, {Int3( 0, -1,  0), 19},
            {Int3( 0,  0,  1), 20}, {Int3( 0,  0, -1), 20}
        },

        // Rank 22
        {
            {Int3( 1,  0,  0), 4}, {Int3(-1,  0,  0), 16},
            {Int3( 0,  1,  0), 18}, {Int3( 0, -1,  0), 20},
            {Int3( 0,  0,  1), 23}, {Int3( 0,  0, -1), 23}
        },

        // Rank 23
        {
            {Int3( 1,  0,  0), 5}, {Int3(-1,  0,  0), 17},
            {Int3( 0,  1,  0), 19}, {Int3( 0, -1,  0), 21},
            {Int3( 0,  0,  1), 22}, {Int3( 0,  0, -1), 22}
        },
    };

    for (int r = 0; r < size; r++)
        for (auto& offset_neighbor : neighbor_matrix[r])
        {
            Int3 offset = offset_neighbor.first;
            int neighbor = offset_neighbor.second;

            EXPECT_EQ(topology.getNeighbor(r, offset), neighbor);
        }
}