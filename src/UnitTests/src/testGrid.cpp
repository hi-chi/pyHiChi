#include "TestingUtility.h"

template <class gridType>
class GridTest : public BaseGridFixture<gridType> {
};

typedef ::testing::Types<YeeGrid, SimpleGrid> types;
TYPED_TEST_CASE(GridTest, types);

/*TYPED_TEST(GridTest, LoadDumpE)
{
    int numCells = grid->numCells.volume();
    FP3 minValue(0.135, -124.234, -3.1234);
    FP3 maxValue(1835.623023, -54.2135, 2.641);
    std::vector<FP3> expectedE = randomVectors(numCells, minValue, maxValue);
    Int3 minIndex(0, 0, 0);
    Int3 maxIndex = grid->numCells;
    grid->loadE(expectedE.data(), &minIndex, &maxIndex);
    std::vector<FP3> actualE(numCells);
    grid->dumpE(actualE.data(), &minIndex, &maxIndex);
    for (int idx = 0; idx < numCells; ++idx)
    {
        ASSERT_EQ_FP3(actualE[idx], expectedE[idx]);
    }
}


TYPED_TEST(GridTest, LoadDumpB)
{
    int numCells = grid->numCells.volume();
    FP3 minValue(-12.324, 84.3872, 23.12378);
    FP3 maxValue(-8.9724, 223.97234, 125.453);
    std::vector<FP3> expectedB = randomVectors(numCells, minValue, maxValue);
    Int3 minIndex(0, 0, 0);
    Int3 maxIndex = grid->numCells;
    grid->loadB(expectedB.data(), &minIndex, &maxIndex);
    std::vector<FP3> actualB(numCells);
    grid->dumpB(actualB.data(), &minIndex, &maxIndex);
    for (int idx = 0; idx < numCells; ++idx)
    {
        ASSERT_EQ_FP3(actualB[idx], expectedB[idx]);
    }
}*/


TYPED_TEST(GridTest, InterpolateB)
{
    // Set Bx(x,y,z) = x, By(x,y,z) = y, Bz(x,y,z) = z
    auto grid = this->grid;
    for (int i = 0; i < grid->numCells.x; i++)
        for (int j = 0; j < grid->numCells.y; j++)
            for (int k = 0; k < grid->numCells.z; k++)
            {
                grid->Bx(i, j, k) = grid->BxPosition(i, j, k).x;
                grid->By(i, j, k) = grid->ByPosition(i, j, k).y;
                grid->Bz(i, j, k) = grid->BzPosition(i, j, k).z;
            }
    // Check interpolation in random points gives B = coords of point.
    for (int testIdx = 0; testIdx < 100; ++testIdx)
    {
        FP3 coords = this->internalPoint(), actualB, tmp;
        actualB = grid->getB(coords);
        ASSERT_NEAR_FP3(coords, actualB);
    }
}


TYPED_TEST(GridTest, InterpolateE)
{
    // Set Ex(x,y,z) = x, Ey(x,y,z) = y, Ez(x,y,z) = z
    auto grid = this->grid;
    for (int i = 0; i < grid->numCells.x; i++)
        for (int j = 0; j < grid->numCells.y; j++)
            for (int k = 0; k < grid->numCells.z; k++)
            {
                grid->Ex(i, j, k) = grid->ExPosition(i, j, k).x;
                grid->Ey(i, j, k) = grid->EyPosition(i, j, k).y;
                grid->Ez(i, j, k) = grid->EzPosition(i, j, k).z;
            }
    // Check interpolation in random points gives E = coords of point.
    for (int testIdx = 0; testIdx < 100; ++testIdx)
    {
        FP3 coords = this->internalPoint(), actualE;
        actualE = grid->getE(coords);
        ASSERT_NEAR_FP3(coords, actualE);
    }
}


TYPED_TEST(GridTest, InterpolateJ)
{
    // Set Jx(x,y,z) = x, Jy(x,y,z) = y, Jz(x,y,z) = z
    auto grid = this->grid;
    for (int i = 0; i < grid->numCells.x; i++)
        for (int j = 0; j < grid->numCells.y; j++)
            for (int k = 0; k < grid->numCells.z; k++)
            {
                grid->Jx(i, j, k) = grid->JxPosition(i, j, k).x;
                grid->Jy(i, j, k) = grid->JyPosition(i, j, k).y;
                grid->Jz(i, j, k) = grid->JzPosition(i, j, k).z;
            }
    // Check interpolation in random points gives J = coords of point.
    for (int testIdx = 0; testIdx < 100; ++testIdx)
    {
        FP3 coords = this->internalPoint(), actualJ;
        actualJ = grid->getJ(coords);
        ASSERT_NEAR_FP3(coords, actualJ);
    }
}

