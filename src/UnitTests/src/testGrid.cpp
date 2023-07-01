#include "TestingUtility.h"

#include "Grid.h"

template <class gridType>
class GridTest : public BaseGridFixture<gridType> {
public:
    // excludes border cells for spectral solvers
    FP3 internalPointNotNearBorders() {
        return this->urandFP3(this->minCoords + this->grid->steps, this->maxCoords - this->grid->steps);
    }
};

typedef ::testing::Types<YeeGrid, SimpleGrid, PSTDGrid, PSATDGrid, PSATDTimeStaggeredGrid> types;
TYPED_TEST_CASE(GridTest, types);


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
        FP3 coords = this->internalPointNotNearBorders(), actualB, tmp;
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
        FP3 coords = this->internalPointNotNearBorders(), actualE;
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
        FP3 coords = this->internalPointNotNearBorders(), actualJ;
        actualJ = grid->getJ(coords);
        ASSERT_NEAR_FP3(coords, actualJ);
    }
}

TYPED_TEST(GridTest, InterpolateBPeriodical)
{
    // Checks whether interpolation is really periodical (for spectral grids)
    const FP value = 10000.0;
    auto grid = this->grid;
    for (int i = 0; i < grid->numCells.x; i++)
        for (int j = 0; j < grid->numCells.y; j++)
            for (int k = 0; k < grid->numCells.z; k++)
            {
                grid->Bx(i, j, k) = value;
                grid->By(i, j, k) = value;
                grid->Bz(i, j, k) = value;
            }
    // Check interpolation in random points
    for (int testIdx = 0; testIdx < 100; ++testIdx)
    {
        FP3 coords = this->internalPoint(), actualB, expectedB;
        actualB = grid->getB(coords);
        expectedB = { value, value, value };
        ASSERT_NEAR_FP3(expectedB, actualB);
    }
}

TYPED_TEST(GridTest, InterpolateEPeriodical)
{
    // Checks whether interpolation is really periodical (for spectral grids)
    const FP value = 10000.0;
    auto grid = this->grid;
    for (int i = 0; i < grid->numCells.x; i++)
        for (int j = 0; j < grid->numCells.y; j++)
            for (int k = 0; k < grid->numCells.z; k++)
            {
                grid->Ex(i, j, k) = value;
                grid->Ey(i, j, k) = value;
                grid->Ez(i, j, k) = value;
            }
    // Check interpolation in random points
    for (int testIdx = 0; testIdx < 100; ++testIdx)
    {
        FP3 coords = this->internalPoint(), actualE, expectedE;
        actualE = grid->getE(coords);
        expectedE = { value, value, value };
        ASSERT_NEAR_FP3(expectedE, actualE);
    }
}

TYPED_TEST(GridTest, InterpolateJPeriodical)
{
    // Checks whether interpolation is really periodical (for spectral grids)
    const FP value = 10000.0;
    auto grid = this->grid;
    for (int i = 0; i < grid->numCells.x; i++)
        for (int j = 0; j < grid->numCells.y; j++)
            for (int k = 0; k < grid->numCells.z; k++)
            {
                grid->Jx(i, j, k) = value;
                grid->Jy(i, j, k) = value;
                grid->Jz(i, j, k) = value;
            }
    // Check interpolation in random points
    for (int testIdx = 0; testIdx < 100; ++testIdx)
    {
        FP3 coords = this->internalPoint(), actualJ, expectedJ;
        actualJ = grid->getJ(coords);
        expectedJ = { value, value, value };
        ASSERT_NEAR_FP3(expectedJ, actualJ);
    }
}

TYPED_TEST(GridTest, GridInitialization)
{
    auto grid = this->grid;
    
    auto fieldFunc = [](FP3 coords) {
        return coords.x + coords.y + coords.z;
    };

    Int3 begin = Int3(0, 0, 0);
    Int3 end = this->grid->numCells;

    for (int i = begin.x; i < end.x; i++)
        for (int j = begin.y; j < end.y; j++)
            for (int k = begin.z; k < end.z; k++) {
                FP3 coords = grid->ExPosition(i, j, k);
                grid->Ex(i, j, k) = fieldFunc(coords);
                coords = grid->EyPosition(i, j, k);
                grid->Ey(i, j, k) = fieldFunc(coords);
                coords = grid->EzPosition(i, j, k);
                grid->Ez(i, j, k) = fieldFunc(coords);
                coords = grid->BxPosition(i, j, k);
                grid->Bx(i, j, k) = fieldFunc(coords);
                coords = grid->ByPosition(i, j, k);
                grid->By(i, j, k) = fieldFunc(coords);
                coords = grid->BzPosition(i, j, k);
                grid->Bz(i, j, k) = fieldFunc(coords);
            }

    for (int i = begin.x; i < end.x; ++i)
        for (int j = begin.y; j < end.y; ++j)
            for (int k = begin.z; k < end.z; ++k)
            {
                FP3 expectedE, actualE;
                FP3 coords = this->grid->ExPosition(i, j, k);
                expectedE.x = fieldFunc(coords);
                coords = this->grid->EyPosition(i, j, k);
                expectedE.y = fieldFunc(coords);
                coords = this->grid->EzPosition(i, j, k);
                expectedE.z = fieldFunc(coords);
                actualE.x = this->grid->Ex(i, j, k);
                actualE.y = this->grid->Ey(i, j, k);
                actualE.z = this->grid->Ez(i, j, k);
                ASSERT_NEAR_FP3(expectedE, actualE);
            }

    for (int i = begin.x; i < end.x; ++i)
        for (int j = begin.y; j < end.y; ++j)
            for (int k = begin.z; k < end.z; ++k)
            {
                FP3 expectedB, actualB;
                FP3 coords = this->grid->BxPosition(i, j, k);
                expectedB.x = fieldFunc(coords);
                coords = this->grid->ByPosition(i, j, k);
                expectedB.y = fieldFunc(coords);
                coords = this->grid->BzPosition(i, j, k);
                expectedB.z = fieldFunc(coords);
                actualB.x = this->grid->Bx(i, j, k);
                actualB.y = this->grid->By(i, j, k);
                actualB.z = this->grid->Bz(i, j, k);
                ASSERT_NEAR_FP3(expectedB, actualB);
            }
}
