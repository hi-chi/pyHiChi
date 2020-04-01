#include "TestingUtility.h"

#include "Psatd.h"

class GridPSATDTimeStraggeredTest : public BaseGridFixture<PSATDGrid> {
public:
    PSATDTimeStraggered * psatd;

    virtual void SetUp() {
        BaseGridFixture<PSATDGrid>::SetUp();
        initializeGrid();
        psatd = new PSATDTimeStraggered(grid);
    }

    void initializeGrid() {
        for (int i = 0; i < grid->numCells.x; i++)
            for (int j = 0; j < grid->numCells.y; j++)
                for (int k = 0; k < grid->numCells.z; k++) {
                    FP3 coords = grid->ExPosition(i, j, k);
                    grid->Ex(i, j, k) = funcE(coords.x, coords.y, coords.z, 0).x;
                    coords = grid->EyPosition(i, j, k);
                    grid->Ey(i, j, k) = funcE(coords.x, coords.y, coords.z, 0).y;
                    coords = grid->EzPosition(i, j, k);
                    grid->Ez(i, j, k) = funcE(coords.x, coords.y, coords.z, 0).z;
                    coords = grid->BxPosition(i, j, k);
                    grid->Bx(i, j, k) = funcB(coords.x, coords.y, coords.z, 0).x;
                    coords = grid->ByPosition(i, j, k);
                    grid->By(i, j, k) = funcB(coords.x, coords.y, coords.z, 0).y;
                    coords = grid->BzPosition(i, j, k);
                    grid->Bz(i, j, k) = funcB(coords.x, coords.y, coords.z, 0).z;
                }
    }

    FP3 funcE(FP x, FP y, FP z, FP t) {
        return FP3(sin(2 * constants::pi / (grid->steps.z * grid->numCells.z)
            *(z - constants::c*t)), 0.0, 0.0);
    }

    FP3 funcB(FP x, FP y, FP z, FP t) {
        return FP3(0.0, sin(2 * constants::pi / (grid->steps.z * grid->numCells.z)
            *(z - constants::c*t)), 0.0);
    }
};

TEST_F(GridPSATDTimeStraggeredTest, ADD_TEST_FFT_PREFIX(PSATDTimeStraggered)) {

    const int numSteps = 100;
    for (int step = 0; step < numSteps; ++step)
        psatd->updateFields();

    FP finalT = grid->dt * numSteps;
    for (int i = 0; i < grid->numCells.x; ++i)
        for (int j = 0; j < grid->numCells.y; ++j)
            for (int k = 0; k < grid->numCells.z; ++k)
            {
                FP3 expectedE, actualE;
                FP3 coords = grid->ExPosition(i, j, k);
                expectedE.x = funcE(coords.x, coords.y, coords.z, finalT).x;
                coords = grid->EyPosition(i, j, k);
                expectedE.y = funcE(coords.x, coords.y, coords.z, finalT).y;
                coords = grid->EzPosition(i, j, k);
                expectedE.z = funcE(coords.x, coords.y, coords.z, finalT).z;
                actualE.x = grid->Ex(i, j, k);
                actualE.y = grid->Ey(i, j, k);
                actualE.z = grid->Ez(i, j, k);
                ASSERT_NEAR_FP3(expectedE, actualE);
            }

    for (int i = 0; i < grid->numCells.x; ++i)
        for (int j = 0; j < grid->numCells.y; ++j)
            for (int k = 0; k < grid->numCells.z; ++k)
            {
                FP3 expectedB, actualB;
                FP3 coords = grid->BxPosition(i, j, k);
                expectedB.x = funcB(coords.x, coords.y, coords.z, finalT).x;
                coords = grid->ByPosition(i, j, k);
                expectedB.y = funcB(coords.x, coords.y, coords.z, finalT).y;
                coords = grid->BzPosition(i, j, k);
                expectedB.z = funcB(coords.x, coords.y, coords.z, finalT).z;
                actualB.x = grid->Bx(i, j, k);
                actualB.y = grid->By(i, j, k);
                actualB.z = grid->Bz(i, j, k);
                ASSERT_NEAR_FP3(expectedB, actualB);
            }
}