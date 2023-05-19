#include "TestingUtility.h"

#include "Pstd.h"
#include "Psatd.h"
#include "PsatdTimeStraggered.h"

template <class FieldSolverType, class GridType>
class PMLTest : public BaseGridFixture<GridType> {
public:
    FieldSolverType * fieldSolver;
    Int3 gridSize;
    Int3 pmlSize;
    FP3 pmlLeftEnd, pmlRightStart;
    FP timeStep = 0.1;

    virtual void SetUp() {
        gridSize = Int3(32, 1, 1);
        pmlSize = Int3(4, 0, 0);
        this->minCoords = FP3(0, 0, 0);
        this->maxCoords = FP3(gridSize.x * constants::c, gridSize.y * constants::c, gridSize.z * constants::c);
        createGrid();
        pmlLeftEnd.x = this->minCoords.x + pmlSize.x * this->grid->steps.x;
        pmlRightStart.x = this->maxCoords.x - pmlSize.x * this->grid->steps.x;
        initializeGrid();
        fieldSolver = new FieldSolverType(this->grid, this->timeStep);
        fieldSolver->setPML(pmlSize.x, pmlSize.y, pmlSize.z);
    }

    void createGrid() {
        FP3 steps((this->maxCoords.x - this->minCoords.x) / gridSize.x,
            (this->maxCoords.y - this->minCoords.y) / gridSize.y,
            (this->maxCoords.z - this->minCoords.z) / gridSize.z);
        this->grid = new GridType(gridSize, this->minCoords, steps, gridSize);
    }

    void initializeGrid() {
        auto grid = this->grid;
        for (int i = pmlSize.x + grid->getNumExternalLeftCells().x;
            i < grid->numCells.x - pmlSize.x - grid->getNumExternalRightCells().x; i++)
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


    void print() {
        auto grid = this->grid;
        for (int i = pmlSize.x + grid->getNumExternalLeftCells().x;
            i < grid->numCells.x - pmlSize.x - grid->getNumExternalRightCells().x; i++)
            std::cout << grid->Ey(i, 0, 0) << " ";
        std::cout<<std::endl;
    }

    FP3 funcE(FP x, FP y, FP z, FP t) {
        return FP3(0.0, sin(2 * constants::pi / (pmlRightStart.x - pmlLeftEnd.x)
            *(x - constants::c*t - pmlLeftEnd.x)), 0.0);
    }

    FP3 funcB(FP x, FP y, FP z, FP t) {
        return FP3(0.0, 0.0, sin(2 * constants::pi / (pmlRightStart.x - pmlLeftEnd.x)
            *(x - constants::c*t - pmlLeftEnd.x)));
    }

    FP computeEnergy() {
        FP energy = 0;
        auto grid = this->grid;
        for (int i = 0; i < gridSize.x; i++)
            for (int j = 0; j < gridSize.y; j++)
                for (int k = 0; k < gridSize.z; k++)
                    energy += pow(grid->Ex(i, j, k), 2) + pow(grid->Ey(i, j, k), 2) + pow(grid->Ez(i, j, k), 2) +
                    pow(grid->Bx(i, j, k), 2) + pow(grid->By(i, j, k), 2) + pow(grid->Bz(i, j, k), 2);
        return energy;
    }

};

typedef PMLTest<PSTD, PSTDGrid> PMLTestPSTD;
typedef PMLTest<PSATD, PSATDGrid> PMLTestPSATD;
typedef PMLTest<PSATDTimeStraggered, PSATDTimeStraggeredGrid> PMLTestPSATDTimeStraggered;

TEST_F(PMLTestPSTD, ADD_TEST_FFT_PREFIX(PmlPstd)) {
    const int numSteps = (int)((pmlRightStart.x - pmlLeftEnd.x) / constants::c / fieldSolver->dt);

    FP startEnergy = computeEnergy();

    for (int step = 0; step < numSteps; ++step)
        fieldSolver->updateFields();

    FP finalEnergy = computeEnergy();

    ASSERT_NEAR(finalEnergy / startEnergy, 0, 0.05);
}

TEST_F(PMLTestPSATD, ADD_TEST_FFT_PREFIX(PmlPsatd)) {
    const int numSteps = (int)((pmlRightStart.x - pmlLeftEnd.x) / constants::c / fieldSolver->dt);

    FP startEnergy = computeEnergy();

    for (int step = 0; step < numSteps; ++step)
        fieldSolver->updateFields();

    FP finalEnergy = computeEnergy();

    ASSERT_NEAR(finalEnergy / startEnergy, 0, 0.05);
}

TEST_F(PMLTestPSATDTimeStraggered, ADD_TEST_FFT_PREFIX(PmlPsatdTimeStraggered)) {
    const int numSteps = (int)((pmlRightStart.x - pmlLeftEnd.x) / constants::c / fieldSolver->dt);

    FP startEnergy = computeEnergy();

    for (int step = 0; step < numSteps; ++step)
        fieldSolver->updateFields();

    FP finalEnergy = computeEnergy();

    ASSERT_NEAR(finalEnergy / startEnergy, 0, 0.05);
}