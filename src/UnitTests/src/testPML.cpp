#include "TestingUtility.h"

#include <type_traits>

#include "Fdtd.h"
#include "Pstd.h"
#include "Psatd.h"
#include "PsatdTimeStaggered.h"

template <class TTypeDefinitionsFieldTest>
class PMLTest : public BaseFixture {
public:

    using FieldSolverType = typename TTypeDefinitionsFieldTest::FieldSolverType;
    using GridType = typename TTypeDefinitionsFieldTest::FieldSolverType::GridType;

    const int dimension = TTypeDefinitionsFieldTest::dimension;
    const CoordinateEnum axis = TTypeDefinitionsFieldTest::axis;

    const int gridSizeLongitudinal = 32;
    const int gridSizeTransverse = 8;
    const int pmlSize1d = 4;

    FP3 pmlLeftEnd;
    FP3 pmlRightStart;

    std::unique_ptr<FieldSolverType> fieldSolver;
    std::unique_ptr<GridType> grid;

    Int3 gridSize;
    Int3 pmlSize;
    FP3 gridStep;
    FP timeStep = 0;
    FP3 minCoords, maxCoords;
    int numSteps = 0;

    FP relatedEnergyThreshold = 0.05;

    virtual void SetUp() {
        gridSize = Int3(1, 1, 1);
        for (int d = 0; d < dimension; d++) {
            gridSize[d] = gridSizeTransverse;
        }
        gridSize[(int)axis] = gridSizeLongitudinal;

        pmlSize = Int3(0, 0, 0);
        pmlSize[(int)axis] = pmlSize1d;

        gridSize += 2 * pmlSize;
        int maxGridSize = gridSizeLongitudinal + 2 * pmlSize1d;

        this->minCoords = FP3(0, 0, 0);
        this->maxCoords = FP3(maxGridSize, maxGridSize, maxGridSize) * constants::c;
        this->gridStep = (this->maxCoords - this->minCoords) / (FP3)this->gridSize;
        this->pmlLeftEnd = this->minCoords + pmlSize * this->gridStep;
        this->pmlRightStart = this->maxCoords - pmlSize * this->gridStep;

        this->grid.reset(new GridType(this->gridSize, this->minCoords, this->gridStep, this->gridSize));

        this->timeStep = 0.5 * FieldSolverType::getCourantConditionTimeStep(this->gridStep);
        this->numSteps = (int)(1.3 * (this->pmlRightStart - this->pmlLeftEnd)[(int)axis] /
            (constants::c * this->timeStep));

        this->fieldSolver.reset(new FieldSolverType(this->grid.get(), this->timeStep));
        this->fieldSolver->setPML(pmlSize.x, pmlSize.y, pmlSize.z);

        initTest();
        initializeGrid();
    }

    virtual void initTest() {}

    void initializeGrid() {
        Int3 begin = pmlSize + grid->getNumExternalLeftCells();
        Int3 end = grid->numCells - pmlSize - grid->getNumExternalRightCells();

        for (int i = begin.x; i < end.x; i++)
            for (int j = begin.y; j < end.y; j++)
                for (int k = begin.z; k < end.z; k++) {
                    FP3 coords = grid->ExPosition(i, j, k);
                    grid->Ex(i, j, k) = funcE(coords.x, coords.y, coords.z).x;
                    coords = grid->EyPosition(i, j, k);
                    grid->Ey(i, j, k) = funcE(coords.x, coords.y, coords.z).y;
                    coords = grid->EzPosition(i, j, k);
                    grid->Ez(i, j, k) = funcE(coords.x, coords.y, coords.z).z;
                    coords = grid->BxPosition(i, j, k);
                    grid->Bx(i, j, k) = funcB(coords.x, coords.y, coords.z).x;
                    coords = grid->ByPosition(i, j, k);
                    grid->By(i, j, k) = funcB(coords.x, coords.y, coords.z).y;
                    coords = grid->BzPosition(i, j, k);
                    grid->Bz(i, j, k) = funcB(coords.x, coords.y, coords.z).z;
                }
    }

    virtual FP fieldFunc(FP x, FP y, FP z) = 0;

    FP3 funcE(FP x, FP y, FP z) {
        CoordinateEnum axisE = CoordinateEnum(((int)axis + 1) % 3);
        FP3 e;
        e[(int)axisE] = fieldFunc(x, y, z);
        return e;
    }

    FP3 funcB(FP x, FP y, FP z) {
        CoordinateEnum axisB = CoordinateEnum(((int)axis + 2) % 3);
        FP3 b;
        b[(int)axisB] = -fieldFunc(x, y, z);
        return b;
    }

    FP computeEnergy() {
        FP energy = 0;
        for (int i = 0; i < grid->numCells.x; i++)
            for (int j = 0; j < grid->numCells.y; j++)
                for (int k = 0; k < grid->numCells.z; k++)
                    energy += pow(grid->Ex(i, j, k), 2) + pow(grid->Ey(i, j, k), 2) + pow(grid->Ez(i, j, k), 2) +
                    pow(grid->Bx(i, j, k), 2) + pow(grid->By(i, j, k), 2) + pow(grid->Bz(i, j, k), 2);
        return energy;
    }

};

#ifndef __USE_FFT__

typedef ::testing::Types <
    TypeDefinitionsFieldTest<FDTD, 1, CoordinateEnum::x>,
    TypeDefinitionsFieldTest<FDTD, 2, CoordinateEnum::x>,
    TypeDefinitionsFieldTest<FDTD, 2, CoordinateEnum::y>,
    TypeDefinitionsFieldTest<FDTD, 3, CoordinateEnum::x>,
    TypeDefinitionsFieldTest<FDTD, 3, CoordinateEnum::y>,
    TypeDefinitionsFieldTest<FDTD, 3, CoordinateEnum::z>
> types;

#else

typedef ::testing::Types <
    TypeDefinitionsFieldTest<FDTD, 1, CoordinateEnum::x>,
    TypeDefinitionsFieldTest<FDTD, 2, CoordinateEnum::x>,
    TypeDefinitionsFieldTest<FDTD, 2, CoordinateEnum::y>,
    TypeDefinitionsFieldTest<FDTD, 3, CoordinateEnum::x>,
    TypeDefinitionsFieldTest<FDTD, 3, CoordinateEnum::y>,
    TypeDefinitionsFieldTest<FDTD, 3, CoordinateEnum::z>,

    TypeDefinitionsFieldTest<PSTD, 1, CoordinateEnum::x>,
    TypeDefinitionsFieldTest<PSTD, 2, CoordinateEnum::x>,
    TypeDefinitionsFieldTest<PSTD, 2, CoordinateEnum::y>,
    TypeDefinitionsFieldTest<PSTD, 3, CoordinateEnum::x>,
    TypeDefinitionsFieldTest<PSTD, 3, CoordinateEnum::y>,
    TypeDefinitionsFieldTest<PSTD, 3, CoordinateEnum::z>,

    TypeDefinitionsFieldTest<PSATD, 1, CoordinateEnum::x>,
    TypeDefinitionsFieldTest<PSATD, 2, CoordinateEnum::x>,
    TypeDefinitionsFieldTest<PSATD, 2, CoordinateEnum::y>,
    TypeDefinitionsFieldTest<PSATD, 3, CoordinateEnum::x>,
    TypeDefinitionsFieldTest<PSATD, 3, CoordinateEnum::y>,
    TypeDefinitionsFieldTest<PSATD, 3, CoordinateEnum::z>,

    TypeDefinitionsFieldTest<PSATDTimeStaggered, 1, CoordinateEnum::x>,
    TypeDefinitionsFieldTest<PSATDTimeStaggered, 2, CoordinateEnum::x>,
    TypeDefinitionsFieldTest<PSATDTimeStaggered, 2, CoordinateEnum::y>,
    TypeDefinitionsFieldTest<PSATDTimeStaggered, 3, CoordinateEnum::x>,
    TypeDefinitionsFieldTest<PSATDTimeStaggered, 3, CoordinateEnum::y>,
    TypeDefinitionsFieldTest<PSATDTimeStaggered, 3, CoordinateEnum::z>
> types;

#endif


template <class TTypeDefinitionsPMLTest>
class PMLTestOnly : public PMLTest<TTypeDefinitionsPMLTest> {
public:

    FP harrisFunction(FP x, FP a, FP b) {
        FP coord = x / (b - a);
        if (coord < 0.0 || coord >= 1.0) return 0.0;
        const FP pi2 = 2.0 * constants::pi;
        return 0.03125 * (10.0 - 15.0 * cos(pi2 * coord) + 6.0 * cos(2.0 * pi2 * coord) - cos(3.0 * pi2 * coord));
    }

    FP waveFunction(FP x, FP a, FP b) {
        const FP omega = 4.0;
        FP coord = x / (b - a);
        return sin(2.0 * constants::pi * omega * coord);
    }

    FP fieldFunc(FP x, FP y, FP z) override {
        int axis0 = (int)this->axis;
        int axis1 = (axis0 + 1) % 3;
        int axis2 = (axis0 + 2) % 3;

        FP3 coord(x, y, z);
        FP3 a = this->pmlLeftEnd, b = this->pmlRightStart;
        FP omega = 4.0;

        FP res = waveFunction(coord[axis0], a[axis0], b[axis0]);
        for (int d = 0; d < this->grid->dimensionality; d++)
            res *= harrisFunction(coord[d], a[d], b[d]);

        return res;
    }

};

TYPED_TEST_SUITE(PMLTestOnly, types);

TYPED_TEST(PMLTestOnly, PmlTest) {

    FP startEnergy = this->computeEnergy();

    for (int step = 0; step < this->numSteps; ++step) {
        this->fieldSolver->updateFields();
    }

    FP finalEnergy = this->computeEnergy();

    ASSERT_NEAR(finalEnergy / startEnergy, 0, this->relatedEnergyThreshold);
}

TYPED_TEST(PMLTestOnly, CanResetPmlWhenChangeTimeStep) {
    FP oldDt = this->fieldSolver->pml->dt;

    this->fieldSolver->setTimeStep(this->fieldSolver->dt * 0.5);

    ASSERT_NE(oldDt, this->fieldSolver->pml->dt);
}


template <class TTypeDefinitionsPMLTest>
class PMLTestPeriodical : public PMLTest<TTypeDefinitionsPMLTest> {
public:

    FP fieldFunc(FP x, FP y, FP z) override {
        FP3 coord(x, y, z);
        int axis0 = (int)this->axis;
        FP3 a = this->pmlLeftEnd, b = this->pmlRightStart;
        return sin((FP)2.0 * constants::pi / (b[axis0] - a[axis0]) * (coord[axis0] - a[axis0]));
    }

    void initTest() override {
        for (int d = 0; d < this->grid->dimensionality; d++)
            if (d != (int)this->axis)
                this->fieldSolver->setPeriodicalBoundaryConditions((CoordinateEnum)d);
    }

};

TYPED_TEST_SUITE(PMLTestPeriodical, types);

TYPED_TEST(PMLTestPeriodical, PmlTestPeriodical) {

    FP startEnergy = this->computeEnergy();

    for (int step = 0; step < this->numSteps; ++step) {
        this->fieldSolver->updateFields();
    }

    FP finalEnergy = this->computeEnergy();

    ASSERT_NEAR(finalEnergy / startEnergy, 0, this->relatedEnergyThreshold);
}
