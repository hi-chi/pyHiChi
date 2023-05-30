#include "TestingUtility.h"

#include "Fdtd.h"
#include "Pstd.h"
#include "Psatd.h"
#include "PsatdTimeStraggered.h"

template <class TTypeDefinitionsFieldTest>
class PMLTest : public BaseFixture {
public:

    using FieldSolverType = typename TTypeDefinitionsFieldTest::FieldSolverType;
    using GridType = typename TTypeDefinitionsFieldTest::GridType;
    using PeriodicalBoundaryConditionType = typename FieldSolverType::PeriodicalBoundaryConditionType;

    const int dimension = TTypeDefinitionsFieldTest::dimension;
    const CoordinateEnum axis = TTypeDefinitionsFieldTest::axis;

    const int gridSizeLongitudinal = 32;
    const int gridSizeTransverse = 8;
    const int pmlSize1d = 4;

    FP3 pmlLeftEnd;
    FP3 pmlRightStart;

    std::unique_ptr<FieldSolverType> fieldSolver;
    std::unique_ptr<GridType> grid;
    std::unique_ptr<PeriodicalBoundaryConditionType> boundaryCondition;

    Int3 gridSize;
    Int3 pmlSize;
    FP3 gridStep;
    FP timeStep = 0;
    FP3 minCoords, maxCoords;

    FP relatedEnergyThreshold = 1e-2;

    virtual void SetUp() {
        gridSize = Int3(1, 1, 1);
        for (int d = 0; d < dimension; d++) {
            gridSize[d] = gridSizeTransverse;
        }
        gridSize[(int)axis] = gridSizeLongitudinal;

        pmlSize = Int3(0, 0, 0);
        pmlSize[(int)axis] = pmlSize1d;

        this->minCoords = FP3(0, 0, 0);
        FP b = gridSizeLongitudinal * constants::c;
        this->maxCoords = FP3(b, b, b);
        this->gridStep = (this->maxCoords - this->minCoords) / (FP3)this->gridSize;
        this->pmlLeftEnd = this->minCoords + pmlSize * this->gridStep;
        this->pmlRightStart = this->maxCoords - pmlSize * this->gridStep;

        this->grid.reset(new GridType(this->gridSize, this->minCoords, this->gridStep, this->gridSize));

        // should satisfy the Courant's condition for all solvers
        this->timeStep = 0.4 * constants::c / grid->steps.norm();

        fieldSolver.reset(new FieldSolverType(this->grid.get(), this->timeStep));
        fieldSolver->setPML(pmlSize.x, pmlSize.y, pmlSize.z);

        initTest();
        initializeGrid();
    }

    virtual void initTest() = 0;

    int getIterationNumber() {
        return (int)(grid->numCells[(int)axis] * grid->steps[(int)axis] /
            (constants::c * fieldSolver->dt));
    }

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
        b[(int)axisB] = fieldFunc(x, y, z);
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

typedef ::testing::Types <
    TypeDefinitionsFieldTest<FDTD, YeeGrid, 1, CoordinateEnum::x>,
    TypeDefinitionsFieldTest<FDTD, YeeGrid, 2, CoordinateEnum::x>,
    TypeDefinitionsFieldTest<FDTD, YeeGrid, 2, CoordinateEnum::y>,
    TypeDefinitionsFieldTest<FDTD, YeeGrid, 3, CoordinateEnum::x>,
    TypeDefinitionsFieldTest<FDTD, YeeGrid, 3, CoordinateEnum::y>,
    TypeDefinitionsFieldTest<FDTD, YeeGrid, 3, CoordinateEnum::z>,
    
    TypeDefinitionsFieldTest<PSTD, PSTDGrid, 1, CoordinateEnum::x>,
    TypeDefinitionsFieldTest<PSTD, PSTDGrid, 2, CoordinateEnum::x>,
    TypeDefinitionsFieldTest<PSTD, PSTDGrid, 2, CoordinateEnum::y>,
    TypeDefinitionsFieldTest<PSTD, PSTDGrid, 3, CoordinateEnum::x>,
    TypeDefinitionsFieldTest<PSTD, PSTDGrid, 3, CoordinateEnum::y>,
    TypeDefinitionsFieldTest<PSTD, PSTDGrid, 3, CoordinateEnum::z>,
    
    TypeDefinitionsFieldTest<PSATD, PSATDGrid, 1, CoordinateEnum::x>,
    TypeDefinitionsFieldTest<PSATD, PSATDGrid, 2, CoordinateEnum::x>,
    TypeDefinitionsFieldTest<PSATD, PSATDGrid, 2, CoordinateEnum::y>,
    TypeDefinitionsFieldTest<PSATD, PSATDGrid, 3, CoordinateEnum::x>,
    TypeDefinitionsFieldTest<PSATD, PSATDGrid, 3, CoordinateEnum::y>,
    TypeDefinitionsFieldTest<PSATD, PSATDGrid, 3, CoordinateEnum::z>,
    
    TypeDefinitionsFieldTest<PSATDTimeStraggered, PSATDTimeStraggeredGrid, 1, CoordinateEnum::x>,
    TypeDefinitionsFieldTest<PSATDTimeStraggered, PSATDTimeStraggeredGrid, 2, CoordinateEnum::x>,
    TypeDefinitionsFieldTest<PSATDTimeStraggered, PSATDTimeStraggeredGrid, 2, CoordinateEnum::y>,
    TypeDefinitionsFieldTest<PSATDTimeStraggered, PSATDTimeStraggeredGrid, 3, CoordinateEnum::x>,
    TypeDefinitionsFieldTest<PSATDTimeStraggered, PSATDTimeStraggeredGrid, 3, CoordinateEnum::y>,
    TypeDefinitionsFieldTest<PSATDTimeStraggered, PSATDTimeStraggeredGrid, 3, CoordinateEnum::z>
> types;


template <class TTypeDefinitionsPMLTest>
class PMLTestOnly : public PMLTest<TTypeDefinitionsPMLTest> {
public:

    FP fieldFunc(FP x, FP y, FP z) override {
        int axis0 = (int)this->axis;
        int axis1 = (axis0 + 1) % 3;
        int axis2 = (axis0 + 2) % 3;
        FP3 coord(x, y, z);
        FP3 a = this->pmlLeftEnd, b = this->pmlRightStart;
        FP omega = 4.0;
        FP wave = sin((FP)2.0 * omega * constants::pi / (b[axis0] - a[axis0]) * (coord[axis0] - a[axis0]));
        FP env1 = cos(constants::pi / (b[axis1] - a[axis1]) * (coord[axis1] - a[axis1] - (a[axis1] + b[axis1]) * 0.5));
        FP env2 = cos(constants::pi / (b[axis2] - a[axis2]) * (coord[axis2] - a[axis2] - (a[axis2] + b[axis2]) * 0.5));
        return wave * env1 * env1 * env2 * env2;
    }

    void initTest() override {
        this->relatedEnergyThreshold = 1e-1;
    }

};

TYPED_TEST_SUITE(PMLTestOnly, types);

TYPED_TEST(PMLTestOnly, PmlTest) {
    // to disable testing of spectral solvers without enabled fftw
#ifndef __USE_FFT__
    SUCCEED();
#else

    const int numSteps = (int)(grid->numCells[(int)axis] * grid->steps[(int)axis] /
        (constants::c * fieldSolver->dt));

    FP startEnergy = computeEnergy();

    for (int step = 0; step < numSteps; ++step) {
        fieldSolver->updateFields();
    }

    FP finalEnergy = computeEnergy();

    ASSERT_NEAR(finalEnergy / startEnergy, 0, relatedEnergyThreshold);

#endif
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
        bool periodicAxis[3] = { true, true, true };
        periodicAxis[(int)axis] = false;

        boundaryCondition.reset(new PeriodicalBoundaryConditionType(
            fieldSolver.get(), periodicAxis[0], periodicAxis[1], periodicAxis[2]));

        fieldSolver->setBoundaryCondition(boundaryCondition.get());
    }

};

TYPED_TEST_SUITE(PMLTestPeriodical, types);

TYPED_TEST(PMLTestPeriodical, PmlTestPeriodical) {
    // to disable testing of spectral solvers without enabled fftw
#ifndef __USE_FFT__
    SUCCEED();
#else

    const int numSteps = (int)(grid->numCells[(int)axis] * grid->steps[(int)axis] /
        (constants::c * fieldSolver->dt));

    FP startEnergy = computeEnergy();

    for (int step = 0; step < numSteps; ++step) {
        fieldSolver->updateFields();
    }

    FP finalEnergy = computeEnergy();

    ASSERT_NEAR(finalEnergy / startEnergy, 0, relatedEnergyThreshold);

#endif
}
