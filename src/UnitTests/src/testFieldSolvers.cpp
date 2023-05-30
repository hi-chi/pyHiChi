#include "TestingUtility.h"

#include "Fdtd.h"
#include "Pstd.h"
#include "Psatd.h"
#include "PsatdTimeStraggered.h"


template <class TFieldSolver, class TGrid, int VDimension, CoordinateEnum VAxis>
struct TypeDefinitionsFieldSolverTest
{
    using FieldSolverType = TFieldSolver;
    using GridType = TGrid;
    static const int dimension = VDimension;
    static const CoordinateEnum axis = VAxis;
};

template <class TTypeDefinitionsFieldSolverTest>
class FieldSolverTest : public BaseFixture {
public:
    using FieldSolverType = typename TTypeDefinitionsFieldSolverTest::FieldSolverType;
    using GridType = typename TTypeDefinitionsFieldSolverTest::GridType;

    const int dimension = TTypeDefinitionsFieldSolverTest::dimension;
    const CoordinateEnum axis = TTypeDefinitionsFieldSolverTest::axis;

    const int gridSizeLongitudinal = 32;
    const int gridSizeTransverse = 8;
    const int pmlSize1d = 4;

    std::unique_ptr<FieldSolverType> fieldSolver;
    std::unique_ptr<GridType> grid;
    Int3 gridSize;
    FP3 gridStep;
    FP timeStep = 0;
    FP3 minCoords, maxCoords;

    std::unique_ptr<typename FieldSolverType::PeriodicalBoundaryConditionType> boundaryCondition;

    virtual void SetUp() {
        this->maxAbsoluteError = 1e-4;
        this->maxRelativeError = 0.5;

        gridSize = Int3(1, 1, 1);
        for (int d = 0; d < dimension; d++) {
            gridSize[d] = gridSizeTransverse;
        }
        gridSize[(int)axis] = gridSizeLongitudinal;

        this->minCoords = FP3(0, 0, 0);
        this->maxCoords = constants::c * (FP3)gridSize;
        this->gridStep = (this->maxCoords - this->minCoords) / (FP3)this->gridSize;

        this->grid.reset(new GridType(this->gridSize, this->minCoords, this->gridStep, this->gridSize));

        // should satisfy the Courant's condition for all solvers
        this->timeStep = 0.4 * constants::c / grid->steps.norm();

        fieldSolver.reset(new FieldSolverType(this->grid.get(), this->timeStep));

        initializeGrid();
    }

    void initializeGrid() {
        Int3 begin = Int3(0, 0, 0);
        Int3 end = grid->numCells;

        for (int i = begin.x; i < end.x; i++)
            for (int j = begin.y; j < end.y; j++)
                for (int k = begin.z; k < end.z; k++) {
                    FP3 coords = grid->ExPosition(i, j, k);
                    grid->Ex(i, j, k) = eTest(coords.x, coords.y, coords.z, 0).x;
                    coords = grid->EyPosition(i, j, k);
                    grid->Ey(i, j, k) = eTest(coords.x, coords.y, coords.z, 0).y;
                    coords = grid->EzPosition(i, j, k);
                    grid->Ez(i, j, k) = eTest(coords.x, coords.y, coords.z, 0).z;
                    coords = grid->BxPosition(i, j, k);
                    grid->Bx(i, j, k) = bTest(coords.x, coords.y, coords.z, 0).x;
                    coords = grid->ByPosition(i, j, k);
                    grid->By(i, j, k) = bTest(coords.x, coords.y, coords.z, 0).y;
                    coords = grid->BzPosition(i, j, k);
                    grid->Bz(i, j, k) = bTest(coords.x, coords.y, coords.z, 0).z;
                }
    }

    FP fieldFunc(FP x, FP y, FP z, FP t) {
        FP3 coord(x, y, z);
        int axis0 = (int)this->axis;
        return sin((FP)2.0 * constants::pi / (maxCoords[axis0] - minCoords[axis0]) *
            (coord[axis0] - constants::c * t - minCoords[axis0]));
    }

    FP3 eTest(FP x, FP y, FP z, FP t) {
        CoordinateEnum axisE = CoordinateEnum(((int)axis + 1) % 3);
        FP3 e;
        e[(int)axisE] = fieldFunc(x, y, z, t);
        return e;
    }

    FP3 bTest(FP x, FP y, FP z, FP t) {
        CoordinateEnum axisB = CoordinateEnum(((int)axis + 2) % 3);
        FP3 b;
        b[(int)axisB] = fieldFunc(x, y, z, t);
        return b;
    }
};

typedef ::testing::Types<
    TypeDefinitionsFieldSolverTest<FDTD, YeeGrid, 1, CoordinateEnum::x>,
    TypeDefinitionsFieldSolverTest<FDTD, YeeGrid, 2, CoordinateEnum::x>,
    TypeDefinitionsFieldSolverTest<FDTD, YeeGrid, 2, CoordinateEnum::y>,
    TypeDefinitionsFieldSolverTest<FDTD, YeeGrid, 3, CoordinateEnum::x>,
    TypeDefinitionsFieldSolverTest<FDTD, YeeGrid, 3, CoordinateEnum::y>,
    TypeDefinitionsFieldSolverTest<FDTD, YeeGrid, 3, CoordinateEnum::z>,

    TypeDefinitionsFieldSolverTest<PSTD, PSTDGrid, 1, CoordinateEnum::x>,
    TypeDefinitionsFieldSolverTest<PSTD, PSTDGrid, 2, CoordinateEnum::x>,
    TypeDefinitionsFieldSolverTest<PSTD, PSTDGrid, 2, CoordinateEnum::y>,
    TypeDefinitionsFieldSolverTest<PSTD, PSTDGrid, 3, CoordinateEnum::x>,
    TypeDefinitionsFieldSolverTest<PSTD, PSTDGrid, 3, CoordinateEnum::y>,
    TypeDefinitionsFieldSolverTest<PSTD, PSTDGrid, 3, CoordinateEnum::z>,

    TypeDefinitionsFieldSolverTest<PSATD, PSATDGrid, 1, CoordinateEnum::x>,
    TypeDefinitionsFieldSolverTest<PSATD, PSATDGrid, 2, CoordinateEnum::x>,
    TypeDefinitionsFieldSolverTest<PSATD, PSATDGrid, 2, CoordinateEnum::y>,
    TypeDefinitionsFieldSolverTest<PSATD, PSATDGrid, 3, CoordinateEnum::x>,
    TypeDefinitionsFieldSolverTest<PSATD, PSATDGrid, 3, CoordinateEnum::y>,
    TypeDefinitionsFieldSolverTest<PSATD, PSATDGrid, 3, CoordinateEnum::z>,

    TypeDefinitionsFieldSolverTest<PSATDTimeStraggered, PSATDTimeStraggeredGrid, 1, CoordinateEnum::x>,
    TypeDefinitionsFieldSolverTest<PSATDTimeStraggered, PSATDTimeStraggeredGrid, 2, CoordinateEnum::x>,
    TypeDefinitionsFieldSolverTest<PSATDTimeStraggered, PSATDTimeStraggeredGrid, 2, CoordinateEnum::y>,
    TypeDefinitionsFieldSolverTest<PSATDTimeStraggered, PSATDTimeStraggeredGrid, 3, CoordinateEnum::x>,
    TypeDefinitionsFieldSolverTest<PSATDTimeStraggered, PSATDTimeStraggeredGrid, 3, CoordinateEnum::y>,
    TypeDefinitionsFieldSolverTest<PSATDTimeStraggered, PSATDTimeStraggeredGrid, 3, CoordinateEnum::z>
> types;
TYPED_TEST_CASE(FieldSolverTest, types);


TYPED_TEST(FieldSolverTest, GridInitializationTest)
{
    FP startT = 0;

    Int3 begin = this->fieldSolver->updateEAreaBegin;
    Int3 end = this->fieldSolver->updateEAreaEnd;

    for (int i = begin.x; i < end.x; ++i)
        for (int j = begin.y; j < end.y; ++j)
            for (int k = begin.z; k < end.z; ++k)
            {
                FP3 expectedE, actualE;
                FP3 coords = this->grid->ExPosition(i, j, k);
                expectedE.x = this->eTest(coords.x, coords.y, coords.z, startT).x;
                coords = this->grid->EyPosition(i, j, k);
                expectedE.y = this->eTest(coords.x, coords.y, coords.z, startT).y;
                coords = this->grid->EzPosition(i, j, k);
                expectedE.z = this->eTest(coords.x, coords.y, coords.z, startT).z;
                actualE.x = this->grid->Ex(i, j, k);
                actualE.y = this->grid->Ey(i, j, k);
                actualE.z = this->grid->Ez(i, j, k);
                ASSERT_NEAR_FP3(expectedE, actualE);
            }

    begin = this->fieldSolver->updateBAreaBegin;
    end = this->fieldSolver->updateBAreaEnd;

    for (int i = begin.x; i < end.x; ++i)
        for (int j = begin.y; j < end.y; ++j)
            for (int k = begin.z; k < end.z; ++k)
            {
                FP3 expectedB, actualB;
                FP3 coords = this->grid->BxPosition(i, j, k);
                expectedB.x = this->bTest(coords.x, coords.y, coords.z, startT).x;
                coords = this->grid->ByPosition(i, j, k);
                expectedB.y = this->bTest(coords.x, coords.y, coords.z, startT).y;
                coords = this->grid->BzPosition(i, j, k);
                expectedB.z = this->bTest(coords.x, coords.y, coords.z, startT).z;
                actualB.x = this->grid->Bx(i, j, k);
                actualB.y = this->grid->By(i, j, k);
                actualB.z = this->grid->Bz(i, j, k);
                ASSERT_NEAR_FP3(expectedB, actualB);
            }
}

TYPED_TEST(FieldSolverTest, PeriodicalFieldSolverTest)
{
    // to disable testing of spectral solvers without enabled fftw
#ifndef __USE_FFT__
    SUCCEED();
#else

    boundaryCondition.reset(new FieldSolverType::PeriodicalBoundaryConditionType(fieldSolver.get()));
    fieldSolver->setBoundaryCondition(boundaryCondition.get());

    const int numSteps = (int)(grid->numInternalCells[(int)axis] * grid->steps[(int)axis] /
        (constants::c * fieldSolver->dt));

    for (int step = 0; step < numSteps; ++step)
    {
        this->fieldSolver->updateFields();
    }

    FP finalT = this->fieldSolver->dt * numSteps;

    Int3 begin = this->fieldSolver->updateEAreaBegin;
    Int3 end = this->fieldSolver->updateEAreaEnd;

    for (int i = begin.x; i < end.x; ++i)
        for (int j = begin.y; j < end.y; ++j)
            for (int k = begin.z; k < end.z; ++k)
            {
                FP3 expectedE, actualE;
                FP3 coords = this->grid->ExPosition(i, j, k);
                expectedE.x = this->eTest(coords.x, coords.y, coords.z, finalT).x;
                coords = this->grid->EyPosition(i, j, k);
                expectedE.y = this->eTest(coords.x, coords.y, coords.z, finalT).y;
                coords = this->grid->EzPosition(i, j, k);
                expectedE.z = this->eTest(coords.x, coords.y, coords.z, finalT).z;
                actualE.x = this->grid->Ex(i, j, k);
                actualE.y = this->grid->Ey(i, j, k);
                actualE.z = this->grid->Ez(i, j, k);
                ASSERT_NEAR_FP3(expectedE, actualE);
            }

    begin = this->fieldSolver->updateBAreaBegin;
    end = this->fieldSolver->updateBAreaEnd;

    for (int i = begin.x; i < end.x; ++i)
        for (int j = begin.y; j < end.y; ++j)
            for (int k = begin.z; k < end.z; ++k)
            {
                FP3 expectedB, actualB;
                FP3 coords = this->grid->BxPosition(i, j, k);
                expectedB.x = this->bTest(coords.x, coords.y, coords.z, finalT).x;
                coords = this->grid->ByPosition(i, j, k);
                expectedB.y = this->bTest(coords.x, coords.y, coords.z, finalT).y;
                coords = this->grid->BzPosition(i, j, k);
                expectedB.z = this->bTest(coords.x, coords.y, coords.z, finalT).z;
                actualB.x = this->grid->Bx(i, j, k);
                actualB.y = this->grid->By(i, j, k);
                actualB.z = this->grid->Bz(i, j, k);
                ASSERT_NEAR_FP3(expectedB, actualB);
            }

#endif
    
}