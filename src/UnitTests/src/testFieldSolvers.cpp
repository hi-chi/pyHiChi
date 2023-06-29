#include "TestingUtility.h"

#include "Fdtd.h"
#include "Pstd.h"
#include "Psatd.h"
#include "PsatdTimeStaggered.h"

template <class TTypeDefinitionsFieldTest>
class FieldSolverTest : public BaseFixture {
public:
    using FieldSolverType = typename TTypeDefinitionsFieldTest::FieldSolverType;
    using GridType = typename TTypeDefinitionsFieldTest::FieldSolverType::GridType;

    const int dimension = TTypeDefinitionsFieldTest::dimension;
    const CoordinateEnum axis = TTypeDefinitionsFieldTest::axis;

    const int gridSizeLongitudinal = 32;
    const int gridSizeTransverse = 8;

    std::unique_ptr<FieldSolverType> fieldSolver;
    std::unique_ptr<GridType> grid;
    Int3 gridSize;
    FP3 gridStep;
    FP timeStep = 0;
    FP3 minCoords, maxCoords;
    int numSteps = 0;

    const FP maxError = 1e-2;

    virtual void SetUp() {
        gridSize = Int3(1, 1, 1);
        for (int d = 0; d < dimension; d++) {
            gridSize[d] = gridSizeTransverse;
        }
        gridSize[(int)axis] = gridSizeLongitudinal;

        this->minCoords = FP3(0, 0, 0);
        this->maxCoords = constants::c * (FP3)gridSize;
        this->gridStep = (this->maxCoords - this->minCoords) / (FP3)this->gridSize;

        this->grid.reset(new GridType(this->gridSize, this->minCoords, this->gridStep, this->gridSize));

        this->timeStep = 0.5 * FieldSolverType::getCourantConditionTimeStep(this->gridStep);
        this->numSteps = (int)((this->maxCoords - this->minCoords)[(int)axis] /
            (constants::c * this->timeStep) * 0.2);

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
                    grid->Ex(i, j, k) = eFunc(coords.x, coords.y, coords.z, 0).x;
                    coords = grid->EyPosition(i, j, k);
                    grid->Ey(i, j, k) = eFunc(coords.x, coords.y, coords.z, 0).y;
                    coords = grid->EzPosition(i, j, k);
                    grid->Ez(i, j, k) = eFunc(coords.x, coords.y, coords.z, 0).z;
                    coords = grid->BxPosition(i, j, k);
                    grid->Bx(i, j, k) = bFunc(coords.x, coords.y, coords.z, 0).x;
                    coords = grid->ByPosition(i, j, k);
                    grid->By(i, j, k) = bFunc(coords.x, coords.y, coords.z, 0).y;
                    coords = grid->BzPosition(i, j, k);
                    grid->Bz(i, j, k) = bFunc(coords.x, coords.y, coords.z, 0).z;
                }
    }

    FP fieldFunc(FP x, FP y, FP z, FP t) {
        FP3 coord(x, y, z);
        int axis0 = (int)this->axis;
        return sin((FP)2.0 * constants::pi / (maxCoords[axis0] - minCoords[axis0]) *
            (coord[axis0] - constants::c * t - minCoords[axis0]));
    }

    FP3 eFunc(FP x, FP y, FP z, FP t) {
        CoordinateEnum axisE = CoordinateEnum(((int)axis + 1) % 3);
        FP3 e;
        e[(int)axisE] = fieldFunc(x, y, z, t);
        return e;
    }

    FP3 bFunc(FP x, FP y, FP z, FP t) {
        CoordinateEnum axisB = CoordinateEnum(((int)axis + 2) % 3);
        FP3 b;
        b[(int)axisB] = fieldFunc(x, y, z, t);
        return b;
    }
};

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
TYPED_TEST_CASE(FieldSolverTest, types);


TYPED_TEST(FieldSolverTest, PeriodicalFieldSolverTest)
{
    // to disable testing of spectral solvers without enabled fftw
#ifndef __USE_FFT__
    SUCCEED();
#else

    this->fieldSolver->setPeriodicalBoundaryCondition();

    for (int step = 0; step < this->numSteps; ++step)
    {
        this->fieldSolver->updateFields();
    }

    FP finalT = this->fieldSolver->dt * this->numSteps;

    Int3 begin = this->fieldSolver->internalIndexBegin;
    Int3 end = this->fieldSolver->internalIndexEnd;

    for (int i = begin.x; i < end.x; ++i)
        for (int j = begin.y; j < end.y; ++j)
            for (int k = begin.z; k < end.z; ++k)
            {
                FP3 expectedE, actualE;
                FP3 coords = this->grid->ExPosition(i, j, k);
                expectedE.x = this->eFunc(coords.x, coords.y, coords.z, finalT).x;
                coords = this->grid->EyPosition(i, j, k);
                expectedE.y = this->eFunc(coords.x, coords.y, coords.z, finalT).y;
                coords = this->grid->EzPosition(i, j, k);
                expectedE.z = this->eFunc(coords.x, coords.y, coords.z, finalT).z;
                actualE.x = this->grid->Ex(i, j, k);
                actualE.y = this->grid->Ey(i, j, k);
                actualE.z = this->grid->Ez(i, j, k);
                ASSERT_NEAR(expectedE.norm(), actualE.norm(), this->maxError);
            }

    for (int i = begin.x; i < end.x; ++i)
        for (int j = begin.y; j < end.y; ++j)
            for (int k = begin.z; k < end.z; ++k)
            {
                FP3 expectedB, actualB;
                FP3 coords = this->grid->BxPosition(i, j, k);
                expectedB.x = this->bFunc(coords.x, coords.y, coords.z, finalT).x;
                coords = this->grid->ByPosition(i, j, k);
                expectedB.y = this->bFunc(coords.x, coords.y, coords.z, finalT).y;
                coords = this->grid->BzPosition(i, j, k);
                expectedB.z = this->bFunc(coords.x, coords.y, coords.z, finalT).z;
                actualB.x = this->grid->Bx(i, j, k);
                actualB.y = this->grid->By(i, j, k);
                actualB.z = this->grid->Bz(i, j, k);
                ASSERT_NEAR(expectedB.norm(), actualB.norm(), this->maxError);
            }

#endif
}

TYPED_TEST(FieldSolverTest, SetTimeStepTest) {

    const int pmlSize = 1, genIndex = 2;
    const FP newTimeStep = this->timeStep * 0.25;
    const Int3 genLeftIndex3d = this->grid->correctNumCellsAccordingToDim(Int3(genIndex, genIndex, genIndex));
    const Int3 genRightIndex3d = this->gridSize - genLeftIndex3d;
    const Int3 pmlSize3d = this->grid->correctNumCellsAccordingToDim(Int3(pmlSize, pmlSize, pmlSize));

    this->fieldSolver->setPeriodicalBoundaryCondition();

    this->fieldSolver->setPML(pmlSize3d.x, pmlSize3d.y, pmlSize3d.z);

    this->fieldSolver->setFieldGenerator(genLeftIndex3d, genRightIndex3d,
        field_generator::defaultFieldFunction, field_generator::defaultFieldFunction,
        field_generator::defaultFieldFunction, field_generator::defaultFieldFunction,
        field_generator::defaultFieldFunction, field_generator::defaultFieldFunction);

    this->fieldSolver->setTimeStep(newTimeStep);

    ASSERT_EQ(newTimeStep, this->fieldSolver->pml->dt);
    ASSERT_EQ(newTimeStep, this->fieldSolver->generator->dt);
}
