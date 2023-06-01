#include "TestingUtility.h"

#include "Fdtd.h"

template <class TTypeDefinitionsFieldTest>
class BoundaryConditionTest : public BaseFixture {
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

    const FP maxError = 0.2;

    BoundaryConditionTest() {
        gridSize = Int3(1, 1, 1);
        for (int d = 0; d < dimension; d++) {
            gridSize[d] = gridSizeTransverse;
        }
        gridSize[(int)axis] = gridSizeLongitudinal;

        FP d = constants::c;
        this->minCoords = FP3(0, 0, 0);
        this->maxCoords = d * (FP3)gridSize;
        this->gridStep = FP3(d, d, d);

        this->grid.reset(new GridType(this->gridSize, this->minCoords, this->gridStep, this->gridSize));

        // should satisfy the Courant's condition for all solvers
        this->timeStep = 0.4;  // = 0.4 * d / constants::c
        this->numSteps = gridSize[(int)axis] * 10 / 4;  // to cross the entire area exactly

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
        int axis0 = (int)this->axis;
        FP L = (maxCoords - minCoords)[axis0];
        FP coord = (FP3(x, y, z)[axis0] - constants::c * t) / L;
        if (coord < 0.0 || coord >= 1.0) coord = 0.0;
        const FP pi2 = 2.0 * constants::pi;
        // H = harris function
        FP H = 0.03125 * (10.0 - 15.0 * cos(pi2 * coord) + 6.0 * cos(2.0 * pi2 * coord) - cos(3.0 * pi2 * coord));
        return H;
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

typedef ::testing::Types <
    TypeDefinitionsFieldTest<FDTD, 1, CoordinateEnum::x>,
    TypeDefinitionsFieldTest<FDTD, 2, CoordinateEnum::x>,
    TypeDefinitionsFieldTest<FDTD, 2, CoordinateEnum::y>,
    TypeDefinitionsFieldTest<FDTD, 3, CoordinateEnum::x>,
    TypeDefinitionsFieldTest<FDTD, 3, CoordinateEnum::y>,
    TypeDefinitionsFieldTest<FDTD, 3, CoordinateEnum::z>
> types;


template <class TTypeDefinitionsFieldTest>
class PeriodicBoundaryConditionTest : public BoundaryConditionTest<TTypeDefinitionsFieldTest> {
public:

    using BoundaryConditionType =
        typename BoundaryConditionTest<TTypeDefinitionsFieldTest>::FieldSolverType::PeriodicalBoundaryConditionType;
    std::unique_ptr<BoundaryConditionType> boundaryCondition;

    PeriodicBoundaryConditionTest() {
        boundaryCondition.reset(new BoundaryConditionType(fieldSolver.get()));
        fieldSolver->setBoundaryCondition(boundaryCondition.get());
    }

};
TYPED_TEST_CASE(PeriodicBoundaryConditionTest, types);

TYPED_TEST(PeriodicBoundaryConditionTest, PeriodicBoundaryConditionTest)
{
    // to disable testing of spectral solvers without enabled fftw
#ifndef __USE_FFT__
    SUCCEED();
#else

    for (int step = 0; step < numSteps; ++step)
    {
        this->fieldSolver->updateFields();
    }

    // signal should be the same as at the beginning
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
                ASSERT_NEAR(expectedE.norm(), actualE.norm(), maxError);
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
                ASSERT_NEAR(expectedB.norm(), actualB.norm(), maxError);
            }

#endif

}


template <class TTypeDefinitionsFieldTest>
class ReflectBoundaryConditionTest : public BoundaryConditionTest<TTypeDefinitionsFieldTest> {
public:

    using BoundaryConditionType =
        typename BoundaryConditionTest<TTypeDefinitionsFieldTest>::FieldSolverType::ReflectBoundaryConditionType;
    std::unique_ptr<BoundaryConditionType> boundaryCondition;

    ReflectBoundaryConditionTest() {
        boundaryCondition.reset(new BoundaryConditionType(fieldSolver.get()));
        fieldSolver->setBoundaryCondition(boundaryCondition.get());
    }

};
TYPED_TEST_CASE(ReflectBoundaryConditionTest, types);

TYPED_TEST(ReflectBoundaryConditionTest, ReflectBoundaryConditionTest)
{
    // to disable testing of spectral solvers without enabled fftw
#ifndef __USE_FFT__
    SUCCEED();
#else

    for (int step = 0; step < numSteps; ++step)
    {
        this->fieldSolver->updateFields();
    }

    // signal should be the same as at the beginning
    // because signal is symmetric
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
                ASSERT_NEAR(expectedE.norm(), actualE.norm(), maxError);
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
                ASSERT_NEAR(expectedB.norm(), actualB.norm(), maxError);
            }

#endif

}
