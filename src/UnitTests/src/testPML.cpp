#include "TestingUtility.h"

#include "Enums.h"

#include "Fdtd.h"
#include "Pstd.h"
#include "Psatd.h"
#include "PsatdTimeStraggered.h"


template <class TFieldSolver, class TGrid, int VDimension, CoordinateEnum VAxis>
struct TypeDefinitionsPMLTest
{
    using FieldSolverType = TFieldSolver;
    using GridType = TGrid;
    static const int dimension = VDimension;
    static const CoordinateEnum axis = VAxis;
};

template <class TTypeDefinitionsPMLTest>
class PMLTest : public BaseFixture {
public:

    using FieldSolverType = typename TTypeDefinitionsPMLTest::FieldSolverType;
    using GridType = typename TTypeDefinitionsPMLTest::GridType;
    using PeriodicalFieldGeneratorType = typename FieldSolverType::PeriodicalFieldGeneratorType;

    const int dimension = TTypeDefinitionsPMLTest::dimension;
    const CoordinateEnum axis = TTypeDefinitionsPMLTest::axis;

    const int gridSize1d = 32;
    const int pmlSize1d = 4;

    FP3 pmlLeftEnd;
    FP3 pmlRightStart;

    std::unique_ptr<FieldSolverType> fieldSolver;
    std::unique_ptr<GridType> grid;
    std::unique_ptr<PeriodicalFieldGeneratorType> generator;

    Int3 gridSize;
    Int3 pmlSize;
    FP3 gridStep;
    FP timeStep = 0;
    FP3 minCoords, maxCoords;

    const FP relatedEnergyThreshold = 1e-2;

    virtual void SetUp() {
        gridSize = Int3(1, 1, 1);
        for (int d = 0; d < dimension; d++) {
            gridSize[d] = gridSize1d;
        }

        pmlSize = Int3(0, 0, 0);
        pmlSize[(int)axis] = pmlSize1d;

        this->minCoords = FP3(0, 0, 0);
        this->maxCoords = FP3(gridSize.x * constants::c, gridSize.y * constants::c, gridSize.z * constants::c);
        this->gridStep = FP3(constants::c, constants::c, constants::c);
        this->pmlLeftEnd = this->minCoords + pmlSize * this->gridStep;
        this->pmlRightStart = this->maxCoords - pmlSize * this->gridStep;

        this->grid.reset(new GridType(this->gridSize, this->minCoords, this->gridStep, this->gridSize));
        initializeGrid();

        // should satisfy the Courant's condition for all solvers
        this->timeStep = 0.1 * constants::c / grid->steps[(int)axis];

        fieldSolver.reset(new FieldSolverType(this->grid.get(), this->timeStep));
        fieldSolver->setPML(pmlSize.x, pmlSize.y, pmlSize.z);

        generator.reset(new PeriodicalFieldGeneratorType(fieldSolver.get()));
        fieldSolver->setFieldGenerator(generator.get());

    }

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
        FP3 coord(x, y, z);
        CoordinateEnum axisE = CoordinateEnum(((int)axis + 1) % 3);
        FP3 e;
        e[(int)axisE] = sin((FP)2.0 * constants::pi / (pmlRightStart[(int)axis] - pmlLeftEnd[(int)axis]) *
            (coord[(int)axis] - constants::c * t - pmlLeftEnd[(int)axis]));
        return e;
    }

    FP3 funcB(FP x, FP y, FP z, FP t) {
        FP3 coord(x, y, z);
        CoordinateEnum axisB = CoordinateEnum(((int)axis + 2) % 3);
        FP3 b;
        b[(int)axisB] = sin((FP)2.0 * constants::pi / (pmlRightStart[(int)axis] - pmlLeftEnd[(int)axis]) *
            (coord[(int)axis] - constants::c * t - pmlLeftEnd[(int)axis]));
        return b;
    }

    FP computeEnergy() {
        FP energy = 0;
        for (int i = 0; i < grid->numCells.x; i++)
            for (int j = 0; j < grid->numCells.y; j++)
                for (int k = 0; k < grid->numCells.z; k++)
                    energy += pow(grid->Ex(i, j, k), 2) + pow(grid->Ey(i, j, k), 2) + pow(grid->Ez(i, j, k), 2) +
                    pow(grid->Bx(i, j, k), 2) + pow(grid->By(i, j, k), 2) + pow(grid->Bz(i, j, k), 2);
        return energy * gridStep.volume();
    }

};

typedef ::testing::Types <
    TypeDefinitionsPMLTest<FDTD, YeeGrid, 1, CoordinateEnum::x>,
    TypeDefinitionsPMLTest<FDTD, YeeGrid, 2, CoordinateEnum::x>,
    TypeDefinitionsPMLTest<FDTD, YeeGrid, 2, CoordinateEnum::y>,
    TypeDefinitionsPMLTest<FDTD, YeeGrid, 3, CoordinateEnum::x>,
    TypeDefinitionsPMLTest<FDTD, YeeGrid, 3, CoordinateEnum::y>,
    TypeDefinitionsPMLTest<FDTD, YeeGrid, 3, CoordinateEnum::z>,

    TypeDefinitionsPMLTest<PSTD, PSTDGrid, 1, CoordinateEnum::x>,
    TypeDefinitionsPMLTest<PSTD, PSTDGrid, 2, CoordinateEnum::x>,
    TypeDefinitionsPMLTest<PSTD, PSTDGrid, 2, CoordinateEnum::y>,
    TypeDefinitionsPMLTest<PSTD, PSTDGrid, 3, CoordinateEnum::x>,
    TypeDefinitionsPMLTest<PSTD, PSTDGrid, 3, CoordinateEnum::y>,
    TypeDefinitionsPMLTest<PSTD, PSTDGrid, 3, CoordinateEnum::z>,
    
    TypeDefinitionsPMLTest<PSATD, PSATDGrid, 1, CoordinateEnum::x>,
    TypeDefinitionsPMLTest<PSATD, PSATDGrid, 2, CoordinateEnum::x>,
    TypeDefinitionsPMLTest<PSATD, PSATDGrid, 2, CoordinateEnum::y>,
    TypeDefinitionsPMLTest<PSATD, PSATDGrid, 3, CoordinateEnum::x>,
    TypeDefinitionsPMLTest<PSATD, PSATDGrid, 3, CoordinateEnum::y>,
    TypeDefinitionsPMLTest<PSATD, PSATDGrid, 3, CoordinateEnum::z>,
    
    TypeDefinitionsPMLTest<PSATDTimeStraggered, PSATDTimeStraggeredGrid, 1, CoordinateEnum::x>,
    TypeDefinitionsPMLTest<PSATDTimeStraggered, PSATDTimeStraggeredGrid, 2, CoordinateEnum::x>,
    TypeDefinitionsPMLTest<PSATDTimeStraggered, PSATDTimeStraggeredGrid, 2, CoordinateEnum::y>,
    TypeDefinitionsPMLTest<PSATDTimeStraggered, PSATDTimeStraggeredGrid, 3, CoordinateEnum::x>,
    TypeDefinitionsPMLTest<PSATDTimeStraggered, PSATDTimeStraggeredGrid, 3, CoordinateEnum::y>,
    TypeDefinitionsPMLTest<PSATDTimeStraggered, PSATDTimeStraggeredGrid, 3, CoordinateEnum::z>
> types;
TYPED_TEST_SUITE(PMLTest, types);


TYPED_TEST(PMLTest, PmlTest) {
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
