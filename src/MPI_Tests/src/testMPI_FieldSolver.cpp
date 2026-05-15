#include "TestingUtility.h"

#include "gtest/gtest.h"

#include "Fdtd.h"
#include "Pstd.h"
#include "Psatd.h"
#include "PsatdTimeStaggered.h"

#include "HiChi_MPI.h"

template <class TTypeDefinitionsFieldTest>
class MPI_FieldSolverTest : public BaseFixture {
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
    FP3 mainMinCoords = FP3(0, 0, 0);
    FP3 mainMaxCoords = FP3(0, 0, 0);
    FP timeStep = 0;
    FP3 minCoords, maxCoords;
    int numSteps = 0;

    const FP maxError = 1e-2;

    // MPI
    std::shared_ptr<mpi::Topology> topology;
    std::shared_ptr<mpi::FieldExchanger> fieldExchanger;

    virtual void SetUp() {
        int mpi_size;
        
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

        Int3 sections = Int3(4, 1, 1);
        int topologySize = sections.x * sections.y * sections.z;

        if (topologySize > mpi_size)
        {
            GTEST_SKIP() << "Not enough processes for FieldSolverTest";
            return;
        }

        topology = std::make_shared<mpi::Topology>(sections, mpi::Topology::LoopType::LoopXYZ, mpi_size);

        if (!topology->isValidOnThisRank())
        {
            GTEST_SKIP() << "This rank is not required";
            return;
        }

        int mpi_rank;
        MPI_Comm_rank(topology->getTopologyCommunicator(), &mpi_rank);

        Int3 mainGridSize = Int3(1, 1, 1);
        for (int d = 0; d < dimension; d++) {
            mainGridSize[d] = gridSizeTransverse;
        }
        mainGridSize[(int)axis] = gridSizeLongitudinal;

        mainMinCoords = FP3(0, 0, 0);
        mainMaxCoords = constants::c * (FP3)mainGridSize;
        this->gridStep = (mainMaxCoords - mainMinCoords) / (FP3)mainGridSize;

        int div = mainGridSize.x / 4;
        std::vector<int> divisions[3] = {
            {div, 2 * div, 3 * div},
            {}, 
            {}
        };
        mpi::GridSlicer::getSubGridParameters(  this->minCoords, this->gridSize, 
                                                mainMinCoords, mainGridSize, 
                                                this->gridStep, divisions, 
                                                mpi_rank, *topology);

        this->grid.reset(new GridType(this->gridSize, this->minCoords, this->gridStep, this->gridSize));

        fieldExchanger = std::make_shared<mpi::FieldExchanger>(grid->numCells, grid->numExternalCells);

        this->timeStep = 0.5 * FieldSolverType::getCourantConditionTimeStep(this->gridStep);
        this->numSteps = (int)((mainMaxCoords - mainMinCoords)[(int)axis] /
            (constants::c * this->timeStep) * 0.2);

        fieldSolver.reset(new FieldSolverType(this->grid.get(), this->timeStep));

        mpi::BoundaryType boundaries[6] = { 
            mpi::BoundaryType::Periodic, // +X
            mpi::BoundaryType::Periodic, // -X
            mpi::BoundaryType::Periodic, // +Y
            mpi::BoundaryType::Periodic, // -Y
            mpi::BoundaryType::Periodic, // +Z
            mpi::BoundaryType::Periodic, // -Z
        };
        using BoundaryManager = mpi::FieldBoundaryManager<FieldSolverType, GridType, TTypeDefinitionsFieldTest::ReflectBoundaryConditionType>;
        BoundaryManager::setupBoundaryConditions(fieldSolver, boundaries, topology, fieldExchanger, mpi_rank);

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
        return sin((FP)2.0 * constants::pi / (mainMaxCoords[axis0] - mainMinCoords[axis0]) *
            (coord[axis0] - constants::c * t - mainMinCoords[axis0]));
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

#ifndef __USE_FFT__

typedef ::testing::Types <
    MPI_TypeDefinitionsFieldTest<FDTD, ReflectBoundaryConditionMonoDirectionFdtd, 3, CoordinateEnum::x>,
    MPI_TypeDefinitionsFieldTest<FDTD, ReflectBoundaryConditionMonoDirectionFdtd, 3, CoordinateEnum::y>,
    MPI_TypeDefinitionsFieldTest<FDTD, ReflectBoundaryConditionMonoDirectionFdtd, 3, CoordinateEnum::z>
> types;

#else

typedef ::testing::Types <
    MPI_TypeDefinitionsFieldTest<FDTD, ReflectBoundaryConditionMonoDirectionFdtd, 3, CoordinateEnum::x>,
    MPI_TypeDefinitionsFieldTest<FDTD, ReflectBoundaryConditionMonoDirectionFdtd, 3, CoordinateEnum::y>,
    MPI_TypeDefinitionsFieldTest<FDTD, ReflectBoundaryConditionMonoDirectionFdtd, 3, CoordinateEnum::z>

    //MPI_TypeDefinitionsFieldTest<PSTD, void, 3, CoordinateEnum::x>
    // MPI_TypeDefinitionsFieldTest<PSTD, void, 3, CoordinateEnum::y>,
    // MPI_TypeDefinitionsFieldTest<PSTD, void, 3, CoordinateEnum::z>

    // MPI_TypeDefinitionsFieldTest<PSATD, void, 3, CoordinateEnum::x>,
    // MPI_TypeDefinitionsFieldTest<PSATD, void, 3, CoordinateEnum::y>,
    // MPI_TypeDefinitionsFieldTest<PSATD, void, 3, CoordinateEnum::z>,

    // MPI_TypeDefinitionsFieldTest<PSATDTimeStaggered, void, 3, CoordinateEnum::x>,
    // MPI_TypeDefinitionsFieldTest<PSATDTimeStaggered, void, 3, CoordinateEnum::y>,
    // MPI_TypeDefinitionsFieldTest<PSATDTimeStaggered, void, 3, CoordinateEnum::z>
> types;

#endif

TYPED_TEST_CASE(MPI_FieldSolverTest, types);


TYPED_TEST(MPI_FieldSolverTest, PeriodicalFieldSolverTest)
{
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
                    
                EXPECT_NEAR(expectedE.norm(), actualE.norm(), this->maxError);
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
                EXPECT_NEAR(expectedB.norm(), actualB.norm(), this->maxError);
            }
}