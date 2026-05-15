#include "TestingUtility.h"

#include "Fdtd.h"
#include "Pstd.h"
#include "Psatd.h"
#include "PsatdTimeStaggered.h"

#include "HiChi_MPI.h"

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
    Int3 gridSize = Int3(0, 0, 0);
    FP3 gridStep = FP3(0, 0, 0);
    FP timeStep = 0;
    FP3 mainMinCoords = FP3(0, 0, 0);
    FP3 mainMaxCoords = FP3(0, 0, 0);
    FP3 minCoords = FP3(0, 0, 0);
    FP3 maxCoords = FP3(0, 0, 0);
    int numSteps = 0;

    const FP maxError = 0.2;

    // MPI
    Int3 sections;
    std::shared_ptr<mpi::Topology> topology;
    std::shared_ptr<mpi::FieldExchanger> fieldExchanger;

    BoundaryConditionTest() {}

    virtual void SetUp() override
    {
        int mpi_size;
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

        sections = Int3(2, 2, 1);

        int topologySize = sections.x * sections.y * sections.z;

        if (topologySize > mpi_size)
            return;

        bool loopMask[3] = {true, true, true};
        loopMask[static_cast<int>(axis)] = false;

        topology = std::make_shared<mpi::Topology>(sections, loopMask, mpi_size);

        if (!topology->isValidOnThisRank())
            return;

        int mpi_rank;
        MPI_Comm_rank(topology->getTopologyCommunicator(), &mpi_rank);

        Int3 mainGridSize = Int3(1, 1, 1);
        for (int d = 0; d < dimension; d++) {
            mainGridSize[d] = gridSizeTransverse;
        }
        mainGridSize[(int)axis] = gridSizeLongitudinal;

        mainMinCoords = FP3(0, 0, 0);
        mainMaxCoords = (FP3)mainGridSize * constants::c;
        this->gridStep = (mainMaxCoords - mainMinCoords) / (FP3)mainGridSize;

        std::vector<int> divisions[3] = {
            {mainGridSize.x / 2},
            {mainGridSize.y / 2}, 
            {}
        };
        mpi::GridSlicer::getSubGridParameters(  this->minCoords, this->gridSize, 
                                                mainMinCoords, mainGridSize, 
                                                this->gridStep, divisions, 
                                                mpi_rank, *topology);

        this->grid.reset(new GridType(this->gridSize, this->minCoords, this->gridStep, this->gridSize));

        fieldExchanger = std::make_shared<mpi::FieldExchanger>(grid->numCells, grid->numExternalCells);

        this->timeStep = 0.5 * FieldSolverType::getCourantConditionTimeStep(this->gridStep);
        // wave goes throught area 2 times to check both left and right boundary conditions
        this->numSteps = 2 * (int)((this->maxCoords - this->minCoords)[(int)axis] / 
            (constants::c * this->timeStep));

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

    FP harrisFunction(FP x, FP t, FP a, FP b) {
        FP coord = (x - constants::c * t) / (b - a);
        if (coord < 0.0 || coord >= 1.0) return 0.0;
        const FP pi2 = 2.0 * constants::pi;
        return 0.03125 * (10.0 - 15.0 * cos(pi2 * coord) + 6.0 * cos(2.0 * pi2 * coord) - cos(3.0 * pi2 * coord));
    }

    FP fieldFunc(FP x, FP y, FP z, FP t) {
        int axis0 = (int)this->axis;
        return harrisFunction(FP3(x, y, z)[axis0], t, mainMinCoords[axis0], mainMaxCoords[axis0]);
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

// NOTE: Periodic Boundary conditions work as long as testMPI_FieldExchanger passes

template <class TTypeDefinitionsFieldTest>
class ReflectBoundaryConditionTest : public BoundaryConditionTest<TTypeDefinitionsFieldTest> {
public:
    ReflectBoundaryConditionTest() {}

    virtual void SetUp() override {
        // this->fieldSolver->setReflectBoundaryConditions(this->axis);

        // for (int d = 0; d < 3; d++) {
        //     int dim = ((int)this->axis + d) % 3;
        //     if (dim < this->grid->dimensionality)
        //         this->fieldSolver->setPeriodicalBoundaryConditions((CoordinateEnum)dim);
        // }

        BoundaryConditionTest<TTypeDefinitionsFieldTest>::SetUp();

        int mpi_size;
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

        int topologySize = this->sections.x * this->sections.y * this->sections.z;

        if (topologySize > mpi_size)
        {    
            GTEST_SKIP() << "Not enough processes";
            return;
        }

        if (!this->topology->isValidOnThisRank())
        {
            GTEST_SKIP() << "This rank is not required";
            return;
        }

        int mpi_rank;
        MPI_Comm_rank(this->topology->getTopologyCommunicator(), &mpi_rank);

        mpi::BoundaryType boundaries[6] = { 
            mpi::BoundaryType::Periodic, // +X
            mpi::BoundaryType::Periodic, // -X
            mpi::BoundaryType::Periodic, // +Y
            mpi::BoundaryType::Periodic, // -Y
            mpi::BoundaryType::Periodic, // +Z
            mpi::BoundaryType::Periodic, // -Z
        };

        int reflectIndex = static_cast<int>(axis) * 2;
        boundaries[reflectIndex] = mpi::BoundaryType::Reflect;
        boundaries[reflectIndex + 1] = mpi::BoundaryType::Reflect;

        using BoundaryManager = mpi::FieldBoundaryManager<FieldSolverType, GridType, TTypeDefinitionsFieldTest::ReflectBoundaryConditionType>;
        BoundaryManager::setupBoundaryConditions(this->fieldSolver, boundaries, this->topology, this->fieldExchanger, mpi_rank);
    }

};

typedef ::testing::Types <
    MPI_TypeDefinitionsFieldTest<FDTD, ReflectBoundaryConditionMonoDirectionFdtd, 3, CoordinateEnum::x>,
    MPI_TypeDefinitionsFieldTest<FDTD, ReflectBoundaryConditionMonoDirectionFdtd, 3, CoordinateEnum::y>,
    MPI_TypeDefinitionsFieldTest<FDTD, ReflectBoundaryConditionMonoDirectionFdtd, 3, CoordinateEnum::z>
> typesReflect;

TYPED_TEST_CASE(ReflectBoundaryConditionTest, typesReflect);

TYPED_TEST(ReflectBoundaryConditionTest, MixedPeriodicAndReflectBoundaryConditionTest)
{
    for (int step = 0; step < this->numSteps; ++step)
    {
        this->fieldSolver->updateFields();
    }

    // signal should be the same as at the beginning
    // because signal is symmetric
    FP startT = 0;

    Int3 begin = this->fieldSolver->internalIndexBegin;
    Int3 end = this->fieldSolver->internalIndexEnd;

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
                ASSERT_NEAR((expectedE - actualE).norm(), 0.0, this->maxError);
            }

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
                ASSERT_NEAR((expectedB - actualB).norm(), 0.0, this->maxError);
            }
}
