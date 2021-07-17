#include "TestingUtility.h"

#include <sstream>
#include <memory>
#include <type_traits>

#include "ScalarField.h"
#include "Grid.h"
#include "Fdtd.h"
#include "Pstd.h"
#include "Psatd.h"

template <class T1, class T2, bool isFirstType>
struct ActivateType {
    using type = T1;
};

template <class T1, class T2>
struct ActivateType<T1, T2, false> {
    using type = T2;
};

template <class TGrid>
class SaveLoadGridTest : public testing::Test {
public:
    Int3 gridSize;
    FP3 minCoords, maxCoords;
    std::unique_ptr<TGrid> grid, grid2;

    using DataType = typename ActivateType<complexFP, FP, TGrid::isComplex>::type;

    SaveLoadGridTest() : gridSize(11, 21, 2), minCoords(0.0, 0.0, 0.0), maxCoords(1.0, 1.0, 1.0) {
        grid.reset(new TGrid(gridSize, minCoords,
            (maxCoords - minCoords) / ((FP3)gridSize), gridSize));
        grid2.reset(new TGrid());

        for (int i = 0; i < this->grid->Ex.getSize().x; i++)
            for (int j = 0; j < this->grid->Ex.getSize().y; j++)
                for (int k = 0; k < this->grid->Ex.getSize().z; k++) {
                    FP3 coords = grid->ExPosition(i, j, k);
                    grid->Ex(i, j, k) = eTest(coords.x, coords.y, coords.z, 0.0).x;
                    coords = grid->EyPosition(i, j, k);
                    grid->Ey(i, j, k) = eTest(coords.x, coords.y, coords.z, 0.0).y;
                    coords = grid->EzPosition(i, j, k);
                    grid->Ez(i, j, k) = eTest(coords.x, coords.y, coords.z, 0.0).z;

                    coords = grid->BxPosition(i, j, k);
                    grid->Bx(i, j, k) = bTest(coords.x, coords.y, coords.z, 0.0).x;
                    coords = grid->ByPosition(i, j, k);
                    grid->By(i, j, k) = bTest(coords.x, coords.y, coords.z, 0.0).y;
                    coords = grid->BzPosition(i, j, k);
                    grid->Bz(i, j, k) = bTest(coords.x, coords.y, coords.z, 0.0).z;
                }
    }

    FP3 eTest(FP x, FP y, FP z, FP t) {
        return FP3(0, sin(2 * constants::pi * (-constants::c * t + x)), 0);
    }

    FP3 bTest(FP x, FP y, FP z, FP t) {
        return FP3(0, 0, sin(2 * constants::pi * (-constants::c * t + x)));
    }

    bool areGridsEqual(const TGrid& sf1, const TGrid& sf2) {
        return sf1.globalGridDims == sf2.globalGridDims &&
            sf1.steps == sf2.steps &&
            sf1.numInternalCells == sf2.numInternalCells &&
            sf1.sizeStorage == sf2.sizeStorage &&
            sf1.numCells == sf2.numCells &&
            sf1.origin == sf2.origin &&
            sf1.dimensionality == sf2.dimensionality &&
            areFieldsEqual(sf1.Ex, sf2.Ex) && areFieldsEqual(sf1.Ey, sf2.Ey) && areFieldsEqual(sf1.Ez, sf2.Ez) &&
            areFieldsEqual(sf1.Bx, sf2.Bx) && areFieldsEqual(sf1.By, sf2.By) && areFieldsEqual(sf1.Bz, sf2.Bz) &&
            areFieldsEqual(sf1.Jx, sf2.Jx) && areFieldsEqual(sf1.Jy, sf2.Jy) && areFieldsEqual(sf1.Jz, sf2.Jz);
    }

    bool areFieldsEqual(const ScalarField<DataType>& sf1, const ScalarField<DataType>& sf2) {
        if (sf1.getSize() != sf2.getSize())
            return false;
        for (int i = 0; i < sf1.getSize().x; i++)
            for (int j = 0; j < sf1.getSize().y; j++)
                for (int k = 0; k < sf1.getSize().z; k++)
                    if (!(sf1(i, j, k) == sf2(i, j, k)))
                        return false;
        return true;
    }
};

typedef ::testing::Types<
    SimpleGrid,
    YeeGrid,
    PSTDGrid,
    PSATDGrid,
    PSATDTimeStraggeredGrid,
    Grid<complexFP, GridTypes::PSTDGridType>,
    Grid<complexFP, GridTypes::PSATDGridType>,
    Grid<complexFP, GridTypes::PSATDTimeStraggeredGridType>
> typesSaveLoadGridTest;
TYPED_TEST_CASE(SaveLoadGridTest, typesSaveLoadGridTest);

TYPED_TEST(SaveLoadGridTest, can_save_and_load_grid)
{
    std::stringstream sstr;
    grid->save(sstr);
    grid2->load(sstr);

    ASSERT_TRUE(areGridsEqual(*grid, *grid2));
}


template <class TGrid, class TSolver>
struct GridSolverTypes {
    using GridType = TGrid;
    using SolverType = TSolver;
};

template <class TGridSolverTypes>
class SaveLoadGridSolverTest : public BaseGridFixture<typename TGridSolverTypes::GridType> {
public:
    using TGrid = typename TGridSolverTypes::GridType;
    using TSolver = typename TGridSolverTypes::SolverType;

    std::unique_ptr<TGrid> grid2;
    std::unique_ptr<TSolver> solver, solver2;

    const int numSteps = 100;

    virtual void SetUp() {
        BaseGridFixture<TGrid>::SetUp();
        solver.reset(new TSolver(grid, this->timeStep));
        initializeGrid();
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
            * (z - constants::c * t)), 0.0, 0.0);
    }

    FP3 funcB(FP x, FP y, FP z, FP t) {
        return FP3(0.0, sin(2 * constants::pi / (grid->steps.z * grid->numCells.z)
            * (z - constants::c * t)), 0.0);
    }

    void saveGridAndSolver(std::ostream& ostr) {
        grid->save(ostr);
        solver->save(ostr);
    }

    void loadGridAndSolver(std::istream& istr) {
        grid2.reset(new TGrid());
        grid2->load(istr);
        solver2.reset(new TSolver(grid2.get(), this->timeStep));
        solver2->load(istr);
    }
};

typedef ::testing::Types<
    GridSolverTypes<YeeGrid, FDTD>,
    GridSolverTypes<PSTDGrid, PSTD>,
    GridSolverTypes<PSATDGrid, PSATD>,
    GridSolverTypes<PSATDTimeStraggeredGrid, PSATDTimeStraggered>,
    GridSolverTypes<PSATDGrid, PSATDPoisson>,
    GridSolverTypes<PSATDTimeStraggeredGrid, PSATDTimeStraggeredPoisson>
> typesSaveLoadGridSolverTest;
TYPED_TEST_CASE(SaveLoadGridSolverTest, typesSaveLoadGridSolverTest);

TYPED_TEST(SaveLoadGridSolverTest, can_save_and_load_grid_and_solver)
{
    for (int step = 0; step < numSteps / 2; ++step)
        solver->updateFields();

    std::stringstream sstr;
    saveGridAndSolver(sstr);

    for (int step = 0; step < numSteps / 2; ++step)
        solver->updateFields();

    loadGridAndSolver(sstr);

    for (int step = 0; step < numSteps / 2; ++step)
        solver2->updateFields();

    for (int i = 0; i < grid->numCells.x; ++i)
        for (int j = 0; j < grid->numCells.y; ++j)
            for (int k = 0; k < grid->numCells.z; ++k)
            {
                FP3 expectedE(grid->Ex(i, j, k), grid->Ey(i, j, k), grid->Ez(i, j, k));
                FP3 expectedB(grid->Bx(i, j, k), grid->By(i, j, k), grid->Bz(i, j, k));
                FP3 actualE(grid2->Ex(i, j, k), grid2->Ey(i, j, k), grid2->Ez(i, j, k));
                FP3 actualB(grid2->Bx(i, j, k), grid2->By(i, j, k), grid2->Bz(i, j, k));
                ASSERT_NEAR_FP3(expectedE, actualE);
                ASSERT_NEAR_FP3(expectedB, actualB);
            }
}