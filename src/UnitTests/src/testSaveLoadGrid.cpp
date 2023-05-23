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
    std::unique_ptr<TGrid> grid;

    using DataType = typename ActivateType<complexFP, FP, TGrid::isComplex>::type;
    using GridType = TGrid;

    SaveLoadGridTest() : gridSize(11, 21, 2), minCoords(0.0, 0.0, 0.0), maxCoords(1.0, 1.0, 1.0) {
        grid.reset(new TGrid(gridSize, minCoords,
            (maxCoords - minCoords) / ((FP3)gridSize), gridSize));

        for (int i = 0; i < this->grid->Ex.getSize().x; i++)
            for (int j = 0; j < this->grid->Ex.getSize().y; j++)
                for (int k = 0; k < this->grid->Ex.getSize().z; k++) {
                    FP3 coords = this->grid->ExPosition(i, j, k);
                    this->grid->Ex(i, j, k) = eTest(coords.x, coords.y, coords.z, 0.0).x;
                    coords = this->grid->EyPosition(i, j, k);
                    this->grid->Ey(i, j, k) = eTest(coords.x, coords.y, coords.z, 0.0).y;
                    coords = this->grid->EzPosition(i, j, k);
                    this->grid->Ez(i, j, k) = eTest(coords.x, coords.y, coords.z, 0.0).z;

                    coords = this->grid->BxPosition(i, j, k);
                    this->grid->Bx(i, j, k) = bTest(coords.x, coords.y, coords.z, 0.0).x;
                    coords = this->grid->ByPosition(i, j, k);
                    this->grid->By(i, j, k) = bTest(coords.x, coords.y, coords.z, 0.0).y;
                    coords = this->grid->BzPosition(i, j, k);
                    this->grid->Bz(i, j, k) = bTest(coords.x, coords.y, coords.z, 0.0).z;
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
    PSATDTimeStraggeredGrid
> typesSaveLoadGridTest;
TYPED_TEST_CASE(SaveLoadGridTest, typesSaveLoadGridTest);

TYPED_TEST(SaveLoadGridTest, can_save_and_load_grid)
{
    using TGrid = typename SaveLoadGridTest<TypeParam>::GridType;

    std::stringstream sstr;
    this->grid->save(sstr);

    std::unique_ptr<TGrid> grid2(new TGrid());
    grid2->load(sstr);

    ASSERT_TRUE(this->areGridsEqual(*(this->grid), *grid2));
}

TYPED_TEST(SaveLoadGridTest, cannot_save_and_load_grid_of_wrong_type)
{
    using TGrid = typename SaveLoadGridTest<TypeParam>::GridType;
    
    std::stringstream sstr;
    this->grid->save(sstr);

    std::unique_ptr<YeeGrid> grid2(new YeeGrid());
    if (!std::is_same<TGrid, YeeGrid>::value)
        ASSERT_ANY_THROW(grid2->load(sstr));
    else ASSERT_NO_THROW(grid2->load(sstr));
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
        this->solver.reset(new TSolver(this->grid, this->timeStep));
        initializeGrid();
    }

    void initializeGrid() {
        for (int i = 0; i < this->grid->numCells.x; i++)
            for (int j = 0; j < this->grid->numCells.y; j++)
                for (int k = 0; k < this->grid->numCells.z; k++) {
                    FP3 coords = this->grid->ExPosition(i, j, k);
                    this->grid->Ex(i, j, k) = funcE(coords.x, coords.y, coords.z, 0).x;
                    coords = this->grid->EyPosition(i, j, k);
                    this->grid->Ey(i, j, k) = funcE(coords.x, coords.y, coords.z, 0).y;
                    coords = this->grid->EzPosition(i, j, k);
                    this->grid->Ez(i, j, k) = funcE(coords.x, coords.y, coords.z, 0).z;
                    coords = this->grid->BxPosition(i, j, k);
                    this->grid->Bx(i, j, k) = funcB(coords.x, coords.y, coords.z, 0).x;
                    coords = this->grid->ByPosition(i, j, k);
                    this->grid->By(i, j, k) = funcB(coords.x, coords.y, coords.z, 0).y;
                    coords = this->grid->BzPosition(i, j, k);
                    this->grid->Bz(i, j, k) = funcB(coords.x, coords.y, coords.z, 0).z;
                }
    }

    FP3 funcE(FP x, FP y, FP z, FP t) {
        return FP3(sin(2 * constants::pi / (this->grid->steps.z * this->grid->numCells.z)
            * (z - constants::c * t)), 0.0, 0.0);
    }

    FP3 funcB(FP x, FP y, FP z, FP t) {
        return FP3(0.0, sin(2 * constants::pi / (this->grid->steps.z * this->grid->numCells.z)
            * (z - constants::c * t)), 0.0);
    }

    void saveGridAndSolver(std::ostream& ostr) {
        this->grid->save(ostr);
        this->solver->save(ostr);
    }

    void loadGridAndSolver(std::istream& istr) {
        this->grid2.reset(new TGrid());
        this->grid2->load(istr);
        this->solver2.reset(new TSolver(this->grid2.get(), this->timeStep));
        this->solver2->load(istr);
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
    for (int step = 0; step < this->numSteps / 2; ++step)
        this->solver->updateFields();

    std::stringstream sstr;
    this->saveGridAndSolver(sstr);

    for (int step = 0; step < this->numSteps / 2; ++step)
        this->solver->updateFields();

    this->loadGridAndSolver(sstr);

    for (int step = 0; step < this->numSteps / 2; ++step)
        this->solver2->updateFields();

    for (int i = 0; i < this->grid->numCells.x; ++i)
        for (int j = 0; j < this->grid->numCells.y; ++j)
            for (int k = 0; k < this->grid->numCells.z; ++k)
            {
                FP3 expectedE(this->grid->Ex(i, j, k), this->grid->Ey(i, j, k), this->grid->Ez(i, j, k));
                FP3 expectedB(this->grid->Bx(i, j, k), this->grid->By(i, j, k), this->grid->Bz(i, j, k));
                FP3 actualE(this->grid2->Ex(i, j, k), this->grid2->Ey(i, j, k), this->grid2->Ez(i, j, k));
                FP3 actualB(this->grid2->Bx(i, j, k), this->grid2->By(i, j, k), this->grid2->Bz(i, j, k));
                ASSERT_NEAR_FP3(expectedE, actualE);
                ASSERT_NEAR_FP3(expectedB, actualB);
            }
}
