#include "TestingUtility.h"

#include <sstream>
#include <memory>
#include <type_traits>
#include <functional>

#include "ScalarField.h"
#include "Grid.h"
#include "Fdtd.h"
#include "Pstd.h"
#include "Psatd.h"
#include "PsatdTimeStaggered.h"
#include "AnalyticalFieldSolver.h"

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
    PSATDTimeStaggeredGrid
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


template <class TSolver>
class SaveLoadGridSolverTest : public BaseGridFixture<typename TSolver::GridType> {
public:
    using TGrid = typename TSolver::GridType;

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

#ifndef __USE_FFT__

typedef ::testing::Types<
    FDTD
> typesSaveLoadSolverModulesTest;

#else

typedef ::testing::Types<
    FDTD,
    PSTD,
    PSATD,
    PSATDTimeStaggered,
    PSATDPoisson,
    PSATDTimeStaggeredPoisson
> typesSaveLoadGridSolverTest;

#endif

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


template <class TSolver>
class SaveLoadSolverModulesTest : public BaseGridFixture<typename TSolver::GridType> {
public:
    using TGrid = typename TSolver::GridType;
    using PeriodicalBCType = typename TSolver::PeriodicalBoundaryConditionType;

    std::unique_ptr<TGrid> grid2;
    std::unique_ptr<TSolver> solver, solver2;

    const int numSteps = 10;

    virtual void SetUp() {
        this->maxAbsoluteError = std::is_base_of<RealFieldSolver<typename TSolver::SchemeParams>, TSolver>::value?
            (FP)1e-1: (FP)1e-5;
        this->maxRelativeError = (FP)1e-1;
        this->timeStep = 1e-13;

        Int3 gridSize(19, 10, 11);
        this->minCoords = FP3(-1.0, 0.0, 0.0);
        this->maxCoords = FP3(1.0, 1.0, 1.0);
        FP3 steps = (this->maxCoords - this->minCoords) / (FP3)gridSize;
        this->grid = new TGrid(gridSize, this->minCoords, steps, gridSize);

        this->solver.reset(new TSolver(this->grid, this->timeStep));
        this->solver->setPML(4, 3, 3);
        this->solver->setFieldGenerator(Int3(6, 5, 5), gridSize - Int3(6, 5, 5),
            field_generator::defaultFieldFunction, field_generator::defaultFieldFunction,
            field_generator::defaultFieldFunction, field_generator::defaultFieldFunction,
            field_generator::defaultFieldFunction, field_generator::defaultFieldFunction);
        this->solver->setPeriodicalBoundaryConditions();

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
        this->solver2.reset(new TSolver(this->grid2.get()));
        this->solver2->load(istr);
    }

    bool compareFPVectors(const std::vector<FP>& v1, const std::vector<FP>& v2, const FP maxAbsError) {
        if (v1.size() != v2.size()) return false;
        for (int i = 0; i < v1.size(); i++)
            if (v1[i] - v2[i] < -maxAbsError || v1[i] - v2[i] >= maxAbsError) return false;
        return true;
    }

    bool comparePmlSplitGrids(FP maxAbsError) {
        auto splitGrid1 = this->solver->pml->splitGrid.get();
        auto splitGrid2 = this->solver2->pml->splitGrid.get();

        bool res = true;
        res = res && (splitGrid1->index == splitGrid2->index);
        res = res && compareFPVectors(splitGrid1->bxy, splitGrid2->bxy, maxAbsError);
        res = res && compareFPVectors(splitGrid1->bxz, splitGrid2->bxz, maxAbsError);
        res = res && compareFPVectors(splitGrid1->byx, splitGrid2->byx, maxAbsError);
        res = res && compareFPVectors(splitGrid1->byz, splitGrid2->byz, maxAbsError);
        res = res && compareFPVectors(splitGrid1->bzx, splitGrid2->bzx, maxAbsError);
        res = res && compareFPVectors(splitGrid1->bzy, splitGrid2->bzy, maxAbsError);
        res = res && compareFPVectors(splitGrid1->exy, splitGrid2->exy, maxAbsError);
        res = res && compareFPVectors(splitGrid1->exz, splitGrid2->exz, maxAbsError);
        res = res && compareFPVectors(splitGrid1->eyx, splitGrid2->eyx, maxAbsError);
        res = res && compareFPVectors(splitGrid1->eyz, splitGrid2->eyz, maxAbsError);
        res = res && compareFPVectors(splitGrid1->ezx, splitGrid2->ezx, maxAbsError);
        res = res && compareFPVectors(splitGrid1->ezy, splitGrid2->ezy, maxAbsError);

        return res;
    }

    bool comparePmlBase() {
        auto pml1 = this->solver->pml.get();
        auto pml2 = this->solver2->pml.get();

        bool res = true;
        res = res && (pml1->dt == pml2->dt);
        res = res && (pml1->domainIndexBegin == pml2->domainIndexBegin);
        res = res && (pml1->domainIndexEnd == pml2->domainIndexEnd);
        res = res && (pml1->sizePML == pml2->sizePML);
        res = res && (pml1->leftPmlBorder == pml2->leftPmlBorder);
        res = res && (pml1->rightPmlBorder == pml2->rightPmlBorder);
        res = res && (pml1->leftPmlBorderCoord == pml2->leftPmlBorderCoord);
        res = res && (pml1->rightPmlBorderCoord == pml2->rightPmlBorderCoord);
        res = res && (pml1->leftGlobalBorderCoord == pml2->leftGlobalBorderCoord);
        res = res && (pml1->rightGlobalBorderCoord == pml2->rightGlobalBorderCoord);

        return res;
    }

    bool comparePmlDerived(FP maxAbsError);

    bool comparePml(FP maxAbsError) {
        if (!this->solver->pml) return false;
        if (!this->solver2->pml) return false;
        return comparePmlBase() && comparePmlSplitGrids(maxAbsError) && comparePmlDerived(maxAbsError);
    }

    bool compareFields(FP maxAbsError) {
        for (int i = 0; i < this->grid->numCells.x; ++i)
            for (int j = 0; j < this->grid->numCells.y; ++j)
                for (int k = 0; k < this->grid->numCells.z; ++k)
                {
                    FP3 expectedE(this->grid->Ex(i, j, k), this->grid->Ey(i, j, k), this->grid->Ez(i, j, k));
                    FP3 expectedB(this->grid->Bx(i, j, k), this->grid->By(i, j, k), this->grid->Bz(i, j, k));
                    FP3 actualE(this->grid2->Ex(i, j, k), this->grid2->Ey(i, j, k), this->grid2->Ez(i, j, k));
                    FP3 actualB(this->grid2->Bx(i, j, k), this->grid2->By(i, j, k), this->grid2->Bz(i, j, k));
                    if ((expectedE - actualE).norm() > maxAbsError) return false;
                    if ((expectedB - actualB).norm() > maxAbsError) return false;
                }
        return true;
    }

    bool compareGenerator() {
        auto gen1 = this->solver->generator.get();
        auto gen2 = this->solver2->generator.get();

        if (!gen1) return false;
        if (!gen2) return false;

        bool res = true;
        res = res && (gen1->dt == gen2->dt);
        res = res && (gen1->domainIndexBegin == gen2->domainIndexBegin);
        res = res && (gen1->domainIndexEnd == gen2->domainIndexEnd);
        res = res && (gen1->leftGeneratorIndex == gen2->leftGeneratorIndex);
        res = res && (gen1->rightGeneratorIndex == gen2->rightGeneratorIndex);
        res = res && (gen1->isLeftBorderEnabled == gen2->isLeftBorderEnabled);
        res = res && (gen1->isRightBorderEnabled == gen2->isRightBorderEnabled);
        // TODO: compare generator functions

        return res;
    }

    bool compareBC() {
        PeriodicalBCType* bc1[3] = {
            dynamic_cast<PeriodicalBCType*>(this->solver->boundaryConditions[0].get()),
            dynamic_cast<PeriodicalBCType*>(this->solver->boundaryConditions[1].get()),
            dynamic_cast<PeriodicalBCType*>(this->solver->boundaryConditions[2].get())
        };

        PeriodicalBCType* bc2[3] = {
            dynamic_cast<PeriodicalBCType*>(this->solver2->boundaryConditions[0].get()),
            dynamic_cast<PeriodicalBCType*>(this->solver2->boundaryConditions[1].get()),
            dynamic_cast<PeriodicalBCType*>(this->solver2->boundaryConditions[2].get())
        };

        if (!bc1[0] || !bc1[1] || !bc1[2]) return false;
        if (!bc2[0] || !bc2[1] || !bc2[2]) return false;

        bool res = true;
        for (int d = 0; d < 3; d++) {
            res = res && (bc1[d]->axis == bc2[d]->axis);
            res = res && (bc1[d]->leftBorderIndex == bc2[d]->leftBorderIndex);
            res = res && (bc1[d]->rightBorderIndex == bc2[d]->rightBorderIndex);
        }

        return res;
    }
};

// compare pml for spectral solvers
template <class TSolver>
bool SaveLoadSolverModulesTest<TSolver>::comparePmlDerived(FP maxAbsError) {
    auto pml1 = this->solver->pml.get();
    auto pml2 = this->solver2->pml.get();

    bool res = true;
    res = res && (pml1->complexDomainIndexBegin == pml2->complexDomainIndexBegin);
    res = res && (pml1->complexDomainIndexEnd == pml2->complexDomainIndexEnd);
    res = res && compareFPVectors(pml1->bCoeffX, pml2->bCoeffX, maxAbsError);
    res = res && compareFPVectors(pml1->bCoeffY, pml2->bCoeffY, maxAbsError);
    res = res && compareFPVectors(pml1->bCoeffZ, pml2->bCoeffZ, maxAbsError);
    res = res && compareFPVectors(pml1->eCoeffX, pml2->eCoeffX, maxAbsError);
    res = res && compareFPVectors(pml1->eCoeffY, pml2->eCoeffY, maxAbsError);
    res = res && compareFPVectors(pml1->eCoeffZ, pml2->eCoeffZ, maxAbsError);

    return res;
}

// compare pml for FDTD
template <>
bool SaveLoadSolverModulesTest<FDTD>::comparePmlDerived(FP maxAbsError) {
    auto pml1 = this->solver->pml.get();
    auto pml2 = this->solver2->pml.get();

    bool res = true;
    res = res && compareFPVectors(pml1->bCoeff1X, pml2->bCoeff1X, maxAbsError);
    res = res && compareFPVectors(pml1->bCoeff1Y, pml2->bCoeff1Y, maxAbsError);
    res = res && compareFPVectors(pml1->bCoeff1Z, pml2->bCoeff1Z, maxAbsError);
    res = res && compareFPVectors(pml1->eCoeff1X, pml2->eCoeff1X, maxAbsError);
    res = res && compareFPVectors(pml1->eCoeff1Y, pml2->eCoeff1Y, maxAbsError);
    res = res && compareFPVectors(pml1->eCoeff1Z, pml2->eCoeff1Z, maxAbsError);
    res = res && compareFPVectors(pml1->bCoeff2X, pml2->bCoeff2X, maxAbsError);
    res = res && compareFPVectors(pml1->bCoeff2Y, pml2->bCoeff2Y, maxAbsError);
    res = res && compareFPVectors(pml1->bCoeff2Z, pml2->bCoeff2Z, maxAbsError);
    res = res && compareFPVectors(pml1->eCoeff2X, pml2->eCoeff2X, maxAbsError);
    res = res && compareFPVectors(pml1->eCoeff2Y, pml2->eCoeff2Y, maxAbsError);
    res = res && compareFPVectors(pml1->eCoeff2Z, pml2->eCoeff2Z, maxAbsError);

    return res;
}

#ifndef __USE_FFT__

typedef ::testing::Types<
    FDTD
> typesSaveLoadSolverModulesTest;

#else

typedef ::testing::Types<
    FDTD,
    PSTD,
    PSATD,
    PSATDTimeStaggered,
    PSATDPoisson,
    PSATDTimeStaggeredPoisson
> typesSaveLoadSolverModulesTest;

#endif

TYPED_TEST_CASE(SaveLoadSolverModulesTest, typesSaveLoadSolverModulesTest);

TYPED_TEST(SaveLoadSolverModulesTest, can_save_and_load_grid_solver_and_modules)
{
    for (int step = 0; step < this->numSteps; ++step)
        this->solver->updateFields();

    std::stringstream sstr;
    this->saveGridAndSolver(sstr);
    this->loadGridAndSolver(sstr);

    ASSERT_EQ(this->solver->dt, this->solver2->dt);
    ASSERT_TRUE(this->compareFields(std::numeric_limits<FP>::epsilon()));
    ASSERT_TRUE(this->comparePml(std::numeric_limits<FP>::epsilon()));
    ASSERT_TRUE(this->compareGenerator());
    ASSERT_TRUE(this->compareBC());
}

TYPED_TEST(SaveLoadSolverModulesTest, old_and_new_grid_solvers_and_modules_are_independent)
{
    for (int step = 0; step < this->numSteps / 2; ++step)
        this->solver->updateFields();

    std::stringstream sstr;
    this->saveGridAndSolver(sstr);

    for (int step = 0; step < this->numSteps / 2; ++step)
        this->solver->updateFields();

    this->loadGridAndSolver(sstr);

    ASSERT_FALSE(this->compareFields(this->maxAbsoluteError));

    for (int step = 0; step < this->numSteps / 2; ++step)
        this->solver2->updateFields();

    ASSERT_EQ(this->solver->dt, this->solver2->dt);
    ASSERT_TRUE(this->compareFields(this->maxAbsoluteError));
    ASSERT_TRUE(this->comparePml(this->maxAbsoluteError));
    ASSERT_TRUE(this->compareGenerator());
    ASSERT_TRUE(this->compareBC());
}


class SaveLoadSolverAnalyticalFieldSolver : public BaseFixture {
public:

    std::unique_ptr<AnalyticalField> field, field2;
    std::unique_ptr<AnalyticalFieldSolver> solver, solver2;

    std::function<FP(FP, FP, FP, FP)> funcField;
    std::function<FP(FP, FP, FP, FP)> zeroFunc;

    const int numSteps = 10;
    FP timeStep = 0.0;
    FP3 minCoords, maxCoords;

    virtual void SetUp() {
        this->maxAbsoluteError = (FP)1e-8;
        this->maxRelativeError = (FP)1e-4;
        this->timeStep = 1e-13;
        this->minCoords = FP3(-1.0, 0.0, 0.0);
        this->maxCoords = FP3(1.0, 1.0, 1.0);

        this->funcField = std::function<FP(FP, FP, FP, FP)>(
            [this](FP x, FP y, FP z, FP t)->FP {
                // TODO: save/load function
                //return sin(2 * constants::pi / (this->maxCoords.x - this->minCoords.x)
                //    * (x - constants::c * t));
                // temporary solution
                return (FP)0.0;
            });

        this->zeroFunc = [this](FP x, FP y, FP z, FP t)->FP {
            return (FP)0.0;
        };

        this->field.reset(new AnalyticalField(this->zeroFunc, this->funcField, this->zeroFunc,
            this->zeroFunc, this->zeroFunc, this->funcField));  // Ey, Bz
        this->solver.reset(new AnalyticalFieldSolver(this->field.get(), this->timeStep));
    }

    void saveGridAndSolver(std::ostream& ostr) {
        this->field->save(ostr);
        this->solver->save(ostr);
    }

    void loadGridAndSolver(std::istream& istr) {
        this->field2.reset(new AnalyticalField());
        this->field2->load(istr);
        this->solver2.reset(new AnalyticalFieldSolver(this->field2.get()));
        this->solver2->load(istr);
    }

    bool compareFields(FP maxAbsError) {
        Int3 gridSize(11, 12, 7);
        FP3 gridStep = (this->maxCoords - this->minCoords) / (FP3)gridSize;

        for (int i = 0; i < gridSize.x; i++)
            for (int j = 0; j < gridSize.y; j++)
                for (int k = 0; k < gridSize.z; k++) {
                    FP3 coord = this->minCoords + gridStep * FP3(i, j, k);
                    FP x = coord.x, y = coord.y, z = coord.z;
                    FP t = this->solver->getTime();

                    FP3 E(this->field->getE(coord));
                    FP3 B(this->field->getB(coord));
                    FP3 E2(this->field2->getE(coord));
                    FP3 B2(this->field2->getB(coord));
                    FP3 expE(this->zeroFunc(x, y, z, t), this->funcField(x, y, z, t), this->zeroFunc(x, y, z, t));
                    FP3 expB(this->zeroFunc(x, y, z, t), this->zeroFunc(x, y, z, t), this->funcField(x, y, z, t));

                    if ((expE - E).norm() > maxAbsError) return false;
                    if ((expB - B).norm() > maxAbsError) return false;
                    if ((expE - E2).norm() > maxAbsError) return false;
                    if ((expB - B2).norm() > maxAbsError) return false;
                    if ((E - E2).norm() > maxAbsError) return false;
                    if ((B - B2).norm() > maxAbsError) return false;
                }

        return true;
    }
};

TEST_F(SaveLoadSolverAnalyticalFieldSolver, can_save_and_load_analytical_field_solver)
{
    for (int step = 0; step < this->numSteps; ++step)
        this->solver->updateFields();

    std::stringstream sstr;
    this->saveGridAndSolver(sstr);
    this->loadGridAndSolver(sstr);

    ASSERT_NEAR_FP(this->field->globalTime, this->field2->globalTime);
    ASSERT_NEAR_FP(this->solver->dt, this->solver2->dt);
    ASSERT_TRUE(this->compareFields(std::numeric_limits<FP>::epsilon()));
}

TEST_F(SaveLoadSolverAnalyticalFieldSolver, old_and_new_analytical_field_solvers_are_independent)
{
    for (int step = 0; step < this->numSteps / 2; ++step)
        this->solver->updateFields();

    std::stringstream sstr;
    this->saveGridAndSolver(sstr);

    for (int step = 0; step < this->numSteps / 2; ++step)
        this->solver->updateFields();

    this->loadGridAndSolver(sstr);

    // TODO: uncomment when function save/load is implemented
    //ASSERT_FALSE(this->compareFields(std::numeric_limits<FP>::epsilon()));

    for (int step = 0; step < this->numSteps / 2; ++step)
        this->solver2->updateFields();

    ASSERT_NEAR_FP(this->field->globalTime, this->field2->globalTime);
    ASSERT_NEAR_FP(this->solver->dt, this->solver2->dt);
    ASSERT_TRUE(this->compareFields(std::numeric_limits<FP>::epsilon()));
}
