#include "TestingUtility.h"

#include "Fdtd.h"

template <class TTypeDefinitionsFieldTest>
class FieldGeneratorTest : public BaseFixture {
public:
    using FieldSolverType = typename TTypeDefinitionsFieldTest::FieldSolverType;
    using GridType = typename TTypeDefinitionsFieldTest::FieldSolverType::GridType;
    using BoundaryConditionType = typename FieldSolverType::PeriodicalBoundaryConditionType;

    const int dimension = TTypeDefinitionsFieldTest::dimension;
    const CoordinateEnum axis = TTypeDefinitionsFieldTest::axis;

    const int gridSizeLongitudinal = 24;
    const int gridSizeTransverse = 12;
    const int generatorIndex = 4;

    std::unique_ptr<FieldSolverType> fieldSolver;
    std::unique_ptr<GridType> grid;
    Int3 gridSize;
    FP3 gridStep;
    FP timeStep = 0;
    FP3 minCoords, maxCoords;
    Int3 generatorStartIndex, generatorEndIndex;
    int numSteps = 0;

    FP relatedAmpThreshold = 0.03;

    virtual void SetUp() {
        this->maxAbsoluteError = 1e-4;
        this->maxRelativeError = 0.5;

        gridSize = Int3(1, 1, 1);
        for (int d = 0; d < dimension; d++) {
            gridSize[d] = gridSizeTransverse;
        }
        gridSize[(int)axis] = gridSizeLongitudinal;

        this->minCoords = FP3(0, 0, 0);
        this->maxCoords = (FP3)gridSize * constants::c;
        this->gridStep = (this->maxCoords - this->minCoords) / (FP3)this->gridSize;

        this->grid.reset(new GridType(this->gridSize, this->minCoords, this->gridStep, this->gridSize));

        this->timeStep = 0.5 * FieldSolverType::getCourantConditionTimeStep(this->gridStep);

        fieldSolver.reset(new FieldSolverType(this->grid.get(), this->timeStep));
        fieldSolver->template setBoundaryCondition<BoundaryConditionType>();
           
        this->numSteps = 2 * (int)((this->maxCoords - this->minCoords)[(int)axis] /
            (constants::c * fieldSolver->dt));

        initializeGenerator();
    }

    void initializeGenerator() {
        generatorStartIndex = Int3(1, 1, 1);
        generatorEndIndex = Int3(1, 1, 1);
        for (int d = 0; d < dimension; d++) {
            generatorStartIndex[d] = generatorIndex;
            generatorEndIndex[d] = gridSize[d] - generatorIndex;
        }

        int axis0 = (int)this->axis;
        int axis1 = ((int)axis + 1) % 3;
        int axis2 = ((int)axis + 2) % 3;

        FP L = (this->maxCoords - this->minCoords)[(int)axis];

        auto fieldFunc = [axis0, L](FP x, FP y, FP z, FP t) {
            FP coord = (FP3(x, y, z)[axis0] - constants::c * t) / L;
            if (coord < -1.0 || coord >= 0.0) coord = 0.0;
            const FP pi2 = 2.0 * constants::pi;
            // H = harris function
            FP H = 0.03125 * (10.0 - 15.0 * cos(pi2 * coord) + 6.0 * cos(2.0 * pi2 * coord) - cos(3.0 * pi2 * coord));
            return H;
        };
        
        auto bFunc = [axis2, fieldFunc](FP x, FP y, FP z, FP t) {
            FP3 b;
            b[axis2] = fieldFunc(x, y, z, t);
            return b;
        };

        auto eFunc = [axis1, fieldFunc](FP x, FP y, FP z, FP t) {
            FP3 e;
            e[axis1] = fieldFunc(x, y, z, t);
            return e;
        };

        auto bxFunc = [bFunc](FP x, FP y, FP z, FP t) { return bFunc(x, y, z, t).x; };
        auto byFunc = [bFunc](FP x, FP y, FP z, FP t) { return bFunc(x, y, z, t).y; };
        auto bzFunc = [bFunc](FP x, FP y, FP z, FP t) { return bFunc(x, y, z, t).z; };

        auto exFunc = [eFunc](FP x, FP y, FP z, FP t) { return eFunc(x, y, z, t).x; };
        auto eyFunc = [eFunc](FP x, FP y, FP z, FP t) { return eFunc(x, y, z, t).y; };
        auto ezFunc = [eFunc](FP x, FP y, FP z, FP t) { return eFunc(x, y, z, t).z; };

        fieldSolver->setFieldGenerator(generatorStartIndex, generatorEndIndex,
            bxFunc, byFunc, bzFunc, exFunc, eyFunc, ezFunc);
    }

    FP computeAmp() {
        FP amp = 0;
        for (int i = 0; i < grid->numCells.x; i++)
            for (int j = 0; j < grid->numCells.y; j++)
                for (int k = 0; k < grid->numCells.z; k++) {
                    FP dEnergy = pow(grid->Ex(i, j, k), 2) + pow(grid->Ey(i, j, k), 2) + pow(grid->Ez(i, j, k), 2) +
                        pow(grid->Bx(i, j, k), 2) + pow(grid->By(i, j, k), 2) + pow(grid->Bz(i, j, k), 2);
                    FP curAmp = sqrt(dEnergy * 0.5);
                    amp = curAmp > amp ? curAmp : amp;
                }
        return amp;
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
TYPED_TEST_CASE(FieldGeneratorTest, types);


TYPED_TEST(FieldGeneratorTest, FieldGeneratorTestPlaneWave) {
    // to disable testing of spectral solvers without enabled fftw
#ifndef __USE_FFT__
    SUCCEED();
#else

    for (int step = 0; step < this->numSteps; ++step) {
        this->fieldSolver->updateFields();
    }

    FP finalAmp = this->computeAmp();

    ASSERT_NEAR(finalAmp, 0, this->relatedAmpThreshold);

#endif
}