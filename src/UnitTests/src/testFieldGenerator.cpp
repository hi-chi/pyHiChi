#include "TestingUtility.h"

#include "Fdtd.h"
#include "BoundaryConditionFdtd.h"
#include "FieldGeneratorFdtd.h"

template <class TTypeDefinitionsFieldTest>
class FieldGeneratorTest : public BaseFixture {
public:
    using FieldSolverType = typename TTypeDefinitionsFieldTest::FieldSolverType;
    using GridType = typename TTypeDefinitionsFieldTest::GridType;
    using BoundaryConditionType = typename FieldSolverType::PeriodicalBoundaryConditionType;

    const int dimension = TTypeDefinitionsFieldTest::dimension;
    const CoordinateEnum axis = TTypeDefinitionsFieldTest::axis;

    const int gridSizeLongitudinal = 32;
    const int gridSizeTransverse = 16;
    const int generatorIndex = 4;

    std::unique_ptr<FieldSolverType> fieldSolver;
    std::unique_ptr<GridType> grid;
    std::unique_ptr<BoundaryConditionType> boundaryCondition;
    Int3 gridSize;
    FP3 gridStep;
    FP timeStep = 0;
    FP3 minCoords, maxCoords;
    Int3 generatorStartIndex, generatorEndIndex;

    FP relatedAmpThreshold = 0.03;

    virtual void SetUp() {
        gridSize = Int3(1, 1, 1);
        for (int d = 0; d < dimension; d++) {
            gridSize[d] = gridSizeTransverse;
        }
        gridSize[(int)axis] = gridSizeLongitudinal;

        this->minCoords = FP3(0, 0, 0);
        FP b = gridSizeLongitudinal * constants::c;
        this->maxCoords = FP3(b, b, b);
        this->gridStep = (this->maxCoords - this->minCoords) / (FP3)this->gridSize;

        this->grid.reset(new GridType(this->gridSize, this->minCoords, this->gridStep, this->gridSize));

        // should satisfy the Courant's condition for all solvers
        this->timeStep = 0.4 * constants::c / grid->steps.norm();

        fieldSolver.reset(new FieldSolverType(this->grid.get(), this->timeStep));

        boundaryCondition.reset(new BoundaryConditionType(fieldSolver.get()));
        fieldSolver->setBoundaryCondition(boundaryCondition.get());

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
    TypeDefinitionsFieldTest<FDTD, YeeGrid, 1, CoordinateEnum::x>,
    TypeDefinitionsFieldTest<FDTD, YeeGrid, 2, CoordinateEnum::x>,
    TypeDefinitionsFieldTest<FDTD, YeeGrid, 2, CoordinateEnum::y>,
    TypeDefinitionsFieldTest<FDTD, YeeGrid, 3, CoordinateEnum::x>,
    TypeDefinitionsFieldTest<FDTD, YeeGrid, 3, CoordinateEnum::y>,
    TypeDefinitionsFieldTest<FDTD, YeeGrid, 3, CoordinateEnum::z>
> types;
TYPED_TEST_CASE(FieldGeneratorTest, types);


TYPED_TEST(FieldGeneratorTest, FieldGeneratorTestPlaneWave) {
    // to disable testing of spectral solvers without enabled fftw
#ifndef __USE_FFT__
    SUCCEED();
#else

    const int numSteps = 2 * (int)((this->maxCoords - this->minCoords)[(int)axis] /
        (constants::c * fieldSolver->dt));

    for (int step = 0; step < numSteps; ++step) {
        fieldSolver->updateFields();
    }

    FP finalAmp = computeAmp();

    ASSERT_NEAR(finalAmp, 0, relatedAmpThreshold);

#endif
}

