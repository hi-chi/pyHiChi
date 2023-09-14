#include "TestingUtility.h"

#include "FourierTransform.h"
#include "Pstd.h"

class FourierTransformTest : public BaseFixture {
public:

    const int nx = 10, ny = 7, nz = 8;
    const FP a = 0.0, b = 10.0;
    const FP dx = (b - a) / nx, dy = (b - a) / ny, dz = (b - a) / nz;
    const int stepX = 1, stepY = 3, stepZ = 6;
    const FP shiftX = stepX * dx, shiftY = stepY * dy, shiftZ = stepZ * dz;

    Int3 size;
    Int3 sizeComplex;

    ScalarField<FP> field;
    ScalarField<complexFP> complexField;

    ArrayFourierTransform3d fourierTransform;

    FourierTransformTest() :
        size(nx, ny, nz),
        sizeComplex(fourier_transform::getSizeOfComplexArray(Int3(nx, ny, nz))),
        field(size),
        complexField(sizeComplex)
    {
        fourierTransform.initialize(field.getData(), complexField.getData(), size);

        setField(&FourierTransformTest::fSin3);

        maxAbsoluteError = (FP)1e-7;
        maxRelativeError = (FP)1e-4;
    }

    ~FourierTransformTest() {
    }

    FP fSin3(int i, int j, int k) {
        const double A = 1000;
        return A * sin(2 * constants::pi * (i * dx + j * dy + k * dz) / (b - a));
    }

    FP fSin(int i, int j, int k) {
        return sin(2 * constants::pi * i * dx / (b - a));
    }

    void setField(FP(FourierTransformTest::* func)(int, int, int)) {
        for (int i = 0; i < size.x; i++)
            for (int j = 0; j < size.y; j++)
                for (int k = 0; k < size.z; k++)
                    field(i, j, k) = (this->*func)(i, j, k);
    }

    FP3 getFrequencyVector(const Int3& index) {
        FP kx(2 * constants::pi * ((index.x <= size.x / 2) ? index.x : index.x - size.x) / (b - a));
        FP ky(2 * constants::pi * ((index.y <= size.y / 2) ? index.y : index.y - size.y) / (b - a));
        FP kz(2 * constants::pi * ((index.z <= size.z / 2) ? index.z : index.z - size.z) / (b - a));
        return FP3(kx, ky, kz);
    }

    void doShift() {
        for (int i = 0; i < sizeComplex.x; i++)
            for (int j = 0; j < sizeComplex.y; j++)
                for (int k = 0; k < sizeComplex.z; k++) {
                    FP3 freq = getFrequencyVector(Int3(i, j, k));
                    complexField(i, j, k) *= complexFP::createInTrig(1.0,
                        -shiftX * freq.x - shiftY * freq.y - shiftZ * freq.z);
                }
    }

};

TEST_F(FourierTransformTest, DirectAndInverseTransform) {

#ifdef __USE_FFT__

    fourierTransform.doFourierTransform(fourier_transform::Direction::RtoC);
    fourierTransform.doFourierTransform(fourier_transform::Direction::CtoR);

    for (int i = 0; i < size.x; i++)
        for (int j = 0; j < size.y; j++)
            for (int k = 0; k < size.z; k++)
                ASSERT_NEAR_FP(fSin3(i, j, k), field(i, j, k));

    // shows that some tests were skipped
#else
    GTEST_SKIP();
#endif
}

TEST_F(FourierTransformTest, TransformSinus) {

#ifdef __USE_FFT__

    setField(&FourierTransformTest::fSin);

    fourierTransform.doFourierTransform(fourier_transform::Direction::RtoC);

    for (int i = 0; i < sizeComplex.x; i++)
        for (int j = 0; j < sizeComplex.y; j++)
            for (int k = 0; k < sizeComplex.z; k++)
                if (!(i == 1 && j == 0 && k == 0 || i == sizeComplex.x - 1 && j == 0 && k == 0))
                    ASSERT_NEAR_MODULE_COMPLEXFP(complexFP(0), complexField(i, j, k));

    // shows that some tests were skipped
#else
    GTEST_SKIP();
#endif
}

TEST_F(FourierTransformTest, Lag) {

#ifdef __USE_FFT__

    fourierTransform.doFourierTransform(fourier_transform::Direction::RtoC);
    doShift();
    fourierTransform.doFourierTransform(fourier_transform::Direction::CtoR);

    for (int i = 0; i < size.x; i++)
        for (int j = 0; j < size.y; j++)
            for (int k = 0; k < size.z; k++)
                ASSERT_NEAR_FP(fSin3(i - stepX, j - stepY, k - stepZ), field(i, j, k));

    // shows that some tests were skipped
#else
    GTEST_SKIP();
#endif
}



class FourierTransformSolverTest : public BaseFixture {
public:

    const int nx = 10, ny = 7, nz = 8;
    const FP a = 0.0, b = 10.0;
    const FP dx = (b - a) / nx, dy = (b - a) / ny, dz = (b - a) / nz;
    const int stepX = 1, stepY = 3, stepZ = 6;
    const FP shiftX = stepX * dx, shiftY = stepY * dy, shiftZ = stepZ * dz;
    const FP dt = 1.0 / (constants::c * 4.0);

    Int3 size;
    Int3 sizeComplex;
    Int3 minCoord;
    FP3 d;

    std::unique_ptr<PSTD> pstd;

    std::unique_ptr<Grid<FP, GridTypes::PSTDGridType>> grid;
    SpectralGrid<FP, complexFP>* complexGrid = 0;

    void SetUp() {
        size = Int3(nx, ny, nz);
        minCoord = Int3(a, a, a);
        d = FP3(dx, dy, dz);

        grid.reset(new Grid<FP, GridTypes::PSTDGridType>(size, minCoord, d, size));
        pstd.reset(new PSTD(grid.get(), dt));

        sizeComplex = fourier_transform::getSizeOfComplexArray(grid->numCells);
        complexGrid = pstd->complexGrid.get();

        setBz(&FourierTransformSolverTest::fSin3);

        maxAbsoluteError = (FP)1e-7;
        maxRelativeError = (FP)1e-4;
    }

    FP fSin3(int i, int j, int k) {
        const double A = 1000;
        return A * sin(2 * constants::pi * (i * dx + j * dy + k * dz) / (b - a));
    }

    FP fSin(int i, int j, int k) {
        return sin(2 * constants::pi * i * dx / (b - a));
    }

    void setBz(FP(FourierTransformSolverTest::* func)(int, int, int)) {
        for (int i = 0; i < grid->numCells.x; i++)
            for (int j = 0; j < grid->numCells.y; j++)
                for (int k = 0; k < grid->numCells.z; k++)
                    grid->Bz(i, j, k) = (this->*func)(i, j, k);
    }

    FP3 getFrequencyVector(const Int3& index) {
        FP kx(2 * constants::pi * ((index.x <= size.x / 2) ? index.x : index.x - size.x) / (b - a));
        FP ky(2 * constants::pi * ((index.y <= size.y / 2) ? index.y : index.y - size.y) / (b - a));
        FP kz(2 * constants::pi * ((index.z <= size.z / 2) ? index.z : index.z - size.z) / (b - a));
        return FP3(kx, ky, kz);
    }

    void doShift() {
        for (int i = 0; i < complexGrid->numCells.x; i++)
            for (int j = 0; j < complexGrid->numCells.y; j++)
                for (int k = 0; k < complexGrid->numCells.z; k++) {
                    FP3 freq = getFrequencyVector(Int3(i, j, k));
                    complexGrid->Bz(i, j, k) *= complexFP::createInTrig(1.0,
                        -shiftX * freq.x - shiftY * freq.y - shiftZ * freq.z);
                }
    }

};


TEST_F(FourierTransformSolverTest, DirectAndInverseTransform) {

#ifdef __USE_FFT__

    pstd->fourierTransform.doFourierTransform(FieldEnum::B, CoordinateEnum::z, fourier_transform::Direction::RtoC);
    pstd->fourierTransform.doFourierTransform(FieldEnum::B, CoordinateEnum::z, fourier_transform::Direction::CtoR);

    for (int i = 0; i < grid->numCells.x; i++)
        for (int j = 0; j < grid->numCells.y; j++)
            for (int k = 0; k < grid->numCells.z; k++)
                ASSERT_NEAR_FP(fSin3(i, j, k), grid->Bz(i, j, k));

    // shows that some tests were skipped
#else
    GTEST_SKIP();
#endif
}

TEST_F(FourierTransformSolverTest, TransformSinus) {

#ifdef __USE_FFT__

    setBz(&FourierTransformSolverTest::fSin);

    pstd->fourierTransform.doFourierTransform(FieldEnum::B, CoordinateEnum::z, fourier_transform::Direction::RtoC);

    for (int i = 0; i < complexGrid->numCells.x; i++)
        for (int j = 0; j < complexGrid->numCells.y; j++)
            for (int k = 0; k < complexGrid->numCells.z; k++)
                if (!(i == 1 && j == 0 && k == 0 || i == complexGrid->numCells.x - 1 && j == 0 && k == 0))
                    ASSERT_NEAR_MODULE_COMPLEXFP(complexFP(0), complexGrid->Bz(i, j, k));

    // shows that some tests were skipped
#else
    GTEST_SKIP();
#endif
}

TEST_F(FourierTransformSolverTest, Lag) {

#ifdef __USE_FFT__

    pstd->fourierTransform.doFourierTransform(FieldEnum::B, CoordinateEnum::z, fourier_transform::Direction::RtoC);
    doShift();
    pstd->fourierTransform.doFourierTransform(FieldEnum::B, CoordinateEnum::z, fourier_transform::Direction::CtoR);

    for (int i = 0; i < grid->numCells.x; i++)
        for (int j = 0; j < grid->numCells.y; j++)
            for (int k = 0; k < grid->numCells.z; k++)
                ASSERT_NEAR_FP(fSin3(i - stepX, j - stepY, k - stepZ), grid->Bz(i, j, k));

    // shows that some tests were skipped
#else
    GTEST_SKIP();
#endif
}
