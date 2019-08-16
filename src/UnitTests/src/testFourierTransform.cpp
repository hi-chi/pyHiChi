#include "TestingUtility.h"

#include "FourierTransform.h"

class FourierTransformTest : public BaseFixture {
public:

    const int nx = 10, ny = 5, nz = 8;
    const FP a = 0.0, b = 10.0;
    const FP dx = (b - a) / nx, dy = (b - a) / ny, dz = (b - a) / nz;
    const int stepX = 1, stepY = 3, stepZ = 6;
    const FP shiftX = stepX * dx, shiftY = stepY * dy, shiftZ = stepZ * dz;

    Int3 size;
    Int3 sizeComplex;

    ScalarField<FP> field;
    ScalarField<complexFP> complexField;

    FourierTransformTest() :
        size(nx, ny, nz),
        sizeComplex(FourierTransform::getSizeOfComplex(Int3(nx, ny, nz))),
        field(size),
        complexField(sizeComplex)
    {
        setField(&FourierTransformTest::fSin3);

        maxAbsoluteError = (FP)1e-7;
        maxRelativeError = (FP)1e-4;
    }

    ~FourierTransformTest() {
    }

    FP fSin3(int i, int j, int k) {
        const double A = 1000;
        return A * sin(2 * constants::pi*(i * dx + j * dy + k * dz) / (b - a));
    }

    FP fSin(int i, int j, int k) {
        return sin(2 * constants::pi * i * dx / (b - a));
    }

    void setField(FP (FourierTransformTest::* func)(int, int, int)) {
        for (int i = 0; i < size.x; i++)
            for (int j = 0; j < size.y; j++)
                for (int k = 0; k < size.z; k++)
                    field(i, j, k) = (this->*func)(i, j, k);
    }

    FP3 getFrequencyVector(const Int3 & index) {
        FP kx(2 * constants::pi*((index.x <= size.x / 2) ? index.x : index.x - size.x) / (b - a));
        FP ky(2 * constants::pi*((index.y <= size.y / 2) ? index.y : index.y - size.y) / (b - a));
        FP kz(2 * constants::pi*((index.z <= size.z / 2) ? index.z : index.z - size.z) / (b - a));
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

TEST_F(FourierTransformTest, ADD_TEST_FFT_PREFIX(DirectAndInverseTransform)) {

    FourierTransform::doFourierTransform(field, complexField, FourierTransformDirection::RtoC);
    FourierTransform::doFourierTransform(field, complexField, FourierTransformDirection::CtoR);

    for (int i = 0; i < size.x; i++)
        for (int j = 0; j < size.y; j++)
            for (int k = 0; k < size.z; k++)
                ASSERT_NEAR_FP(fSin3(i, j, k), field(i, j, k));
}

TEST_F(FourierTransformTest, ADD_TEST_FFT_PREFIX(TransformSinus)) {

    setField(&FourierTransformTest::fSin);

    FourierTransform::doFourierTransform(field, complexField, FourierTransformDirection::RtoC);

    for (int i = 0; i < sizeComplex.x; i++)
        for (int j = 0; j < sizeComplex.y; j++)
            for (int k = 0; k < sizeComplex.z; k++)
                if (!(i == 1 && j == 0 && k == 0 || i == sizeComplex.x - 1 && j == 0 && k == 0))
                    ASSERT_NEAR_MODULE_COMPLEXFP(complexFP(0), complexField(i, j, k));
}

TEST_F(FourierTransformTest, ADD_TEST_FFT_PREFIX(Lag)) {

    FourierTransform::doFourierTransform(field, complexField, FourierTransformDirection::RtoC);
    doShift();
    FourierTransform::doFourierTransform(field, complexField, FourierTransformDirection::CtoR);

    for (int i = 0; i < size.x; i++)
        for (int j = 0; j < size.y; j++)
            for (int k = 0; k < size.z; k++)
                ASSERT_NEAR_FP(fSin3(i - stepX, j - stepY, k - stepZ), field(i, j, k));
}