#include "TestingUtility.h"

#include "sycl/DeviceSYCL.h"
#include "sycl/FdtdSYCL.h"

using namespace pfc;
enum axis
{
    X = 0,
    Y = 1,
    Z = 2
};

template<axis a>
class ax
{};

typedef ax<X> axisX;
typedef ax<Y> axisY;
typedef ax<Z> axisZ;

template <class axis>
class GridFDTDSYCLTest : public BaseGridFixture<YeeGrid> {
public:
    virtual void SetUp() {
        BaseGridFixture<YeeGrid>::SetUp();

        fdtd = new sycl_pfc::FDTD(grid, this->timeStep);
    }

    sycl_pfc::FDTD* fdtd;

    ~GridFDTDSYCLTest() {
        //delete fdtd;
    }

    FP3 eTest(FP x, FP y, FP z, FP t);
    FP3 bTest(FP x, FP y, FP z, FP t);
};

template<>
FP3 GridFDTDSYCLTest<axisX>::eTest(FP x, FP y, FP z, FP t)
{
    return FP3(0, sin(2 * constants::pi * (-constants::c * t + x)), 0);
}

template<>
FP3 GridFDTDSYCLTest<axisX>::bTest(FP x, FP y, FP z, FP t)
{
    return FP3(0, 0, sin(2 * constants::pi * (-constants::c * t + x)));
}

template<>
FP3 GridFDTDSYCLTest<axisY>::eTest(FP x, FP y, FP z, FP t)
{
    return FP3(0, 0, sin(2 * constants::pi * (-constants::c * t + y)));
}

template<>
FP3 GridFDTDSYCLTest<axisY>::bTest(FP x, FP y, FP z, FP t)
{
    return FP3(sin(2 * constants::pi * (-constants::c * t + y)), 0, 0);
}

template<>
FP3 GridFDTDSYCLTest<axisZ>::eTest(FP x, FP y, FP z, FP t)
{
    return FP3(sin(2 * constants::pi * (-constants::c * t + z)), 0, 0);
}

template<>
FP3 GridFDTDSYCLTest<axisZ>::bTest(FP x, FP y, FP z, FP t)
{
    return FP3(0, sin(2 * constants::pi * (-constants::c * t + z)), 0);
}

typedef ::testing::Types<
    axisX,
    axisY,
    axisZ
> types;
TYPED_TEST_CASE(GridFDTDSYCLTest, types);

class GridFDTDSYCLTestZ : public GridFDTDSYCLTest<axisZ>
{};

TYPED_TEST(GridFDTDSYCLTest, FDTDSYCL_CPU_periodical)
//TEST_F(GridFDTDTestZ, FDTD_periodical)
{
    PeriodicalFieldGeneratorYee periodicalBC(this->fdtd);
    this->fdtd->setFieldGenerator(&periodicalBC);
    this->fdtd->changeDevice(sycl_pfc::Devices::CPU);

    int numCells = this->grid->numCells.volume();
    // e0(x, y, z) = eTest(x, y, z, 0), b0(x, y, z) = bTest(x, y, z, dt/2)
    for (int i = 0; i < this->grid->numCells.x; ++i)
        for (int j = 0; j < this->grid->numCells.y; ++j)
            for (int k = 0; k < this->grid->numCells.z; ++k)
            {
                FP3 coords = this->grid->ExPosition(i, j, k);
                this->grid->Ex(i, j, k) = this->eTest(coords.x, coords.y, coords.z, 0).x;
                coords = this->grid->EyPosition(i, j, k);
                this->grid->Ey(i, j, k) = this->eTest(coords.x, coords.y, coords.z, 0).y;
                coords = this->grid->EzPosition(i, j, k);
                this->grid->Ez(i, j, k) = this->eTest(coords.x, coords.y, coords.z, 0).z;
                coords = this->grid->BxPosition(i, j, k);
                this->grid->Bx(i, j, k) = this->bTest(coords.x, coords.y, coords.z, 0).x;
                coords = this->grid->ByPosition(i, j, k);
                this->grid->By(i, j, k) = this->bTest(coords.x, coords.y, coords.z, 0).y;
                coords = this->grid->BzPosition(i, j, k);
                this->grid->Bz(i, j, k) = this->bTest(coords.x, coords.y, coords.z, 0).z;
            }

    FP startT = 0;
    for (int i = 0; i < this->grid->numCells.x; ++i)
        for (int j = 0; j < this->grid->numCells.y; ++j)
            for (int k = 0; k < this->grid->numCells.z; ++k)
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
                ASSERT_NEAR_FP3(expectedE, actualE);
            }

    const int numSteps = 512;
    for (int step = 0; step < numSteps; ++step)
    {
        this->fdtd->updateFields();
    }

    FP finalT = this->fdtd->dt * numSteps;
    for (int i = 1; i < this->grid->numCells.x - 1; ++i)
        for (int j = 1; j < this->grid->numCells.y - 1; ++j)
            for (int k = 1; k < this->grid->numCells.z - 1; ++k)
            {
                FP3 expectedE, actualE;
                FP3 coords = this->grid->ExPosition(i, j, k);
                expectedE.x = this->eTest(coords.x, coords.y, coords.z, finalT).x;
                coords = this->grid->EyPosition(i, j, k);
                expectedE.y = this->eTest(coords.x, coords.y, coords.z, finalT).y;
                coords = this->grid->EzPosition(i, j, k);
                expectedE.z = this->eTest(coords.x, coords.y, coords.z, finalT).z;
                actualE.x = this->grid->Ex(i, j, k);
                actualE.y = this->grid->Ey(i, j, k);
                actualE.z = this->grid->Ez(i, j, k);
                ASSERT_NEAR_FP3(expectedE, actualE);
            }

    for (int i = 1; i < this->grid->numCells.x - 1; ++i)
        for (int j = 1; j < this->grid->numCells.y - 1; ++j)
            for (int k = 1; k < this->grid->numCells.z - 1; ++k)
            {
                FP3 expectedB, actualB;
                FP3 coords = this->grid->BxPosition(i, j, k);
                expectedB.x = this->bTest(coords.x, coords.y, coords.z, finalT).x;
                coords = this->grid->ByPosition(i, j, k);
                expectedB.y = this->bTest(coords.x, coords.y, coords.z, finalT).y;
                coords = this->grid->BzPosition(i, j, k);
                expectedB.z = this->bTest(coords.x, coords.y, coords.z, finalT).z;
                actualB.x = this->grid->Bx(i, j, k);
                actualB.y = this->grid->By(i, j, k);
                actualB.z = this->grid->Bz(i, j, k);
                ASSERT_NEAR_FP3(expectedB, actualB);
            }
}

TYPED_TEST(GridFDTDSYCLTest, FDTDSYCL_GPU_periodical)
//TEST_F(GridFDTDTestZ, FDTD_periodical)
{
    PeriodicalFieldGeneratorYee periodicalBC(this->fdtd);
    this->fdtd->setFieldGenerator(&periodicalBC);
    this->fdtd->changeDevice(sycl_pfc::Devices::GPU);

    int numCells = this->grid->numCells.volume();
    // e0(x, y, z) = eTest(x, y, z, 0), b0(x, y, z) = bTest(x, y, z, dt/2)
    for (int i = 0; i < this->grid->numCells.x; ++i)
        for (int j = 0; j < this->grid->numCells.y; ++j)
            for (int k = 0; k < this->grid->numCells.z; ++k)
            {
                FP3 coords = this->grid->ExPosition(i, j, k);
                this->grid->Ex(i, j, k) = this->eTest(coords.x, coords.y, coords.z, 0).x;
                coords = this->grid->EyPosition(i, j, k);
                this->grid->Ey(i, j, k) = this->eTest(coords.x, coords.y, coords.z, 0).y;
                coords = this->grid->EzPosition(i, j, k);
                this->grid->Ez(i, j, k) = this->eTest(coords.x, coords.y, coords.z, 0).z;
                coords = this->grid->BxPosition(i, j, k);
                this->grid->Bx(i, j, k) = this->bTest(coords.x, coords.y, coords.z, 0).x;
                coords = this->grid->ByPosition(i, j, k);
                this->grid->By(i, j, k) = this->bTest(coords.x, coords.y, coords.z, 0).y;
                coords = this->grid->BzPosition(i, j, k);
                this->grid->Bz(i, j, k) = this->bTest(coords.x, coords.y, coords.z, 0).z;
            }

    FP startT = 0;
    for (int i = 0; i < this->grid->numCells.x; ++i)
        for (int j = 0; j < this->grid->numCells.y; ++j)
            for (int k = 0; k < this->grid->numCells.z; ++k)
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
                ASSERT_NEAR_FP3(expectedE, actualE);
            }

    const int numSteps = 512;
    for (int step = 0; step < numSteps; ++step)
    {
        this->fdtd->updateFields();
    }

    FP finalT = this->fdtd->dt * numSteps;
    for (int i = 1; i < this->grid->numCells.x - 1; ++i)
        for (int j = 1; j < this->grid->numCells.y - 1; ++j)
            for (int k = 1; k < this->grid->numCells.z - 1; ++k)
            {
                FP3 expectedE, actualE;
                FP3 coords = this->grid->ExPosition(i, j, k);
                expectedE.x = this->eTest(coords.x, coords.y, coords.z, finalT).x;
                coords = this->grid->EyPosition(i, j, k);
                expectedE.y = this->eTest(coords.x, coords.y, coords.z, finalT).y;
                coords = this->grid->EzPosition(i, j, k);
                expectedE.z = this->eTest(coords.x, coords.y, coords.z, finalT).z;
                actualE.x = this->grid->Ex(i, j, k);
                actualE.y = this->grid->Ey(i, j, k);
                actualE.z = this->grid->Ez(i, j, k);
                ASSERT_NEAR_FP3(expectedE, actualE);
            }

    for (int i = 1; i < this->grid->numCells.x - 1; ++i)
        for (int j = 1; j < this->grid->numCells.y - 1; ++j)
            for (int k = 1; k < this->grid->numCells.z - 1; ++k)
            {
                FP3 expectedB, actualB;
                FP3 coords = this->grid->BxPosition(i, j, k);
                expectedB.x = this->bTest(coords.x, coords.y, coords.z, finalT).x;
                coords = this->grid->ByPosition(i, j, k);
                expectedB.y = this->bTest(coords.x, coords.y, coords.z, finalT).y;
                coords = this->grid->BzPosition(i, j, k);
                expectedB.z = this->bTest(coords.x, coords.y, coords.z, finalT).z;
                actualB.x = this->grid->Bx(i, j, k);
                actualB.y = this->grid->By(i, j, k);
                actualB.z = this->grid->Bz(i, j, k);
                ASSERT_NEAR_FP3(expectedB, actualB);
            }
}