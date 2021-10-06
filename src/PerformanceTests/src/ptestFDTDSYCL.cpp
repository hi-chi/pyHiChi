#include "TestingUtility.h"

#include "sycl/DeviceSYCL.h"
#include "sycl/FdtdSYCL.h"

using namespace pfc;

template <class axis>
class GridFDTDSYCLTest : public BaseGridFixture<YeeGrid> {
public:
    virtual void SetUp(const ::benchmark::State& st)
    {
        BaseGridFixture<YeeGrid>::SetUp(st);

        fdtd = new sycl_pfc::FDTD(grid, this->timeStep);
    }

    sycl_pfc::FDTD* fdtd;

    ~GridFDTDSYCLTest() {
        delete fdtd;
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

static void CustomArguments(benchmark::internal::Benchmark* b) {
    b->Args({ 16 });
    b->Iterations(5);
}

using GridFDTDSYCLTestZ = GridFDTDSYCLTest<axisZ>;
BENCHMARK_DEFINE_F(GridFDTDSYCLTestZ, FDTDSYCL_CPU_periodical_vtune)(benchmark::State& state) {
    PeriodicalFieldGeneratorYee periodicalBC(this->fdtd);
    this->fdtd->setFieldGenerator(&periodicalBC);
    const int numSteps = state.range_x();
    this->fdtd->changeDevice(sycl_pfc::Devices::CPU);

    while (state.KeepRunning()) {
        for (int step = 0; step < numSteps; ++step)
        {
            this->fdtd->updateFields();
        }
    }
}
BENCHMARK_REGISTER_F(GridFDTDSYCLTestZ, FDTDSYCL_CPU_periodical_vtune)->Apply(CustomArguments)->Unit(benchmark::kSecond);

BENCHMARK_DEFINE_F(GridFDTDSYCLTestZ, FDTDSYCL_GPU_periodical_vtune)(benchmark::State& state) {
    PeriodicalFieldGeneratorYee periodicalBC(this->fdtd);
    this->fdtd->setFieldGenerator(&periodicalBC);
    const int numSteps = state.range_x();
    this->fdtd->changeDevice(sycl_pfc::Devices::GPU);

    while (state.KeepRunning()) {
        for (int step = 0; step < numSteps; ++step)
        {
            this->fdtd->updateFields();
        }
    }
}
BENCHMARK_REGISTER_F(GridFDTDSYCLTestZ, FDTDSYCL_GPU_periodical_vtune)->Apply(CustomArguments)->Unit(benchmark::kSecond);
