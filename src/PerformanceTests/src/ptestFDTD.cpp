#include "TestingUtility.h"

#include "Fdtd.h"

template <class axis>
class GridFDTDTest : public BaseGridFixture<YeeGrid> {
public:
    virtual void SetUp(const ::benchmark::State& st)
    {
        BaseGridFixture<YeeGrid>::SetUp(st);

        fdtd = new FDTD(grid, this->timeStep);
    }

    FDTD* fdtd;

    ~GridFDTDTest() {
        delete fdtd;
    }

    FP3 eTest(FP x, FP y, FP z, FP t);
    FP3 bTest(FP x, FP y, FP z, FP t);
};

template<>
FP3 GridFDTDTest<axisX>::eTest(FP x, FP y, FP z, FP t)
{
    return FP3(0, sin(2 * constants::pi * (-constants::c * t + x)), 0);
}

template<>
FP3 GridFDTDTest<axisX>::bTest(FP x, FP y, FP z, FP t)
{
    return FP3(0, 0, sin(2 * constants::pi * (-constants::c * t + x)));
}

template<>
FP3 GridFDTDTest<axisY>::eTest(FP x, FP y, FP z, FP t)
{
    return FP3(0, 0, sin(2 * constants::pi * (-constants::c * t + y)));
}

template<>
FP3 GridFDTDTest<axisY>::bTest(FP x, FP y, FP z, FP t)
{
    return FP3(sin(2 * constants::pi * (-constants::c * t + y)), 0, 0);
}

template<>
FP3 GridFDTDTest<axisZ>::eTest(FP x, FP y, FP z, FP t)
{
    return FP3(sin(2 * constants::pi * (-constants::c * t + z)), 0, 0);
}

template<>
FP3 GridFDTDTest<axisZ>::bTest(FP x, FP y, FP z, FP t)
{
    return FP3(0, sin(2 * constants::pi * (-constants::c * t + z)), 0);
}

static void CustomArguments(benchmark::internal::Benchmark* b) {
    b->Args({ 16 });
    b->Iterations(5);
}

using GridFDTDTestZ = GridFDTDTest<axisZ>;
BENCHMARK_DEFINE_F(GridFDTDTestZ, FDTD_periodical_vtune)(benchmark::State& state) {
    PeriodicalFieldGeneratorYee periodicalBC(this->fdtd);
    this->fdtd->setFieldGenerator(&periodicalBC);
    const int numSteps = state.range_x();

    while (state.KeepRunning()) {
        for (int step = 0; step < numSteps; ++step)
        {
            this->fdtd->updateFields();
        }
    }
}
BENCHMARK_REGISTER_F(GridFDTDTestZ, FDTD_periodical_vtune)->Apply(CustomArguments)->Unit(benchmark::kSecond);
