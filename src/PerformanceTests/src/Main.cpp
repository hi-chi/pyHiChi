#include <benchmark/benchmark.h>

#include "Particle.h"
#include "sycl/DeviceSYCL.h"

#include <sycl/CL/sycl.hpp>

namespace pfc {
    namespace sycl_pfc {
        AllDevices node;
    }
    sycl_pfc::sycl_vector<ParticleType> ParticleInfo::typesVector({ {Constants<FP>::electronMass(), Constants<FP>::electronCharge()},//electron
                                    {Constants<FP>::electronMass(), -Constants<FP>::electronCharge()},//positron
                                    {Constants<FP>::protonMass(), 0.0},//proton
                                    {Constants<FP>::electronMass(), 0.0 } }, sycl_pfc::node.default_device);
    const ParticleType* ParticleInfo::types = ParticleInfo::typesVector.data();
    short ParticleInfo::numTypes = sizeParticleTypes;
}

int main(int argc, char** argv) {
    try
    {
        ::benchmark::Initialize(&argc, argv);
        if (::benchmark::ReportUnrecognizedArguments(argc, argv)) return 1;
            ::benchmark::RunSpecifiedBenchmarks();
    }
    catch (std::exception& e)
    {
        std::cout << e.what() << std::endl;
    }
}