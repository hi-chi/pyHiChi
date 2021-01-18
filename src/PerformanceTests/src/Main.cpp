#include <benchmark/benchmark.h>

#include "Particle.h"

namespace pfc {

    namespace ParticleInfo {
        std::vector<ParticleType> typesVector;
        const ParticleType* types;
        short numTypes;
    } // namespace ParticleInfo

} // namespace pica

BENCHMARK_MAIN();