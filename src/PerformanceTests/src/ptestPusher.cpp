#include "TestingUtility.h"

#include "ParticleArray.h"
#include "Pusher.h"

#include <chrono>

static void CustomArguments(benchmark::internal::Benchmark* b) {
    b->Args({ 1000000, 100 });
    b->Iterations(3);
}

using particleArrayAoS = PusherTest<ParticleArrayAoS3d>;
BENCHMARK_DEFINE_F(particleArrayAoS, pusher)(benchmark::State& state) {
    BorisPusher pusher;
    while (state.KeepRunning()) {
        for(size_t iter = 0; iter < state.range_y(); iter ++)
            pusher(particles, fields, dt);
    }
}
BENCHMARK_REGISTER_F(particleArrayAoS, pusher)->Apply(CustomArguments)->Unit(benchmark::kSecond);

using particleArraySoA = PusherTest<ParticleArray3d>;
BENCHMARK_DEFINE_F(particleArraySoA, pusher)(benchmark::State& state) {
    BorisPusher pusher;
    while (state.KeepRunning()) {
        for (size_t iter = 0; iter < state.range_y(); iter++)
            pusher(particles, fields, dt);
    }
}
BENCHMARK_REGISTER_F(particleArraySoA, pusher)->Apply(CustomArguments)->Unit(benchmark::kSecond);
