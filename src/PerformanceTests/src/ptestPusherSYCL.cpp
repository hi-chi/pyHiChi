#include "TestingUtility.h"

#include "sycl/ParticleArraySYCL.h"
#include "sycl/PusherSYCL.h"

#include <chrono>
#include <iostream>

static void CustomArgumentsDRAM(benchmark::internal::Benchmark* b) { //600MB field(240MB) + particles(360MB)
    b->Args({ 1000000, 100 });
    b->Iterations(3);
}

static void CustomArgumentsCache(benchmark::internal::Benchmark* b) { //15MB field(6MB) + particles(9MB)
    b->Args({ 250000, 40000 });
    b->Iterations(3);
}

using particleArrayAoS = PusherTestSYCL<sycl_pfc::ParticleArrayAoS3d>;
BENCHMARK_DEFINE_F(particleArrayAoS, pusherSYCL_CPU)(benchmark::State& state) {
    sycl_pfc::BorisPusherSYCL pusher;
    while (state.KeepRunning()) {
        for(size_t iter = 0; iter < state.range_y(); iter ++)
            pusher(sycl_pfc::Devices::CPU, particles, fields, dt);
    }
}
BENCHMARK_REGISTER_F(particleArrayAoS, pusherSYCL_CPU)->Apply(CustomArgumentsDRAM)->Unit(benchmark::kSecond);
//BENCHMARK_REGISTER_F(particleArrayAoS, pusherSYCL_CPU)->Apply(CustomArgumentsCache)->Unit(benchmark::kSecond);

BENCHMARK_DEFINE_F(particleArrayAoS, pusherSYCL_GPU)(benchmark::State& state) {
    sycl_pfc::BorisPusherSYCL pusher;
    while (state.KeepRunning()) {
        for (size_t iter = 0; iter < state.range_y(); iter++)
            pusher(sycl_pfc::Devices::GPU, particles, fields, dt);
    }
}
BENCHMARK_REGISTER_F(particleArrayAoS, pusherSYCL_GPU)->Apply(CustomArgumentsDRAM)->Unit(benchmark::kSecond);

using particleArraySoA = PusherTestSYCL<sycl_pfc::ParticleArray3d>;
BENCHMARK_DEFINE_F(particleArraySoA, pusherSYCL_CPU)(benchmark::State& state) {
    sycl_pfc::BorisPusherSYCL pusher;
    while (state.KeepRunning()) {
        for (size_t iter = 0; iter < state.range_y(); iter++)
           pusher(sycl_pfc::Devices::CPU, particles, fields, dt);
    }
}
BENCHMARK_REGISTER_F(particleArraySoA, pusherSYCL_CPU)->Apply(CustomArgumentsDRAM)->Unit(benchmark::kSecond);
//BENCHMARK_REGISTER_F(particleArraySoA, pusherSYCL_CPU)->Apply(CustomArgumentsCache)->Unit(benchmark::kSecond);

BENCHMARK_DEFINE_F(particleArraySoA, pusherSYCL_GPU)(benchmark::State& state) {
    sycl_pfc::BorisPusherSYCL pusher;
    while (state.KeepRunning()) {
        for (size_t iter = 0; iter < state.range_y(); iter++)
            pusher(sycl_pfc::Devices::GPU, particles, fields, dt);
    }
}
BENCHMARK_REGISTER_F(particleArraySoA, pusherSYCL_GPU)->Apply(CustomArgumentsDRAM)->Unit(benchmark::kSecond);