#include "TestingUtility.h"

#include "Handler.h"

template <class SpeciesArrayType>
class PusherTest : public SpeciesTest<SpeciesArrayType> {
};


typedef ::testing::Types<
    Species<Three, Electron, ParticleRepresentation_AoS>,
    Species<Three, Positron, ParticleRepresentation_AoS>,
    Species<Three, Proton, ParticleRepresentation_AoS>,
    Species<Three, Electron, ParticleRepresentation_SoA>,
    Species<Three, Positron, ParticleRepresentation_SoA>,
    Species<Three, Proton, ParticleRepresentation_SoA>
> types;
TYPED_TEST_CASE(PusherTest, types);

TYPED_TEST(PusherTest, ScalarBorisPusherSaveEnergyForOneParticle)
{
    typedef typename SpeciesTest<TypeParam>::SpeciesArray SpeciesArray;
    typedef typename SpeciesTest<TypeParam>::MomentumType MomentumType;

    SpeciesArray speciesParticles;
    FP energy;
    speciesParticles.pushBack(this->randomParticle(speciesParticles.getType()));
    MomentumType p = speciesParticles.back().getP();
    energy = p.norm2();

    ScalarBorisPusher scalarPusher;
    FP Ex = 0.0, Ey = 0.0, Ez = 0.0;
    FP Bx = 1.0, By = 1.0, Bz = 1.0;
    FP timeStep = 0.01;

    auto particle = speciesParticles.back();
    auto field = ValueField(Ex, Ey, Ez, Bx, By, Bz);
    handleChunk(scalarPusher, particle, field, timeStep);

    p = speciesParticles.back().getP();
    ASSERT_NEAR_FP(energy, p.norm2());
}

TYPED_TEST(PusherTest, ScalarBorisPusherSaveEnergyForChunk)
{
    typedef typename SpeciesTest<TypeParam>::SpeciesArray SpeciesArray;
    typedef typename SpeciesTest<TypeParam>::MomentumType MomentumType;

    SpeciesArray speciesParticles;
    std::vector<FP> energy;
    int numParticles = 12;
    for (int i = 0; i < numParticles; i++)
    {
        speciesParticles.pushBack(this->randomParticle(speciesParticles.getType()));
        MomentumType p = speciesParticles[i].getP();
        energy.push_back(p.norm2());
    }

    ScalarBorisPusher scalarPusher;
    FP Ex = 0.0, Ey = 0.0, Ez = 0.0;
    FP Bx = 1.0, By = 1.0, Bz = 1.0;
    FP timeStep = 0.01;
    const int chunkSize = 4;
    const int numChunks = numParticles / chunkSize;
    for (int chunkIdx = 0; chunkIdx < numChunks; chunkIdx++)
    {
        const int chunkStartIdx = chunkIdx * chunkSize;
        auto chunk = Chunk<typename SpeciesArray::ParticleProxyType, chunkSize>{
            speciesParticles[chunkStartIdx],
            speciesParticles[chunkStartIdx + 1],
            speciesParticles[chunkStartIdx + 2],
            speciesParticles[chunkStartIdx + 3]};
        auto fields = Chunk<ValueField, chunkSize>{
            ValueField(Ex, Ey, Ez, Bx, By, Bz),
            ValueField(Ex, Ey, Ez, Bx, By, Bz),
            ValueField(Ex, Ey, Ez, Bx, By, Bz),
            ValueField(Ex, Ey, Ez, Bx, By, Bz) };
        handleChunk(scalarPusher, chunk, fields, timeStep);
    }

    for (int i = 0; i < numParticles; i++)
    {
        MomentumType p = speciesParticles[i].getP();
        ASSERT_NEAR_FP(energy[i], p.norm2());
    }
}


TYPED_TEST(PusherTest, VectorizedBorisPusherSaveEnergyForOneParticle)
{
    typedef typename SpeciesTest<TypeParam>::SpeciesArray SpeciesArray;
    typedef typename SpeciesTest<TypeParam>::MomentumType MomentumType;

    SpeciesArray speciesParticles;
    FP energy;
    speciesParticles.pushBack(this->randomParticle(speciesParticles.getType()));
    MomentumType p = speciesParticles.back().getP();
    energy = p.norm2();

    VectorizedBorisPusher vectorizedPusher;
    FP Ex = 0.0, Ey = 0.0, Ez = 0.0;
    FP Bx = 1.0, By = 1.0, Bz = 1.0;
    FP timeStep = 0.01;

    auto particle = speciesParticles.back();
    auto field = ValueField(Ex, Ey, Ez, Bx, By, Bz);
    handleChunk(vectorizedPusher, particle, field, timeStep);

    p = speciesParticles.back().getP();
    ASSERT_NEAR_FP(energy, p.norm2());
}

TYPED_TEST(PusherTest, VectorizedBorisPusherSaveEnergyForChunk)
{
    typedef typename SpeciesTest<TypeParam>::SpeciesArray SpeciesArray;
    typedef typename SpeciesTest<TypeParam>::MomentumType MomentumType;

    SpeciesArray speciesParticles;
    std::vector<FP> energy;
    int numParticles = 12;
    for (int i = 0; i < numParticles; i++)
    {
        speciesParticles.pushBack(this->randomParticle(speciesParticles.getType()));
        MomentumType p = speciesParticles[i].getP();
        energy.push_back(p.norm2());
    }

    VectorizedBorisPusher vectorizedPusher;
    FP Ex = 0.0, Ey = 0.0, Ez = 0.0;
    FP Bx = 1.0, By = 1.0, Bz = 1.0;
    FP timeStep = 0.01;
    const int chunkSize = 4;
    const int numChunks = numParticles / chunkSize;
    for (int chunkIdx = 0; chunkIdx < numChunks; chunkIdx++)
    {
        const int chunkStartIdx = chunkIdx * chunkSize;
        auto chunk = Chunk<typename SpeciesArray::ParticleProxyType, chunkSize>{
            speciesParticles[chunkStartIdx],
            speciesParticles[chunkStartIdx + 1],
            speciesParticles[chunkStartIdx + 2],
            speciesParticles[chunkStartIdx + 3] };
        auto fields = Chunk<ValueField, chunkSize>{
            ValueField(Ex, Ey, Ez, Bx, By, Bz),
            ValueField(Ex, Ey, Ez, Bx, By, Bz),
            ValueField(Ex, Ey, Ez, Bx, By, Bz),
            ValueField(Ex, Ey, Ez, Bx, By, Bz) };
        handleChunk(vectorizedPusher, chunk, fields, timeStep);
    }

    for (int i = 0; i < numParticles; i++)
    {
        MomentumType p = speciesParticles[i].getP();
        ASSERT_NEAR_FP(energy[i], p.norm2());
    }
}