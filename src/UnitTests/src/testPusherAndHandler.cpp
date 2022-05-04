#include "TestingUtility.h"

#include "Pusher.h"

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

TYPED_TEST(PusherTest, BorisPusherSaveEnergyForOneParticle)
{
    typedef typename SpeciesTest<TypeParam>::SpeciesArray SpeciesArray;
    typedef typename SpeciesTest<TypeParam>::MomentumType MomentumType;

    SpeciesArray speciesParticles;
    FP energy;
    speciesParticles.pushBack(this->randomParticle(speciesParticles.getType()));
    MomentumType p = speciesParticles.back().getP();
    energy = p.norm2();

    BorisPusher scalarPusher;
    FP Ex = 0.0, Ey = 0.0, Ez = 0.0;
    FP Bx = 1.0, By = 1.0, Bz = 1.0;
    FP timeStep = 0.01;

    auto particle = speciesParticles.back();
    auto field = ValueField(Ex, Ey, Ez, Bx, By, Bz);
    scalarPusher(&particle, field, timeStep);

    p = speciesParticles.back().getP();
    ASSERT_NEAR_FP(energy, p.norm2());
}

TYPED_TEST(PusherTest, BorisPusherSaveEnergyForChunk)
{
    typedef typename SpeciesTest<TypeParam>::SpeciesArray SpeciesArray;
    typedef typename SpeciesTest<TypeParam>::MomentumType MomentumType;

    SpeciesArray speciesParticles;
    std::vector<FP> energy;
    std::vector<ValueField> fields;
    FP Ex = 0.0, Ey = 0.0, Ez = 0.0;
    FP Bx = 1.0, By = 1.0, Bz = 1.0;
    int numParticles = 12;
    for (int i = 0; i < numParticles; i++)
    {
        speciesParticles.pushBack(this->randomParticle(speciesParticles.getType()));
        MomentumType p = speciesParticles[i].getP();
        energy.push_back(p.norm2());
        fields.push_back(ValueField(Ex, Ey, Ez, Bx, By, Bz));
    }

    BorisPusher scalarPusher;
    FP timeStep = 0.01;
    scalarPusher(&speciesParticles, fields, timeStep);
    
    for (int i = 0; i < numParticles; i++)
    {
        MomentumType p = speciesParticles[i].getP();
        ASSERT_NEAR_FP(energy[i], p.norm2());
    }
}

TYPED_TEST(PusherTest, VayPusherSaveEnergyForOneParticle)
{
    typedef typename SpeciesTest<TypeParam>::SpeciesArray SpeciesArray;
    typedef typename SpeciesTest<TypeParam>::MomentumType MomentumType;

    SpeciesArray speciesParticles;
    FP energy;
    speciesParticles.pushBack(this->randomParticle(speciesParticles.getType()));
    MomentumType p = speciesParticles.back().getP();
    energy = p.norm2();

    VayPusher scalarPusher;
    FP Ex = 0.0, Ey = 0.0, Ez = 0.0;
    FP Bx = 1.0, By = 1.0, Bz = 1.0;
    FP timeStep = 0.01;

    auto particle = speciesParticles.back();
    auto field = ValueField(Ex, Ey, Ez, Bx, By, Bz);
    scalarPusher(&particle, field, timeStep);

    p = speciesParticles.back().getP();
    ASSERT_NEAR_FP(energy, p.norm2());
}

TYPED_TEST(PusherTest, VayPusherSaveEnergyForChunk)
{
    typedef typename SpeciesTest<TypeParam>::SpeciesArray SpeciesArray;
    typedef typename SpeciesTest<TypeParam>::MomentumType MomentumType;

    SpeciesArray speciesParticles;
    std::vector<FP> energy;
    std::vector<ValueField> fields;
    FP Ex = 0.0, Ey = 0.0, Ez = 0.0;
    FP Bx = 1.0, By = 1.0, Bz = 1.0;
    int numParticles = 12;
    for (int i = 0; i < numParticles; i++)
    {
        speciesParticles.pushBack(this->randomParticle(speciesParticles.getType()));
        MomentumType p = speciesParticles[i].getP();
        energy.push_back(p.norm2());
        fields.push_back(ValueField(Ex, Ey, Ez, Bx, By, Bz));
    }

    VayPusher scalarPusher;
    FP timeStep = 0.01;
    scalarPusher(&speciesParticles, fields, timeStep);

    for (int i = 0; i < numParticles; i++)
    {
        MomentumType p = speciesParticles[i].getP();
        ASSERT_NEAR_FP(energy[i], p.norm2());
    }
}