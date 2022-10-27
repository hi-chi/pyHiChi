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

TYPED_TEST(PusherTest, VayPusherRelativisticAccelerationInStaticField)
{
    typedef typename SpeciesTest<TypeParam>::SpeciesArray SpeciesArray;
    typedef typename SpeciesTest<TypeParam>::MomentumType MomentumType;
    typedef typename SpeciesTest<TypeParam>::PositionType PositionType;

    SpeciesArray speciesParticles;
    speciesParticles.pushBack(this->randomParticle(speciesParticles.getType()));
    auto particle = speciesParticles.back();

    PositionType startPosition = { (FP)0, (FP)0, (FP)0 };
    particle.setPosition(startPosition);
    MomentumType startVelocity = { (FP)0, (FP)0, (FP)0 };
    particle.setVelocity(startVelocity);

    VayPusher scalarPusher;
    FP E0 = 1.0;
    FP Ex = E0, Ey = 0.0, Ez = 0.0;
    FP Bx = 0.0, By = 0.0, Bz = 0.0;
    auto field = ValueField(Ex, Ey, Ez, Bx, By, Bz);
   
    FP N = 1000;
    FP timeStep = particle.getMass() * Constants<FP>::lightVelocity() / (particle.getCharge() * E0 * N);

    for (int i = 0; i < N; i++) {
        scalarPusher(&particle, field, timeStep);
    }

    FP p_final = particle.getMass() * Constants<FP>::lightVelocity();
    MomentumType p = particle.getMomentum();
    MomentumType finalMomentum = { p_final, (FP)0, (FP)0 };
    ASSERT_NEAR_FP3(finalMomentum, p);

    FP r_final= particle.getMass() * Constants<FP>::lightVelocity() * Constants<FP>::lightVelocity() * (sqrt((FP)2) - (FP)1) / (particle.getCharge() * E0);
    PositionType r = particle.getPosition();
    PositionType finalPosition = { r_final, (FP)0, (FP)0 };
    FP err = (fabs(r_final - r[0]) / fabs(r_final));
    ASSERT_NEAR(err, 0, 0.001);
}

TYPED_TEST(PusherTest, VayPusherOscillationInStaticMagneticField)
{
    typedef typename SpeciesTest<TypeParam>::SpeciesArray SpeciesArray;
    typedef typename SpeciesTest<TypeParam>::MomentumType MomentumType;
    typedef typename SpeciesTest<TypeParam>::PositionType PositionType;
    typedef typename SpeciesTest<TypeParam>::GammaType GammaType;

    SpeciesArray speciesParticles;
    speciesParticles.pushBack(this->randomParticle(speciesParticles.getType()));
    auto particle = speciesParticles.back();

    PositionType startPosition = { (FP)0, (FP)0, (FP)0 };
    particle.setPosition(startPosition);
    MomentumType startVelocity = { (FP)10E-5 * Constants<FP>::lightVelocity(), (FP)0, (FP)0 };
    particle.setVelocity(startVelocity);
    MomentumType p0 = particle.getMomentum();
    GammaType gamma = particle.getGamma();

    VayPusher scalarPusher;
    FP B0 = 10.0;
    FP Ex = 0.0, Ey = 0.0, Ez = 0.0;
    FP Bx = 0.0, By = 0.0, Bz = B0;
    auto field = ValueField(Ex, Ey, Ez, Bx, By, Bz);

    FP N = 1000;
    FP timeStep = Constants<FP>::pi() * particle.getMass() * Constants<FP>::lightVelocity() * gamma / (fabs(particle.getCharge()) * B0 * N);

    for (int i = 0; i < N; i++) {
        scalarPusher(&particle, field, timeStep);
    }

    MomentumType p = particle.getMomentum();
    MomentumType finalMomentum = { -p0[0], (FP)0, (FP)0 };
    ASSERT_NEAR_FP3(finalMomentum, p);

    FP r_final = -(FP)2 * p0[0] * Constants<FP>::lightVelocity() / (particle.getCharge() * B0);
    PositionType r = particle.getPosition();
    PositionType finalPosition = { (FP)0, r_final , (FP)0 };
    ASSERT_NEAR_FP(finalPosition[1], r[1]);
}
