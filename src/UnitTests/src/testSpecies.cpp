#include "TestingUtility.h"

#include "Species.h"

typedef ::testing::Types<
    Species<One, Electron, ParticleRepresentation_AoS>,
    Species<Two, Electron, ParticleRepresentation_AoS>,
    Species<Three, Electron, ParticleRepresentation_AoS>,
    Species<One, Positron, ParticleRepresentation_AoS>,
    Species<Two, Positron, ParticleRepresentation_AoS>,
    Species<Three, Positron, ParticleRepresentation_AoS>,
    Species<One, Proton, ParticleRepresentation_AoS>,
    Species<Two, Proton, ParticleRepresentation_AoS>,
    Species<Three, Proton, ParticleRepresentation_AoS>,
    Species<One, Electron, ParticleRepresentation_SoA>,
    Species<Two, Electron, ParticleRepresentation_SoA>,
    Species<Three, Electron, ParticleRepresentation_SoA>,
    Species<One, Positron, ParticleRepresentation_SoA>,
    Species<Two, Positron, ParticleRepresentation_SoA>,
    Species<Three, Positron, ParticleRepresentation_SoA>,
    Species<One, Proton, ParticleRepresentation_SoA>,
    Species<Two, Proton, ParticleRepresentation_SoA>,
    Species<Three, Proton, ParticleRepresentation_SoA>
> types;
TYPED_TEST_CASE(SpeciesTest, types);

TYPED_TEST(SpeciesTest, DefaultConstructor)
{
    typedef typename SpeciesTest<TypeParam>::SpeciesArray SpeciesArray;

    SpeciesArray speciesParticles;

    ASSERT_EQ(0, speciesParticles.size());
}

TYPED_TEST(SpeciesTest, Size)
{
    typedef typename SpeciesTest<TypeParam>::SpeciesArray SpeciesArray;

    SpeciesArray speciesParticles;
    int numParticles = 12;
    for (int i = 0; i < numParticles; i++)
        speciesParticles.pushBack(this->randomParticle(speciesParticles.getType()));

    ASSERT_EQ(numParticles, speciesParticles.size());
}

TYPED_TEST(SpeciesTest, PushBackSameType)
{
    typedef typename SpeciesTest<TypeParam>::SpeciesArray SpeciesArray;

    SpeciesArray speciesParticles;
    for (int i = 0; i < 17; i++)
        speciesParticles.pushBack(this->randomParticle(speciesParticles.getType()));
    ASSERT_EQ(17, speciesParticles.size());
}

TYPED_TEST(SpeciesTest, NoPushBackOtherType)
{
    typedef typename SpeciesTest<TypeParam>::SpeciesArray SpeciesArray;

    SpeciesArray speciesParticles;
    for (int i = 0; i < 17; i++)
        speciesParticles.pushBack(this->randomParticle(static_cast<ParticleTypes>(speciesParticles.getType() + 1)));
    ASSERT_EQ(0, speciesParticles.size());
}
