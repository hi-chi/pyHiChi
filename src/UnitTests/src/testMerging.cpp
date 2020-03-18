#include "TestingUtility.h"

#include "Particle.h"
#include "ParticleArray.h"
#include "Merging.h"


using namespace pfc;

typedef ::testing::Types<
    ParticleArray<Three, ParticleRepresentation_AoS>::Type,
    ParticleArray<Three, ParticleRepresentation_SoA>::Type
> types;
TYPED_TEST_CASE(ThinningTest, types);

TYPED_TEST(ThinningTest, kmeansMerging)
{
    typedef typename ThinningTest<TypeParam>::ParticleArray ParticleArray;

    ParticleArray particles;
    int numberParticles = 100;
    this->addRandomParticlesWithSameWeight(particles, numberParticles);
    FP originalTotalWeight = this->totalWeight(particles);

    Merging<ParticleArray>::merge_with_kmeans(particles, numberParticles/2, 30);

    FP modifiedTotalWeight = this->totalWeight(particles);
    ASSERT_NEAR_FP(originalTotalWeight, modifiedTotalWeight);
    ASSERT_TRUE(particles.size() <= numberParticles / 2);
}
