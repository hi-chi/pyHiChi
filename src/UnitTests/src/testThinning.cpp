#include "TestingUtility.h"

#include "Particle.h"
#include "ParticleArray.h"
#include "Thinning.h"


using namespace pfc;

typedef ::testing::Types<
    ParticleArray<Three, ParticleRepresentation_AoS>::Type,
    ParticleArray<Three, ParticleRepresentation_SoA>::Type
> types;
TYPED_TEST_CASE(ThinningTest, types);

TYPED_TEST(ThinningTest, simpleThinning)
{
    typedef typename ThinningTest<TypeParam>::ParticleArray ParticleArray;
    
    Thinning<ParticleArray> thin;
    
    ParticleArray particles;
    int numberParticles = 100;
    this->addRandomParticlesWithSameWeight(particles, numberParticles);
    FP originalTotalWeight = this->totalWeight(particles);

    thin.simple(particles, numberParticles/2);

    FP modifiedTotalWeight = this->totalWeight(particles);
    ASSERT_NEAR_FP(originalTotalWeight, modifiedTotalWeight);
    ASSERT_TRUE(particles.size() == numberParticles / 2);
}

TYPED_TEST(ThinningTest, levelingThinning)
{
    typedef typename ThinningTest<TypeParam>::ParticleArray ParticleArray;

    Thinning<ParticleArray> thin;
    
    ParticleArray particles;
    int numberParticles = 100;
    this->addRandomParticles(particles, numberParticles);
    FP originalTotalWeight = this->totalWeight(particles);

    thin.leveling(particles);

    FP modifiedTotalWeight = this->totalWeight(particles);
    ASSERT_TRUE(particles.size() < numberParticles);
}

TYPED_TEST(ThinningTest, numberConservativeThinning)
{
    typedef typename ThinningTest<TypeParam>::ParticleArray ParticleArray;

    Thinning<ParticleArray> thin;
    
    ParticleArray particles;
    int numberParticles = 100;
    this->addRandomParticles(particles, numberParticles);
    FP originalTotalWeight = this->totalWeight(particles);

    thin.numberConservative(particles, numberParticles/2);

    FP modifiedTotalWeight = this->totalWeight(particles);
    ASSERT_NEAR_FP(originalTotalWeight, modifiedTotalWeight);
    ASSERT_TRUE(particles.size() <= numberParticles / 2);
}

TYPED_TEST(ThinningTest, energyConservativeThinning)
{
    typedef typename ThinningTest<TypeParam>::ParticleArray ParticleArray;

    Thinning<ParticleArray> thin;
    
    ParticleArray particles;
    int numberParticles = 100;
    this->addRandomParticles(particles, numberParticles);
    FP originalTotalEnergy = this->totalEnergy(particles);

    thin.energyConservative(particles, numberParticles / 2);

    FP modifiedTotalEnergy = this->totalEnergy(particles);
    ASSERT_NEAR_FP(originalTotalEnergy, modifiedTotalEnergy);
    ASSERT_TRUE(particles.size() <= numberParticles / 2);
}