#include "TestingUtility.h"

#include "Particle.h"
#include "ParticleArray.h"


using namespace pfc;

typedef ::testing::Types<
	ParticleArray<One, ParticleRepresentation_AoS>::Type,
	ParticleArray<Two, ParticleRepresentation_AoS>::Type,
	ParticleArray<Three, ParticleRepresentation_AoS>::Type,
	ParticleArray<One, ParticleRepresentation_SoA>::Type,
	ParticleArray<Two, ParticleRepresentation_SoA>::Type,
	ParticleArray<Three, ParticleRepresentation_SoA>::Type
> types;
TYPED_TEST_CASE(ParticleArrayTest, types);

TYPED_TEST(ParticleArrayTest, DefaultConstructor)
{
	typedef typename ParticleArrayTest<TypeParam>::ParticleArray ParticleArray;

	ParticleArray particles;

	ASSERT_EQ(0, particles.size());
}

TYPED_TEST(ParticleArrayTest, Assignment)
{
	typedef typename ParticleArrayTest<TypeParam>::ParticleArray ParticleArray;

	ParticleArray particles, particlesCopy, particlesAnotherCopy;
	for (int i = 0; i < 9; i++)
		particles.pushBack(this->randomParticle());
	for (int i = 0; i < 11; i++)
		particlesCopy.pushBack(this->randomParticle());
	particlesCopy = particles;
	particlesAnotherCopy = particles;
	ASSERT_TRUE(this->eqParticleArrays(particles, particlesCopy));
	ASSERT_TRUE(this->eqParticleArrays(particles, particlesAnotherCopy));
}

TYPED_TEST(ParticleArrayTest, Size)
{
	typedef typename ParticleArrayTest<TypeParam>::ParticleArray ParticleArray;

	ParticleArray particles;
	int numParticles = 12;
	for (int i = 0; i < numParticles; i++)
		particles.pushBack(this->randomParticle());

	ASSERT_EQ(numParticles, particles.size());
}

TYPED_TEST(ParticleArrayTest, IndexAccess)
{
	typedef typename ParticleArrayTest<TypeParam>::ParticleArray ParticleArray;
	typedef typename ParticleArrayTest<TypeParam>::Particle ParticleType;
	typedef typename ParticleArray::ParticleProxyType ParticleProxyType;

	ParticleArray particles;
	const int numParticles = 15;
	ParticleType particleArray[numParticles];
	for (int i = 0; i < numParticles; i++) {
		ParticleType particle = this->randomParticle();
		particleArray[i] = particle;
		particles.pushBack(particle);
	}
	ASSERT_EQ(numParticles, particles.size());
	for (int i = 0; i < numParticles; i++) {
		ParticleProxyType proxyParticleFromArray(particleArray[i]);
		ParticleProxyType proxyParticle(particles[i]);
	   	EXPECT_TRUE(this->eqParticles_(proxyParticleFromArray, proxyParticle));
	}
}

TYPED_TEST(ParticleArrayTest, IteratorAccess)
{
	typedef typename ParticleArrayTest<TypeParam>::ParticleArray ParticleArray;
	typedef typename ParticleArrayTest<TypeParam>::Particle ParticleType;
	typedef typename ParticleArray::ParticleProxyType ParticleProxyType;

	ParticleArray particles;
	const int numParticles = 15;
	ParticleType particleArray[numParticles];
	for (int i = 0; i < numParticles; i++) {
		ParticleType particle = this->randomParticle();
		particleArray[i] = particle;
		particles.pushBack(particle);
	}
	ASSERT_EQ(numParticles, particles.size());
	int i = 0;
	for (auto elem = particles.begin(); elem != particles.end(); elem++, i++) {
		ParticleProxyType proxyParticleFromArray(particleArray[i]);
		ParticleProxyType proxyParticle(*elem);
		EXPECT_TRUE(this->eqParticles_(proxyParticleFromArray, proxyParticle));
	}
}

TYPED_TEST(ParticleArrayTest, Back)
{
	typedef typename ParticleArrayTest<TypeParam>::ParticleArray ParticleArray;
	typedef typename ParticleArrayTest<TypeParam>::Particle ParticleType;
	typedef typename ParticleArray::ParticleProxyType ParticleProxyType;

	ParticleArray particles;
	for (int i = 0; i < 13; i++) {
		ParticleType particle = this->randomParticle();
		particles.pushBack(particle);
		ParticleProxyType proxyParticle(particle), proxyBack(particles.back());
		EXPECT_TRUE(this->eqParticles_(proxyParticle, proxyBack));
	}
}

TYPED_TEST(ParticleArrayTest, PushBack)
{
	typedef typename ParticleArrayTest<TypeParam>::ParticleArray ParticleArray;

	ParticleArray particles;
	for (int i = 0; i < 17; i++)
		particles.pushBack(this->randomParticle());
	ParticleArray particlesCopy(particles);
	ASSERT_TRUE(this->eqParticleArrays(particles, particlesCopy));
}

TYPED_TEST(ParticleArrayTest, PopBack)
{
	typedef typename ParticleArrayTest<TypeParam>::ParticleArray ParticleArray;
	typedef typename ParticleArray::ParticleProxyType ParticleProxyType;

	ParticleArray particles;
	for (int i = 0; i < 33; i++)
		particles.pushBack(this->randomParticle());
	for (int i = particles.size(); i > 0; i--) {
		ParticleArray particlesCopy = particles;
		particles.popBack();
		ASSERT_EQ(particlesCopy.size(), particles.size() + 1);
		for (int j = 0; j < particles.size(); j++)
		{
		    ParticleProxyType particle(particles[j]), particleCopy(particlesCopy[j]);
			EXPECT_TRUE(this->eqParticles_(particle, particleCopy));
		}
	}
}