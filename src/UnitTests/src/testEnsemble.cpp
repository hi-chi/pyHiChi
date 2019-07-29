#include "TestingUtility.h"

#include "Ensemble.h"

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

TYPED_TEST(ParticleArrayTest, EnsembleDefaultConstructor)
{
	typedef typename ParticleArrayTest<TypeParam>::ParticleArray ParticleArray;

	Ensemble<ParticleArray> particles;

	ASSERT_EQ(0, particles.size());
}

TYPED_TEST(ParticleArrayTest, EnsembleAssignment)
{
	typedef typename ParticleArrayTest<TypeParam>::ParticleArray ParticleArray;

	Ensemble<ParticleArray> particles, particlesCopy, particlesAnotherCopy;
	for (int i = 0; i < 9; i++)
		particles.addParticle(this->randomParticle());
	for (int i = 0; i < 11; i++)
		particlesCopy.addParticle(this->randomParticle());
	particlesCopy = particles;
	particlesAnotherCopy = particles;
	for (int t = 0; t < sizeParticleTypes; t++)
	{
		ASSERT_TRUE(this->eqParticleArrays(particles[t], particlesCopy[t]));
		ASSERT_TRUE(this->eqParticleArrays(particles[t], particlesAnotherCopy[t]));
		ASSERT_TRUE(this->eqParticleArrays(particles[particleNames[t]], particlesCopy[t]));
	}
}

TYPED_TEST(ParticleArrayTest, EnsembleSize)
{
	typedef typename ParticleArrayTest<TypeParam>::ParticleArray ParticleArray;

	Ensemble<ParticleArray> particles;
	int numParticles = 12;
	for (int i = 0; i < numParticles; i++)
		particles.addParticle(this->randomParticle(Electron));
	ASSERT_EQ(numParticles, particles[Electron].size());
	numParticles = 15;
	for (int i = 0; i < numParticles; i++)
		particles.addParticle(this->randomParticle(Positron));
	ASSERT_EQ(numParticles, particles[Positron].size());
	numParticles = 13;
	for (int i = 0; i < numParticles; i++)
		particles.addParticle(this->randomParticle(Proton));
	ASSERT_EQ(numParticles, particles[Proton].size());
	numParticles = 40;
	ASSERT_EQ(numParticles, particles.size());
}

TYPED_TEST(ParticleArrayTest, EnsembleAccess)
{
	typedef typename ParticleArrayTest<TypeParam>::ParticleArray ParticleArray;
	typedef typename ParticleArrayTest<TypeParam>::Particle ParticleType;
	typedef typename ParticleArray::ParticleProxyType ParticleProxyType;

	Ensemble<ParticleArray> particles;
	const int numParticles = 15;
	ParticleType particleArray[numParticles];
	for (int i = 0; i < numParticles; i++) {
		ParticleType particle = this->randomParticle(Electron);
		particleArray[i] = particle;
		particles.addParticle(particle);
	}
	ASSERT_EQ(numParticles, particles.size());
	for (int i = 0; i < numParticles; i++) {
		ParticleProxyType proxyParticleFromArray(particleArray[i]);
		ParticleProxyType proxyParticle(particles[Electron][i]);
		EXPECT_TRUE(this->eqParticles_(proxyParticleFromArray, proxyParticle));
	}
}
