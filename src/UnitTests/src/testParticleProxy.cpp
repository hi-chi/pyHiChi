#include "TestingUtility.h"

#include "Constants.h"
#include "Particle.h"

using namespace pfc;

template <class ParticleType>
class ParticleTest : public BaseParticleFixture<ParticleType> {
public:
	using typename BaseParticleFixture<ParticleType>::Particle;
	using typename BaseParticleFixture<ParticleType>::PositionType;
	using typename BaseParticleFixture<ParticleType>::MomentumType;
	using typename BaseParticleFixture<ParticleType>::GammaType;
	using typename BaseParticleFixture<ParticleType>::WeightType;
	using typename BaseParticleFixture<ParticleType>::TypeIndexType;
};

typedef ::testing::Types<Particle1d, Particle2d, Particle3d> types;
TYPED_TEST_CASE(ParticleTest, types);

TYPED_TEST(ParticleTest, DefaultConstructorProxy)
{
	typedef typename ParticleTest<TypeParam>::Particle ParticleType;
	typedef typename ParticleTest<TypeParam>::PositionType PositionType;
	typedef typename ParticleTest<TypeParam>::MomentumType MomentumType;
	typedef typename ParticleTest<TypeParam>::GammaType GammaType;
	typedef typename ParticleTest<TypeParam>::WeightType WeightType;
	static const Dimension dimension = Dimension(VectorDimensionHelper<PositionType>::dimension);

	ParticleType particle;
	ParticleProxy<dimension> particleProxy(particle);
	

	ASSERT_EQ_VECTOR(PositionType(), particleProxy.getPosition(), this->dimension);
	ASSERT_EQ_VECTOR(MomentumType(), particleProxy.getMomentum(), this->momentumDimension);
	ASSERT_EQ_VECTOR(MomentumType(), particleProxy.getVelocity(), this->momentumDimension);
	ASSERT_EQ(static_cast<WeightType>(1.0), particleProxy.getWeight());
	ASSERT_EQ(static_cast<GammaType>(1.0), particleProxy.getGamma());
}

TYPED_TEST(ParticleTest, ConstructorProxy)
{
	typedef typename ParticleTest<TypeParam>::Particle ParticleType;
	typedef typename ParticleTest<TypeParam>::PositionType PositionType;
	typedef typename ParticleTest<TypeParam>::MomentumType MomentumType;
	typedef typename ParticleTest<TypeParam>::GammaType GammaType;
	typedef typename ParticleTest<TypeParam>::WeightType WeightType;
	typedef typename ParticleTest<TypeParam>::TypeIndexType TypeIndexType;

	static const Dimension dimension = Dimension(VectorDimensionHelper<PositionType>::dimension);

	PositionType position = this->getPosition(3.1, -32.1, 4.3e-5);
	MomentumType momentum(-231.3e9, 0.0, 1.23e-5);
	MassType mass = Constants<MassType>::electronMass();
	ChargeType charge = Constants<ChargeType>::electronCharge();
	WeightType weight = static_cast<WeightType>(1.4e2);
	GammaType expectedGamma = sqrt((FP)1 + momentum.norm2() / sqr(mass * Constants<GammaType >::c()));
	this->maxRelativeError = 1e-12;

	ParticleType particle(position, momentum, weight, ParticleTypes::Electron);
	ParticleProxy<dimension> particleProxy(particle);

	ASSERT_EQ_VECTOR(position, particleProxy.getPosition(), this->dimension);
	ASSERT_EQ_VECTOR(momentum, particleProxy.getMomentum(), this->momentumDimension);
	ASSERT_EQ(mass, particleProxy.getMass());
	ASSERT_EQ(charge, particleProxy.getCharge());
	ASSERT_EQ(weight, particleProxy.getWeight());
	ASSERT_NEAR_FP(expectedGamma, particleProxy.getGamma());
}

TYPED_TEST(ParticleTest, PiecemealConstructorProxy)
{
	typedef typename ParticleTest<TypeParam>::Particle ParticleType;
	typedef typename ParticleTest<TypeParam>::PositionType PositionType;
	typedef typename ParticleTest<TypeParam>::MomentumType MomentumType;
	typedef typename ParticleTest<TypeParam>::GammaType GammaType;
	typedef typename ParticleTest<TypeParam>::WeightType WeightType;
	typedef typename ParticleTest<TypeParam>::TypeIndexType TypeIndexType;
	static const Dimension dimension = Dimension(VectorDimensionHelper<PositionType>::dimension);

	PositionType position = this->getPosition(-12.34, 0.2, 423.12e-2);
	MomentumType momentum(3254.23, -123.324, 1.23e5);
	MassType mass = Constants<MassType>::electronMass();
	ChargeType charge = Constants<ChargeType>::electronCharge();
	WeightType weight = static_cast<WeightType>(1.4e2);
	GammaType gamma = sqrt((FP)1 + momentum.norm2() / sqr(mass * Constants<GammaType >::c()));
	TypeIndexType type = ParticleTypes::Electron;
	MomentumType p = momentum / (Constants<GammaType>::c() * mass);


	ParticleProxy<dimension> particleProxy(position, p, weight, type, gamma);

	ASSERT_EQ_VECTOR(position, particleProxy.getPosition(), this->dimension);
	ASSERT_EQ_VECTOR(momentum, particleProxy.getMomentum(), this->momentumDimension);
	ASSERT_EQ(mass, particleProxy.getMass());
	ASSERT_EQ(charge, particleProxy.getCharge());
	ASSERT_EQ(weight, particleProxy.getWeight());
	ASSERT_EQ(gamma, particleProxy.getGamma());
}

TYPED_TEST(ParticleTest, CopyConstructorProxy)
{
	typedef typename ParticleTest<TypeParam>::Particle ParticleType;
	typedef typename ParticleTest<TypeParam>::PositionType PositionType;
	typedef typename ParticleTest<TypeParam>::MomentumType MomentumType;
	typedef typename ParticleTest<TypeParam>::GammaType GammaType;
	typedef typename ParticleTest<TypeParam>::WeightType WeightType;
	typedef typename ParticleTest<TypeParam>::TypeIndexType TypeIndexType;

	static const Dimension dimension = Dimension(VectorDimensionHelper<PositionType>::dimension);

	PositionType position = this->getPosition(-134.12, 412.6342, 2346.562);
	MomentumType momentum(-4531.23e5, 6534.123e3, 12.32);
	WeightType weight = 213.51;
	ParticleType particle(position, momentum, weight, ParticleTypes::Electron);
	ParticleProxy<dimension> particleProxy(particle);

	ParticleProxy<dimension> copyParticleProxy(particleProxy);

	ASSERT_TRUE(this->eqParticles_(particleProxy, copyParticleProxy));
}

TYPED_TEST(ParticleTest, GetSetPositionProxy)
{
	typedef typename ParticleTest<TypeParam>::Particle ParticleType;
	typedef typename ParticleTest<TypeParam>::PositionType PositionType;
	static const Dimension dimension = Dimension(VectorDimensionHelper<PositionType>::dimension);

	ParticleType particle = this->randomParticle();
	ParticleProxy<dimension> particleProxy(particle);
	PositionType newPosition = this->getPosition(54.126, -431.35, 35.65);
	PositionType newPosition2 = this->getPosition(12.34, -4567.12, 34.56);

	particle.setPosition(newPosition);

	ASSERT_EQ_VECTOR(newPosition, particleProxy.getPosition(), this->dimension);

	particleProxy.setPosition(newPosition2);

	ASSERT_EQ_VECTOR(newPosition2, particleProxy.getPosition(), this->dimension);
	ASSERT_EQ_VECTOR(newPosition2, particle.getPosition(), this->dimension);
}

TYPED_TEST(ParticleTest, GetSetProxyPositionProxy)
{
	typedef typename ParticleTest<TypeParam>::Particle ParticleType;
	typedef typename ParticleTest<TypeParam>::PositionType PositionType;
	static const Dimension dimension = Dimension(VectorDimensionHelper<PositionType>::dimension);
	typedef typename ParticleProxy<dimension>::PositionTypeProxy PositionTypeProxy;

	ParticleType particle = this->randomParticle();
	ParticleProxy<dimension> particleProxy(particle);
	PositionType newPosition = this->getPosition(54.126, -431.35, 35.65);

	PositionTypeProxy particlePosition = particleProxy.getProxyPosition();
	particlePosition = newPosition;

	ASSERT_EQ_VECTOR(newPosition, particle.getPosition(), this->dimension);
	ASSERT_EQ_VECTOR(newPosition, particleProxy.getPosition(), this->dimension);
}

TYPED_TEST(ParticleTest, GetSetMomentumProxy)
{
	typedef typename ParticleTest<TypeParam>::Particle ParticleType;
	typedef typename ParticleTest<TypeParam>::PositionType PositionType;
	typedef typename ParticleTest<TypeParam>::MomentumType MomentumType;
	typedef typename ParticleTest<TypeParam>::GammaType GammaType;
	static const Dimension dimension = Dimension(VectorDimensionHelper<PositionType>::dimension);

	ParticleType particle = this->randomParticle();
	ParticleProxy<dimension> particleProxy(particle);
	MomentumType newMomentum(54.12e+4, -543.63e-2, 643.165e5);
	MomentumType newMomentum2(12.34e+5, -567.89e-3, 123.456e7);
	this->maxRelativeError = 1e-12;

	particle.setMomentum(newMomentum);
	GammaType expectedGamma = sqrt((FP)1 + newMomentum.norm2() / sqr(particle.getMass() * Constants<GammaType >::c()));

	ASSERT_NEAR_VECTOR(newMomentum, particleProxy.getMomentum());
	ASSERT_NEAR_FP(expectedGamma, particleProxy.getGamma());


	particleProxy.setMomentum(newMomentum2);
	expectedGamma = sqrt((FP)1 + newMomentum2.norm2() / sqr(particle.getMass() * Constants<GammaType >::c()));

	ASSERT_NEAR_VECTOR(newMomentum2, particleProxy.getMomentum());
	ASSERT_NEAR_FP(expectedGamma, particleProxy.getGamma());
	ASSERT_NEAR_VECTOR(newMomentum2, particle.getMomentum());
	ASSERT_NEAR_FP(expectedGamma, particle.getGamma());
}

TYPED_TEST(ParticleTest, GetProxyPAndGetProxyGammaProxy)
{
	typedef typename ParticleTest<TypeParam>::Particle ParticleType;
	typedef typename ParticleTest<TypeParam>::PositionType PositionType;
	typedef typename ParticleTest<TypeParam>::MomentumType MomentumType;
	typedef typename ParticleTest<TypeParam>::GammaType GammaType;
	static const Dimension dimension = Dimension(VectorDimensionHelper<PositionType>::dimension);
	typedef typename ParticleProxy<dimension>::MomentumTypeProxy MomentumTypeProxy;
	typedef typename ParticleProxy<dimension>::GammaTypeProxy GammaTypeProxy;


	ParticleType particle = this->randomParticle();
	ParticleProxy<dimension> particleProxy(particle);
	MomentumType newP(54.12e+4, -543.63e-2, 643.165e5);
	this->maxRelativeError = 1e-12;

	MomentumTypeProxy particleP = particleProxy.getProxyP();
	GammaTypeProxy particleGamma = particleProxy.getProxyGamma();
	GammaType newGamma = sqrt((FP)1 + newP.norm2());

	particleP = newP;
	particleGamma.get() = newGamma;

	ASSERT_NEAR_VECTOR(newP, particleProxy.getP());
	ASSERT_NEAR_FP(newGamma, particleProxy.getGamma());
	ASSERT_NEAR_VECTOR(newP, particle.getP());
	ASSERT_NEAR_FP(newGamma, particle.getGamma());
}

TYPED_TEST(ParticleTest, GetSetVelocityProxy)
{
	typedef typename ParticleTest<TypeParam>::Particle ParticleType;
	typedef typename ParticleTest<TypeParam>::PositionType PositionType;
	typedef typename ParticleTest<TypeParam>::MomentumType MomentumType;
	static const Dimension dimension = Dimension(VectorDimensionHelper<PositionType>::dimension);

	ParticleType particle = this->randomParticle();
	ParticleProxy<dimension> particleProxy(particle);
	MomentumType newVelocity(5243.1654, -56.23e5, -65.237e-4);
	MomentumType newVelocity2(1234.5678, -90.12e4, -56.789e-1);

	particle.setVelocity(newVelocity);
	this->maxRelativeError = 1e-12;

	ASSERT_NEAR_VECTOR(newVelocity, particleProxy.getVelocity());

	particleProxy.setVelocity(newVelocity2);

	ASSERT_NEAR_VECTOR(newVelocity2, particleProxy.getVelocity());
	ASSERT_NEAR_VECTOR(newVelocity2, particle.getVelocity());
}

TYPED_TEST(ParticleTest, GetGammaProxy)
{
	typedef typename ParticleTest<TypeParam>::Particle ParticleType;
	typedef typename ParticleTest<TypeParam>::PositionType PositionType;
	typedef typename ParticleTest<TypeParam>::GammaType GammaType;
	static const Dimension dimension = Dimension(VectorDimensionHelper<PositionType>::dimension);

	ParticleType particle = this->randomParticle();
	ParticleProxy<dimension> particleProxy(particle);
	GammaType expectedGamma = sqrt(static_cast<GammaType>(1.0) + particle.getMomentum().norm2() /
		sqr(particle.getMass() * Constants<GammaType>::c()));
	this->maxRelativeError = 1e-12;

	ASSERT_NEAR_FP(expectedGamma, particleProxy.getGamma());
}

TYPED_TEST(ParticleTest, GetMassProxy)
{
	typedef typename ParticleTest<TypeParam>::Particle ParticleType;
	typedef typename ParticleTest<TypeParam>::PositionType PositionType;
	static const Dimension dimension = Dimension(VectorDimensionHelper<PositionType>::dimension);

	ParticleType particle = this->randomParticle();
	ParticleProxy<dimension> particleProxy(particle);

	ASSERT_EQ(Constants<MassType>::electronMass(), particleProxy.getMass());
}

TYPED_TEST(ParticleTest, GetChargeProxy)
{
	typedef typename ParticleTest<TypeParam>::Particle ParticleType;
	typedef typename ParticleTest<TypeParam>::PositionType PositionType;
	static const Dimension dimension = Dimension(VectorDimensionHelper<PositionType>::dimension);

	ParticleType particle = this->randomParticle();
	ParticleProxy<dimension> particleProxy(particle);

	ASSERT_EQ(Constants<ChargeType>::electronCharge(), particleProxy.getCharge());
}

TYPED_TEST(ParticleTest, GetSetWeightProxy)
{
	typedef typename ParticleTest<TypeParam>::Particle ParticleType;
	typedef typename ParticleTest<TypeParam>::PositionType PositionType;
	typedef typename ParticleTest<TypeParam>::WeightType WeightType;
	static const Dimension dimension = Dimension(VectorDimensionHelper<PositionType>::dimension);

	ParticleType particle = this->randomParticle();
	ParticleProxy<dimension> particleProxy(particle);
	WeightType newWeight = 6523.54;
	WeightType newWeight2 = 123.45;

	particle.setWeight(newWeight);

	ASSERT_EQ(newWeight, particleProxy.getWeight());

	particleProxy.setWeight(newWeight2);

	ASSERT_NEAR_FP(newWeight2, particleProxy.getWeight());
	ASSERT_NEAR_FP(newWeight2, particle.getWeight());
}

TYPED_TEST(ParticleTest, GetProxyWeightProxy)
{
	typedef typename ParticleTest<TypeParam>::Particle ParticleType;
	typedef typename ParticleTest<TypeParam>::PositionType PositionType;
	typedef typename ParticleTest<TypeParam>::WeightType WeightType;
	static const Dimension dimension = Dimension(VectorDimensionHelper<PositionType>::dimension);
	typedef typename ParticleProxy<dimension>::WeightTypeProxy WeightTypeProxy;

	ParticleType particle = this->randomParticle();
	ParticleProxy<dimension> particleProxy(particle);
	WeightType newWeight = 6523.54;

	WeightTypeProxy particleWeight = particleProxy.getProxyWeight();
	particleWeight.get() = newWeight;

	ASSERT_NEAR_FP(newWeight, particleProxy.getWeight());
	ASSERT_NEAR_FP(newWeight, particle.getWeight());
}