#include "TestingUtility.h"
#include "Constants.h"
#include "Particle.h"

using namespace pfc;
template <class ParticleType>
class SaveLoadParticle : public BaseParticleFixture<ParticleType> {
public:
    using typename BaseParticleFixture<ParticleType>::Particle;
    using typename BaseParticleFixture<ParticleType>::PositionType;
    using typename BaseParticleFixture<ParticleType>::MomentumType;
    using typename BaseParticleFixture<ParticleType>::GammaType;
    using typename BaseParticleFixture<ParticleType>::WeightType;

    stringstream stream;

    MomentumType momentum;
    PositionType position;
    WeightType weight;
    ParticleTypes type;
    GammaType expectedGamma;
    ParticleType particle;

    virtual void SetUp()
    {
        maxRelativeError = 1e-12;
    }
};

typedef ::testing::Types<Particle1d, Particle2d, Particle3d> types;
TYPED_TEST_CASE(SaveLoadParticle, types);

TYPED_TEST(SaveLoadParticle, testParticle)
{
    using ParticleType = typename SaveLoadParticle<TypeParam>::Particle;

    momentum = MomentumType(-231.3e9, 0.0, 1.23e-5);
    position = getPosition(-0.06, 13e-5, -2e3);
    weight = static_cast<WeightType>(1.4e2);
    type = ParticleTypes::Proton;
    particle = ParticleType(position, momentum, weight, type);
    
    particle.save(stream);
    particle = ParticleType(); // reset particle
    ParticleType particleRes;
    particleRes.load(stream);
    
    ASSERT_EQ_VECTOR(position, particleRes.getPosition(), this->dimension);
    ASSERT_EQ_VECTOR(momentum, particleRes.getMomentum(), this->momentumDimension);
    ASSERT_EQ(Constants<MassType>::protonMass(), particleRes.getMass());
    ASSERT_EQ(-Constants<ChargeType>::electronCharge(), particleRes.getCharge());
    ASSERT_EQ(weight, particleRes.getWeight());
}
