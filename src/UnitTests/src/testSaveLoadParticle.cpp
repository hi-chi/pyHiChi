#include "TestingUtility.h"
#include "Constants.h"
#include "Particle.h"
#include "ParticleArray.h"
#include "Ensemble.h"

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
        BaseParticleFixture<ParticleType>::maxRelativeError = 1e-12;
    }
};

typedef ::testing::Types<Particle1d, Particle2d, Particle3d> types;
TYPED_TEST_CASE(SaveLoadParticle, types);

TYPED_TEST(SaveLoadParticle, testParticle)
{
    using ParticleType = typename SaveLoadParticle<TypeParam>::Particle;
    using MomentumType = typename SaveLoadParticle<TypeParam>::MomentumType;
    using WeightType = typename SaveLoadParticle<TypeParam>::WeightType;

    this->momentum = MomentumType(-231.3e9, 0.0, 1.23e-5);
    this->position = this->getPosition(-0.06, 13e-5, -2e3);
    this->weight = static_cast<WeightType>(1.4e2);
    this->type = ParticleTypes::Proton;
    this->particle = ParticleType(this->position, this->momentum, this->weight, this->type);
    
    this->particle.save(this->stream);
    this->particle = ParticleType(); // reset particle
    ParticleType particleRes;
    particleRes.load(this->stream);
    
    ASSERT_EQ_VECTOR(this->position, particleRes.getPosition(), this->dimension);
    ASSERT_EQ_VECTOR(this->momentum, particleRes.getMomentum(), this->momentumDimension);
    ASSERT_EQ(Constants<MassType>::protonMass(), particleRes.getMass());
    ASSERT_EQ(-Constants<ChargeType>::electronCharge(), particleRes.getCharge());
    ASSERT_EQ(this->weight, particleRes.getWeight());
}



typedef ::testing::Types<
    ParticleArray<One, ParticleRepresentation_AoS>::Type,
    ParticleArray<Two, ParticleRepresentation_AoS>::Type,
    ParticleArray<Three, ParticleRepresentation_AoS>::Type,
    ParticleArray<One, ParticleRepresentation_SoA>::Type,
    ParticleArray<Two, ParticleRepresentation_SoA>::Type,
    ParticleArray<Three, ParticleRepresentation_SoA>::Type
> typesArray;
TYPED_TEST_CASE(ParticleArrayTest, typesArray);

TYPED_TEST(ParticleArrayTest, testParticlesSaveLoad)
{
    typedef typename ParticleArrayTest<TypeParam>::ParticleArray ParticleArray;

    ParticleArray particles, particlesCopy, res;
    stringstream stream;
    for (int i = 0; i < 51; i++)
    {
        auto p = this->randomParticle();
        particles.pushBack(p);
        particlesCopy.pushBack(p);
    }
    particles.save(stream);
    particles = ParticleArray(); // reset ParticleArray
    res.load(stream);
    ASSERT_TRUE(this->eqParticleArrays(res, particlesCopy));
}

TYPED_TEST(ParticleArrayTest, testParticlesSaveLoadProcessingWrongType)
{
    using WrongParticleArrayType = ParticleArrayAoS<One>;
    using CorrectParticleArrayType = typename ParticleArrayTest<TypeParam>::ParticleArray;

    CorrectParticleArrayType particles;
    WrongParticleArrayType res;
    stringstream stream;
    for (int i = 0; i < 51; i++)
    {
        auto p = this->randomParticle();
        particles.pushBack(p);
    }
    particles.save(stream);
    particles = CorrectParticleArrayType(); // reset ParticleArray
    if (!std::is_same<CorrectParticleArrayType, WrongParticleArrayType>::value)
        ASSERT_ANY_THROW(res.load(stream));
    else ASSERT_NO_THROW(res.load(stream));
}

TYPED_TEST(ParticleArrayTest, testEnsembleSaveLoad)
{
    typedef typename ParticleArrayTest<TypeParam>::ParticleArray ParticleArray;

    Ensemble<ParticleArray> particles, particlesCopy, res;
    stringstream stream;
    for (int i = 0; i < 8; i++)
    {
        auto p = this->randomParticle(Electron);
        particles.addParticle(p);
        particlesCopy.addParticle(p);
    }
    for (int i = 0; i < 17; i++)
    {
        auto p = this->randomParticle(Proton);
        particles.addParticle(p);
        particlesCopy.addParticle(p);
    }
    for (int i = 0; i < 3; i++)
    {
        auto p = this->randomParticle(Photon);
        particles.addParticle(p);
        particlesCopy.addParticle(p);
    }

    particles.save(stream);
    particles = Ensemble<ParticleArray>(); // reset ParticleEnsemble
    res.load(stream);

    for (int t = 0; t < sizeParticleTypes; t++)
    {
        ASSERT_TRUE(this->eqParticleArrays(res[t], particlesCopy[t]));
        ASSERT_TRUE(this->eqParticleArrays(res[particleNames[t]], particlesCopy[t]));
    }
}

TYPED_TEST(ParticleArrayTest, testEnsembleSaveLoadProcessingWrongType)
{
    using WrongParticleArrayType = ParticleArrayAoS<One>;
    using CorrectParticleArrayType = typename ParticleArrayTest<TypeParam>::ParticleArray;

    Ensemble<CorrectParticleArrayType> particles;
    Ensemble<WrongParticleArrayType> res;
    stringstream stream;
    for (int i = 0; i < 8; i++)
    {
        auto p = this->randomParticle(Electron);
        particles.addParticle(p);
    }
    for (int i = 0; i < 17; i++)
    {
        auto p = this->randomParticle(Proton);
        particles.addParticle(p);
    }
    for (int i = 0; i < 3; i++)
    {
        auto p = this->randomParticle(Photon);
        particles.addParticle(p);
    }

    particles.save(stream);
    particles = Ensemble<CorrectParticleArrayType>(); // reset ParticleArray

    if (!std::is_same<CorrectParticleArrayType, WrongParticleArrayType>::value)
        ASSERT_ANY_THROW(res.load(stream));
    else ASSERT_NO_THROW(res.load(stream));
}

