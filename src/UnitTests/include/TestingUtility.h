#pragma once

#include "Constants.h"
#include "Enums.h"
#include "Grid.h"
#include "Particle.h"
#include "ParticleTraits.h"
#include "ParticleTypes.h"
#include "Vectors.h"
#include "VectorsProxy.h"

#include "gtest/gtest.h"

#include <vector>

using namespace pfc;

// if project was built without FFT then don't run tests with FFT
#ifdef __USE_FFT__
#define ADD_TEST_FFT_PREFIX(name) name
#else
#define ADD_TEST_FFT_PREFIX(name) DISABLED_ ## name
#endif

#define ASSERT_EQ_COMPLEXFP(expected, actual) \
    ASSERT_EQ((expected).real, (actual).real); \
    ASSERT_EQ((expected).imag, (actual).imag); \

#define ASSERT_EQ_FP3(expected, actual) \
    ASSERT_EQ((expected).x, (actual).x); \
    ASSERT_EQ((expected).y, (actual).y); \
    ASSERT_EQ((expected).z, (actual).z); \

#define ASSERT_EQ_INT3(expected, actual) \
    ASSERT_EQ((expected).x, (actual).x); \
    ASSERT_EQ((expected).y, (actual).y); \
    ASSERT_EQ((expected).z, (actual).z); \

#define ASSERT_EQ_COMPLEXFP3(expected, actual) \
    ASSERT_EQ((expected).x, (actual).x); \
    ASSERT_EQ((expected).y, (actual).y); \
    ASSERT_EQ((expected).z, (actual).z); \

#define ASSERT_EQ_VECTOR(expected, actual, dimension) \
    for (int d = 0; d < dimension; d++) \
        ASSERT_EQ(((expected))[d], ((actual))[d]); \

// Assert two FP3s are nearly equal:
// if expected value is not near zero, expect relative error is smaller than
// m_maxRelativeError, else expect absolute error is smaller than
// m_maxAbsoluteError.
#define ASSERT_NEAR_VECTOR(expected, actual) \
    if ((expected).norm() > this->maxAbsoluteError) \
        ASSERT_LE(dist(expected, actual) / (expected).norm(), this->maxRelativeError); \
    else \
        ASSERT_LE(dist(expected, actual), this->maxAbsoluteError);

#define ASSERT_NEAR_FP3(expected, actual) \
    if ((expected).norm() > this->maxAbsoluteError) \
        ASSERT_LE(dist(expected, actual) / (expected).norm(), this->maxRelativeError); \
    else \
        ASSERT_LE(dist(expected, actual), this->maxAbsoluteError); \

#define ASSERT_NEAR_FP(expected, actual) \
    if (fabs(expected) > this->maxAbsoluteError) \
        ASSERT_LE(fabs(expected - actual) / fabs(expected), this->maxRelativeError); \
    else \
        ASSERT_LE(fabs(expected - actual), this->maxAbsoluteError); \

#define ASSERT_NEAR_COMPLEXFP(expected, actual) \
    ASSERT_NEAR_FP((expected).real, (actual).real); \
    ASSERT_NEAR_FP((expected).imag, (actual).imag); \

#define ASSERT_NEAR_MODULE_COMPLEXFP(expected, actual) \
    if ((expected).getModule() > this->maxAbsoluteError) \
        ASSERT_LE((expected - actual).getModule() / (expected).getModule(), this->maxRelativeError); \
    else \
        ASSERT_LE((expected - actual).getModule(), this->maxAbsoluteError); \

class BaseFixture : public testing::Test {
protected:

    virtual void SetUp();
    virtual void TearDown();

    // Return whether two FP3s have coords differ by not larger than eps each.
    bool nearFP3(const FP3 & a, const FP3 & b, const FP eps);

    // Get uniformly distributed in [a, b) pseudo-random number.
    FP urand(FP a, FP b) const;
    int urandInt(int a, int b);
    // Get distributed in [a, b) pseudo-random vector.
    FP3 urandFP3(FP3 a, FP3 b);
    Int3 urandInt3(Int3 a, Int3 b);
    // Get _n_ random vectors between _minValue_ and _maxValue_.
    std::vector<FP3> randomVectors(int n, FP3 & minValue, FP3 & maxValue);

    FP maxAbsoluteError; // max absolute error that counts as "passed"
    FP maxRelativeError; // max relative error that counts as "passed"
};


template<class ParticleType>
class BaseParticleFixture : public BaseFixture {
public:
    typedef ParticleType Particle;
    typedef typename ParticleTraits<Particle>::PositionType PositionType;
    typedef typename ParticleTraits<Particle>::MomentumType MomentumType;
    typedef typename ParticleTraits<Particle>::GammaType GammaType;
    typedef typename ParticleTraits<Particle>::WeightType WeightType;
    typedef typename ScalarType<MomentumType>::Type Real;
    typedef typename ParticleTraits<Particle>::TypeIndexType TypeIndexType;
    static const int dimension = VectorDimensionHelper<PositionType>::dimension;
    static const int momentumDimension = VectorDimensionHelper<MomentumType>::dimension;

    virtual void SetUp()
    {
        BaseFixture::SetUp();

        ParticleInfo::typesVector = { {constants::electronMass, constants::electronCharge},//electron
                                    {constants::electronMass, -constants::electronCharge},//positron
                                    {constants::protonMass, -constants::electronCharge},//proton
                                    {constants::electronMass, 0.0}//photon
                                    };
        ParticleInfo::types = &ParticleInfo::typesVector[0];
        ParticleInfo::numTypes = sizeParticleTypes;
    }

    // Helper function to unify initialization of positions for 1d, 2d and 3d
    // In 1d y, z are ignored, in 2d z is ignored
    PositionType getPosition(Real x, Real y, Real z) const
    {
        Real positionArray[] = { x, y, z };
        PositionType position;
        for (int d = 0; d < dimension; d++)
            position[d] = positionArray[d];
        return position;
    }

    Particle randomParticle(TypeIndexType type = Electron)
    {
        Real minPosition = -10;
        Real maxPosition = 10;
        return randomParticle(getPosition(minPosition, minPosition, minPosition),
            getPosition(maxPosition, maxPosition, maxPosition), type);
    }

    Particle randomParticle(PositionType minPosition, PositionType maxPosition, TypeIndexType type)
    {
        PositionType position;
        for (int d = 0; d < dimension; d++)
            position[d] = urand(minPosition[d], maxPosition[d]);
        Real minMomentum = -10;
        Real maxMomentum = 10;
        MomentumType momentum(urand(minMomentum, maxMomentum),
            urand(minMomentum, maxMomentum), urand(minMomentum, maxMomentum));
        WeightType weight = static_cast<WeightType>(urand(1e-5, 1e5));
        return Particle(position, momentum, weight, type);
    }

    template<class ConstParticleRef1, class ConstParticleRef2>
    bool eqParticles_(ConstParticleRef1& a, ConstParticleRef2& b) const
    {
        return (a.getPosition() == b.getPosition()) &&
            (a.getMomentum() == b.getMomentum()) &&
            (a.getMass() == b.getMass()) &&
            (a.getCharge() == b.getCharge()) &&
            (a.getWeight() == b.getWeight());
    }
};

template <class ParticleArrayType>
class ParticleArrayTest : public BaseParticleFixture<typename ParticleArrayType::ParticleType> {
public:
    typedef ParticleArrayType ParticleArray;
    typedef typename ParticleArrayType::ParticleType Particle;

    bool eqParticleArrays(ParticleArray& a, ParticleArray& b) const
    {
        if (a.size() != b.size())
            return false;
        for (int i = 0; i < a.size(); i++)
        {
            auto fromA = a[i], fromB = b[i];
            if (!this->eqParticles_(fromA, fromB))
                return false;
        }
        return true;
    }
};


template <class SpeciesArrayType>
class SpeciesTest : public ParticleArrayTest<SpeciesArrayType> {
public:
    typedef SpeciesArrayType SpeciesArray;
    typedef typename SpeciesArrayType::ParticleType Particle;

    bool eqParticleArrays(SpeciesArray& a, SpeciesArray& b) const
    {
        if (a.size() != b.size())
            return false;
        for (int i = 0; i < a.size(); i++)
        {
            if (!this->eqParticles_(a[i], b[i]))
                return false;
        }
        return true;
    }
};

template <class ParticleArrayType>
class ThinningTest : public BaseParticleFixture<typename ParticleArrayType::ParticleType> {
public:
    typedef ParticleArrayType ParticleArray;
    typedef typename ParticleArrayType::ParticleType Particle;
    typedef typename ParticleTraits<Particle>::PositionType PositionType;
    typedef typename ParticleTraits<Particle>::MomentumType MomentumType;
    typedef typename ParticleTraits<Particle>::TypeIndexType TypeIndexType;
    static const int dimension = VectorDimensionHelper<PositionType>::dimension;
    static const int momentumDimension = VectorDimensionHelper<MomentumType>::dimension;

    Particle randomParticleWithWeight(FP weight, TypeIndexType type = Electron)
    {
        FP limitPosition = 10;
        PositionType position;
        for (int d = 0; d < dimension; d++)
            position[d] = this->urand(-limitPosition, limitPosition);
        FP minMomentum = -10;
        FP maxMomentum = 10;
        MomentumType momentum(this->urand(minMomentum, maxMomentum),
            this->urand(minMomentum, maxMomentum), this->urand(minMomentum, maxMomentum));
        return Particle(position, momentum, weight, type);
    }

    void addRandomParticlesWithSameWeight(ParticleArray& particles, size_t size, TypeIndexType type = Electron)
    {
        FP newWeight = this->urand(1e-5, 1e5);
        for (size_t idx = 0; idx < size; idx++)
        {
            Particle newParticle = randomParticleWithWeight(newWeight, type);
            particles.pushBack(newParticle);
        }
    }

    void addRandomParticles(ParticleArray& particles, size_t size, TypeIndexType type = Electron)
    {
        for (size_t idx = 0; idx < size; idx++)
        {
            Particle newParticle = this->randomParticle(type);
            particles.pushBack(newParticle);
        }
    }

    void addRandomParticles(ParticleArray& particles, size_t size, PositionType minPosition,
        PositionType maxPosition, TypeIndexType type)
    {
        for (size_t idx = 0; idx < size; idx++)
        {
            Particle newParticle = randomParticle(minPosition, maxPosition, type);
            particles.pushBack(newParticle);
        }
    }

    FP totalWeight(ParticleArray& a) const
    {
        FP sumWeights = 0.0;
        for (int idx = 0; idx < a.size(); idx++)
        {
            sumWeights += a[idx].getWeight();
        }
        return sumWeights;
    }

    double totalEnergy(ParticleArray& a) const
    {
        FP sumEnergy = 0.0;
        for (int idx = 0; idx < a.size(); idx++)
        {
            FP mass = a[idx].getMass();
            FP c = Constants<FP>::lightVelocity();
            FP mc = mass * c;
            FP energyParticle = sqrt(mc * mc + a[idx].getMomentum().norm2());
            sumEnergy += energyParticle * a[idx].getWeight();
        }
        return sumEnergy;
    }
};

template<class gridType>
class BaseGridFixture : public BaseFixture {
public:
    gridType * grid;
    FP3 minCoords;
    FP3 maxCoords;
    FP timeStep;
protected:
    virtual void SetUp() {
        BaseFixture::SetUp();
        maxAbsoluteError = (FP)1e-2;
        maxRelativeError = (FP)1e-1;
        timeStep = 1e-15;

        Int3 gridSize(11, 5, 6);
        minCoords = FP3(-1.0, 0.0, 0.0);
        maxCoords = FP3(1.0, 1.0, 1.0);
        FP3 steps((maxCoords.x - minCoords.x) / gridSize.x,
            (maxCoords.y - minCoords.y) / gridSize.y,
            (maxCoords.z - minCoords.z) / gridSize.z);
        grid = new gridType(gridSize, minCoords, steps, gridSize);
    }

    FP3 internalPoint() {
        return urandFP3(minCoords, maxCoords);
    }

    ~BaseGridFixture()
    {
        delete(grid);
    }
};

// necessary to run test with different fields and dimensions (1d, 2d, 3d)
// wave proparates along VAxis (x, y, z)
template <class TFieldSolver, class TGrid, int VDimension, CoordinateEnum VAxis>
struct TypeDefinitionsFieldTest
{
    using FieldSolverType = TFieldSolver;
    using GridType = TGrid;
    static const int dimension = VDimension;
    static const CoordinateEnum axis = VAxis;
};
