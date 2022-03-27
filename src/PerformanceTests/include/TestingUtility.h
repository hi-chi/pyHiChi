#pragma once

#include "Constants.h"
#include "Grid.h"
#include "Particle.h"
#include "ParticleTraits.h"
#include "ParticleTypes.h"
#include "FieldValue.h"
#include "Vectors.h"
#include "VectorsProxy.h"

#include "benchmark/benchmark.h"

#include <vector>

#include "sycl/DeviceSYCL.h"

class BaseFixture : public benchmark::Fixture {
public:

    virtual void SetUp(const ::benchmark::State& st)
    {
        srand(1);
    }
    virtual void TearDown(const ::benchmark::State&) {}
    
    // Get uniformly distributed in [a, b) pseudo-random number.
    FP urand(FP a, FP b) const {
        return a + (b - a) * ((FP)rand()) / RAND_MAX;
    }
    int urandInt(int a, int b) {
        return a + rand() % (b - a + 1);
    }
    // Get distributed in [a, b) pseudo-random vector.
    FP3 urandFP3(FP3 a, FP3 b) {
        FP3 result;
        result.x = urand(a.x, b.x);
        result.y = urand(a.y, b.y);
        result.z = urand(a.z, b.z);
        return result;
    }
    Int3 urandInt3(Int3 a, Int3 b) {
        Int3 result;
        result.x = urandInt(a.x, b.x);
        result.y = urandInt(a.y, b.y);
        result.z = urandInt(a.z, b.z);
        return result;
    }
    // Get _n_ random vectors between _minValue_ and _maxValue_.
    std::vector<FP3> randomVectors(int n, FP3 & minValue, FP3 & maxValue) {
        std::vector<FP3> result(n);
        for (int i = 0; i < n; ++i)
            result[i] = urandFP3(minValue, maxValue);
        return result;
    }
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

    using BaseFixture::urand;

    virtual void SetUp(const ::benchmark::State& st)
    {
        BaseFixture::SetUp(st);
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
        Real minPosition = -1e-4;
        Real maxPosition = 1e-4;
        return randomParticle(getPosition(minPosition, minPosition, minPosition),
            getPosition(maxPosition, maxPosition, maxPosition), type);
    }

    Particle randomParticle(PositionType minPosition, PositionType maxPosition, TypeIndexType type)
    {
        PositionType position;
        for (int d = 0; d < dimension; d++)
            position[d] = urand(minPosition[d], maxPosition[d]);
        Real minMomentum = -1e-10;
        Real maxMomentum = 1-10;
        MomentumType momentum(urand(minMomentum, maxMomentum),
            urand(minMomentum, maxMomentum), urand(minMomentum, maxMomentum));
        WeightType weight = static_cast<WeightType>(urand((FP)1e-5, (FP)1e5));
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
class ParticleArrayFixture : public BaseParticleFixture<typename ParticleArrayType::ParticleType> {
public:

    typedef ParticleArrayType ParticleArray;
    typedef typename ParticleArrayType::ParticleType Particle;

    using BaseParticleFixture<Particle>::randomParticle;

    virtual void SetUp(const ::benchmark::State& st)
    {
        BaseParticleFixture<Particle>::SetUp(st);
        particles = new ParticleArray();
        for (size_t index = 0; index < st.range_x(); index++)
            particles->pushBack(randomParticle());
    }

    virtual void TearDown(const ::benchmark::State&)
    {
        delete (particles);
    }

    ParticleArray * particles;
};

template <class ParticleArrayType>
class PusherTest : public ParticleArrayFixture<ParticleArrayType> {
public:
    typedef ParticleArrayType ParticleArray;
    typedef typename ParticleArrayType::ParticleType Particle;

    using BaseFixture::urandFP3;


    virtual void SetUp(const ::benchmark::State& st)
    {
        ParticleArrayFixture<ParticleArray>::SetUp(st);
        dt = 1e-17;
        FP3 minField(-1e14, -1e14, -1e14), maxField(1e14, 1e14, 1e14);

        for (size_t index = 0; index < st.range_x(); index++)
            fields.push_back(ValueField(urandFP3(minField, maxField), urandFP3(minField, maxField)));
    }

    virtual void TearDown(const ::benchmark::State& st)
    {
        fields.clear();
        ParticleArrayFixture<ParticleArray>::TearDown(st);
    }

    std::vector<ValueField> fields;
    FP dt;
};


template <class ParticleArrayType>
class PusherTestSYCL : public ParticleArrayFixture<ParticleArrayType> {
public:
    typedef ParticleArrayType ParticleArray;
    typedef typename ParticleArrayType::ParticleType Particle;

    using BaseFixture::urandFP3;
    
    virtual void SetUp(const ::benchmark::State& st)
    {
        fields = new sycl_pfc::sycl_vector<ValueField>{ sycl_pfc::node.default_device };
        ParticleArrayFixture<ParticleArray>::SetUp(st);
        dt = 1e-17;
        FP3 minField(-1e14, -1e14, -1e14), maxField(1e14, 1e14, 1e14);

        for (size_t index = 0; index < st.range_x(); index++)
            fields->push_back(ValueField(urandFP3(minField, maxField), urandFP3(minField, maxField)));
    }

    virtual void TearDown(const ::benchmark::State& st)
    {
        fields->clear();
        ParticleArrayFixture<ParticleArray>::TearDown(st);
    }

    sycl_pfc::sycl_vector<ValueField> * fields;
    FP dt;
};


template<class gridType>
class BaseGridFixture : public BaseFixture {
public:
    gridType* grid;
    FP3 minCoords;
    FP3 maxCoords;
    FP timeStep;
protected:
    virtual void SetUp(const ::benchmark::State& st)
    {
        BaseFixture::SetUp(st);
        timeStep = 1e-15;

        Int3 gridSize(768, 768, 768);
        minCoords = FP3(-1.0, -1.0, -1.0);
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

enum axis
{
    X = 0,
    Y = 1,
    Z = 2
};

template<axis a>
class ax
{};

typedef ax<X> axisX;
typedef ax<Y> axisY;
typedef ax<Z> axisZ;