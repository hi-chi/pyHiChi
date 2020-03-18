#pragma once
#include "Constants.h"
#include "Species.h"
#include "FieldValue.h"

#include <array>
#include <vector>

// Example of scalar particle handler,
// only provide implementation for a single particle
// So this is for prototyping and beginner users

template<class Object, size_t T_size>
using Chunk = std::array<Object, T_size>;


class ScalarParticlePusher
{
    template<class T_Particle>
    void operator()(T_Particle& particle, ValueField field, FP timeStep) {};
};

class VectorizedParticlePusher
{
    template<class T_Particle>
    void operator()(T_Particle& particle, ValueField field, FP timeStep) {};

    template<class T_Particle, int T_size>
    void operator()(Chunk<T_Particle, T_size> chunk, Chunk<ValueField, T_size> fields, FP timeStep) {};
};

class ScalarBorisPusher : public ScalarParticlePusher
{
public:

    void operator()(ParticleProxy3d& particle, ValueField field, FP timeStep)
    {
        FP3 e = field.getE();
        FP3 b = field.getB();
        FP eCoeff = timeStep * particle.getCharge() / (2 * particle.getMass()*Constants<FP>::lightVelocity());
        FP3 eMomentum = e * eCoeff;
        FP3 um = particle.getP() + eMomentum;
        FP3 t = b * eCoeff / sqrt((FP)1 + um.norm2());
        FP3 uprime = um + cross(um, t);
        FP3 s = t * (FP)2 / ((FP)1 + t.norm2());
        particle.setP(eMomentum + um + cross(uprime, s));
        particle.setPosition(particle.getPosition() + timeStep * particle.getVelocity());

    }
    
    void operator()(Particle3d& particle, ValueField field, FP timeStep)
    {
        FP3 e = field.getE();
        FP3 b = field.getB();
        FP eCoeff = timeStep * particle.getCharge() / (2 * particle.getMass()*Constants<FP>::lightVelocity());
        FP3 eMomentum = e * eCoeff;
        FP3 um = particle.getP() + eMomentum;
        FP3 t = b * eCoeff / sqrt((FP)1 + um.norm2());
        FP3 uprime = um + cross(um, t);
        FP3 s = t * (FP)2 / ((FP)1 + t.norm2());
        particle.setP(eMomentum + um + cross(uprime, s));
        particle.setPosition(particle.getPosition() + timeStep * particle.getVelocity());

    }
private:

    FP3 u, u_, velocity;
    FP b, d, type_a;
};

class VectorizedBorisPusher : public VectorizedParticlePusher
{
public:
    template<class T_Particle>
    void operator()(T_Particle& particle, ValueField field, FP timeStep)
    {
        BorisPusher.operator()(particle, field, timeStep);
    }

    template<class T_Particle, size_t T_size>
    void operator()(Chunk<T_Particle, T_size> chunk, Chunk<ValueField, T_size> fields, FP timeStep)
    {
        for (int i = 0; i < T_size; i++)
        {
            BorisPusher.operator()(chunk[i], fields[i], timeStep);
        }
    }

private:
    ScalarBorisPusher BorisPusher;
};

class ScalarRadiationReaction : public ScalarParticlePusher
{
public:

    void operator()(ParticleProxy3d& particle, ValueField field, FP timeStep)
    {
        if (particle.getType() == Electron || particle.getType() == Positron)
        {
            
            FP3 e = field.getE();
            FP3 b = field.getB();
            FP3 v = particle.getVelocity();
            FP gamma = particle.getGamma();
            FP c = Constants<FP>::lightVelocity();
            FP electronCharge = Constants<FP>::electronCharge();
            FP electronMass = Constants<FP>::electronMass();
            FP3 dp = timeStep * (2.0 / 3.0) * sqr(sqr(electronCharge) / (electronMass * sqr(c))) *
                (VP(e, b) + (1 / c) * (VP(b, VP(b, v)) + SP(v, e) * e) -
                (1 / c) * sqr(gamma) * (sqr(e + (1 / c) * VP(v, b)) - sqr(SP(e, v) / c)) * v);
            
            particle.setMomentum(particle.getMomentum() + dp);
        }
    }

    void operator()(Particle3d& particle, ValueField field, FP timeStep)
    {
        if (particle.getType() == Electron || particle.getType() == Positron)
        {
            FP3 e = field.getE();
            FP3 b = field.getB();
            FP3 v = particle.getVelocity();
            FP gamma = particle.getGamma();
            FP c = Constants<FP>::lightVelocity();
            FP electronCharge = Constants<FP>::electronCharge();
            FP electronMass = Constants<FP>::electronMass();
            FP3 dp = timeStep * (2.0 / 3.0) * sqr(sqr(electronCharge) / (electronMass * sqr(c))) *
                (VP(e, b) + (1 / c) * (VP(b, VP(b, v)) + SP(v, e) * e) -
                (1 / c) * sqr(gamma) * (sqr(e + (1 / c) * VP(v, b)) - sqr(SP(e, v) / c)) * v);
            particle.setMomentum(particle.getMomentum() + dp);
        }
    }
};

class VectorizedRadiationReaction : public VectorizedParticlePusher
{
public:
    template<class T_Particle>
    void operator()(T_Particle& particle, ValueField field, FP timeStep)
    {
        RadiationReaction.operator()(particle, field, timeStep);
    }

    template<class T_Particle, size_t T_size>
    void operator()(Chunk<T_Particle, T_size> chunk, Chunk<ValueField, T_size> fields, FP timeStep)
    {
        for (int i = 0; i < T_size; i++)
        {
            RadiationReaction.operator()(chunk[i], fields[i], timeStep);
        }
    }

private:
    ScalarRadiationReaction RadiationReaction;
};