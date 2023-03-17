#pragma once

#include "globals.hh"
#include "G4VPrimaryGenerator.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryVertex.hh"
#include "G4ParticleMomentum.hh"

class G4Event;

// The class is the G4ParticleGun class without the G3GunParticleMessenger instance
class MyParticleGun : public G4VPrimaryGenerator
{
public:

    MyParticleGun();
    MyParticleGun(G4int numberofparticles);
    MyParticleGun(G4ParticleDefinition* particleDef,
        G4int numberofparticles = 1);
    // Costructors. "numberofparticles" is the number of particles to be
    // shot at one invokation of GeneratePrimaryVertex() method.
    // All particles are shot with the same physical quantities.

    virtual ~MyParticleGun();

    MyParticleGun(const MyParticleGun&) = delete;
    const MyParticleGun& operator=(const MyParticleGun&) = delete;
    G4bool operator==(const MyParticleGun&) const = delete;
    G4bool operator!=(const MyParticleGun&) const = delete;

    virtual void GeneratePrimaryVertex(G4Event* evt);
    // Creates a primary vertex at the given point
    // and put primary particles to it.

  // Followings are the Set methods for the particle properties.
  // SetParticleDefinition() should be called first.  
  // By using SetParticleMomentum(), both particle_momentum_direction and
  // particle_energy(Kinetic Energy) are set.
  //
    void SetParticleDefinition(G4ParticleDefinition* aParticleDefinition);
    void SetParticleEnergy(G4double aKineticEnergy);
    void SetParticleMomentum(G4double aMomentum);
    void SetParticleMomentum(G4ParticleMomentum aMomentum);
    inline void SetParticleMomentumDirection(G4ParticleMomentum aMomDirection)
    {
        particle_momentum_direction = aMomDirection.unit();
    }
    inline void SetParticleCharge(G4double aCharge)
    {
        particle_charge = aCharge;
    }
    inline void SetParticlePolarization(G4ThreeVector aVal)
    {
        particle_polarization = aVal;
    }
    inline void SetNumberOfParticles(G4int i)
    {
        NumberOfParticlesToBeGenerated = i;
    }

    inline G4ParticleDefinition* GetParticleDefinition() const
    {
        return particle_definition;
    }
    inline G4ParticleMomentum GetParticleMomentumDirection() const
    {
        return particle_momentum_direction;
    }
    inline G4double GetParticleEnergy() const
    {
        return particle_energy;
    }
    inline G4double GetParticleMomentum() const
    {
        return particle_momentum;
    }
    inline G4double GetParticleCharge() const
    {
        return particle_charge;
    }
    inline G4ThreeVector GetParticlePolarization() const
    {
        return particle_polarization;
    }
    inline G4int GetNumberOfParticles() const
    {
        return NumberOfParticlesToBeGenerated;
    }

protected:

    virtual void SetInitialValues();

    G4int                 NumberOfParticlesToBeGenerated = 0;
    G4ParticleDefinition* particle_definition = nullptr;
    G4ParticleMomentum    particle_momentum_direction;
    G4double              particle_energy = 0.0;
    G4double              particle_momentum = 0.0;
    G4double              particle_charge = 0.0;
    G4ThreeVector         particle_polarization;

};



