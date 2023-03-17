#include "MyParticleGun.h"

#include "G4SystemOfUnits.hh"
#include "G4PrimaryParticle.hh"
#include "G4Event.hh"
#include "G4ios.hh"

MyParticleGun::MyParticleGun()
{
    SetInitialValues();
}

MyParticleGun::MyParticleGun(G4int numberofparticles)
{
    SetInitialValues();
    NumberOfParticlesToBeGenerated = numberofparticles;
}

MyParticleGun::MyParticleGun(G4ParticleDefinition* particleDef,
    G4int numberofparticles)
{
    SetInitialValues();
    NumberOfParticlesToBeGenerated = numberofparticles;
    SetParticleDefinition(particleDef);
}

void MyParticleGun::SetInitialValues()
{
    NumberOfParticlesToBeGenerated = 1;
    particle_definition = nullptr;
    G4ThreeVector zero;
    particle_momentum_direction = (G4ParticleMomentum)zero;
    particle_energy = 0.0;
    particle_momentum = 0.0;
    particle_position = zero;
    particle_time = 0.0;
    particle_polarization = zero;
    particle_charge = 0.0;
}

MyParticleGun::~MyParticleGun()
{
}

void MyParticleGun::
SetParticleDefinition(G4ParticleDefinition* aParticleDefinition)
{
    if (aParticleDefinition == nullptr)
    {
        G4Exception("MyParticleGun::SetParticleDefinition()", "Event0101",
            FatalException, "Null pointer is given.");
    }
    if (aParticleDefinition->IsShortLived())
    {
        if (aParticleDefinition->GetDecayTable() == nullptr)
        {
            G4ExceptionDescription ED;
            ED << "MyParticleGun does not support shooting a short-lived "
                << "particle without a valid decay table." << G4endl;
            ED << "MyParticleGun::SetParticleDefinition for "
                << aParticleDefinition->GetParticleName() << " is ignored." << G4endl;
            G4Exception("MyParticleGun::SetParticleDefinition()", "Event0102",
                JustWarning, ED);
            return;
        }
    }
    particle_definition = aParticleDefinition;
    particle_charge = particle_definition->GetPDGCharge();
    if (particle_momentum > 0.0)
    {
        G4double mass = particle_definition->GetPDGMass();
        particle_energy =
            std::sqrt(particle_momentum * particle_momentum + mass * mass) - mass;
    }
}

void MyParticleGun::SetParticleEnergy(G4double aKineticEnergy)
{
    particle_energy = aKineticEnergy;
    if (particle_momentum > 0.0)
    {
        particle_momentum = 0.0;
    }
}

void MyParticleGun::SetParticleMomentum(G4double aMomentum)
{
    if (particle_definition == nullptr)
    {
        particle_momentum = aMomentum;
        particle_energy = aMomentum;
    }
    else
    {
        G4double mass = particle_definition->GetPDGMass();
        particle_momentum = aMomentum;
        particle_energy =
            std::sqrt(particle_momentum * particle_momentum + mass * mass) - mass;
    }
}

void MyParticleGun::SetParticleMomentum(G4ParticleMomentum aMomentum)
{
    if (particle_definition == nullptr)
    {
        particle_momentum_direction = aMomentum.unit();
        particle_momentum = aMomentum.mag();
        particle_energy = aMomentum.mag();
    }
    else
    {
        G4double mass = particle_definition->GetPDGMass();
        particle_momentum = aMomentum.mag();
        particle_momentum_direction = aMomentum.unit();
        particle_energy =
            std::sqrt(particle_momentum * particle_momentum + mass * mass) - mass;
    }
}

void MyParticleGun::GeneratePrimaryVertex(G4Event* evt)
{
    if (particle_definition == nullptr)
    {
        G4ExceptionDescription ED;
        ED << "Particle definition is not defined." << G4endl;
        ED << "MyParticleGun::SetParticleDefinition() has to be invoked beforehand."
            << G4endl;
        G4Exception("MyParticleGun::GeneratePrimaryVertex()", "Event0109",
            FatalException, ED);
        return;
    }

    // Create a new vertex
    //
    G4PrimaryVertex* vertex =
        new G4PrimaryVertex(particle_position, particle_time);

    // Create new primaries and set them to the vertex
    //
    G4double mass = particle_definition->GetPDGMass();
    for (G4int i = 0; i < NumberOfParticlesToBeGenerated; ++i)
    {
        G4PrimaryParticle* particle =
            new G4PrimaryParticle(particle_definition);
        particle->SetKineticEnergy(particle_energy);
        particle->SetMass(mass);
        particle->SetMomentumDirection(particle_momentum_direction);
        particle->SetCharge(particle_charge);
        particle->SetPolarization(particle_polarization.x(),
            particle_polarization.y(),
            particle_polarization.z());
        vertex->SetPrimary(particle);
    }
    evt->AddPrimaryVertex(vertex);
}
