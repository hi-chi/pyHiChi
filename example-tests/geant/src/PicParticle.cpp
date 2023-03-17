#include "PicParticle.h"

#include "G4ParticleTable.hh"

using namespace pfc;

void PicParticle::SetType(std::string stype)
{
    if (stype == "Electron") type = ParticleType::Electron;
    if (stype == "Positron") type = ParticleType::Positron;
    if (stype == "Photon") type = ParticleType::Photon;
}

FP PicParticle::GetEnergy() const {
    if (type == ParticleType::Photon) {
        return momentum.norm() * constants::c;
    }
    if (type == ParticleType::Electron || type == ParticleType::Positron) {
        FP p0 = mass * constants::c;
        return sqrt(1 + sqr(momentum.norm() / p0)) * p0 * constants::c;
    }
    return momentum.norm() * constants::c * constants::c / velocity.norm();
}

FP3 PicParticle::GetDirection() const
{
    return momentum / momentum.norm();
}

G4ParticleDefinition* PicParticle::GetParticleDefinition() const {
    G4String particleName;
    if (type == ParticleType::Electron) {
        particleName = "e-";
    }
    else if (type == ParticleType::Positron) {
        particleName = "e+";
    }
    else if (type == ParticleType::Photon) {
        particleName = "gamma";
    }
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    return particleTable->FindParticle(particleName);
}
