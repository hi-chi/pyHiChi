#pragma once
#include <vector>
#include <string>
#include <iostream>

#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

#include "Vectors.h"
#include "Constants.h"

class G4ParticleDefinition;

const pfc::FP erg = MeV / 1.602176634e-6;

enum class ParticleType {
    Electron,
    Positron,
    Photon
};

// интегрированная из Picador частица
struct PicParticle {
    size_t id = 0;  // идентификатор
    ParticleType type;  // тип частицы
    pfc::FP emissionTime = 0.;  // время испускания = время вылета за границу расчетной области Picador
    // физические параметры частицы
    pfc::FP mass = 0., charge = 0., factor = 0.;
    pfc::FP3 position, velocity, momentum;

    void SetType(std::string stype);

    pfc::FP GetEnergy() const;
    pfc::FP3 GetDirection() const;

    G4ParticleDefinition* GetParticleDefinition() const;
};
