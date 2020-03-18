#pragma once
#include <vector>
#include <string>

#include "Dimension.h"
#include "Particle.h"

namespace pfc {
    enum ParticleTypes {
        Electron = 0, 
        Positron = 1, 
        Proton = 2,
        Photon = 3
    };
    const int sizeParticleTypes = 4;
    const vector<std::string> particleNames = { "Electron", "Positron", "Proton", "Photon" };


    //void createTypes();

} // namespace pfc


