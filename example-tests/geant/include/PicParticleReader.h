#pragma once
#include <vector>
#include "PicParticleContainer.h"

class PicParticleReader
{
public:

    static void ReadFromFile(const std::vector<std::string>& fileNames,
        PicParticleContainer& particles);

};

