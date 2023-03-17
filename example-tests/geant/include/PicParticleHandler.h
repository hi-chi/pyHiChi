#pragma once
#include "PicParticle.h"
#include "PicParticleContainer.h"


// класс, отвечающий за предварительную обработку частиц, приходящих из Picador
class PicParticleHandler {
public:

    // отбрасывает частицы, которые не попадут в область расчета Geant
    // оставшиеся частицы двигает в эту расчетную область
    void FilterAndMoveParticles(PicParticleContainer& particles,
        pfc::FP3 minGeantBox, pfc::FP3 maxGeantBox, pfc::FP timeStep) {
        FilterParticles(particles, minGeantBox, maxGeantBox, timeStep);
    }

private:

    void FilterParticles(PicParticleContainer& particles,
        pfc::FP3 minBox, pfc::FP3 maxBox, pfc::FP timeStep);
    void MoveParticle(PicParticle* particle, pfc::FP timeStep);
};
