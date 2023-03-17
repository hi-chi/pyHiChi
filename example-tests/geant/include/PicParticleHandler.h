#pragma once
#include "PicParticle.h"
#include "PicParticleContainer.h"


// �����, ���������� �� ��������������� ��������� ������, ���������� �� Picador
class PicParticleHandler {
public:

    // ����������� �������, ������� �� ������� � ������� ������� Geant
    // ���������� ������� ������� � ��� ��������� �������
    void FilterAndMoveParticles(PicParticleContainer& particles,
        pfc::FP3 minGeantBox, pfc::FP3 maxGeantBox, pfc::FP timeStep) {
        FilterParticles(particles, minGeantBox, maxGeantBox, timeStep);
    }

private:

    void FilterParticles(PicParticleContainer& particles,
        pfc::FP3 minBox, pfc::FP3 maxBox, pfc::FP timeStep);
    void MoveParticle(PicParticle* particle, pfc::FP timeStep);
};
