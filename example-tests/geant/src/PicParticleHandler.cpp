#include "PicParticleHandler.h"

#include <string>
#include <vector>
#include <fstream>
#include <algorithm>

using namespace pfc;


void PicParticleHandler::MoveParticle(PicParticle* particle, FP timeStep)
{
    // TODO: implement Boris pusher
    particle->position = particle->position + particle->velocity * timeStep;
    particle->emissionTime += timeStep;
}

void PicParticleHandler::FilterParticles(PicParticleContainer& particles,
    FP3 minBox, FP3 maxBox, FP timeStep)
{
    int index = 0;
    while (index < particles.size()) {
        PicParticle* p = particles[index];

        FP3 minParticleSearchBox(std::min(minBox.x, p->position.x),
            std::min(minBox.y, p->position.y), std::min(minBox.z, p->position.z));
        FP3 maxParticleSearchBox(std::max(maxBox.x, p->position.x),
            std::max(maxBox.y, p->position.y), std::max(maxBox.z, p->position.z));
        minParticleSearchBox -= minParticleSearchBox * 0.05;
        maxParticleSearchBox += maxParticleSearchBox * 0.05;

        bool filtered = false;
        // ���� �� �� ������, ��� ������� �������������� �� ������� � (minBox, maxBox)
        while (minParticleSearchBox <= p->position && p->position < maxParticleSearchBox) {
            // ���� ���-���� ������ � (minBox, maxBox), �� �������� ��
            if (minBox <= p->position && p->position < maxBox) {
                filtered = true;
                break;
            }
            // ����� ������ ��� pusher'�
            MoveParticle(p, timeStep);
        }

        if (!filtered) {  // ���� �� ��������, �� ������� �� �� �������
            particles.Pop(index);
        }
        else {
            index++;
        }
    }
}
