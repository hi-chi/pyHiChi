#include "Simulation.h"
#include "Fdtd.h"
#include <mpi.h>

namespace pfc {
    namespace ParticleInfo {
        std::vector<ParticleType> typesVector = { {constants::electronMass, constants::electronCharge},//electron
                                    {constants::electronMass, -constants::electronCharge},//positron
                                    {constants::protonMass, -constants::electronCharge},//proton
                                    {constants::electronMass, 0.0 } };//photon
        const ParticleType* types = &ParticleInfo::typesVector[0];
        short numTypes = sizeParticleTypes;
    }
}
    int main(int arc, char* argv[])
    {
        std::shared_ptr<FieldEntity<YeeGrid, FDTD>> ptrField = std::make_shared<FieldEntity<YeeGrid, FDTD>>(Int3(8, 9, 10), FP3(0.0, 0.0, 0.0), FP3(1.0, 1.0, 1.0), 2e-12);
        std::shared_ptr<Ensemble<ParticleArray3d>> ptrEnsemble = std::make_shared<Ensemble<ParticleArray3d>>();
        std::shared_ptr<BorisPusher> ptrPusher = std::make_shared<BorisPusher>();

        //for (int i = 0; i < 9; i++)
        //    ensemble->addParticle(this->randomParticle());
        //ensemble->addParticle(this->randomParticle(Photon));
        //ensemble->addParticle(this->randomParticle(Positron));
        //ensemble->addParticle(this->randomParticle(Proton));
        BaseSimulation* simulation = new Simulation<YeeGrid, FDTD, ParticleArray3d>(ptrField, ptrEnsemble, ptrPusher);
        return 0;
    }