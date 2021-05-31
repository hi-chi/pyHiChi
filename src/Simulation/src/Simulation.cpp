#include "Simulation.h"
#include "Fdtd.h"
#include <mpi.h>

class MPIDomain : public BaseDomain
{
public:
    MPIDomain(int MPI_rank, int MPI_size, Int3 domainSize)
    {
        this->MPI_rank = MPI_rank;
        this->MPI_size = MPI_size;
        this->domainSize = domainSize;
        this->domainIndx = this->getInt3Rank();;
    }
    MPIDomain(Int3 domainRank, int MPI_size, Int3 domainSize)
    {
        this->MPI_size = MPI_size;
        this->domainSize = domainSize;
        this->domainIndx = domainRank;
        this->MPI_rank = this->getLinearRank();
    }
    void init() override
    {
        for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        for (int k = 0; k < 3; k++)
        {
            if (i == 1 && j == 1 && k == 1) neighbors[i][j][k] = this;
            else 
            {
                Int3 domainRank(domainIndx.x + i - 1, domainIndx.y + j - 1, domainIndx.z + k - 1);
                if (isInside(domainRank)) neighbors[i][j][k] = new MPIDomain(domainRank, MPI_size, domainSize);
                else neighbors[i][j][k] = new NullDomain();
            }
        }
    }
};
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
        BaseSimulation* simulation = new Simulation<YeeGrid, FDTD, ParticleArray3d>(ptrField, ptrEnsemble, ptrPusher); //BaseSimulation* simulation = new Simulation<YeeGrid, FDTD>(ptrField);
        return 0;
    }