#include "Simulation.h"
#include "Fdtd.h"
#include <mpi.h>

class MPIDomain : public BaseDomain
{
public:
    MPI_Request req;
    MPI_Status status;
    MPIDomain(int MPI_rank, int MPI_size, Int3 domainSize)
    {
        this->MPI_rank = MPI_rank;
        this->MPI_size = MPI_size;
        this->domainSize = domainSize;
        this->domainIndx = this->getInt3Rank();
        this->init();
    }
    MPIDomain(Int3 domainRank, int MPI_size, Int3 domainSize)
    {
        this->MPI_size = MPI_size;
        this->domainSize = domainSize;
        this->domainIndx = domainRank;
        this->MPI_rank = this->getLinearRank();
        this->init();
    }
    void init()
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
    void iSend(char* ar, int size, MPI_TAG tag) override
    {
        MPI_Isend(ar, size, MPI_CHAR, this->MPI_rank, tag, MPI_COMM_WORLD, &req);
    }

    void waitRecive(MPI_TAG tag) override {} // need add tag
    bool isMessage(MPI_TAG tag) override // need add tag
    {
        MPI_Iprobe(MPI_rank, tag, MPI_COMM_WORLD, &isMPImessage, &status);
        if (isMPImessage == 1)
        {
            isMPImessage = 0;

            return true;
        }
        return false;
    }
    void setMessage(char* ar, int& size, MPI_TAG tag) override 
    {
        MPI_Get_count(&status, MPI_CHAR, &size);
        ar = new char[size];
        MPI_Recv(ar, size, MPI_CHAR, MPI_rank, tag, MPI_COMM_WORLD, &status);
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
    FP urand(FP a, FP b) {
        return a + (b - a) * ((FP)rand()) / RAND_MAX;
    }
    Particle3d randomParticle(FP3 minPosition = { 0.0, 0.0, 0.0 }, FP3 maxPosition = { 1.0, 1.0, 1.0 }, ParticleTypes type = ParticleTypes::Electron)
    {
        FP3 position;
        for (int d = 0; d < 3; d++)
            position[d] = urand(minPosition[d], maxPosition[d]);
        Real minMomentum = -10.0;
        Real maxMomentum =  10.0;
        FP3 momentum(urand(minMomentum, maxMomentum),
            urand(minMomentum, maxMomentum), urand(minMomentum, maxMomentum));
        Real weight = static_cast<Real>(urand(1e-5, 1e5));
        return Particle3d(position, momentum, weight, type);
    }
}
    int main(int argc, char* argv[])
    {
        int size_x = 1, size_y = 1, mpi_size = 1, mpi_rank = 0;
        if (argc > 1)
        {
            size_x = std::stoi(argv[1]);
            size_y = std::stoi(argv[2]);
            MPI_Init(&argc, &argv);
            MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
            MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
            std::cout << mpi_size << " " << mpi_rank << "\n";

        }
        std::shared_ptr<FieldEntity<YeeGrid, FDTD>> ptrField = std::make_shared<FieldEntity<YeeGrid, FDTD>>(Int3(10, 10, 10), FP3(0.0, 0.0, 0.0), FP3(0.1, 0.1, 0.1), 2e-14);
        std::shared_ptr<Ensemble<ParticleArray3d>> ptrEnsemble = std::make_shared<Ensemble<ParticleArray3d>>();
        std::shared_ptr<BorisPusher> ptrPusher = std::make_shared<BorisPusher>();



        for (int i = 0; i < 10000; i++)
            ptrEnsemble->addParticle(randomParticle());
        BaseSimulation* simulation = new Simulation<YeeGrid, FDTD, ParticleArray3d>(ptrField, ptrEnsemble, ptrPusher); //BaseSimulation* simulation = new Simulation<YeeGrid, FDTD>(ptrField);
        for (int i = 0; i < 1; i++)
            simulation->runIteration();
        if (argc > 1) MPI_Finalize();
        return 0;
    }