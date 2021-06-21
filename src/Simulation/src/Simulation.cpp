#include "Simulation.h"
#include "Fdtd.h"
#include <mpi.h>
#include <ctime>
class MPIDomain : public BaseDomain
{
public:
    static int realRank;
    MPI_Request req;
    MPI_Status status;
    MPIDomain(int MPI_rank, int MPI_size, Int3 domainSize)
    {
        this->MPI_rank = MPI_rank;
        this->MPI_size = MPI_size;
        this->domainSize = domainSize;
        this->domainIndx = this->getInt3Rank();
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
        int neib = 0;
        for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        for (int k = 0; k < 3; k++)
        {
            if (i == 1 && j == 1 && k == 1) neighbors[i][j][k] = new SingleDomain();
            else 
            {
                Int3 domainRank(domainIndx.x + i - 1, domainIndx.y + j - 1, domainIndx.z + k - 1);
                if (isInside(domainRank))
                {
                    neighbors[i][j][k] = new MPIDomain(domainRank, MPI_size, domainSize);
                    neib++;
                }
                else neighbors[i][j][k] = new NullDomain();
            }
        }
    }
    void iSend(char* ar, int size, MPI_TAG tag) override
    {
        //std::cout << "send from " << MPIDomain::realRank << " to " << MPI_rank << " size " << size << endl;
        MPI_Isend(ar, size, MPI_CHAR, MPI_rank, tag, MPI_COMM_WORLD, &req);
    }
    void barrier() override
    {
        MPI_Barrier(MPI_COMM_WORLD);
    }
    void Recive(char*& ar, int& size, MPI_TAG tag) override
    {
        int tmp = 0;
        //std::cout << "try_prob from " << MPIDomain::realRank << " to " << MPI_rank << " size " << size << endl;
        while (tmp != 1) { MPI_Iprobe(MPI_rank, tag, MPI_COMM_WORLD, &tmp, &status);}
        //std::cout << "suc_prob from " << MPIDomain::realRank << " to " << MPI_rank << " size " << size << endl;
        MPI_Get_count(&status, MPI_CHAR, &size);
        //std::cout << "recive from " << MPIDomain::realRank << " to " << MPI_rank << " size " << size << endl;
        ar = new char[size];
        MPI_Recv(ar, size, MPI_CHAR, MPI_rank, tag, MPI_COMM_WORLD, &status);
    }
    void waitRecive(MPI_TAG tag) override
    {
        isMPImessage = 0;
    } // need add tag
    bool isMessage(MPI_TAG tag) override // need add tag
    {
        if (isMPImessage)
            return false;
        MPI_Iprobe(MPI_rank, tag, MPI_COMM_WORLD, &isMPImessage, &status);
        if (isMPImessage == 1)
            return true;
        return false;
    }
    int setMessage(char*& ar, MPI_TAG tag) override
    {
        int size = 0;
        MPI_Get_count(&status, MPI_CHAR, &size);
        ar = new char[size];
        MPI_Recv(ar, size, MPI_CHAR, MPI_rank, tag, MPI_COMM_WORLD, &status);
        return size;
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
        Real minMomentum = -1e-18; // -10
        Real maxMomentum =  1e-18; // 10
        FP3 momentum(urand(minMomentum, maxMomentum),
            urand(minMomentum, maxMomentum), urand(minMomentum, maxMomentum));
        Real weight = static_cast<Real>(urand(1e-5, 1e5));
        return Particle3d(position, momentum, weight, type);
    }
}
int MPIDomain::realRank = 0;
    int main(int argc, char* argv[])
    {
        int size_x = 1, size_y = 1, mpi_size = 1, mpi_rank = 0;
        FP3 start = FP3(0.0, 0.0, 0.0);
        FP3 end = FP3(1.0, 1.0, 1.0);
        FP3 step = FP3(0.1, 0.1, 0.1);
        Int3 cells = (end - start) / step, size_domain = cells, mpi_indx = Int3(0, 0, 0);
        if (argc > 1)
        {
            size_x = std::stoi(argv[1]);
            size_y = std::stoi(argv[2]);
            size_domain.x /= size_x;
            size_domain.y /= size_y;
            MPI_Init(&argc, &argv);
            MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
            MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
            std::cout << "rank " << mpi_rank << " size " << mpi_size << endl;
            mpi_indx.x = mpi_rank % size_x;
            mpi_indx.y = mpi_rank / size_x;
            start = size_domain * mpi_indx * step;
            end = start + size_domain * step;
            MPIDomain::realRank = mpi_rank;
        }
        std::shared_ptr<FieldEntity<YeeGrid, FDTD>> ptrField = std::make_shared<FieldEntity<YeeGrid, FDTD>>(size_domain, start, step, 2e-14);
        std::shared_ptr<Ensemble<ParticleArray3d>> ptrEnsemble = std::make_shared<Ensemble<ParticleArray3d>>();
        std::shared_ptr<BorisPusher> ptrPusher = std::make_shared<BorisPusher>();
        for (int i = 0; i < size_domain.x*size_domain.y * 1000; i++)
            ptrEnsemble->addParticle(randomParticle(start, end));

        BaseSimulation* simulation = new Simulation<YeeGrid, FDTD, ParticleArray3d>(ptrField, ptrEnsemble, ptrPusher); //BaseSimulation* simulation = new Simulation<YeeGrid, FDTD>(ptrField);
        simulation->domain.reset(new MPIDomain(mpi_rank, mpi_size, Int3(size_x, size_y, 1)));
        //simulation->domain.reset(new SingleDomain());
        simulation->domain->init();

        int t = clock();
        for (int i = 0; i < 1000; i++)
            simulation->runIteration();
        std::cout << "time " << (clock() - t) << endl;

        if (argc > 1) MPI_Finalize();
        return 0;
    }