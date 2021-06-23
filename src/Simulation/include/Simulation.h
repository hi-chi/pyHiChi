#pragma once
#include "Grid.h"
#include "FieldEntity.h"
#include "Ensemble.h"
#include "Pusher.h"
#include <memory>
#include <fstream>

namespace pfc {
    enum MPI_TAG
    {
        TAG_PARTICLE,
        TAG_FIELD
    };
    struct Logger
    {
        std::unique_ptr<std::ofstream> pLogger;
        Logger()
        {
            pLogger.reset(new std::ofstream("log.txt"));
        }
        std::ostream& operator*() { return *pLogger; }
    };

    class BaseDomain
    {
    public:
        Int3 domainIndx = Int3(0, 0, 0);
        Int3 domainSize = Int3(1, 1, 1);
        int MPI_rank = 0;
        int MPI_size = 1;
        BaseDomain* neighbors[3][3][3];

        std::vector<BaseDomain*> jaggedNeighbors;
        int getLinearRank()
        {
            return domainSize.x * domainSize.y * domainIndx.z + domainSize.x * domainIndx.y + domainIndx.x;
        }
        Int3 getInt3Rank()
        {
            return Int3((MPI_rank % (domainSize.x * domainSize.y)) % domainSize.x, (MPI_rank % (domainSize.x * domainSize.y)) / domainSize.x, MPI_rank / (domainSize.x * domainSize.y));
        }
        Int3 getInt3Rank(int rank)
        {
            return Int3(rank % );
        }
        bool isInside(Int3 indx) { return indx >= Int3(0, 0, 0) && indx < domainSize; }
        virtual void init() {};
        virtual void setLogger(Logger& logger) {}

        virtual bool isMessage(MPI_TAG tag) // need add tag to bool
        {
            if (isNewMessage == 1)
            {
                isNewMessage = 0;
                return true;
            }
            return false;
        }
        virtual void iSend(char* ar, int size, MPI_TAG tag){ isNewMessage = 1; } // need add tag
        virtual int getMessage(char*& ar, MPI_TAG tag) { return 0; }
        virtual void barrier(){}
        virtual void finalize(){}
    protected:
        int isNewMessage = 0;
    };
    class NullDomain : public BaseDomain
    {
    public:
    };
    class SingleDomain : public BaseDomain
    {
    public:
        SingleDomain()
        {
            init();
        }
        void init()
        {
            for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                neighbors[i][j][k] = new NullDomain();
        }
        ~SingleDomain()
        {
            for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                delete neighbors[i][j][k];
        }
    };

    class BaseSimulation {
    public:
        size_t curIteration = 0;
        size_t numIteration = 0;
        std::unique_ptr<BaseDomain> domain = 0;
        Logger logger = Logger();

        virtual void save(std::ostream& ostr) = 0;
        virtual void load(std::istream& ostr) = 0;
        virtual void runIteration() = 0;
        virtual void run() = 0;
    };

    template <class TGrid, class TFieldSolver, class TParticleArray = NoParticleArray>
    class Simulation : public BaseSimulation {
    public:
        std::shared_ptr<FieldEntity<TGrid, TFieldSolver>> field;
        std::shared_ptr<Ensemble<TParticleArray>> ensemble;
        std::shared_ptr<BorisPusher> particlePusher;

        Simulation(const std::shared_ptr<FieldEntity<TGrid, TFieldSolver>>& field) : field(field) {
            domain.reset(new SingleDomain());
        }
        Simulation(const std::shared_ptr<FieldEntity<TGrid, TFieldSolver>>& field,
                   const std::shared_ptr<Ensemble<TParticleArray>>& ensemble,
                   const std::shared_ptr<BorisPusher>& particlePusher): field(field), ensemble(ensemble), particlePusher(particlePusher) {
            domain.reset(new SingleDomain());
        }

        void save(std::ostream& ostr) override
        {
            field->save(ostr);
            if (ensemble) ensemble->save(ostr);
        }
        void load(std::istream& istr) override
        {
            field->load(istr);
            if (ensemble) ensemble->load(istr);
        }
        int sizeParts[3][3][3] = {0};
        void runIteration() override
        {
            if (domain->MPI_rank == 0)
                *logger << curIteration << endl;
            *logger << "particle in " << domain->MPI_rank << " " << ensemble->size() << endl;
            field->getFieldSolver()->updateFields();
            if (particlePusher)
            {
                for (int i = 0; i < pfc::sizeParticleTypes; i++)
                {
                    TParticleArray* particleArr = &(*ensemble)[i];
                    particlePusher->operator()(particleArr, field->getGrid(), field->getFieldSolver()->dt);
                    std::vector<Particle3d> migratedParticles[3][3][3];

                    particleArr->cutMigratedParticles(migratedParticles, field->getGrid()->getMinInternalCoord(), field->getGrid()->getMaxInternalCoord());
                    for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                    for (int k = 0; k < 3; k++)
                    {
                        char* ar = (char*)&migratedParticles[i][j][k];
                        const int size = migratedParticles[i][j][k].size() * sizeof(Particle3d);
                        domain->neighbors[i][j][k]->iSend(ar, size, MPI_TAG::TAG_PARTICLE);
                        sizeParts[i][j][k] += migratedParticles[i][j][k].size();
                    }
                    BaseDomain** tmpNeighbors = (BaseDomain**)domain->neighbors;
                    for (int count = 0; count < 27;)
                    for (int i = 0; i < 3*3*3; i++)
                    {
                        if (tmpNeighbors[i]->isMessage(MPI_TAG::TAG_PARTICLE))
                        {
                            count++;
                            char* ar = 0;
                            int byte_size = tmpNeighbors[i]->getMessage(ar, MPI_TAG::TAG_PARTICLE);
                            Particle3d* parts = (Particle3d*)ar;
                            *logger << "Start recieve " << domain->getInt3Rank(i) << endl;
                            if (byte_size > 0)
                            {
                                *logger << "Recieve from: MPI_rank " << domain->MPI_rank << " byte_size " << byte_size << endl;
                                for (int j = 0; j < byte_size / sizeof(Particle3d); j++)
                                    particleArr->pushBack(parts[j]);
                                delete ar;
                            }
                            *logger << "End recieve" << endl;
                        }
                    }
                    domain->barrier();
                }
            }
            curIteration++;
        }
        void run()
        {
            for (size_t i = 0; i < numIteration; i++)
            {
                runIteration();
            }
        }
    protected:
    };
}