#pragma once
#include "Grid.h"
#include "FieldEntity.h"
#include "Ensemble.h"
#include "Pusher.h"
#include <memory>

namespace pfc {
    enum MPI_TAG
    {
        TAG_PARTICLE,
        TAG_FIELD
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
        int getLinearRank(Int3 indx)
        {
            return domainSize.x * domainSize.y * indx.z + domainSize.x * indx.y + indx.x;
        }
        int getLinearRank()
        {
            return domainSize.x * domainSize.y * domainIndx.z + domainSize.x * domainIndx.y + domainIndx.x;
        }
        Int3 getInt3Rank(int indx)
        {
            return Int3((indx % (domainSize.x * domainSize.y)) % domainSize.x, (indx % (domainSize.x * domainSize.y)) / domainSize.x, indx / (domainSize.x * domainSize.y));
        }
        Int3 getInt3Rank()
        {
            return Int3((MPI_rank % (domainSize.x * domainSize.y)) % domainSize.x, (MPI_rank % (domainSize.x * domainSize.y)) / domainSize.x, MPI_rank / (domainSize.x * domainSize.y));
        }
        bool isInside(Int3 indx) { return indx >= Int3(0, 0, 0) && indx < domainSize; }
        void init() {};

        virtual void waitRecive(MPI_TAG tag) { isMPImessage = 1; } // need add tag
        virtual bool isMessage(MPI_TAG tag) // need add tag
        {
            if (isMPImessage == 1)
            {
                isMPImessage = 0;
                return true;
            }
            return false;
        }

        virtual void iSend(char* ar, int size, MPI_TAG tag){}
        virtual void iRecive(char* ar, int size, MPI_TAG tag) {}
        virtual void setMessage(char* ar, int &size, MPI_TAG tag) {}
        virtual void finalize(){}
    protected:
        int isMPImessage = 0;
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
            this->init();
        }
        void init()
        {
            for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
            {
                if (i == 1 && j == 1 && k == 1) neighbors[i][j][k] = this;
                else neighbors[i][j][k] = new NullDomain();
            }
        }
        ~SingleDomain()
        {
            for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                if (!(i == 1 && j == 1 && k == 1)) delete neighbors[i][j][k];
        }
    };

    class BaseSimulation {
    public:
        size_t curIteration = 0;
        size_t numIteration = 0;
        std::unique_ptr<BaseDomain> domain = 0;

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
                        domain->neighbors[i][j][k]->waitRecive(MPI_TAG::TAG_PARTICLE);
                        sizeParts[i][j][k] += migratedParticles[i][j][k].size();
                    }
                    BaseDomain** tmpNeighbors = (BaseDomain**)domain->neighbors;
                    int count = 0;
                    while (count != 27)
                    for (int i = 0; i < 3*3*3; i++)
                    {
                        if (tmpNeighbors[i]->isMessage(MPI_TAG::TAG_PARTICLE))
                        {
                            count++;
                            int size = 0;
                            char* ar = 0;
                            tmpNeighbors[i]->setMessage(ar, size, MPI_TAG::TAG_PARTICLE);
                            Particle3d* parts = (Particle3d*)ar;
                            // copy ar to Ensemble
                            for (int i = 0; i < size / sizeof(Particle3d); i++)
                                particleArr->pushBack(*parts);
                            if (ar) delete ar;
                        }
                    }
                }

            }
            std::cout << curIteration << " ";
            for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                if (sizeParts[i][j][k] > 0) std::cout << Int3(i, j, k) << " " << sizeParts[i][j][k] << " ";
            std::cout << "\n";
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