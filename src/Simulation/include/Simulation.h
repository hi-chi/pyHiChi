#pragma once
#include "Grid.h"
#include "FieldEntity.h"
#include "Ensemble.h"
#include "Pusher.h"
#include <memory>

namespace pfc {

    class SimulationBase {
    public:
        virtual void save(std::ostream& ostr) = 0;
        virtual void load(std::istream& ostr) = 0;
        virtual void run() = 0;
    };

    template <class TGrid, class TFieldSolver, class TParticleArray = NoParticleArray>
    class Simulation : public SimulationBase {
    public:
        std::shared_ptr<FieldEntity<TGrid, TFieldSolver>> field;
        std::shared_ptr<Ensemble<TParticleArray>> ensemble;
        std::shared_ptr<BorisPusher> particlePusher;

        Simulation(const std::shared_ptr<FieldEntity<TGrid, TFieldSolver>>& field): field(field) {}
        Simulation(const std::shared_ptr<FieldEntity<TGrid, TFieldSolver>>& field,
                   const std::shared_ptr<Ensemble<TParticleArray>>& ensemble,
                   const std::shared_ptr<BorisPusher>& particlePusher): field(field), ensemble(ensemble), particlePusher(particlePusher) {}

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
        void run() override
        {
            if (particlePusher) (*particlePusher)(ensemble.get(), field->getGrid(), field->getFieldSolver()->dt);
        }
    protected:
    };
}