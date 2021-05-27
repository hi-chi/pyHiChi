#pragma once
#include "Grid.h"
#include "FieldEntity.h"

#include <memory>

namespace pfc {

    template <class TGrid, class TFieldSolver>
    class Simulation {
    public:
        std::shared_ptr<FieldEntity<TGrid, TFieldSolver>> field;
        Simulation(const std::shared_ptr<FieldEntity<TGrid, TFieldSolver>>& field);

        void save(std::ostream& ostr);
        void load(std::istream& istr);

    protected:
    };

    template<class TGrid, class TFieldSolver>
    inline Simulation<TGrid, TFieldSolver>::Simulation(const std::shared_ptr<FieldEntity<TGrid, TFieldSolver>>& field)
    {
        this->field = field;
    }

    template<class TGrid, class TFieldSolver>
    inline void Simulation<TGrid, TFieldSolver>::save(std::ostream& ostr)
    {
        field->save(ostr);
    }

    template<class TGrid, class TFieldSolver>
    inline void Simulation<TGrid, TFieldSolver>::load(std::istream& istr)
    {
        field->load(istr);
    }
}