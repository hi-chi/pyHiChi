#pragma once

#include "FieldBoundaryCondition.h"
#include "MPI_FieldExchanger.h"
#include "MPI_Topology.h"

namespace pfc
{
    template <class GridType>
    class MPI_FieldBoundaryCondition : public FieldBoundaryCondition<GridType>
    {
    protected:
        std::shared_ptr<mpi::FieldExchanger> exchanger;
        std::shared_ptr<mpi::Topology> topology;
        int rank;

    void checkCommunicator()
    {
        if (!topology->isValidOnThisRank())
            throw std::runtime_error("Attempted to create MPI_FieldBoundaryCondition on rank, that does not participate in topology!");
    }

    public:

        MPI_FieldBoundaryCondition(GridType* grid,
            Int3 leftBorderIndex, Int3 rightBorderIndex, CoordinateEnum axis, 
            std::shared_ptr<mpi::FieldExchanger> fieldExchanger,
            std::shared_ptr<mpi::Topology> topology, int mpi_rank) :

            FieldBoundaryCondition(grid, leftBorderIndex, rightBorderIndex, axis), 
            exchanger(fieldExchanger), topology(topology), rank(mpi_rank)
        { 
            checkCommunicator(); 
        }

        // constructor for loading
        explicit MPI_FieldBoundaryCondition(GridType* grid,
            Int3 leftBorderIndex, Int3 rightBorderIndex, 
            std::shared_ptr<mpi::FieldExchanger> fieldExchanger,
            std::shared_ptr<mpi::Topology> topology, int mpi_rank) :

            FieldBoundaryCondition(grid, leftBorderIndex, rightBorderIndex), 
            exchanger(fieldExchanger), topology(topology), rank(mpi_rank)
        {
            checkCommunicator();
        }

        void generateB(FP time) override;
        void generateE(FP time) override;

        FieldBoundaryCondition<GridType>* createInstance(
            GridType* grid, Int3 leftBorderIndex, Int3 rightBorderIndex, CoordinateEnum axis) override {
            return new MPI_FieldBoundaryCondition(grid, leftBorderIndex, rightBorderIndex, axis, exchanger, topology, rank);
        }
    };

    template <class GridType>
    inline void MPI_FieldBoundaryCondition<GridType>::generateE(FP time)
    {
        exchanger->performExchangeOverAxis(grid->Ex.getData(), axis, *topology, rank, topology->getTopologyCommunicator());
        exchanger->performExchangeOverAxis(grid->Ey.getData(), axis, *topology, rank, topology->getTopologyCommunicator());
        exchanger->performExchangeOverAxis(grid->Ez.getData(), axis, *topology, rank, topology->getTopologyCommunicator());
    }

    template <class GridType>
    inline void MPI_FieldBoundaryCondition<GridType>::generateB(FP time)
    {
        exchanger->performExchangeOverAxis(grid->Bx.getData(), axis, *topology, rank, topology->getTopologyCommunicator());
        exchanger->performExchangeOverAxis(grid->By.getData(), axis, *topology, rank, topology->getTopologyCommunicator());
        exchanger->performExchangeOverAxis(grid->Bz.getData(), axis, *topology, rank, topology->getTopologyCommunicator());
    }
}