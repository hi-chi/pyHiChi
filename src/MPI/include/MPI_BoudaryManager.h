#pragma once

#include "MPI_Topology.h"
#include "MPI_BoundaryCondition.h"
#include "MPI_FieldExchanger.h"

#include "FieldSolver.h"

#include <type_traits>

namespace mpi
{

enum class BoundaryType
{
    Periodic, // Does nothing, Periodic conditions are handled automatically through MPI_BoundaryCondition
    Reflect,
    // ...
};

// Sets up field boundary conditions for subgrids
// Template parameters specify FieldBoundaryCondition decendants that will be used
// IMPORTANT: specified (in template parameters) FieldBoundaryCondition classes must work for a SINGLE DIRECTIOn over a single axis!
// Such classes are usually have "MonoDirection" in their name
template<   class FieldSolverType, class GridType,
            class MonoDirectionReflectBoundaryConditionType>
class FieldBoundaryManager
{
private:
    FieldBoundaryManager() {}

public:

    // Setup field boundary conditions for subgrids, based on global grid boundary requirements
    // global_bounaries are specified in the following order: +X, -X, +Y, -Y, +Z, -Z
    static void setupBoundaryConditions(std::unique_ptr<FieldSolverType>& fieldSolver, BoundaryType global_bounaries[6], 
                                        std::shared_ptr<mpi::Topology> topology, std::shared_ptr<mpi::FieldExchanger> fieldExchanger, int rank)
    {
        for (int dir = 0; dir < 6; dir++)
        {
            Direction direction = static_cast<Direction>(dir);
            int neighbor = topology->getNeighbor(rank, direction);

            int axis = dir / 2;
            bool first = !static_cast<bool>( dir % 2 );
                
            std::unique_ptr<pfc::FieldBoundaryCondition<GridType>>& boundary_ptr = 
                (first) ? (fieldSolver->boundaryConditions[axis].first) : (fieldSolver->boundaryConditions[axis].second);

            // MPI_FieldBoundaryConditions
            if (neighbor != MPI_INVALID_RANK)
            {
                Direction inverted_direction = invertDirection(direction);
                int inverted_neighbor = topology->getNeighbor(rank, inverted_direction);

                std::unique_ptr<pfc::FieldBoundaryCondition<GridType>>& inverted_boundary_ptr = 
                    (!first) ? (fieldSolver->boundaryConditions[axis].first) : (fieldSolver->boundaryConditions[axis].second);

                // Only create a new one when the opposite direction does not have one
                // or in cases when the opposite direction has a diferent boundary condition
                if (inverted_boundary_ptr.get() == nullptr || inverted_neighbor == MPI_INVALID_RANK)
                {    
                    boundary_ptr.reset(new pfc::MPI_FieldBoundaryCondition<GridType>(fieldSolver->grid, fieldSolver->domainIndexBegin, 
                                                                                     fieldSolver->domainIndexEnd, static_cast<pfc::CoordinateEnum>(axis), 
                                                                                     fieldExchanger, topology, rank));
                }
            }

            // Other boundary condition
            else
            {
                pfc::SideEnum direction = (first) ? pfc::SideEnum::RIGHT : pfc::SideEnum::LEFT;

                switch (global_bounaries[dir])
                {
                case BoundaryType::Periodic:
                    throw(std::runtime_error("MPI boundary condition does not match topoplogy looping configuration!"));
                    break;
                case BoundaryType::Reflect:
                    if constexpr (!std::is_same_v<MonoDirectionReflectBoundaryConditionType, void>) // For spectral solvers that do not have reflective boundary conditions
                        boundary_ptr.reset(new MonoDirectionReflectBoundaryConditionType(fieldSolver->grid, fieldSolver->domainIndexBegin, 
                                                                fieldSolver->domainIndexEnd, static_cast<pfc::CoordinateEnum>(axis), direction));
                default: break;
                }
            }
        }
    }
};

}