#include "MPI_Utilities.h"

namespace mpi
{

// Returns the inverse of the specified direction
Direction invertDirection(const Direction& direction)
{
    switch (direction)
    {
    case Direction::pX: return Direction::nX;
    case Direction::nX: return Direction::pX;
    case Direction::pY: return Direction::nY;
    case Direction::nY: return Direction::pY;
    case Direction::pZ: return Direction::nZ;
    case Direction::nZ: return Direction::pZ;

    default: return Direction::pX;
    }
}

}