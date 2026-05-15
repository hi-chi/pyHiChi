#pragma once

#include <mpi.h>
#include <vector>

#include "FP.h"
#include "Vectors.h"
#include "MPI_Utilities.h"


namespace mpi
{
// Wrapper class for mpi field utility functions, CAN NOT BE INSTANCED!
// All methods are static
class FieldUtils final
{
private:
    // Private constructor to prevent instancing
    FieldUtils() {}

private:
    // Defines an mpi data type for transmitting sub arrays of the grid
    static void defineSubArrayType(MPI_Datatype &out_type, int sizes[3], int sub_sizes[3], int starts[3]);

    // Resolving "send" operation parameters based on direction
    static void resolveSendParameters(int out_sizes[3], int out_sub_sizes[3], int out_starts[3], 
        const pfc::Int3& num_cells, int num_external_cells, Direction direction);

    // Resolving "recv" operation parameters based on direction
    static void resolveRecvParameters(int out_sizes[3], int out_sub_sizes[3], int out_starts[3], 
        const pfc::Int3& num_cells, int num_external_cells, Direction direction);

public:

    // Determines mpi data type required to send field data in the required direction
    static void defineTransmission_Send(MPI_Datatype& out_type, pfc::Int3 grid_num_cells, int num_external_cells, Direction direction);

     // Determines mpi data type required to received field data from the required direction (direction relative to the sender)
    static void defineTransmission_Recv(MPI_Datatype& out_type, pfc::Int3 grid_num_cells, int num_external_cells, Direction direction);
};
}