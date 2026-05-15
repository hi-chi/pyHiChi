#pragma once
#include "Vectors.h"
#include <vector>

namespace mpi
{

// Wrapper class for mpi utility functions related to dividing a grid between nodes
// Can not be instanced
// All methods are static
class GridSlicer final
{
private:
    GridSlicer() {}

private:
    // ...

public:

    // Slices a min gird into sections and outputs grid parameters to a sub-grid associated with the specified rank
    // Grid divisions are specified in an array of std::vectors, divisions[0] - divisions along the X axis, [1] - Y, [2] - Z
    static void getSubGridParameters(pfc::FP3& out_origin, pfc::Int3& out_size, 
                                     const pfc::FP3& main_origin, const pfc::Int3& main_size, const pfc::FP3& step,
                                     std::vector<int> divisions[3],
                                     int rank, const class Topology& topology);
};

}