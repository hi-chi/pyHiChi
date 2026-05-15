#include "MPI_GridSlicer.h"
#include "MPI_Topology.h"

namespace mpi
{

// Slices a min gird into sections and outputs grid parameters to a sub-grid associated with the specified rank
// Grid divisions are specified in an array of std::vectors, divisions[0] - divisions along the X axis, [1] - Y, [2] - Z
void GridSlicer::getSubGridParameters(pfc::FP3& out_origin, pfc::Int3& out_size, 
                                      const pfc::FP3& main_origin, const pfc::Int3& main_size, const pfc::FP3& step,
                                      std::vector<int> divisions[3],
                                      int rank, const Topology& topology)
{
    pfc::Int3 node = topology.getNodeSection(rank);
    const pfc::Int3& sections = topology.getSections();

    pfc::Int3 lb;
    pfc::Int3 rb;

    lb.x = (node.x - 1 < 0) ? 0 : (divisions[0][node.x - 1]);
    rb.x = (node.x == sections.x - 1) ? main_size.x : (divisions[0][node.x]);

    lb.y = (node.y - 1 < 0) ? 0 : (divisions[1][node.y - 1]);
    rb.y = (node.y == sections.y - 1) ? main_size.y : (divisions[1][node.y]);

    lb.z = (node.z - 1 < 0) ? 0 : (divisions[2][node.z - 1]);
    rb.z = (node.z == sections.z - 1) ? main_size.z : (divisions[2][node.z]);

    // Outputs
    out_origin = main_origin + step * lb;
    out_size = rb - lb;
}

}