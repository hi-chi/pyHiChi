#pragma once
#include "Vectors.h"
#include "macros.h"

#include <vector>

namespace pfc {

    class PmlSplitGrid
    {
    public:

        PmlSplitGrid() {}
        PmlSplitGrid(Int3 leftInnerCornerIndex, Int3 rightInnerCornerIndex,
            Int3 leftOuterCornerIndex, Int3 rightOuterCornerIndex);

        forceinline int getNumPmlNodes() const { return index.size(); }
        forceinline Int3 getIndex3d(int idx) { return index[idx]; }
        forceinline std::vector<Int3>& getIndices() { return index; }

        std::vector<FP> bxy, bxz, byx, byz, bzx, bzy;  // split magnetic field
        std::vector<FP> exy, exz, eyx, eyz, ezx, ezy;  // split electric field
        // first index (x, y, z) is component, second one is propagation direction

        std::vector<Int3> index;  // natural 3d indexes of nodes in PML
        // indices go through area from left corner to right corner

    };

    inline PmlSplitGrid::PmlSplitGrid(
        Int3 leftInnerCornerIndex, Int3 rightInnerCornerIndex,
        Int3 leftOuterCornerIndex, Int3 rightOuterCornerIndex)
    {
        const Int3 begin = leftOuterCornerIndex;
        const Int3 end = rightOuterCornerIndex;

        for (int i = begin.x; i < end.x; i++)
            for (int j = begin.y; j < end.y; j++)
                for (int k = begin.z; k < end.z; k++)
                {
                    bool xBoundaryPml = (i < leftInnerCornerIndex.x) || (i >= rightInnerCornerIndex.x);
                    bool yBoundaryPml = (j < leftInnerCornerIndex.y) || (j >= rightInnerCornerIndex.y);
                    bool zBoundaryPml = (k < leftInnerCornerIndex.z) || (k >= rightInnerCornerIndex.z);
                    if (xBoundaryPml || yBoundaryPml || zBoundaryPml)
                    {
                        index.push_back(Int3(i, j, k));
                    }
                }

        const int indexSize = index.size();

        bxy.resize(indexSize, 0);
        bxz.resize(indexSize, 0);
        byx.resize(indexSize, 0);
        byz.resize(indexSize, 0);
        bzx.resize(indexSize, 0);
        bzy.resize(indexSize, 0);

        exy.resize(indexSize, 0);
        exz.resize(indexSize, 0);
        eyx.resize(indexSize, 0);
        eyz.resize(indexSize, 0);
        ezx.resize(indexSize, 0);
        ezy.resize(indexSize, 0);
    }
}
