#pragma once
#include "Vectors.h"
#include "macros.h"

#include <vector>

namespace pfc {

    class PmlSplitGrid
    {
    public:

        PmlSplitGrid(Int3 leftInnerCornerIndex, Int3 rightInnerCornerIndex,
            Int3 leftOuterCornerIndex, Int3 rightOuterCornerIndex);

        forceinline int getNumPmlNodes() const { return index.size(); }
        forceinline Int3 getIndex3d(int idx) { return index[idx]; }
        forceinline std::vector<Int3>& getIndices() { return index; }

        void save(std::ostream& ostr);
        void load(std::istream& istr);

        void resizeFields(int size);

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
        resizeFields(indexSize);
    }

    inline void PmlSplitGrid::resizeFields(int size) {
        bxy.resize(size, 0);
        bxz.resize(size, 0);
        byx.resize(size, 0);
        byz.resize(size, 0);
        bzx.resize(size, 0);
        bzy.resize(size, 0);
        
        exy.resize(size, 0);
        exz.resize(size, 0);
        eyx.resize(size, 0);
        eyz.resize(size, 0);
        ezx.resize(size, 0);
        ezy.resize(size, 0);
    }

    inline void PmlSplitGrid::save(std::ostream& ostr) {
        const int size = bxy.size();
        ostr.write((char*)&size, sizeof(size));

        ostr.write((char*)bxy.data(), sizeof(FP) * size);
        ostr.write((char*)bxz.data(), sizeof(FP) * size);
        ostr.write((char*)byx.data(), sizeof(FP) * size);
        ostr.write((char*)byz.data(), sizeof(FP) * size);
        ostr.write((char*)bzx.data(), sizeof(FP) * size);
        ostr.write((char*)bzy.data(), sizeof(FP) * size);

        ostr.write((char*)exy.data(), sizeof(FP) * size);
        ostr.write((char*)exz.data(), sizeof(FP) * size);
        ostr.write((char*)eyx.data(), sizeof(FP) * size);
        ostr.write((char*)eyz.data(), sizeof(FP) * size);
        ostr.write((char*)ezx.data(), sizeof(FP) * size);
        ostr.write((char*)ezy.data(), sizeof(FP) * size);

        ostr.write((char*)index.data(), sizeof(Int3) * size);
    }

    inline void PmlSplitGrid::load(std::istream& istr) {
        int size = 0;
        istr.read((char*)&size, sizeof(size));

        resizeFields(size);
        index.resize(size);

        istr.read((char*)bxy.data(), sizeof(FP) * size);
        istr.read((char*)bxz.data(), sizeof(FP) * size);
        istr.read((char*)byx.data(), sizeof(FP) * size);
        istr.read((char*)byz.data(), sizeof(FP) * size);
        istr.read((char*)bzx.data(), sizeof(FP) * size);
        istr.read((char*)bzy.data(), sizeof(FP) * size);

        istr.read((char*)exy.data(), sizeof(FP) * size);
        istr.read((char*)exz.data(), sizeof(FP) * size);
        istr.read((char*)eyx.data(), sizeof(FP) * size);
        istr.read((char*)eyz.data(), sizeof(FP) * size);
        istr.read((char*)ezx.data(), sizeof(FP) * size);
        istr.read((char*)ezy.data(), sizeof(FP) * size);

        istr.read((char*)index.data(), sizeof(Int3) * size);
    }
}
