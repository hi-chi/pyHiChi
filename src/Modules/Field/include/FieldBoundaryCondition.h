#pragma once
#include "Enums.h"
#include "Vectors.h"
#include "Grid.h"

namespace pfc
{

    template<class TGrid>
    class FieldBoundaryCondition
    {
    public:

        FieldBoundaryCondition(TGrid* grid, Int3 leftBorderIndex,
            Int3 rightBorderIndex, CoordinateEnum axis) :
            grid(grid), axis(axis),
            leftBorderIndex(leftBorderIndex), rightBorderIndex(rightBorderIndex)
        {}

        // constructor for loading
        FieldBoundaryCondition(TGrid* grid, Int3 leftBorderIndex,
            Int3 rightBorderIndex) : grid(grid),
            leftBorderIndex(leftBorderIndex), rightBorderIndex(rightBorderIndex),
            axis(CoordinateEnum::x)
        {}

        // polymorfic class
        virtual ~FieldBoundaryCondition() {}
        virtual FieldBoundaryCondition<TGrid>* createInstance(
            TGrid* grid, Int3 leftBorderIndex, Int3 rightBorderIndex, CoordinateEnum axis) = 0;

        virtual void generateB(FP time) = 0;
        virtual void generateE(FP time) = 0;

        virtual void save(std::ostream& ostr);
        virtual void load(std::istream& istr);

        TGrid* grid = nullptr;
        Int3 leftBorderIndex, rightBorderIndex;

        CoordinateEnum axis;
    };

    template<class TGrid>
    inline void FieldBoundaryCondition<TGrid>::save(std::ostream& ostr)
    {
        ostr.write((char*)&axis, sizeof(axis));
    }

    template<class TGrid>
    inline void FieldBoundaryCondition<TGrid>::load(std::istream& istr)
    {
        istr.read((char*)&axis, sizeof(axis));
    }
}
