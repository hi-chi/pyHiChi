#pragma once

namespace pfc {
    enum GridTypes {
        YeeGridType = 0,
        StraightGridType = 1,
        PSTDGridType = 2,
        PSATDGridType = 3, 
        PSATDTimeStraggeredGridType = 4
    };

    /* Yee grid for FDTD method. Implementation is based on Computational
    Electrodynamics: The Finite-Difference Time-Domain Method by Allen Taflove.
    The class provides index-based access and trilinear interpolation from 8
    nearest grid values.
    Components of the electric and magnetic fields correspond to different points.
    Index (i, j, k) refers to cell (i, j, k) being the cell between nodes (i, j, k)
    and (i + 1, j + 1, k + 1). It corresponds to Ex(i, j+1/2, k+1/2),
    Ey(i+1/2, j, k+1/2), Ez(i+1/2, j+1/2, k), same for the components of J,
    Bx(i+1/2, j, k), By(i, j+1/2, k), Bz(i, j, k+1/2).
    There is also EB_timeShift = dt/2 time difference between E(J) and B as needed
    for FDTD.
    For field exchanges and interpolation we need to store some grid values
    that correspond to points outside of the main grid area, size of internal and
    external area are given as constant members. */

    // auxiliary structs to describe some properties of a grid type

    template <GridTypes gridTypes>
    struct LabelFieldsSpatialStraggered {
        static const bool ifFieldsSpatialStraggered = gridTypes == GridTypes::YeeGridType;
    };

    template <GridTypes gridTypes>
    struct LabelFieldsTimeStraggered {
        static const bool ifFieldsTimeStraggered = gridTypes == GridTypes::YeeGridType ||
            gridTypes == GridTypes::PSTDGridType ||
            gridTypes == GridTypes::PSATDTimeStraggeredGridType;
    };

    template <GridTypes gridTypes>
    struct LabelMethodRequiredNumberOfExternalCells {
        static const int numExternalCells = (gridTypes == GridTypes::YeeGridType ||
            gridTypes == GridTypes::StraightGridType) ? 2 : 0;
    };

} // namespace pfc