#pragma once

#include "macros.h"

#include "GridTypes.h"
#include "ScalarField.h"
#include "Vectors.h"
#include "Constants.h"

namespace pfc {

    enum InterpolationType {
        Interpolation_CIC, Interpolation_TSC,
        Interpolation_SecondOrder, Interpolation_FourthOrder, Interpolation_PCS
    };

    template<typename Data, GridTypes gridType_>
    class Grid :
        // next labels define some properties of grid
        public LabelFieldsSpatialStraggered<gridType_>,
        public LabelFieldsTimeStraggered<gridType_>
    {

    public:

        static const GridTypes gridType = gridType_;

        Grid(const Int3 & _numInternalCells,
            const FP3 & minCoords, const FP3 & _steps,
            const Int3 & globalGridDims);
        Grid(const Int3 & _numAllCells,
            const Int3 & globalGridDims);  // for complex grids only
        Grid(const Int3 & _numAllCells, const Int3 & globalGridDims, 
            Grid<FP, gridType_>* grid);  // 'grid' and 'this' will have common memory

        // copy constructor, can make shallow copies
        Grid(const Grid& grid, bool ifShallowCopy = false);

        forceinline const FP3 BxPosition(int x, int y, int z) const
        {
            return baseCoords(x, y, z) + shiftBx;
        }
        forceinline const FP3 ByPosition(int x, int y, int z) const
        {
            return baseCoords(x, y, z) + shiftBy;
        }
        forceinline const FP3 BzPosition(int x, int y, int z) const
        {
            return baseCoords(x, y, z) + shiftBz;
        }
        forceinline const FP3 ExPosition(int x, int y, int z) const
        {
            return baseCoords(x, y, z) + shiftEJx;
        }
        forceinline const FP3 EyPosition(int x, int y, int z) const
        {
            return baseCoords(x, y, z) + shiftEJy;
        }
        forceinline const FP3 EzPosition(int x, int y, int z) const
        {
            return baseCoords(x, y, z) + shiftEJz;
        }
        forceinline const FP3 JxPosition(int x, int y, int z) const
        {
            return baseCoords(x, y, z) + shiftEJx;
        }
        forceinline const FP3 JyPosition(int x, int y, int z) const
        {
            return baseCoords(x, y, z) + shiftEJy;
        }
        forceinline const FP3 JzPosition(int x, int y, int z) const
        {
            return baseCoords(x, y, z) + shiftEJz;
        }

        void getFieldsXYZ(FP x, FP y, FP z, FP3 & e, FP3 & b) const
        {
            FP3 coords(x, y, z);
            getFields(coords, e, b);
        }
        void getFields(const FP3& coords, FP3 & e, FP3 & b) const
        {
            (this->*interpolationFields)(coords, e, b);
        }
        
        virtual FP3 getJ(const FP3& coords) const;
        virtual FP3 getE(const FP3& coords) const;
        virtual FP3 getB(const FP3& coords) const;

        void getFieldsCIC(const FP3& coords, FP3 & e, FP3 & b) const;
        void getFieldsTSC(const FP3& coords, FP3 & e, FP3 & b) const;
        void getFieldsSecondOrder(const FP3& coords, FP3 & e, FP3 & b) const;
        void getFieldsFourthOrder(const FP3& coords, FP3 & e, FP3 & b) const;
        void getFieldsPCS(const FP3& coords, FP3 & e, FP3 & b) const;

        FP getEx(const FP3& coords) const
        {
            return (this->*interpolationEx)(coords);
        }
        FP getEy(const FP3& coords) const
        {
            return (this->*interpolationEy)(coords);
        }
        FP getEz(const FP3& coords) const
        {
            return (this->*interpolationEz)(coords);
        }
        FP getBx(const FP3& coords) const
        {
            return (this->*interpolationBx)(coords);
        }
        FP getBy(const FP3& coords) const
        {
            return (this->*interpolationBy)(coords);
        }
        FP getBz(const FP3& coords) const
        {
            return (this->*interpolationBz)(coords);
        }
        FP getJx(const FP3& coords) const
        {
            return (this->*interpolationJx)(coords);
        }
        FP getJy(const FP3& coords) const
        {
            return (this->*interpolationJy)(coords);
        }
        FP getJz(const FP3& coords) const
        {
            return (this->*interpolationJz)(coords);
        }

        FP getExCIC(const FP3& coords) const {
            return getFieldCIC(coords, Ex, shiftEJx);
        }
        FP getEyCIC(const FP3& coords) const {
            return getFieldCIC(coords, Ey, shiftEJy);
        }
        FP getEzCIC(const FP3& coords) const {
            return getFieldCIC(coords, Ez, shiftEJz);
        }
        FP getBxCIC(const FP3& coords) const {
            return getFieldCIC(coords, Bx, shiftBx);
        }
        FP getByCIC(const FP3& coords) const {
            return getFieldCIC(coords, By, shiftBy);
        }
        FP getBzCIC(const FP3& coords) const {
            return getFieldCIC(coords, Bz, shiftBz);
        }
        FP getJxCIC(const FP3& coords) const {
            return getFieldCIC(coords, Jx, shiftEJx);
        }
        FP getJyCIC(const FP3& coords) const {
            return getFieldCIC(coords, Jy, shiftEJy);
        }
        FP getJzCIC(const FP3& coords) const {
            return getFieldCIC(coords, Jz, shiftEJz);
        }
        FP getExTSC(const FP3& coords) const {
            return getFieldTSC(coords, Ex, shiftEJx);
        }
        FP getEyTSC(const FP3& coords) const {
            return getFieldTSC(coords, Ey, shiftEJy);
        }
        FP getEzTSC(const FP3& coords) const {
            return getFieldTSC(coords, Ez, shiftEJz);
        }
        FP getBxTSC(const FP3& coords) const {
            return getFieldTSC(coords, Bx, shiftBx);
        }
        FP getByTSC(const FP3& coords) const {
            return getFieldTSC(coords, By, shiftBy);
        }
        FP getBzTSC(const FP3& coords) const {
            return getFieldTSC(coords, Bz, shiftBz);
        }
        FP getJxTSC(const FP3& coords) const {
            return getFieldTSC(coords, Jx, shiftEJx);
        }
        FP getJyTSC(const FP3& coords) const {
            return getFieldTSC(coords, Jy, shiftEJy);
        }
        FP getJzTSC(const FP3& coords) const {
            return getFieldTSC(coords, Jz, shiftEJz);
        }
        FP getExSecondOrder(const FP3& coords) const {
            return getFieldSecondOrder(coords, Ex, shiftEJx);
        }
        FP getEySecondOrder(const FP3& coords) const {
            return getFieldSecondOrder(coords, Ey, shiftEJy);
        }
        FP getEzSecondOrder(const FP3& coords) const {
            return getFieldSecondOrder(coords, Ez, shiftEJz);
        }
        FP getBxSecondOrder(const FP3& coords) const {
            return getFieldSecondOrder(coords, Bx, shiftBx);
        }
        FP getBySecondOrder(const FP3& coords) const {
            return getFieldSecondOrder(coords, By, shiftBy);
        }
        FP getBzSecondOrder(const FP3& coords) const {
            return getFieldSecondOrder(coords, Bz, shiftBz);
        }
        FP getJxSecondOrder(const FP3& coords) const {
            return getFieldSecondOrder(coords, Jx, shiftEJx);
        }
        FP getJySecondOrder(const FP3& coords) const {
            return getFieldSecondOrder(coords, Jy, shiftEJy);
        }
        FP getJzSecondOrder(const FP3& coords) const {
            return getFieldSecondOrder(coords, Jz, shiftEJz);
        }
        FP getExFourthOrder(const FP3& coords) const {
            return getFieldFourthOrder(coords, Ex, shiftEJx);
        }
        FP getEyFourthOrder(const FP3& coords) const {
            return getFieldFourthOrder(coords, Ey, shiftEJy);
        }
        FP getEzFourthOrder(const FP3& coords) const {
            return getFieldFourthOrder(coords, Ez, shiftEJz);
        }
        FP getBxFourthOrder(const FP3& coords) const {
            return getFieldFourthOrder(coords, Bx, shiftBx);
        }
        FP getByFourthOrder(const FP3& coords) const {
            return getFieldFourthOrder(coords, By, shiftBy);
        }
        FP getBzFourthOrder(const FP3& coords) const {
            return getFieldFourthOrder(coords, Bz, shiftBz);
        }
        FP getJxFourthOrder(const FP3& coords) const {
            return getFieldFourthOrder(coords, Jx, shiftEJx);
        }
        FP getJyFourthOrder(const FP3& coords) const {
            return getFieldFourthOrder(coords, Jy, shiftEJy);
        }
        FP getJzFourthOrder(const FP3& coords) const {
            return getFieldFourthOrder(coords, Jz, shiftEJz);
        }
        FP getExPCS(const FP3& coords) const {
            return getFieldPCS(coords, Ex, shiftEJx);
        }
        FP getEyPCS(const FP3& coords) const {
            return getFieldPCS(coords, Ey, shiftEJy);
        }
        FP getEzPCS(const FP3& coords) const {
            return getFieldPCS(coords, Ez, shiftEJz);
        }
        FP getBxPCS(const FP3& coords) const {
            return getFieldPCS(coords, Bx, shiftBx);
        }
        FP getByPCS(const FP3& coords) const {
            return getFieldPCS(coords, By, shiftBy);
        }
        FP getBzPCS(const FP3& coords) const {
            return getFieldPCS(coords, Bz, shiftBz);
        }
        FP getJxPCS(const FP3& coords) const {
            return getFieldPCS(coords, Jx, shiftEJx);
        }
        FP getJyPCS(const FP3& coords) const {
            return getFieldPCS(coords, Jy, shiftEJy);
        }
        FP getJzPCS(const FP3& coords) const {
            return getFieldPCS(coords, Jz, shiftEJz);
        }

        /*void dumpE(FP3 * e, const Int3 * minCellIdx, const Int3 * maxCellIdx);
        void dumpB(FP3 * b, const Int3 * minCellIdx, const Int3 * maxCellIdx);
        void dumpCurrents(FP3 * currents, const Int3 * minCellIdx, const Int3 * maxCellIdx);

        void loadE(const FP3 * e, const Int3 * minCellIdx, const Int3 * maxCellIdx);
        void loadB(const FP3 * b, const Int3 * minCellIdx, const Int3 * maxCellIdx);
        void loadCurrents(const FP3 * currents, const Int3 * minCellIdx, const Int3 * maxCellIdx);*/

        /* Make all current density values zero. */
        void zeroizeJ();

        const Int3 getNumExternalLeftCells() const
        {
            Int3 result(2, 2, 2);
            for (int d = 0; d < 3; d++)
                if (globalGridDims[d] == 1)
                    result[d] = 0;
            return result;
        }

        const Int3 getNumExternalRightCells() const
        {
            return getNumExternalLeftCells();
        }

        void setInterpolationType(InterpolationType type);
        InterpolationType getInterpolationType() const;

        const Int3 globalGridDims;  // important to initialize it first
        const FP3 steps;
        const Int3 numInternalCells;
        const Int3 numCells;
        const Int3 sizeStorage;  // sometimes can be larger than numCells
        const FP3 origin;
        const int dimensionality;

        ScalarField<Data> Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz;
        
    private:

        // 3d shifts of the field in the cell
        const FP3 shiftEJx, shiftEJy, shiftEJz,
            shiftBx, shiftBy, shiftBz;

        /* Get grid index and normalized internal coords in [0, 0, 0]..(1, 1, 1) for
        given physical coords and shift. */
        void getGridCoords(const FP3 & coords, const FP3 & shift, Int3 & idx,
            FP3 & internalCoords) const
        {
            idx.x = (int)((coords.x - origin.x - shift.x) / steps.x);
            idx.y = (int)((coords.y - origin.y - shift.y) / steps.y);
            idx.z = (int)((coords.z - origin.z - shift.z) / steps.z);
            internalCoords = (coords - baseCoords(idx.x, idx.y, idx.z) - shift) / steps;
        }

        void getClosestGridCoords(const FP3 & coords, const FP3 & shift, Int3 & idx,
            FP3 & internalCoords) const
        {
            idx.x = (int)((coords.x - origin.x - shift.x) / steps.x + 0.5);
            idx.y = (int)((coords.y - origin.y - shift.y) / steps.y + 0.5);
            idx.z = (int)((coords.z - origin.z - shift.z) / steps.z + 0.5);
            internalCoords = (coords - baseCoords(idx.x, idx.y, idx.z) - shift) / steps;
        }

        /* Get base coords of element (i, j, k) so that its real coords are
        base coords + corresponding shift. */
        forceinline const FP3 baseCoords(int i, int j, int k) const
        {
            return origin + FP3(i, j, k) * steps;
        }

        // if coords is inside of the area that grid defines
        forceinline bool isInside(const FP3 & coords, const FP3 & shift) const
        {
            FP3 minCoords = origin + shift * steps;
            FP3 maxCoords = minCoords + (numCells - Int3(1, 1, 1)) * steps;
            return coords >= minCoords && coords <= maxCoords;
        }

        FP getFieldCIC(const FP3& coords, const ScalarField<Data>& field, const FP3 & shift) const;
        FP getFieldTSC(const FP3& coords, const ScalarField<Data>& field, const FP3 & shift) const;
        FP getFieldSecondOrder(const FP3& coords, const ScalarField<Data>& field, const FP3 & shift) const;
        FP getFieldFourthOrder(const FP3& coords, const ScalarField<Data>& field, const FP3 & shift) const;
        FP getFieldPCS(const FP3& coords, const ScalarField<Data>& field, const FP3 & shift) const;

        InterpolationType interpolationType;
        void (Grid::*interpolationFields)(const FP3&, FP3&, FP3&) const;
        FP (Grid::*interpolationEx)(const FP3&) const;
        FP(Grid::*interpolationEy)(const FP3&) const;
        FP(Grid::*interpolationEz)(const FP3&) const;
        FP(Grid::*interpolationBx)(const FP3&) const;
        FP(Grid::*interpolationBy)(const FP3&) const;
        FP(Grid::*interpolationBz)(const FP3&) const;
        FP(Grid::*interpolationJx)(const FP3&) const;
        FP(Grid::*interpolationJy)(const FP3&) const;
        FP(Grid::*interpolationJz)(const FP3&) const;
    };

    typedef Grid<FP, GridTypes::YeeGridType> YeeGrid;
    typedef Grid<FP, GridTypes::StraightGridType> SimpleGrid;
    typedef Grid<FP, GridTypes::PSTDGridType> PSTDGrid;
    typedef Grid<FP, GridTypes::PSATDGridType> PSATDGrid;
    typedef Grid<FP, GridTypes::PSATDTimeStraggeredGridType> PSATDTimeStraggeredGrid;

    // create deep or shallow copy
    template<typename Data, GridTypes gridType_>
    inline Grid<Data, gridType_>::Grid(const Grid<Data, gridType_>& grid, bool ifShallowCopy) :
        globalGridDims(grid.globalGridDims),
        steps(grid.steps),
        numInternalCells(grid.numInternalCells),
        numCells(grid.numCells),
        sizeStorage(grid.sizeStorage),
        shiftEJx(grid.shiftEJx), shiftEJy(grid.shiftEJy), shiftEJz(grid.shiftEJz),
        shiftBx(grid.shiftBx), shiftBy(grid.shiftBy), shiftBz(grid.shiftBz),
        origin(grid.origin),
        dimensionality(grid.dimensionality),
        Ex(grid.Ex, ifShallowCopy),
        Ey(grid.Ey, ifShallowCopy),
        Ez(grid.Ez, ifShallowCopy),
        Bx(grid.Bx, ifShallowCopy),
        By(grid.By, ifShallowCopy),
        Bz(grid.Bz, ifShallowCopy),
        Jx(grid.Jx, ifShallowCopy),
        Jy(grid.Jy, ifShallowCopy),
        Jz(grid.Jz, ifShallowCopy)
    {
        setInterpolationType(grid.interpolationType);
    }

    template <>
    inline Grid<FP, GridTypes::YeeGridType>::Grid(const Int3 & _numCells, const FP3 & minCoords,
        const FP3 & _steps, const Int3 & _globalGridDims) :
        globalGridDims(_globalGridDims),
        steps(_steps),
        numInternalCells(_numCells),
        numCells(numInternalCells + getNumExternalLeftCells() + getNumExternalRightCells()),
        sizeStorage(numCells),
        Ex(sizeStorage), Ey(sizeStorage), Ez(sizeStorage),
        Bx(sizeStorage), By(sizeStorage), Bz(sizeStorage),
        Jx(sizeStorage), Jy(sizeStorage), Jz(sizeStorage),
        shiftEJx(FP3(0, 0.5, 0.5) * steps),
        shiftEJy(FP3(0.5, 0, 0.5) * steps),
        shiftEJz(FP3(0.5, 0.5, 0) * steps),
        shiftBx(FP3(0.5, 0, 0) * steps),
        shiftBy(FP3(0, 0.5, 0) * steps),
        shiftBz(FP3(0, 0, 0.5) * steps),
        origin(minCoords.x - steps.x * getNumExternalLeftCells().x,
            minCoords.y - steps.y * getNumExternalLeftCells().y,
            minCoords.z - steps.z * getNumExternalLeftCells().z),
        dimensionality((_globalGridDims.x != 1) + (_globalGridDims.y != 1) + (_globalGridDims.z != 1))
    {
        setInterpolationType(Interpolation_CIC);
        setInterpolationType(Interpolation_CIC);
    }

    template<>
    inline Grid<FP, GridTypes::StraightGridType>::Grid(const Int3 & _numInternalCells, 
        const FP3 & minCoords, const FP3 & _steps, const Int3 & _globalGridDims) :
        globalGridDims(_globalGridDims),
        steps(_steps),
        numInternalCells(_numInternalCells),
        numCells(numInternalCells + getNumExternalLeftCells() + getNumExternalRightCells()),
        sizeStorage(numCells),
        Ex(sizeStorage), Ey(sizeStorage), Ez(sizeStorage),
        Bx(sizeStorage), By(sizeStorage), Bz(sizeStorage),
        Jx(sizeStorage), Jy(sizeStorage), Jz(sizeStorage),
        shiftEJx(FP3(0, 0, 0) * steps),
        shiftEJy(FP3(0, 0, 0) * steps),
        shiftEJz(FP3(0, 0, 0) * steps),
        shiftBx(FP3(0, 0, 0) * steps),
        shiftBy(FP3(0, 0, 0) * steps),
        shiftBz(FP3(0, 0, 0) * steps),
        origin(minCoords.x - steps.x * getNumExternalLeftCells().x,
            minCoords.y - steps.y * getNumExternalLeftCells().y,
            minCoords.z - steps.z * getNumExternalLeftCells().z),
        dimensionality((_globalGridDims.x != 1) + (_globalGridDims.y != 1) + (_globalGridDims.z != 1))
    {
        setInterpolationType(Interpolation_CIC);
    }


    // SPECTRAL GRIDS

    // PSTD

    template<>
    inline Grid<complexFP, GridTypes::PSTDGridType>::Grid(const Int3 & _numInternalCells,
        const Int3 & _globalGridDims) :
        globalGridDims(_globalGridDims),
        numInternalCells(_numInternalCells),
        numCells(numInternalCells),
        sizeStorage(numCells),
        Ex(sizeStorage), Ey(sizeStorage), Ez(sizeStorage),
        Bx(sizeStorage), By(sizeStorage), Bz(sizeStorage),
        Jx(sizeStorage), Jy(sizeStorage), Jz(sizeStorage),
        shiftEJx(FP3(0, 0, 0) * steps),
        shiftEJy(FP3(0, 0, 0) * steps),
        shiftEJz(FP3(0, 0, 0) * steps),
        shiftBx(FP3(0, 0, 0) * steps),
        shiftBy(FP3(0, 0, 0) * steps),
        shiftBz(FP3(0, 0, 0) * steps),
        dimensionality((_globalGridDims.x != 1) + (_globalGridDims.y != 1) + (_globalGridDims.z != 1))
    {
        setInterpolationType(Interpolation_CIC);
    }

    template<>
    inline Grid<complexFP, GridTypes::PSTDGridType>::Grid(const Int3 & _numInternalCells,
        const Int3 & _globalGridDims, Grid<FP, GridTypes::PSTDGridType>* grid) :
        globalGridDims(_globalGridDims),
        numInternalCells(_numInternalCells),
        numCells(numInternalCells),
        sizeStorage(numCells),
        Ex(reinterpret_cast<complexFP*>(grid->Ex.getData()), sizeStorage),
        Ey(reinterpret_cast<complexFP*>(grid->Ey.getData()), sizeStorage),
        Ez(reinterpret_cast<complexFP*>(grid->Ez.getData()), sizeStorage),
        Bx(reinterpret_cast<complexFP*>(grid->Bx.getData()), sizeStorage),
        By(reinterpret_cast<complexFP*>(grid->By.getData()), sizeStorage),
        Bz(reinterpret_cast<complexFP*>(grid->Bz.getData()), sizeStorage),
        Jx(reinterpret_cast<complexFP*>(grid->Jx.getData()), sizeStorage),
        Jy(reinterpret_cast<complexFP*>(grid->Jy.getData()), sizeStorage),
        Jz(reinterpret_cast<complexFP*>(grid->Jz.getData()), sizeStorage),
        shiftEJx(FP3(0, 0, 0) * steps),
        shiftEJy(FP3(0, 0, 0) * steps),
        shiftEJz(FP3(0, 0, 0) * steps),
        shiftBx(FP3(0, 0, 0) * steps),
        shiftBy(FP3(0, 0, 0) * steps),
        shiftBz(FP3(0, 0, 0) * steps),
        dimensionality((_globalGridDims.x != 1) + (_globalGridDims.y != 1) + (_globalGridDims.z != 1))
    {
        setInterpolationType(Interpolation_CIC);
    }

    template<>
    inline Grid<FP, GridTypes::PSTDGridType>::Grid(const Int3 & _numInternalCells,
        const FP3 & minCoords, const FP3 & _steps, const Int3 & _globalGridDims) :
        globalGridDims(_globalGridDims),
        steps(_steps),
        numInternalCells(_numInternalCells),
        numCells(numInternalCells),
        sizeStorage(Int3(numCells.x, numCells.y, 2 * (numCells.z / 2 + 1))),
        Ex(sizeStorage), Ey(sizeStorage), Ez(sizeStorage),
        Bx(sizeStorage), By(sizeStorage), Bz(sizeStorage),
        Jx(sizeStorage), Jy(sizeStorage), Jz(sizeStorage),
        shiftEJx(FP3(0, 0, 0) * steps),
        shiftEJy(FP3(0, 0, 0) * steps),
        shiftEJz(FP3(0, 0, 0) * steps),
        shiftBx(FP3(0, 0, 0) * steps),
        shiftBy(FP3(0, 0, 0) * steps),
        shiftBz(FP3(0, 0, 0) * steps),
        origin(minCoords),
        dimensionality((_globalGridDims.x != 1) + (_globalGridDims.y != 1) + (_globalGridDims.z != 1))
    {
        setInterpolationType(Interpolation_CIC);
    }

    // PSATD

    template<>
    inline Grid<complexFP, GridTypes::PSATDGridType>::Grid(const Int3 & _numInternalCells,
        const Int3 & _globalGridDims) :
        globalGridDims(_globalGridDims),
        numInternalCells(_numInternalCells),
        numCells(numInternalCells),
        sizeStorage(numCells),
        Ex(sizeStorage), Ey(sizeStorage), Ez(sizeStorage),
        Bx(sizeStorage), By(sizeStorage), Bz(sizeStorage),
        Jx(sizeStorage), Jy(sizeStorage), Jz(sizeStorage),
        shiftEJx(FP3(0, 0, 0) * steps),
        shiftEJy(FP3(0, 0, 0) * steps),
        shiftEJz(FP3(0, 0, 0) * steps),
        shiftBx(FP3(0, 0, 0) * steps),
        shiftBy(FP3(0, 0, 0) * steps),
        shiftBz(FP3(0, 0, 0) * steps),
        dimensionality((_globalGridDims.x != 1) + (_globalGridDims.y != 1) + (_globalGridDims.z != 1))
    {
        setInterpolationType(Interpolation_CIC);
    }

    template<>
    inline Grid<complexFP, GridTypes::PSATDGridType>::Grid(const Int3 & _numInternalCells,
        const Int3 & _globalGridDims, Grid<FP, GridTypes::PSATDGridType>* grid) :
        globalGridDims(_globalGridDims),
        numInternalCells(_numInternalCells),
        numCells(numInternalCells),
        sizeStorage(numCells),
        Ex(reinterpret_cast<complexFP*>(grid->Ex.getData()), sizeStorage),
        Ey(reinterpret_cast<complexFP*>(grid->Ey.getData()), sizeStorage),
        Ez(reinterpret_cast<complexFP*>(grid->Ez.getData()), sizeStorage),
        Bx(reinterpret_cast<complexFP*>(grid->Bx.getData()), sizeStorage),
        By(reinterpret_cast<complexFP*>(grid->By.getData()), sizeStorage),
        Bz(reinterpret_cast<complexFP*>(grid->Bz.getData()), sizeStorage),
        Jx(reinterpret_cast<complexFP*>(grid->Jx.getData()), sizeStorage),
        Jy(reinterpret_cast<complexFP*>(grid->Jy.getData()), sizeStorage),
        Jz(reinterpret_cast<complexFP*>(grid->Jz.getData()), sizeStorage),
        shiftEJx(FP3(0, 0, 0) * steps),
        shiftEJy(FP3(0, 0, 0) * steps),
        shiftEJz(FP3(0, 0, 0) * steps),
        shiftBx(FP3(0, 0, 0) * steps),
        shiftBy(FP3(0, 0, 0) * steps),
        shiftBz(FP3(0, 0, 0) * steps),
        dimensionality((_globalGridDims.x != 1) + (_globalGridDims.y != 1) + (_globalGridDims.z != 1))
    {
        setInterpolationType(Interpolation_CIC);
    }

    template<>
    inline Grid<FP, GridTypes::PSATDGridType>::Grid(const Int3 & _numInternalCells,
        const FP3 & minCoords, const FP3 & _steps, const Int3 & _globalGridDims) :
        globalGridDims(_globalGridDims),
        steps(_steps),
        numInternalCells(_numInternalCells),
        numCells(numInternalCells),
        sizeStorage(Int3(numCells.x, numCells.y, 2 * (numCells.z / 2 + 1))),
        Ex(sizeStorage), Ey(sizeStorage), Ez(sizeStorage),
        Bx(sizeStorage), By(sizeStorage), Bz(sizeStorage),
        Jx(sizeStorage), Jy(sizeStorage), Jz(sizeStorage),
        shiftEJx(FP3(0, 0, 0) * steps),
        shiftEJy(FP3(0, 0, 0) * steps),
        shiftEJz(FP3(0, 0, 0) * steps),
        shiftBx(FP3(0, 0, 0) * steps),
        shiftBy(FP3(0, 0, 0) * steps),
        shiftBz(FP3(0, 0, 0) * steps),
        origin(minCoords),
        dimensionality((_globalGridDims.x != 1) + (_globalGridDims.y != 1) + (_globalGridDims.z != 1))
    {
        setInterpolationType(Interpolation_CIC);
    }

    // PSATDTimeStraggered

    template<>
    inline Grid<complexFP, GridTypes::PSATDTimeStraggeredGridType>::Grid(const Int3 & _numInternalCells,
        const Int3 & _globalGridDims) :
        globalGridDims(_globalGridDims),
        numInternalCells(_numInternalCells),
        numCells(numInternalCells),
        sizeStorage(numCells),
        Ex(sizeStorage), Ey(sizeStorage), Ez(sizeStorage),
        Bx(sizeStorage), By(sizeStorage), Bz(sizeStorage),
        Jx(sizeStorage), Jy(sizeStorage), Jz(sizeStorage),
        shiftEJx(FP3(0, 0, 0) * steps),
        shiftEJy(FP3(0, 0, 0) * steps),
        shiftEJz(FP3(0, 0, 0) * steps),
        shiftBx(FP3(0, 0, 0) * steps),
        shiftBy(FP3(0, 0, 0) * steps),
        shiftBz(FP3(0, 0, 0) * steps),
        dimensionality((_globalGridDims.x != 1) + (_globalGridDims.y != 1) + (_globalGridDims.z != 1))
    {
        setInterpolationType(Interpolation_CIC);       
    }

    // PSATDTimeStraggered

    template<>
    inline Grid<complexFP, GridTypes::PSATDTimeStraggeredGridType>::Grid(const Int3 & _numInternalCells,
        const Int3 & _globalGridDims, Grid<FP, GridTypes::PSATDTimeStraggeredGridType>* grid) :
        globalGridDims(_globalGridDims),
        numInternalCells(_numInternalCells),
        numCells(numInternalCells),
        sizeStorage(numCells),
        Ex(reinterpret_cast<complexFP*>(grid->Ex.getData()), sizeStorage),
        Ey(reinterpret_cast<complexFP*>(grid->Ey.getData()), sizeStorage),
        Ez(reinterpret_cast<complexFP*>(grid->Ez.getData()), sizeStorage),
        Bx(reinterpret_cast<complexFP*>(grid->Bx.getData()), sizeStorage),
        By(reinterpret_cast<complexFP*>(grid->By.getData()), sizeStorage),
        Bz(reinterpret_cast<complexFP*>(grid->Bz.getData()), sizeStorage),
        Jx(reinterpret_cast<complexFP*>(grid->Jx.getData()), sizeStorage),
        Jy(reinterpret_cast<complexFP*>(grid->Jy.getData()), sizeStorage),
        Jz(reinterpret_cast<complexFP*>(grid->Jz.getData()), sizeStorage),
        shiftEJx(FP3(0, 0, 0) * steps),
        shiftEJy(FP3(0, 0, 0) * steps),
        shiftEJz(FP3(0, 0, 0) * steps),
        shiftBx(FP3(0, 0, 0) * steps),
        shiftBy(FP3(0, 0, 0) * steps),
        shiftBz(FP3(0, 0, 0) * steps),
        dimensionality((_globalGridDims.x != 1) + (_globalGridDims.y != 1) + (_globalGridDims.z != 1))
    {
        setInterpolationType(Interpolation_CIC);
    }

    template<>
    inline Grid<FP, GridTypes::PSATDTimeStraggeredGridType>::Grid(const Int3 & _numInternalCells,
        const FP3 & minCoords, const FP3 & _steps, const Int3 & _globalGridDims) :
        globalGridDims(_globalGridDims),
        steps(_steps),
        numInternalCells(_numInternalCells),
        numCells(numInternalCells),
        sizeStorage(Int3(numCells.x, numCells.y, 2 * (numCells.z / 2 + 1))),
        Ex(sizeStorage), Ey(sizeStorage), Ez(sizeStorage),
        Bx(sizeStorage), By(sizeStorage), Bz(sizeStorage),
        Jx(sizeStorage), Jy(sizeStorage), Jz(sizeStorage),
        shiftEJx(FP3(0, 0, 0) * steps),
        shiftEJy(FP3(0, 0, 0) * steps),
        shiftEJz(FP3(0, 0, 0) * steps),
        shiftBx(FP3(0, 0, 0) * steps),
        shiftBy(FP3(0, 0, 0) * steps),
        shiftBz(FP3(0, 0, 0) * steps),
        origin(minCoords),
        dimensionality((_globalGridDims.x != 1) + (_globalGridDims.y != 1) + (_globalGridDims.z != 1))
    {
        setInterpolationType(Interpolation_CIC);
    }


    // end SPECTRAL GRIDS
    
    template< typename Data, GridTypes gT>
    inline FP Grid<Data, gT>::getFieldCIC(const FP3& coords, const ScalarField<Data>& field, const FP3 & shift) const
    {
        Int3 idx;
        FP3 internalCoords;
        getGridCoords(coords, shift, idx, internalCoords);
        return field.interpolateCIC(idx, internalCoords);
    }

    template< typename Data, GridTypes gT>
    inline FP Grid<Data, gT>::getFieldTSC(const FP3& coords, const ScalarField<Data>& field, const FP3 & shift) const
    {
        Int3 idx;
        FP3 internalCoords;
        getClosestGridCoords(coords, shift, idx, internalCoords);
        return field.interpolateTSC(idx, internalCoords);
    }

    template< typename Data, GridTypes gT>
    inline FP Grid<Data, gT>::getFieldSecondOrder(const FP3& coords, const ScalarField<Data>& field, const FP3 & shift) const
    {
        Int3 idx;
        FP3 internalCoords;
        getClosestGridCoords(coords, shift, idx, internalCoords);
        return field.interpolateSecondOrder(idx, internalCoords);
    }

    template< typename Data, GridTypes gT>
    inline FP Grid<Data, gT>::getFieldFourthOrder(const FP3& coords, const ScalarField<Data>& field, const FP3 & shift) const
    {
        Int3 idx;
        FP3 internalCoords;
        getClosestGridCoords(coords, shift, idx, internalCoords);
        return field.interpolateFourthOrder(idx, internalCoords);
    }

    template< typename Data, GridTypes gT>
    inline FP Grid<Data, gT>::getFieldPCS(const FP3& coords, const ScalarField<Data>& field, const FP3 & shift) const
    {
        Int3 idx;
        FP3 internalCoords;
        getGridCoords(coords, shift, idx, internalCoords);
        return field.interpolatePCS(idx, internalCoords);
    }

    template< typename Data, GridTypes gT>
    inline void Grid<Data, gT>::getFieldsCIC(const FP3& coords, FP3 & e, FP3 & b) const
    {
        /* For each component of E and B get grid index and internal coords,
        use it as base index and coefficients of interpolation. */
        Int3 idx;
        FP3 internalCoords;

        getGridCoords(coords, shiftEJx, idx, internalCoords);
        e.x = Ex.interpolateCIC(idx, internalCoords);
        getGridCoords(coords, shiftEJy, idx, internalCoords);
        e.y = Ey.interpolateCIC(idx, internalCoords);
        getGridCoords(coords, shiftEJz, idx, internalCoords);
        e.z = Ez.interpolateCIC(idx, internalCoords);

        getGridCoords(coords, shiftBx, idx, internalCoords);
        b.x = Bx.interpolateCIC(idx, internalCoords);
        getGridCoords(coords, shiftBy, idx, internalCoords);
        b.y = By.interpolateCIC(idx, internalCoords);
        getGridCoords(coords, shiftBz, idx, internalCoords);
        b.z = Bz.interpolateCIC(idx, internalCoords);
    }

    template< typename Data, GridTypes gT>
    inline void Grid<Data, gT>::getFieldsTSC(const FP3& coords, FP3 & e, FP3 & b) const
    {
        Int3 idx;
        FP3 internalCoords;

        getClosestGridCoords(coords, shiftEJx, idx, internalCoords);
        e.x = Ex.interpolateTSC(idx, internalCoords);
        getClosestGridCoords(coords, shiftEJy, idx, internalCoords);
        e.y = Ey.interpolateTSC(idx, internalCoords);
        getClosestGridCoords(coords, shiftEJz, idx, internalCoords);
        e.z = Ez.interpolateTSC(idx, internalCoords);

        getClosestGridCoords(coords, shiftBx, idx, internalCoords);
        b.x = Bx.interpolateTSC(idx, internalCoords);
        getClosestGridCoords(coords, shiftBy, idx, internalCoords);
        b.y = By.interpolateTSC(idx, internalCoords);
        getClosestGridCoords(coords, shiftBz, idx, internalCoords);
        b.z = Bz.interpolateTSC(idx, internalCoords);
    }

    template< typename Data, GridTypes gT>
    inline void Grid<Data, gT>::getFieldsSecondOrder(const FP3& coords, FP3 & e, FP3 & b) const
    {
        Int3 idx;
        FP3 internalCoords;

        getClosestGridCoords(coords, shiftEJx, idx, internalCoords);
        e.x = Ex.interpolateSecondOrder(idx, internalCoords);
        getClosestGridCoords(coords, shiftEJy, idx, internalCoords);
        e.y = Ey.interpolateSecondOrder(idx, internalCoords);
        getClosestGridCoords(coords, shiftEJz, idx, internalCoords);
        e.z = Ez.interpolateSecondOrder(idx, internalCoords);

        getClosestGridCoords(coords, shiftBx, idx, internalCoords);
        b.x = Bx.interpolateSecondOrder(idx, internalCoords);
        getClosestGridCoords(coords, shiftBy, idx, internalCoords);
        b.y = By.interpolateSecondOrder(idx, internalCoords);
        getClosestGridCoords(coords, shiftBz, idx, internalCoords);
        b.z = Bz.interpolateSecondOrder(idx, internalCoords);
    }

    template< typename Data, GridTypes gT>
    inline void Grid<Data, gT>::getFieldsFourthOrder(const FP3& coords, FP3 & e, FP3 & b) const
    {
        Int3 idx;
        FP3 internalCoords;

        getClosestGridCoords(coords, shiftEJx, idx, internalCoords);
        e.x = Ex.interpolateFourthOrder(idx, internalCoords);
        getClosestGridCoords(coords, shiftEJy, idx, internalCoords);
        e.y = Ey.interpolateFourthOrder(idx, internalCoords);
        getClosestGridCoords(coords, shiftEJz, idx, internalCoords);
        e.z = Ez.interpolateFourthOrder(idx, internalCoords);

        getClosestGridCoords(coords, shiftBx, idx, internalCoords);
        b.x = Bx.interpolateFourthOrder(idx, internalCoords);
        getClosestGridCoords(coords, shiftBy, idx, internalCoords);
        b.y = By.interpolateFourthOrder(idx, internalCoords);
        getClosestGridCoords(coords, shiftBz, idx, internalCoords);
        b.z = Bz.interpolateFourthOrder(idx, internalCoords);
    }

    template< typename Data, GridTypes gT>
    inline void Grid<Data, gT>::getFieldsPCS(const FP3& coords, FP3 & e, FP3 & b) const
    {
        Int3 idx;
        FP3 internalCoords;

        getGridCoords(coords, shiftEJx, idx, internalCoords);
        e.x = Ex.interpolatePCS(idx, internalCoords);
        getGridCoords(coords, shiftEJy, idx, internalCoords);
        e.y = Ey.interpolatePCS(idx, internalCoords);
        getGridCoords(coords, shiftEJz, idx, internalCoords);
        e.z = Ez.interpolatePCS(idx, internalCoords);

        getGridCoords(coords, shiftBx, idx, internalCoords);
        b.x = Bx.interpolatePCS(idx, internalCoords);
        getGridCoords(coords, shiftBy, idx, internalCoords);
        b.y = By.interpolatePCS(idx, internalCoords);
        getGridCoords(coords, shiftBz, idx, internalCoords);
        b.z = Bz.interpolatePCS(idx, internalCoords);
    }

    template< typename Data, GridTypes gT>
    inline FP3 Grid<Data, gT>::getJ(const FP3& coords) const
    {
        // zero fields are outside of area that grid defines
        //if (!isInside(coords, shiftEJx) || !isInside(coords, shiftEJy) || !isInside(coords, shiftEJz))
        //    return FP3(0, 0, 0);

        /* For each component of J get grid index and internal coords,
        use it as base index and coefficients of interpolation. */
        Int3 idx;
        FP3 internalCoords;
        FP3 j;

        getGridCoords(coords, shiftEJx, idx, internalCoords);
        j.x = Jx.interpolateCIC(idx, internalCoords);
        getGridCoords(coords, shiftEJy, idx, internalCoords);
        j.y = Jy.interpolateCIC(idx, internalCoords);
        getGridCoords(coords, shiftEJz, idx, internalCoords);
        j.z = Jz.interpolateCIC(idx, internalCoords);

        return j;
    }

    template< typename Data, GridTypes gT>
    inline FP3 Grid<Data, gT>::getE(const FP3& coords) const
    {
        // zero fields are outside of area that grid defines
        //if (!isInside(coords, shiftEJx) || !isInside(coords, shiftEJy) || !isInside(coords, shiftEJz))
        //    return FP3(0, 0, 0);

        /* For each component of J get grid index and internal coords,
        use it as base index and coefficients of interpolation. */
        Int3 idx;
        FP3 internalCoords;
        FP3 e;

        getGridCoords(coords, shiftEJx, idx, internalCoords);
        e.x = Ex.interpolateCIC(idx, internalCoords);
        getGridCoords(coords, shiftEJy, idx, internalCoords);
        e.y = Ey.interpolateCIC(idx, internalCoords);
        getGridCoords(coords, shiftEJz, idx, internalCoords);
        e.z = Ez.interpolateCIC(idx, internalCoords);

        return e;
    }

    template< typename Data, GridTypes gT>
    inline FP3 Grid<Data, gT>::getB(const FP3& coords) const
    {
        // zero fields are outside of area that grid defines
        //if (!isInside(coords, shiftBx) || !isInside(coords, shiftBy) || !isInside(coords, shiftBz))
        //    return FP3(0, 0, 0);

        /* For each component of J get grid index and internal coords,
        use it as base index and coefficients of interpolation. */
        Int3 idx;
        FP3 internalCoords;
        FP3 b;

        getGridCoords(coords, shiftBx, idx, internalCoords);
        b.x = Bx.interpolateCIC(idx, internalCoords);
        getGridCoords(coords, shiftBy, idx, internalCoords);
        b.y = By.interpolateCIC(idx, internalCoords);
        getGridCoords(coords, shiftBz, idx, internalCoords);
        b.z = Bz.interpolateCIC(idx, internalCoords);

        return b;
    }
    
    template< typename Data, GridTypes gT>
    inline void Grid<Data, gT>::zeroizeJ()
    {
        Jx.zeroize();
        Jy.zeroize();
        Jz.zeroize();
    }

    template< typename Data, GridTypes gT>
    inline void Grid<Data, gT>::setInterpolationType(InterpolationType type)
    {
        interpolationType = type;
        switch (interpolationType)
        {
        case Interpolation_CIC:
            interpolationFields = &Grid<Data, gT>::getFieldsCIC;
            interpolationEx = &Grid<Data, gT>::getExCIC;
            interpolationEy = &Grid<Data, gT>::getEyCIC;
            interpolationEz = &Grid<Data, gT>::getEzCIC;
            interpolationBx = &Grid<Data, gT>::getBxCIC;
            interpolationBy = &Grid<Data, gT>::getByCIC;
            interpolationBz = &Grid<Data, gT>::getBzCIC;
            interpolationJx = &Grid<Data, gT>::getJxCIC;
            interpolationJy = &Grid<Data, gT>::getJyCIC;
            interpolationJz = &Grid<Data, gT>::getJzCIC; break;
        case Interpolation_TSC:
            interpolationFields = &Grid<Data, gT>::getFieldsTSC;
            interpolationEx = &Grid<Data, gT>::getExTSC;
            interpolationEy = &Grid<Data, gT>::getEyTSC;
            interpolationEz = &Grid<Data, gT>::getEzTSC;
            interpolationBx = &Grid<Data, gT>::getBxTSC;
            interpolationBy = &Grid<Data, gT>::getByTSC;
            interpolationBz = &Grid<Data, gT>::getBzTSC;
            interpolationJx = &Grid<Data, gT>::getJxTSC;
            interpolationJy = &Grid<Data, gT>::getJyTSC;
            interpolationJz = &Grid<Data, gT>::getJzTSC; break;
        case Interpolation_PCS:
            interpolationFields = &Grid<Data, gT>::getFieldsPCS;
            interpolationEx = &Grid<Data, gT>::getExPCS;
            interpolationEy = &Grid<Data, gT>::getEyPCS;
            interpolationEz = &Grid<Data, gT>::getEzPCS;
            interpolationBx = &Grid<Data, gT>::getBxPCS;
            interpolationBy = &Grid<Data, gT>::getByPCS;
            interpolationBz = &Grid<Data, gT>::getBzPCS;
            interpolationJx = &Grid<Data, gT>::getJxPCS;
            interpolationJy = &Grid<Data, gT>::getJyPCS;
            interpolationJz = &Grid<Data, gT>::getJzPCS; break;
        case Interpolation_SecondOrder:
            interpolationFields = &Grid<Data, gT>::getFieldsSecondOrder;
            interpolationEx = &Grid<Data, gT>::getExSecondOrder;
            interpolationEy = &Grid<Data, gT>::getEySecondOrder;
            interpolationEz = &Grid<Data, gT>::getEzSecondOrder;
            interpolationBx = &Grid<Data, gT>::getBxSecondOrder;
            interpolationBy = &Grid<Data, gT>::getBySecondOrder;
            interpolationBz = &Grid<Data, gT>::getBzSecondOrder;
            interpolationJx = &Grid<Data, gT>::getJxSecondOrder;
            interpolationJy = &Grid<Data, gT>::getJySecondOrder;
            interpolationJz = &Grid<Data, gT>::getJzSecondOrder; break;
        case Interpolation_FourthOrder:
            interpolationFields = &Grid<Data, gT>::getFieldsFourthOrder;
            interpolationEx = &Grid<Data, gT>::getExFourthOrder;
            interpolationEy = &Grid<Data, gT>::getEyFourthOrder;
            interpolationEz = &Grid<Data, gT>::getEzFourthOrder;
            interpolationBx = &Grid<Data, gT>::getBxFourthOrder;
            interpolationBy = &Grid<Data, gT>::getByFourthOrder;
            interpolationBz = &Grid<Data, gT>::getBzFourthOrder;
            interpolationJx = &Grid<Data, gT>::getJxFourthOrder;
            interpolationJy = &Grid<Data, gT>::getJyFourthOrder;
            interpolationJz = &Grid<Data, gT>::getJzFourthOrder; break;
        }
    }

    template<typename Data, GridTypes gT>
    inline InterpolationType Grid<Data, gT>::getInterpolationType() const
    {
        return interpolationType;
    }

    /*template<>
    inline void Grid<FP, YeeGridType>::dumpB(FP3 * b, const Int3 * minCellIdx, const Int3 * maxCellIdx)
    {
        Int3 numCells = *maxCellIdx - *minCellIdx;
#pragma omp parallel for collapse(3)
        for (int i = 0; i < numCells.x; ++i)
            for (int j = 0; j < numCells.y; ++j)
                for (int k = 0; k < numCells.z; ++k)
                {
                    int idx = numCells.y * numCells.z * i + numCells.z * j + k;
                    Int3 nodeIdx = *minCellIdx + Int3(i, j, k);
                    b[idx].x = Bx(nodeIdx);
                    b[idx].y = By(nodeIdx);
                    b[idx].z = Bz(nodeIdx);
                }
    }

    template<>
    inline void Grid<FP, YeeGridType>::dumpE(FP3 * e, const Int3 * minCellIdx, const Int3 * maxCellIdx)
    {
        Int3 numCells = *maxCellIdx - *minCellIdx;
#pragma omp parallel for collapse(3)
        for (int i = 0; i < numCells.x; ++i)
            for (int j = 0; j < numCells.y; ++j)
                for (int k = 0; k < numCells.z; ++k)
                {
                    int idx = numCells.y * numCells.z * i + numCells.z * j + k;
                    Int3 nodeIdx = *minCellIdx + Int3(i, j, k);
                    e[idx].x = Ex(nodeIdx);
                    e[idx].y = Ey(nodeIdx);
                    e[idx].z = Ez(nodeIdx);
                }
    }

    template<>
    inline void Grid<FP, YeeGridType>::dumpCurrents(FP3 * currents, const Int3 * minCellIdx,
        const Int3 * maxCellIdx)
    {
        Int3 numCells = *maxCellIdx - *minCellIdx;
#pragma omp parallel for collapse(3)
        for (int i = 0; i < numCells.x; ++i)
            for (int j = 0; j < numCells.y; ++j)
                for (int k = 0; k < numCells.z; ++k)
                {
                    int idx = numCells.y * numCells.z * i + numCells.z * j + k;
                    Int3 nodeIdx = *minCellIdx + Int3(i, j, k);
                    currents[idx].x = Jx(nodeIdx);
                    currents[idx].y = Jy(nodeIdx);
                    currents[idx].z = Jz(nodeIdx);
                    idx++;
                }
    }

    template<>
    inline void Grid<FP, YeeGridType>::loadE(const FP3 * e, const Int3 * minCellIdx, const Int3 * maxCellIdx)
    {
        Int3 numCells = *maxCellIdx - *minCellIdx;
#pragma omp parallel for collapse(3)
        for (int i = 0; i < numCells.x; i++)
            for (int j = 0; j < numCells.y; j++)
                for (int k = 0; k < numCells.z; k++)
                {
                    int idx = numCells.y * numCells.z * i + numCells.z * j + k;
                    Int3 nodeIdx = *minCellIdx + Int3(i, j, k);
                    Ex(nodeIdx) = e[idx].x;
                    Ey(nodeIdx) = e[idx].y;
                    Ez(nodeIdx) = e[idx].z;
                }
    }
    
    template<>
    inline void Grid<FP, YeeGridType>::loadB(const FP3 * b, const Int3 * minCellIdx, const Int3 * maxCellIdx)
    {
        Int3 numCells = *maxCellIdx - *minCellIdx;
#pragma omp parallel for collapse(3)
        for (int i = 0; i < numCells.x; ++i)
            for (int j = 0; j < numCells.y; ++j)
                for (int k = 0; k < numCells.z; ++k)
                {
                    int idx = numCells.y * numCells.z * i + numCells.z * j + k;
                    Int3 nodeIdx = *minCellIdx + Int3(i, j, k);
                    Bx(nodeIdx) = b[idx].x;
                    By(nodeIdx) = b[idx].y;
                    Bz(nodeIdx) = b[idx].z;
                }
    }

    template<>
    inline void Grid<FP, YeeGridType>::loadCurrents(const FP3 * currents, const Int3 * minCellIdx, const Int3 * maxCellIdx)
    {
        Int3 numCells = *maxCellIdx - *minCellIdx;
#pragma omp parallel for collapse(3)
        for (int i = 0; i < numCells.x; i++)
            for (int j = 0; j < numCells.y; j++)
                for (int k = 0; k < numCells.z; k++)
                {
                    int idx = numCells.y * numCells.z * i + numCells.z * j + k;
                    Int3 nodeIdx = *minCellIdx + Int3(i, j, k);
                    Jx(nodeIdx) = currents[idx].x;
                    Jy(nodeIdx) = currents[idx].y;
                    Jz(nodeIdx) = currents[idx].z;
                }
    }*/
}