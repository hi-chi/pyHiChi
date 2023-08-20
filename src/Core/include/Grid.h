#pragma once

#include "macros.h"

#include "GridTypes.h"
#include "ScalarField.h"
#include "Vectors.h"
#include "Constants.h"

#include "GridMacros.h"

#include <exception>


namespace pfc {

    enum class InterpolationType {
        Interpolation_CIC, Interpolation_TSC, Interpolation_PCS,
        Interpolation_SecondOrder, Interpolation_FourthOrder
    };

    template<typename Data, GridTypes gridType_>
    class Grid :
        // next labels define some properties of grid
        public LabelFieldsSpatialStaggered<gridType_>,
        /* defines 'numExternalCells' constant depending on method */
        public LabelMethodRequiredNumberOfExternalCells<gridType_>
    {

    public:

        static const GridTypes gridType = gridType_;
        static const bool isComplex = std::is_same<complexFP, Data>::value;

        Grid() {}
        Grid(const Int3& _numInternalCells,
            const FP3& minCoords, const FP3& _steps,
            const Int3& globalGridDims) {}
        Grid(const Grid<Data, gridType_>& grid);

        /* copy values from *this to *grid */
        template <class TGrid>
        void copyValues(TGrid* grid) const;

        // most of the following methods are declared with macros because has almost the same short implementation
        // as a rule, they call one private grid method with appreciate parameters
        // the first argument of each macro ('funcname') is a method name

        // the next methods map grid indices and coords and vice versa

        /* returns coords by index, uses 'baseCoords' method inside */
        /* signature: 'forceinline Int3 funcname(const FP3& coords) const' */
        GRID_GET_POSITION_IMPL(getBaseCoords, ZERO_SHIFT);
        GRID_GET_POSITION_IMPL(BxPosition, shiftBx);
        GRID_GET_POSITION_IMPL(ByPosition, shiftBy);
        GRID_GET_POSITION_IMPL(BzPosition, shiftBz);
        GRID_GET_POSITION_IMPL(ExPosition, shiftEJx);
        GRID_GET_POSITION_IMPL(EyPosition, shiftEJy);
        GRID_GET_POSITION_IMPL(EzPosition, shiftEJz);
        GRID_GET_POSITION_IMPL(JxPosition, shiftEJx);
        GRID_GET_POSITION_IMPL(JyPosition, shiftEJy);
        GRID_GET_POSITION_IMPL(JzPosition, shiftEJz);

        /* returns the closest left grid index for given coords, uses 'getGridCoords' method inside */
        /* signature: 'forceinline Int3 funcname(const FP3& coords) const' */
        GRID_GET_INDEX_IMPL(getBaseIndex, ZERO_SHIFT);
        GRID_GET_INDEX_IMPL(getIndexBx, shiftBx);
        GRID_GET_INDEX_IMPL(getIndexBy, shiftBy);
        GRID_GET_INDEX_IMPL(getIndexBz, shiftBz);
        GRID_GET_INDEX_IMPL(getIndexEx, shiftEJx);
        GRID_GET_INDEX_IMPL(getIndexEy, shiftEJy);
        GRID_GET_INDEX_IMPL(getIndexEz, shiftEJz);
        GRID_GET_INDEX_IMPL(getIndexJx, shiftEJx);
        GRID_GET_INDEX_IMPL(getIndexJy, shiftEJy);
        GRID_GET_INDEX_IMPL(getIndexJz, shiftEJz);     

        /* returns the closest (left or right) grid index for given coords, uses 'getClosestGridCoords' method */
        /* signature: 'forceinline Int3 funcname(const FP3& coords) const' */
        GRID_GET_INDEX_IMPL(getClosestBaseIndex, ZERO_SHIFT);
        GRID_GET_INDEX_IMPL(getClosestIndexBx, shiftBx);
        GRID_GET_INDEX_IMPL(getClosestIndexBy, shiftBy);
        GRID_GET_INDEX_IMPL(getClosestIndexBz, shiftBz);
        GRID_GET_INDEX_IMPL(getClosestIndexEx, shiftEJx);
        GRID_GET_INDEX_IMPL(getClosestIndexEy, shiftEJy);
        GRID_GET_INDEX_IMPL(getClosestIndexEz, shiftEJz);
        GRID_GET_INDEX_IMPL(getClosestIndexJx, shiftEJx);
        GRID_GET_INDEX_IMPL(getClosestIndexJy, shiftEJy);
        GRID_GET_INDEX_IMPL(getClosestIndexJz, shiftEJz);

        /* returns interpolated field value in arbitrary coords, uses 'isInside' method */
        /* signature: 'forceinline bool funcname(const FP3& coords) const' */
        GRID_IS_INSIDE_IMPL(isInsideBaseCoords, ZERO_SHIFT);
        GRID_IS_INSIDE_IMPL(isInsideCoordsBx, shiftBx);
        GRID_IS_INSIDE_IMPL(isInsideCoordsBy, shiftBy);
        GRID_IS_INSIDE_IMPL(isInsideCoordsBz, shiftBz);
        GRID_IS_INSIDE_IMPL(isInsideCoordsEx, shiftEJx);
        GRID_IS_INSIDE_IMPL(isInsideCoordsEy, shiftEJy);
        GRID_IS_INSIDE_IMPL(isInsideCoordsEz, shiftEJz);
        GRID_IS_INSIDE_IMPL(isInsideCoordsJx, shiftEJx);
        GRID_IS_INSIDE_IMPL(isInsideCoordsJy, shiftEJy);
        GRID_IS_INSIDE_IMPL(isInsideCoordsJz, shiftEJz);

        // the next methods interpolate and return field values

        /* returns CIC-interpolated field value in given coords, uses 'getFieldCIC' method */
        /* signature: 'forceinline FP funcname(const FP3& coords) const' */
        GRID_GET_FIELD_CIC_IMPL(getBxCIC, Bx, shiftBx);
        GRID_GET_FIELD_CIC_IMPL(getByCIC, By, shiftBy);
        GRID_GET_FIELD_CIC_IMPL(getBzCIC, Bz, shiftBz);
        GRID_GET_FIELD_CIC_IMPL(getExCIC, Ex, shiftEJx);
        GRID_GET_FIELD_CIC_IMPL(getEyCIC, Ey, shiftEJy);
        GRID_GET_FIELD_CIC_IMPL(getEzCIC, Ez, shiftEJz);
        GRID_GET_FIELD_CIC_IMPL(getJxCIC, Jx, shiftEJx);
        GRID_GET_FIELD_CIC_IMPL(getJyCIC, Jy, shiftEJy);
        GRID_GET_FIELD_CIC_IMPL(getJzCIC, Jz, shiftEJz);

        /* returns TSC-interpolated field value in given coords, uses 'getFieldTSC' method */
        /* signature: 'forceinline FP funcname(const FP3& coords) const' */
        GRID_GET_FIELD_TSC_IMPL(getBxTSC, Bx, shiftBx);
        GRID_GET_FIELD_TSC_IMPL(getByTSC, By, shiftBy);
        GRID_GET_FIELD_TSC_IMPL(getBzTSC, Bz, shiftBz);
        GRID_GET_FIELD_TSC_IMPL(getExTSC, Ex, shiftEJx);
        GRID_GET_FIELD_TSC_IMPL(getEyTSC, Ey, shiftEJy);
        GRID_GET_FIELD_TSC_IMPL(getEzTSC, Ez, shiftEJz);
        GRID_GET_FIELD_TSC_IMPL(getJxTSC, Jx, shiftEJx);
        GRID_GET_FIELD_TSC_IMPL(getJyTSC, Jy, shiftEJy);
        GRID_GET_FIELD_TSC_IMPL(getJzTSC, Jz, shiftEJz);

        /* returns PCS-interpolated field value in given coords, uses 'getFieldPCS' method */
        /* signature: 'forceinline FP funcname(const FP3& coords) const' */
        GRID_GET_FIELD_PCS_IMPL(getBxPCS, Bx, shiftBx);
        GRID_GET_FIELD_PCS_IMPL(getByPCS, By, shiftBy);
        GRID_GET_FIELD_PCS_IMPL(getBzPCS, Bz, shiftBz);
        GRID_GET_FIELD_PCS_IMPL(getExPCS, Ex, shiftEJx);
        GRID_GET_FIELD_PCS_IMPL(getEyPCS, Ey, shiftEJy);
        GRID_GET_FIELD_PCS_IMPL(getEzPCS, Ez, shiftEJz);
        GRID_GET_FIELD_PCS_IMPL(getJxPCS, Jx, shiftEJx);
        GRID_GET_FIELD_PCS_IMPL(getJyPCS, Jy, shiftEJy);
        GRID_GET_FIELD_PCS_IMPL(getJzPCS, Jz, shiftEJz);

        /* returns second order interpolated field value in given coords, uses 'getFieldSecondOrder' method */
        /* signature: 'forceinline FP funcname(const FP3& coords) const' */
        GRID_GET_FIELD_SECOND_ORDER_IMPL(getBxSecondOrder, Bx, shiftBx);
        GRID_GET_FIELD_SECOND_ORDER_IMPL(getBySecondOrder, By, shiftBy);
        GRID_GET_FIELD_SECOND_ORDER_IMPL(getBzSecondOrder, Bz, shiftBz);
        GRID_GET_FIELD_SECOND_ORDER_IMPL(getExSecondOrder, Ex, shiftEJx);
        GRID_GET_FIELD_SECOND_ORDER_IMPL(getEySecondOrder, Ey, shiftEJy);
        GRID_GET_FIELD_SECOND_ORDER_IMPL(getEzSecondOrder, Ez, shiftEJz);
        GRID_GET_FIELD_SECOND_ORDER_IMPL(getJxSecondOrder, Jx, shiftEJx);
        GRID_GET_FIELD_SECOND_ORDER_IMPL(getJySecondOrder, Jy, shiftEJy);
        GRID_GET_FIELD_SECOND_ORDER_IMPL(getJzSecondOrder, Jz, shiftEJz);

        /* returns fourth order interpolated field value in given coords, uses 'getFieldFourthOrder' method */
        /* signature: 'forceinline FP funcname(const FP3& coords) const' */
        GRID_GET_FIELD_FOURTH_ORDER_IMPL(getBxFourthOrder, Bx, shiftBx);
        GRID_GET_FIELD_FOURTH_ORDER_IMPL(getByFourthOrder, By, shiftBy);
        GRID_GET_FIELD_FOURTH_ORDER_IMPL(getBzFourthOrder, Bz, shiftBz);
        GRID_GET_FIELD_FOURTH_ORDER_IMPL(getExFourthOrder, Ex, shiftEJx);
        GRID_GET_FIELD_FOURTH_ORDER_IMPL(getEyFourthOrder, Ey, shiftEJy);
        GRID_GET_FIELD_FOURTH_ORDER_IMPL(getEzFourthOrder, Ez, shiftEJz);
        GRID_GET_FIELD_FOURTH_ORDER_IMPL(getJxFourthOrder, Jx, shiftEJx);
        GRID_GET_FIELD_FOURTH_ORDER_IMPL(getJyFourthOrder, Jy, shiftEJy);
        GRID_GET_FIELD_FOURTH_ORDER_IMPL(getJzFourthOrder, Jz, shiftEJz);

        /* return interpolated fields e, b in given coords */
        void getFieldsCIC(const FP3& coords, FP3& e, FP3& b) const {
            e = FP3(getExCIC(coords), getEyCIC(coords), getEzCIC(coords));
            b = FP3(getBxCIC(coords), getByCIC(coords), getBzCIC(coords));
        }
        void getFieldsTSC(const FP3& coords, FP3& e, FP3& b) const {
            e = FP3(getExTSC(coords), getEyTSC(coords), getEzTSC(coords));
            b = FP3(getBxTSC(coords), getByTSC(coords), getBzTSC(coords));
        }
        void getFieldsPCS(const FP3& coords, FP3& e, FP3& b) const {
            e = FP3(getExPCS(coords), getEyPCS(coords), getEzPCS(coords));
            b = FP3(getBxPCS(coords), getByPCS(coords), getBzPCS(coords));
        }
        void getFieldsSecondOrder(const FP3& coords, FP3& e, FP3& b) const {
            e = FP3(getExSecondOrder(coords), getEySecondOrder(coords), getEzSecondOrder(coords));
            b = FP3(getBxSecondOrder(coords), getBySecondOrder(coords), getBzSecondOrder(coords));
        }
        void getFieldsFourthOrder(const FP3& coords, FP3& e, FP3& b) const {
            e = FP3(getExFourthOrder(coords), getEyFourthOrder(coords), getEzFourthOrder(coords));
            b = FP3(getBxFourthOrder(coords), getByFourthOrder(coords), getBzFourthOrder(coords));
        }

        // the next methods interpolate and return field values, using default interpolation
        // this interpolation is written in the 'interpolationType' variable and can be changed
        // the start default interpolation method is CIC

        void setInterpolationType(InterpolationType type);
        InterpolationType getInterpolationType() const;

        /* returns interpolated field value in arbitrary coords */
        /* signature: 'forceinline FP funcname(const FP3& coords) const' */
        GRID_GET_FIELD_IMPL(getBx, interpolationBx);
        GRID_GET_FIELD_IMPL(getBy, interpolationBy);
        GRID_GET_FIELD_IMPL(getBz, interpolationBz);
        GRID_GET_FIELD_IMPL(getEx, interpolationEx);
        GRID_GET_FIELD_IMPL(getEy, interpolationEy);
        GRID_GET_FIELD_IMPL(getEz, interpolationEz);
        GRID_GET_FIELD_IMPL(getJx, interpolationJx);
        GRID_GET_FIELD_IMPL(getJy, interpolationJy);
        GRID_GET_FIELD_IMPL(getJz, interpolationJz);

        FP3 getB(const FP3& coords) const {
            return FP3(getBx(coords), getBy(coords), getBz(coords));
        }
        FP3 getE(const FP3& coords) const {
            return FP3(getEx(coords), getEy(coords), getEz(coords));
        }
        FP3 getJ(const FP3& coords) const {
            return FP3(getJx(coords), getJy(coords), getJz(coords));
        }

        void getFields(const FP3& coords, FP3& e, FP3& b) const
        {
            (this->*interpolationFields)(coords, e, b);
        }
        void getFields(FP x, FP y, FP z, FP3& e, FP3& b) const
        {
            getFields(FP3(x, y, z), e, b);
        }

        /* Make all current density values zero. */
        void zeroizeJ();

        /* Returns numCells with zeros where globalGridDims == 1 */
        const Int3 correctNumCellsAccordingToDim(const Int3& numCells) const {
            Int3 result = numCells;
            for (int d = 0; d < 3; d++)
                if (globalGridDims[d] == 1)
                    result[d] = 0;
            return result;
        }
        const Int3 getNumExternalLeftCells() const
        {
            const int nCells = LabelMethodRequiredNumberOfExternalCells<gridType_>::numExternalCells;
            return correctNumCellsAccordingToDim(Int3(nCells, nCells, nCells));
        }
        const Int3 getNumExternalRightCells() const
        {
            return getNumExternalLeftCells();
        }

        void save(std::ostream& ostr);
        void load(std::istream& istr);

        Int3 globalGridDims;  // global grid size, important to initialize it first
        FP3 steps;
        Int3 numInternalCells;
        Int3 numCells;
        Int3 sizeStorage;  // sometimes can be larger than numCells
        FP3 origin;
        int dimensionality = 0;

        ScalarField<Data> Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz;

        // 3d shifts of the field in the cell
        FP3 shiftBx, shiftBy, shiftBz, shiftEJx, shiftEJy, shiftEJz;

    private:

        /* returns grid index and normalized internal coords in [0, 0, 0]..(1, 1, 1) for
        given physical coords and shift. */
        void getGridCoords(const FP3& coords, const FP3& shift, Int3& idx,
            FP3& internalCoords) const
        {
            idx.x = (int)((coords.x - origin.x - shift.x) / steps.x);
            idx.y = (int)((coords.y - origin.y - shift.y) / steps.y);
            idx.z = (int)((coords.z - origin.z - shift.z) / steps.z);
            internalCoords = (coords - baseCoords(idx.x, idx.y, idx.z) - shift) / steps;
        }

        void getClosestGridCoords(const FP3& coords, const FP3& shift, Int3& idx,
            FP3& internalCoords) const
        {
            idx.x = (int)((coords.x - origin.x - shift.x) / steps.x + 0.5);
            idx.y = (int)((coords.y - origin.y - shift.y) / steps.y + 0.5);
            idx.z = (int)((coords.z - origin.z - shift.z) / steps.z + 0.5);
            internalCoords = (coords - baseCoords(idx.x, idx.y, idx.z) - shift) / steps;
        }

        /* returns base coords of element (i, j, k) so that its real coords are
        base coords + corresponding shift. */
        forceinline const FP3 baseCoords(int i, int j, int k) const
        {
            return origin + FP3(i, j, k) * steps;
        }

        /* returns true if coords is inside of the area that grid defines (without external cells)*/
        forceinline bool isInside(const FP3& coords, const FP3& shift) const
        {
            FP3 minCoords = origin + shift * steps;
            FP3 maxCoords = minCoords + (numCells - Int3(1, 1, 1)) * steps;
            return coords >= minCoords && coords < maxCoords;
        }

        void checkGridSizeAndOverlaps() {
            if (this->numInternalCells < this->getNumExternalLeftCells() + this->getNumExternalRightCells()) {
                std::string exc = "ERROR: grid size should be larger than both overlaps";
                std::cout << exc << std::endl;
                throw std::logic_error(exc);
            }
        }

        FP getFieldCIC(const FP3& coords, const ScalarField<Data>& field, const FP3& shift) const;
        FP getFieldTSC(const FP3& coords, const ScalarField<Data>& field, const FP3& shift) const;
        FP getFieldPCS(const FP3& coords, const ScalarField<Data>& field, const FP3& shift) const;
        FP getFieldSecondOrder(const FP3& coords, const ScalarField<Data>& field, const FP3& shift) const;
        FP getFieldFourthOrder(const FP3& coords, const ScalarField<Data>& field, const FP3& shift) const;

        InterpolationType interpolationType;
        void (Grid::* interpolationFields)(const FP3&, FP3&, FP3&) const;
        FP (Grid::* interpolationBx)(const FP3&) const;
        FP (Grid::* interpolationBy)(const FP3&) const;
        FP (Grid::* interpolationBz)(const FP3&) const;
        FP (Grid::* interpolationEx)(const FP3&) const;
        FP (Grid::* interpolationEy)(const FP3&) const;
        FP (Grid::* interpolationEz)(const FP3&) const;
        FP (Grid::* interpolationJx)(const FP3&) const;
        FP (Grid::* interpolationJy)(const FP3&) const;
        FP (Grid::* interpolationJz)(const FP3&) const;
    };

    typedef Grid<FP, GridTypes::YeeGridType> YeeGrid;
    typedef Grid<FP, GridTypes::StraightGridType> SimpleGrid;
    typedef Grid<FP, GridTypes::PSTDGridType> PSTDGrid;
    typedef Grid<FP, GridTypes::PSATDGridType> PSATDGrid;
    typedef Grid<FP, GridTypes::PSATDTimeStaggeredGridType> PSATDTimeStaggeredGrid;

    template<typename Data, GridTypes gridType_>
    inline Grid<Data, gridType_>::Grid(const Grid<Data, gridType_>& grid) :
        globalGridDims(grid.globalGridDims),
        steps(grid.steps),
        numInternalCells(grid.numInternalCells),
        numCells(grid.numCells),
        sizeStorage(grid.sizeStorage),
        shiftBx(grid.shiftBx), shiftBy(grid.shiftBy), shiftBz(grid.shiftBz),
        shiftEJx(grid.shiftEJx), shiftEJy(grid.shiftEJy), shiftEJz(grid.shiftEJz),
        origin(grid.origin),
        dimensionality(grid.dimensionality),
        Bx(grid.Bx), By(grid.By), Bz(grid.Bz),
        Ex(grid.Ex), Ey(grid.Ey), Ez(grid.Ez),
        Jx(grid.Jx), Jy(grid.Jy), Jz(grid.Jz)
    {
        checkGridSizeAndOverlaps();
        setInterpolationType(grid.interpolationType);
    }

    template <>
    inline Grid<FP, GridTypes::YeeGridType>::Grid(const Int3& _numCells, const FP3& minCoords,
        const FP3& _steps, const Int3& _globalGridDims) :
        globalGridDims(_globalGridDims),
        steps(_steps),
        numInternalCells(_numCells),
        numCells(numInternalCells + getNumExternalLeftCells() + getNumExternalRightCells()),
        sizeStorage(numCells),
        Ex(numCells, sizeStorage), Ey(numCells, sizeStorage), Ez(numCells, sizeStorage),
        Bx(numCells, sizeStorage), By(numCells, sizeStorage), Bz(numCells, sizeStorage),
        Jx(numCells, sizeStorage), Jy(numCells, sizeStorage), Jz(numCells, sizeStorage),
        shiftBx(FP3(0.5, 0, 0)* steps),
        shiftBy(FP3(0, 0.5, 0)* steps),
        shiftBz(FP3(0, 0, 0.5)* steps),
        shiftEJx(FP3(0, 0.5, 0.5)* steps),
        shiftEJy(FP3(0.5, 0, 0.5)* steps),
        shiftEJz(FP3(0.5, 0.5, 0)* steps),
        origin(minCoords.x - steps.x * getNumExternalLeftCells().x,
            minCoords.y - steps.y * getNumExternalLeftCells().y,
            minCoords.z - steps.z * getNumExternalLeftCells().z),
        dimensionality((_globalGridDims.x != 1) + (_globalGridDims.y != 1) + (_globalGridDims.z != 1))
    {
        checkGridSizeAndOverlaps();
        setInterpolationType(InterpolationType::Interpolation_CIC);
    }

    template<>
    inline Grid<FP, GridTypes::StraightGridType>::Grid(const Int3& _numInternalCells,
        const FP3& minCoords, const FP3& _steps, const Int3& _globalGridDims) :
        globalGridDims(_globalGridDims),
        steps(_steps),
        numInternalCells(_numInternalCells),
        numCells(numInternalCells + getNumExternalLeftCells() + getNumExternalRightCells()),
        sizeStorage(numCells),
        Bx(numCells, sizeStorage), By(numCells, sizeStorage), Bz(numCells, sizeStorage),
        Ex(numCells, sizeStorage), Ey(numCells, sizeStorage), Ez(numCells, sizeStorage),
        Jx(numCells, sizeStorage), Jy(numCells, sizeStorage), Jz(numCells, sizeStorage),
        shiftBx(FP3(0, 0, 0)* steps),
        shiftBy(FP3(0, 0, 0)* steps),
        shiftBz(FP3(0, 0, 0)* steps),
        shiftEJx(FP3(0, 0, 0)* steps),
        shiftEJy(FP3(0, 0, 0)* steps),
        shiftEJz(FP3(0, 0, 0)* steps),
        origin(minCoords.x - steps.x * getNumExternalLeftCells().x,
            minCoords.y - steps.y * getNumExternalLeftCells().y,
            minCoords.z - steps.z * getNumExternalLeftCells().z),
        dimensionality((_globalGridDims.x != 1) + (_globalGridDims.y != 1) + (_globalGridDims.z != 1))
    {
        checkGridSizeAndOverlaps();
        setInterpolationType(InterpolationType::Interpolation_CIC);
    }


    template<>
    inline Grid<FP, GridTypes::PSTDGridType>::Grid(const Int3& _numInternalCells,
        const FP3& minCoords, const FP3& _steps, const Int3& _globalGridDims) :
        globalGridDims(_globalGridDims),
        steps(_steps),
        numInternalCells(_numInternalCells),
        numCells(numInternalCells),
        sizeStorage(Int3(numCells.x, numCells.y, 2 * (numCells.z / 2 + 1))),
        Bx(numCells, sizeStorage), By(numCells, sizeStorage), Bz(numCells, sizeStorage),
        Ex(numCells, sizeStorage), Ey(numCells, sizeStorage), Ez(numCells, sizeStorage),
        Jx(numCells, sizeStorage), Jy(numCells, sizeStorage), Jz(numCells, sizeStorage),
        shiftBx(FP3(0, 0, 0)* steps),
        shiftBy(FP3(0, 0, 0)* steps),
        shiftBz(FP3(0, 0, 0)* steps),
        shiftEJx(FP3(0, 0, 0)* steps),
        shiftEJy(FP3(0, 0, 0)* steps),
        shiftEJz(FP3(0, 0, 0)* steps),
        origin(minCoords),
        dimensionality((_globalGridDims.x != 1) + (_globalGridDims.y != 1) + (_globalGridDims.z != 1))
    {
        checkGridSizeAndOverlaps();
        setInterpolationType(InterpolationType::Interpolation_CIC);
    }

    template<>
    inline Grid<FP, GridTypes::PSATDGridType>::Grid(const Int3& _numInternalCells,
        const FP3& minCoords, const FP3& _steps, const Int3& _globalGridDims) :
        globalGridDims(_globalGridDims),
        steps(_steps),
        numInternalCells(_numInternalCells),
        numCells(numInternalCells),
        sizeStorage(Int3(numCells.x, numCells.y, 2 * (numCells.z / 2 + 1))),
        Bx(numCells, sizeStorage), By(numCells, sizeStorage), Bz(numCells, sizeStorage),
        Ex(numCells, sizeStorage), Ey(numCells, sizeStorage), Ez(numCells, sizeStorage),
        Jx(numCells, sizeStorage), Jy(numCells, sizeStorage), Jz(numCells, sizeStorage),
        shiftBx(FP3(0, 0, 0)* steps),
        shiftBy(FP3(0, 0, 0)* steps),
        shiftBz(FP3(0, 0, 0)* steps),
        shiftEJx(FP3(0, 0, 0)* steps),
        shiftEJy(FP3(0, 0, 0)* steps),
        shiftEJz(FP3(0, 0, 0)* steps),
        origin(minCoords),
        dimensionality((_globalGridDims.x != 1) + (_globalGridDims.y != 1) + (_globalGridDims.z != 1))
    {
        checkGridSizeAndOverlaps();
        setInterpolationType(InterpolationType::Interpolation_CIC);
    }

    template<>
    inline Grid<FP, GridTypes::PSATDTimeStaggeredGridType>::Grid(const Int3& _numInternalCells,
        const FP3& minCoords, const FP3& _steps, const Int3& _globalGridDims) :
        globalGridDims(_globalGridDims),
        steps(_steps),
        numInternalCells(_numInternalCells),
        numCells(numInternalCells),
        sizeStorage(Int3(numCells.x, numCells.y, 2 * (numCells.z / 2 + 1))),
        Bx(numCells, sizeStorage), By(numCells, sizeStorage), Bz(numCells, sizeStorage),
        Ex(numCells, sizeStorage), Ey(numCells, sizeStorage), Ez(numCells, sizeStorage),
        Jx(numCells, sizeStorage), Jy(numCells, sizeStorage), Jz(numCells, sizeStorage),
        shiftBx(FP3(0, 0, 0)* steps),
        shiftBy(FP3(0, 0, 0)* steps),
        shiftBz(FP3(0, 0, 0)* steps),
        shiftEJx(FP3(0, 0, 0)* steps),
        shiftEJy(FP3(0, 0, 0)* steps),
        shiftEJz(FP3(0, 0, 0)* steps),
        origin(minCoords),
        dimensionality((_globalGridDims.x != 1) + (_globalGridDims.y != 1) + (_globalGridDims.z != 1))
    {
        checkGridSizeAndOverlaps();
        setInterpolationType(InterpolationType::Interpolation_CIC);
    }

    template <typename Data, GridTypes gT>
    template <class TGrid>
    inline void Grid<Data, gT>::copyValues(TGrid* grid) const {
        const int nx = grid->numCells.x, ny = grid->numCells.y, nz = grid->numCells.z;

        OMP_FOR_COLLAPSE()
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                OMP_SIMD()
                for (int k = 0; k < nz; k++) {
                    grid->Bx(i, j, k) = this->getBx(grid->BxPosition(i, j, k));
                    grid->By(i, j, k) = this->getBy(grid->ByPosition(i, j, k));
                    grid->Bz(i, j, k) = this->getBz(grid->BzPosition(i, j, k));

                    grid->Ex(i, j, k) = this->getEx(grid->ExPosition(i, j, k));
                    grid->Ey(i, j, k) = this->getEy(grid->EyPosition(i, j, k));
                    grid->Ez(i, j, k) = this->getEz(grid->EzPosition(i, j, k));

                    grid->Jx(i, j, k) = this->getJx(grid->JxPosition(i, j, k));
                    grid->Jy(i, j, k) = this->getJy(grid->JyPosition(i, j, k));
                    grid->Jz(i, j, k) = this->getJz(grid->JzPosition(i, j, k));
                }
    }

    template< typename Data, GridTypes gT>
    inline FP Grid<Data, gT>::getFieldCIC(const FP3& coords, const ScalarField<Data>& field, const FP3& shift) const
    {
        Int3 idx;
        FP3 internalCoords;
        getGridCoords(coords, shift, idx, internalCoords);
        return field.interpolateCIC(idx, internalCoords);
    }

    template< typename Data, GridTypes gT>
    inline FP Grid<Data, gT>::getFieldTSC(const FP3& coords, const ScalarField<Data>& field, const FP3& shift) const
    {
        Int3 idx;
        FP3 internalCoords;
        getClosestGridCoords(coords, shift, idx, internalCoords);
        return field.interpolateTSC(idx, internalCoords);
    }

    template< typename Data, GridTypes gT>
    inline FP Grid<Data, gT>::getFieldSecondOrder(const FP3& coords, const ScalarField<Data>& field, const FP3& shift) const
    {
        Int3 idx;
        FP3 internalCoords;
        getClosestGridCoords(coords, shift, idx, internalCoords);
        return field.interpolateSecondOrder(idx, internalCoords);
    }

    template< typename Data, GridTypes gT>
    inline FP Grid<Data, gT>::getFieldFourthOrder(const FP3& coords, const ScalarField<Data>& field, const FP3& shift) const
    {
        Int3 idx;
        FP3 internalCoords;
        getClosestGridCoords(coords, shift, idx, internalCoords);
        return field.interpolateFourthOrder(idx, internalCoords);
    }

    template< typename Data, GridTypes gT>
    inline FP Grid<Data, gT>::getFieldPCS(const FP3& coords, const ScalarField<Data>& field, const FP3& shift) const
    {
        Int3 idx;
        FP3 internalCoords;
        getGridCoords(coords, shift, idx, internalCoords);
        return field.interpolatePCS(idx, internalCoords);
    }

    template<typename Data, GridTypes gT>
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
        case InterpolationType::Interpolation_CIC:
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
        case InterpolationType::Interpolation_TSC:
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
        case InterpolationType::Interpolation_PCS:
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
        case InterpolationType::Interpolation_SecondOrder:
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
        case InterpolationType::Interpolation_FourthOrder:
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

    template<typename Data, GridTypes gridType_>
    inline void Grid<Data, gridType_>::save(std::ostream& ostr)
    {
        int savedGridType = (int)this->gridType;
        ostr.write((char*)&savedGridType, sizeof(savedGridType));

        ostr.write((char*)&globalGridDims, sizeof(globalGridDims));
        ostr.write((char*)&steps, sizeof(steps));
        ostr.write((char*)&numInternalCells, sizeof(numInternalCells));
        ostr.write((char*)&numCells, sizeof(numCells));
        ostr.write((char*)&sizeStorage, sizeof(sizeStorage));
        ostr.write((char*)&origin, sizeof(origin));
        ostr.write((char*)&dimensionality, sizeof(dimensionality));

        ostr.write((char*)&shiftBx, sizeof(shiftBx));
        ostr.write((char*)&shiftBy, sizeof(shiftBy));
        ostr.write((char*)&shiftBz, sizeof(shiftBz));
        ostr.write((char*)&shiftEJx, sizeof(shiftEJx));
        ostr.write((char*)&shiftEJy, sizeof(shiftEJy));
        ostr.write((char*)&shiftEJz, sizeof(shiftEJz));

        Bx.save(ostr);
        By.save(ostr);
        Bz.save(ostr);
        Ex.save(ostr);
        Ey.save(ostr);
        Ez.save(ostr);
        Jx.save(ostr);
        Jy.save(ostr);
        Jz.save(ostr);
    }

    template<typename Data, GridTypes gridType_>
    inline void Grid<Data, gridType_>::load(std::istream& istr)
    {
        int loadedGridType = -1;
        istr.read((char*)&loadedGridType, sizeof(loadedGridType));
        if (loadedGridType != (int)this->gridType)
            throw "ERROR: types of loaded grids do not match";

        istr.read((char*)&globalGridDims, sizeof(globalGridDims));
        istr.read((char*)&steps, sizeof(steps));
        istr.read((char*)&numInternalCells, sizeof(numInternalCells));
        istr.read((char*)&numCells, sizeof(numCells));
        istr.read((char*)&sizeStorage, sizeof(sizeStorage));
        istr.read((char*)&origin, sizeof(origin));
        istr.read((char*)&dimensionality, sizeof(dimensionality));

        istr.read((char*)&shiftBx, sizeof(shiftBx));
        istr.read((char*)&shiftBy, sizeof(shiftBy));
        istr.read((char*)&shiftBz, sizeof(shiftBz));
        istr.read((char*)&shiftEJx, sizeof(shiftEJx));
        istr.read((char*)&shiftEJy, sizeof(shiftEJy));
        istr.read((char*)&shiftEJz, sizeof(shiftEJz));

        Bx.load(istr);
        By.load(istr);
        Bz.load(istr);
        Ex.load(istr);
        Ey.load(istr);
        Ez.load(istr);
        Jx.load(istr);
        Jy.load(istr);
        Jz.load(istr);
    }
}
