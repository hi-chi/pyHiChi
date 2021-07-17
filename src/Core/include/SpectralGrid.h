#pragma once

#include "macros.h"

#include "GridTypes.h"
#include "ScalarField.h"
#include "Vectors.h"
#include "Constants.h"
#include "Grid.h"

namespace pfc {

    // wrapper of ScalarField, has common memory
    template <class RealData, class SpectralData>
    class SpectralScalarField {
    public:

        SpectralScalarField(ScalarField<RealData>* scalarField, Int3 size) :
            scalarField(scalarField), size(size) {
            this->raw = reinterpret_cast<SpectralData*>(scalarField->getData());
        }

        SpectralData* getData() {
            return raw;
        }

        Int3 getSize() const {
            return size;
        }

        /* Read-only access by scalar indexes */
        SpectralData operator()(int i, int j, int k) const
        {
            return raw[k + (j + i * size.y) * size.z];
        }

        /* Read-write access by scalar indexes */
        SpectralData& operator()(int i, int j, int k)
        {
            return raw[k + (j + i * size.y) * size.z];
        }

        /* Read-only access by vector index */
        SpectralData operator()(const Int3& index) const
        {
            return raw[index.z + (index.y + index.x * size.y) * size.z];
        }

        /* Read-write access by vector index */
        SpectralData& operator()(const Int3& index)
        {
            return raw[index.z + (index.y + index.x * size.y) * size.z];
        }

    private:

        ScalarField<RealData>* scalarField;
        SpectralData* raw;
        Int3 size;
    };


    // grid that shares memory with main grid
    template <class RealData, class SpectralData>
    class SpectralGrid {
    public:

        template <class GridType>
        SpectralGrid(const Int3& _numCells, const Int3& globalGridDims,
            GridType* grid) {}

        Int3 globalGridDims;
        Int3 numCells;
        Int3 sizeStorage;
        int dimensionality;

        SpectralScalarField<RealData, SpectralData> Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz;
    };

    template<>
    template <class GridType>
    inline SpectralGrid<FP, complexFP>::SpectralGrid(const Int3& _numCells,
        const Int3& _globalGridDims, GridType* grid) :
        globalGridDims(_globalGridDims),
        numCells(_numCells),
        sizeStorage(_numCells),
        Ex(&grid->Ex, sizeStorage), Ey(&grid->Ey, sizeStorage), Ez(&grid->Ez, sizeStorage),
        Bx(&grid->Bx, sizeStorage), By(&grid->By, sizeStorage), Bz(&grid->Bz, sizeStorage),
        Jx(&grid->Jx, sizeStorage), Jy(&grid->Jy, sizeStorage), Jz(&grid->Jz, sizeStorage),
        dimensionality((_globalGridDims.x != 1) + (_globalGridDims.y != 1) + (_globalGridDims.z != 1))
    {}

}