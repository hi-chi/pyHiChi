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
            SpectralScalarField(scalarField, size,
                { size.x, size.y, size.z * 2 })
        {}

        SpectralScalarField(ScalarField<RealData>* scalarField, Int3 size, Int3 sizeStorage) :
            scalarField(scalarField), size(size), sizeStorage(sizeStorage)
        {
            this->raw = reinterpret_cast<SpectralData*>(scalarField->getData());
        }

        SpectralData* getData()
        {
            return raw;
        }

        Int3 getSize() const
        {
            return size;
        }

        Int3 getMemSize() const  // returns memory size in FP
        {
            return sizeStorage;
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
        Int3 size;  // spectral size
        Int3 sizeStorage;   // at the moment FP memory size = sizeStorage.volume() == size.volume()*2
    };


    // grid that shares memory with main grid
    template <class RealData, class SpectralData>
    class SpectralGrid {
    public:

        template <class GridType>
        SpectralGrid(GridType* grid, const Int3& _numSpectralCells) {}

        Int3 numCells;  // spectral numCells
        Int3 sizeStorage;  // in real type memory cells

        SpectralScalarField<RealData, SpectralData> Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz;
    };

    template<>
    template <class GridType>
    inline SpectralGrid<FP, complexFP>::SpectralGrid(GridType* grid,
        const Int3& _numCells) :
        numCells(_numCells),
        sizeStorage(grid->sizeStorage),
        Ex(&grid->Ex, numCells, sizeStorage),
        Ey(&grid->Ey, numCells, sizeStorage),
        Ez(&grid->Ez, numCells, sizeStorage),
        Bx(&grid->Bx, numCells, sizeStorage),
        By(&grid->By, numCells, sizeStorage),
        Bz(&grid->Bz, numCells, sizeStorage),
        Jx(&grid->Jx, numCells, sizeStorage),
        Jy(&grid->Jy, numCells, sizeStorage),
        Jz(&grid->Jz, numCells, sizeStorage)
    {}

}