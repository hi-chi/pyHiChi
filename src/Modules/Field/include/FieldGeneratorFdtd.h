#pragma once
#include "FieldGenerator.h"

namespace pfc
{
    class FieldGeneratorFdtd : public FieldGenerator<YeeGrid>
    {
    public:

        using FunctionType = typename FieldGenerator<YeeGrid>::FunctionType;

        FieldGeneratorFdtd(YeeGrid* grid, FP dt,
            const Int3& domainIndexBegin, const Int3& domainIndexEnd,
            const Int3& leftGenIndex, const Int3& rightGenIndex,
            FunctionType bxFunc, FunctionType byFunc, FunctionType bzFunc,
            FunctionType exFunc, FunctionType eyFunc, FunctionType ezFunc,
            const Int3& isLeftBorderEnabled = Int3(1, 1, 1),
            const Int3& isRightBorderEnabled = Int3(1, 1, 1)) :
            FieldGenerator(grid, dt, domainIndexBegin, domainIndexEnd, 
                leftGenIndex, rightGenIndex,
                bxFunc, byFunc, bzFunc, exFunc, eyFunc, ezFunc,
                isLeftBorderEnabled, isRightBorderEnabled)
        {}

        FieldGeneratorFdtd(YeeGrid* grid, FP dt,
            const Int3& domainIndexBegin, const Int3& domainIndexEnd,
            const Int3& leftGenIndex, const Int3& rightGenIndex,
            const std::array<FunctionType, 3>& xLeftBFunc,  // { bx, by, bz }
            const std::array<FunctionType, 3>& xRightBFunc,
            const std::array<FunctionType, 3>& yLeftBFunc,
            const std::array<FunctionType, 3>& yRightBFunc,
            const std::array<FunctionType, 3>& zLeftBFunc,
            const std::array<FunctionType, 3>& zRightBFunc,
            const std::array<FunctionType, 3>& xLeftEFunc,  // { ex, ey, ez }
            const std::array<FunctionType, 3>& xRightEFunc,
            const std::array<FunctionType, 3>& yLeftEFunc,
            const std::array<FunctionType, 3>& yRightEFunc,
            const std::array<FunctionType, 3>& zLeftEFunc,
            const std::array<FunctionType, 3>& zRightEFunc,
            const Int3& isLeftBorderEnabled = Int3(1, 1, 1),
            const Int3& isRightBorderEnabled = Int3(1, 1, 1)) :
            FieldGenerator(grid, dt, domainIndexBegin, domainIndexEnd, 
                leftGenIndex, rightGenIndex,
                xLeftBFunc, xRightBFunc, yLeftBFunc, yRightBFunc, zLeftBFunc, zRightBFunc,
                xLeftEFunc, xRightEFunc, yLeftEFunc, yRightEFunc, zLeftEFunc, zRightEFunc,
                isLeftBorderEnabled, isRightBorderEnabled)
        {}

        FieldGeneratorFdtd(YeeGrid* grid, FP dt,
            const Int3& domainIndexBegin, const Int3& domainIndexEnd,
            const Int3& leftGenIndex, const Int3& rightGenIndex,
            const std::array<std::array<std::array<FunctionType, 3>, 3>, 2>& bFunc,
            const std::array<std::array<std::array<FunctionType, 3>, 3>, 2>& eFunc,
            const Int3& isLeftBorderEnabled = Int3(1, 1, 1),
            const Int3& isRightBorderEnabled = Int3(1, 1, 1)) :
            FieldGenerator(grid, dt, domainIndexBegin, domainIndexEnd,
                leftGenIndex, rightGenIndex, bFunc, eFunc,
                isLeftBorderEnabled, isRightBorderEnabled)
        {}

        FieldGeneratorFdtd(YeeGrid* grid, FP dt, const Int3& domainIndexBegin,
            const Int3& domainIndexEnd, const FieldGeneratorFdtd& gen) :
            FieldGenerator(grid, dt, domainIndexBegin, domainIndexEnd, gen)
        {}

        // constructor for loading
        FieldGeneratorFdtd(YeeGrid* grid, FP dt,
            const Int3& domainIndexBegin, const Int3& domainIndexEnd) :
            FieldGenerator(grid, dt, domainIndexBegin, domainIndexEnd)
        {}

        void generateB(FP time);
        void generateE(FP time);

    protected:

        FP3 getBCurrents(int side, int dim0,
            const FP3& exCoords, const FP3& eyCoords, const FP3& ezCoords, FP time) const;
        FP3 getECurrents(int side, int dim0,
            const FP3& bxCoords, const FP3& byCoords, const FP3& bzCoords, FP time) const;

        std::vector<Int3> getBGridIndices(int i, int j, int k,
            int dim0, int dim1, int dim2) const;
        std::vector<Int3> getEGridIndices(int i, int j, int k,
            int dim0, int dim1, int dim2) const;

    };

    inline void FieldGeneratorFdtd::generateB(FP time)
    {
        const FP cdt = constants::c * this->dt;
        const FP3 norm_coeffs = FP3(cdt, cdt, cdt) / this->grid->steps;

        for (int dim0 = 0; dim0 < this->grid->dimensionality; dim0++)
        {
            int dim1 = (dim0 + 1) % 3;
            int dim2 = (dim0 + 2) % 3;

            // TODO: compute local index inside a domain
            int genIndexDim0[2] = {
                leftGeneratorIndex[dim0],
                rightGeneratorIndex[dim0]
            };

            int begin1 = isLeftBorderEnabled[dim1] ?
                leftGeneratorIndex[dim1] : this->domainIndexBegin[dim1];
            int begin2 = isLeftBorderEnabled[dim2] ?
                leftGeneratorIndex[dim2] : this->domainIndexBegin[dim2];
            int end1 = isRightBorderEnabled[dim1] ?
                rightGeneratorIndex[dim1] : this->domainIndexEnd[dim1];
            int end2 = isRightBorderEnabled[dim2] ?
                rightGeneratorIndex[dim2] : this->domainIndexEnd[dim2];

            Int3 isBorderEnabled[2] = { isLeftBorderEnabled, isRightBorderEnabled };

            for (int side = 0; side < 2; side++) {
                if (!isBorderEnabled[side][dim0]) continue;

                OMP_FOR_COLLAPSE()
                for (int j = begin1; j < end1; j++)
                    for (int k = begin2; k < end2; k++)
                    {
                        std::vector<Int3> eIndices = getEGridIndices(genIndexDim0[side], j, k, dim0, dim1, dim2);
                        std::vector<Int3> bIndices = getBGridIndices(genIndexDim0[side], j, k, dim0, dim1, dim2);

                        FP3 exCoords = this->grid->ExPosition(eIndices[0].x, eIndices[0].y, eIndices[0].z);
                        FP3 eyCoords = this->grid->EyPosition(eIndices[1].x, eIndices[1].y, eIndices[1].z);
                        FP3 ezCoords = this->grid->EzPosition(eIndices[2].x, eIndices[2].y, eIndices[2].z);

                        FP3 current = getBCurrents(side, dim0, exCoords, eyCoords, ezCoords, time);

                        this->grid->Bx(bIndices[0]) += norm_coeffs[dim0] * current.x;
                        this->grid->By(bIndices[1]) += norm_coeffs[dim0] * current.y;
                        this->grid->Bz(bIndices[2]) += norm_coeffs[dim0] * current.z;
                    }
            }
        }
    }

    inline void FieldGeneratorFdtd::generateE(FP time)
    {
        const FP cdt = constants::c * this->dt;
        const FP3 norm_coeffs = FP3(cdt, cdt, cdt) / this->grid->steps;

        for (int dim0 = 0; dim0 < this->grid->dimensionality; dim0++)
        {
            int dim1 = (dim0 + 1) % 3;
            int dim2 = (dim0 + 2) % 3;

            // TODO: compute local index inside a domain
            int genIndexDim0[2] = {
                leftGeneratorIndex[dim0],
                rightGeneratorIndex[dim0]
            };

            int begin1 = isLeftBorderEnabled[dim1] ?
                leftGeneratorIndex[dim1] : this->domainIndexBegin[dim1];
            int begin2 = isLeftBorderEnabled[dim2] ?
                leftGeneratorIndex[dim2] : this->domainIndexBegin[dim2];
            int end1 = isRightBorderEnabled[dim1] ?
                rightGeneratorIndex[dim1] : this->domainIndexEnd[dim1];
            int end2 = isRightBorderEnabled[dim2] ?
                rightGeneratorIndex[dim2] : this->domainIndexEnd[dim2];

            Int3 isBorderEnabled[2] = { isLeftBorderEnabled, isRightBorderEnabled };

            for (int side = 0; side < 2; side++) {
                if (!isBorderEnabled[side][dim0]) continue;

                OMP_FOR_COLLAPSE()
                for (int j = begin1; j < end1; j++)
                    for (int k = begin2; k < end2; k++)
                    {
                        std::vector<Int3> eIndices = getEGridIndices(genIndexDim0[side], j, k, dim0, dim1, dim2);
                        std::vector<Int3> bIndices = getBGridIndices(genIndexDim0[side], j, k, dim0, dim1, dim2);

                        FP3 bxCoords = this->grid->BxPosition(bIndices[0].x, bIndices[0].y, bIndices[0].z);
                        FP3 byCoords = this->grid->ByPosition(bIndices[1].x, bIndices[1].y, bIndices[1].z);
                        FP3 bzCoords = this->grid->BzPosition(bIndices[2].x, bIndices[2].y, bIndices[2].z);

                        FP3 current = getECurrents(side, dim0, bxCoords, byCoords, bzCoords, time);

                        this->grid->Ex(eIndices[0]) += norm_coeffs[dim0] * current.x;
                        this->grid->Ey(eIndices[1]) += norm_coeffs[dim0] * current.y;
                        this->grid->Ez(eIndices[2]) += norm_coeffs[dim0] * current.z;
                    }
            }
        }
    }

    inline FP3 FieldGeneratorFdtd::getBCurrents(int side, int dim0,
        const FP3& exCoords, const FP3& eyCoords, const FP3& ezCoords, FP time) const
    {
        FP3 current, normal, genField;
        normal[dim0] = (side == 0) ? 1.0 : -1.0;
        genField = FP3(
            eFunc[side][dim0][0](exCoords.x, exCoords.y, exCoords.z, time),
            eFunc[side][dim0][1](eyCoords.x, eyCoords.y, eyCoords.z, time),
            eFunc[side][dim0][2](ezCoords.x, ezCoords.y, ezCoords.z, time)
        );
        current = cross(normal, genField);
        return current;
    }

    inline FP3 FieldGeneratorFdtd::getECurrents(int side, int dim0,
        const FP3& bxCoords, const FP3& byCoords, const FP3& bzCoords, FP time) const
    {
        FP3 current, normal, genField;
        normal[dim0] = (side == 0) ? 1.0 : -1.0;
        genField = FP3(
            bFunc[side][dim0][0](bxCoords.x, bxCoords.y, bxCoords.z, time),
            bFunc[side][dim0][1](byCoords.x, byCoords.y, byCoords.z, time),
            bFunc[side][dim0][2](bzCoords.x, bzCoords.y, bzCoords.z, time)
        );
        current = cross(normal, genField) * (-1.0);
        return current;
    }

    inline std::vector<Int3> FieldGeneratorFdtd::getBGridIndices(int i, int j, int k,
        int dim0, int dim1, int dim2) const {
        Int3 index;
        index[dim0] = i;
        index[dim1] = j;
        index[dim2] = k;

        std::vector<Int3> indices = { index, index, index };

        return indices;
    }

    inline std::vector<Int3> FieldGeneratorFdtd::getEGridIndices(int i, int j, int k,
        int dim0, int dim1, int dim2) const {
        Int3 index;
        index[dim0] = i;
        index[dim1] = j;
        index[dim2] = k;

        std::vector<Int3> indices = { index, index, index };
        indices[dim1][dim0]--;
        indices[dim2][dim0]--;

        return indices;
    }
}
