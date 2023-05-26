#pragma once
#include "FieldGenerator.h"

namespace pfc
{
    class FieldGeneratorFdtd : public FieldGenerator<GridTypes::YeeGridType>
    {
    public:

        FieldGeneratorFdtd(FieldSolver<GridTypes::YeeGridType>* fieldSolver) :
            FieldGenerator(fieldSolver) {}

        void generateB() override;
        void generateE() override;

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

    inline void FieldGeneratorFdtd::generateB()
    {
        Grid<FP, GridTypes::YeeGridType>* grid = fieldSolver->grid;

        const FP time = this->fieldSolver->globalTime;

        const FP cdt = constants::c * fieldSolver->dt;
        const FP3 norm_coeffs = FP3(cdt, cdt, cdt) / grid->steps;

        for (int dim0 = 0; dim0 < grid->dimensionality; dim0++)
        {
            int dim1 = (dim0 + 1) % 3;
            int dim2 = (dim0 + 2) % 3;

            // TODO: compute local index inside a domain
            int generatorIndex[2] = {
                generatorGridIndex[(int)SideEnum::LEFT][dim0],
                generatorGridIndex[(int)SideEnum::RIGHT][dim0]
            };

            int begin1 = isBorderEnabled[(int)SideEnum::LEFT][dim1] ?
                generatorGridIndex[(int)SideEnum::LEFT][dim1] : fieldSolver->internalBAreaBegin[dim1];
            int begin2 = isBorderEnabled[(int)SideEnum::LEFT][dim2] ?
                generatorGridIndex[(int)SideEnum::LEFT][dim2] : fieldSolver->internalBAreaBegin[dim2];
            int end1 = isBorderEnabled[(int)SideEnum::RIGHT][dim1] ?
                generatorGridIndex[(int)SideEnum::RIGHT][dim1] : fieldSolver->internalBAreaEnd[dim1];
            int end2 = isBorderEnabled[(int)SideEnum::RIGHT][dim2] ?
                generatorGridIndex[(int)SideEnum::RIGHT][dim2] : fieldSolver->internalBAreaEnd[dim2];

            for (int side = 0; side < 2; side++) {
                if (!isBorderEnabled[side][dim0]) continue;

                OMP_FOR_COLLAPSE()
                for (int j = begin1; j < end1; j++)
                    for (int k = begin2; k < end2; k++)
                    {
                        std::vector<Int3> eIndices = getEGridIndices(generatorIndex[side], j, k, dim0, dim1, dim2);
                        std::vector<Int3> bIndices = getBGridIndices(generatorIndex[side], j, k, dim0, dim1, dim2);

                        FP3 exCoords = grid->ExPosition(eIndices[0].x, eIndices[0].y, eIndices[0].z);
                        FP3 eyCoords = grid->EyPosition(eIndices[1].x, eIndices[1].y, eIndices[1].z);
                        FP3 ezCoords = grid->EzPosition(eIndices[2].x, eIndices[2].y, eIndices[2].z);

                        FP3 current = getBCurrents(side, dim0, exCoords, eyCoords, ezCoords, time);

                        grid->Bx(bIndices[0]) += norm_coeffs[dim0] * current.x;
                        grid->By(bIndices[1]) += norm_coeffs[dim0] * current.y;
                        grid->Bz(bIndices[2]) += norm_coeffs[dim0] * current.z;
                    }
            }
        }
    }

    inline void FieldGeneratorFdtd::generateE()
    {
        Grid<FP, GridTypes::YeeGridType>* grid = fieldSolver->grid;

        const FP time = this->fieldSolver->globalTime +
            (grid->ifFieldsTimeStraggered ? this->fieldSolver->dt * 0.5 : 0.0);

        const FP cdt = constants::c * fieldSolver->dt;
        const FP3 norm_coeffs = FP3(cdt, cdt, cdt) / grid->steps;

        for (int dim0 = 0; dim0 < grid->dimensionality; dim0++)
        {
            int dim1 = (dim0 + 1) % 3;
            int dim2 = (dim0 + 2) % 3;

            // TODO: compute local index inside a domain
            int generatorIndex[2] = {
                generatorGridIndex[(int)SideEnum::LEFT][dim0],
                generatorGridIndex[(int)SideEnum::RIGHT][dim0]
            };

            int begin1 = isBorderEnabled[(int)SideEnum::LEFT][dim1] ?
                generatorGridIndex[(int)SideEnum::LEFT][dim1] : fieldSolver->internalEAreaBegin[dim1];
            int begin2 = isBorderEnabled[(int)SideEnum::LEFT][dim2] ?
                generatorGridIndex[(int)SideEnum::LEFT][dim2] : fieldSolver->internalEAreaBegin[dim2];
            int end1 = isBorderEnabled[(int)SideEnum::RIGHT][dim1] ?
                generatorGridIndex[(int)SideEnum::RIGHT][dim1] : fieldSolver->internalEAreaEnd[dim1];
            int end2 = isBorderEnabled[(int)SideEnum::RIGHT][dim2] ?
                generatorGridIndex[(int)SideEnum::RIGHT][dim2] : fieldSolver->internalEAreaEnd[dim2];

            for (int side = 0; side < 2; side++) {
                if (!isBorderEnabled[side][dim0]) continue;

                OMP_FOR_COLLAPSE()
                for (int j = begin1; j < end1; j++)
                    for (int k = begin2; k < end2; k++)
                    {
                        std::vector<Int3> eIndices = getEGridIndices(generatorIndex[side], j, k, dim0, dim1, dim2);
                        std::vector<Int3> bIndices = getBGridIndices(generatorIndex[side], j, k, dim0, dim1, dim2);

                        FP3 bxCoords = grid->BxPosition(bIndices[0].x, bIndices[0].y, bIndices[0].z);
                        FP3 byCoords = grid->ByPosition(bIndices[1].x, bIndices[1].y, bIndices[1].z);
                        FP3 bzCoords = grid->BzPosition(bIndices[2].x, bIndices[2].y, bIndices[2].z);

                        FP3 current = getECurrents(side, dim0, bxCoords, byCoords, bzCoords, time);

                        grid->Ex(eIndices[0]) += norm_coeffs[dim0] * current.x;
                        grid->Ey(eIndices[1]) += norm_coeffs[dim0] * current.y;
                        grid->Ez(eIndices[2]) += norm_coeffs[dim0] * current.z;
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
