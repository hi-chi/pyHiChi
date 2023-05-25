#pragma once
#include <vector>

#include "Grid.h"
#include "FieldSolver.h"
#include "Vectors.h"

namespace pfc {
    template<GridTypes gridTypes>
    class FieldSolver;
    template<GridTypes gridTypes>
    class RealFieldSolver;
    template<GridTypes gridTypes>
    class SpectralFieldSolver;

    template<GridTypes gridTypes>
    class Pml
    {
    public:
        Pml(FieldSolver<gridTypes>* _fieldSolver, Int3 _sizePml,
            FP nPmlParam = (FP)4.0, FP r0PmlParam = (FP)1e-8);

        virtual void updateB() {};
        virtual void updateE() {};

        Int3 sizePml;
        Int3 leftPmlBorder, rightPmlBorder;
        std::vector<FP> exy, exz, eyx, eyz, ezx, ezy;  // split electric field
        std::vector<FP> bxy, bxz, byx, byz, bzx, bzy;  // split magnetic field

        std::vector<Int3> eIndex, bIndex;  // natural 3d indexes of nodes in PML

        forceinline FP computeSigma(FP coord, CoordinateEnum axis) const;

    protected:

        FieldSolver<gridTypes>* fieldSolver;

        // absorbing function parameters
        FP3 maxSigma;
        FP n;

        void initializePmlParams(FP _n, FP _r0, const FP3& gridStep);

        void initializeSplitFieldsE(int sizeE);
        void initializeSplitFieldsB(int sizeB);

        void fillPmlIndices();

        void checkGridAndPmlSize(const Int3& globalGridDims) {
            if (2 * sizePml > globalGridDims) {
                std::string exc = "ERROR: grid size should be larger than 2x pml size";
                std::cout << exc << std::endl;
                throw std::exception(exc.c_str());
            }
        }

    };

    template<GridTypes gridTypes>
    void Pml<gridTypes>::initializePmlParams(FP _n, FP _r0, const FP3& gridStep) {
        this->n = _n;

        for (int d = 0; d < 3; d++) {
            if (sizePml[d])
                this->maxSigma[d] = -log(_r0) * (_n + 1) / (2 * sizePml[d] * gridStep[d]);
            else this->maxSigma[d] = 0;
        }
    }

    template<GridTypes gridTypes>
    inline Pml<gridTypes>::Pml(FieldSolver<gridTypes>* _fieldSolver, Int3 _sizePml,
        FP nPmlParam, FP r0PmlParam) :
        sizePml(_sizePml), fieldSolver(_fieldSolver)
    {
        Grid<FP, gridTypes>* grid = fieldSolver->grid;

        checkGridAndPmlSize(grid->globalGridDims);

        // TODO: consider that pml is only in edge domains
        leftPmlBorder = sizePml + grid->getNumExternalLeftCells();
        rightPmlBorder = leftPmlBorder + (grid->numInternalCells - 2 * sizePml);

        initializePmlParams(nPmlParam, r0PmlParam, grid->steps);
        fillPmlIndices();
    }

    template<GridTypes gridTypes>
    void Pml<gridTypes>::fillPmlIndices()
    {
        Grid<FP, gridTypes>* grid = this->fieldSolver->grid;

        Int3 begin = this->fieldSolver->updateEAreaBegin;
        Int3 end = this->fieldSolver->updateEAreaEnd;

        for (int i = begin.x; i < end.x; i++)
            for (int j = begin.y; j < end.y; j++)
                for (int k = begin.z; k < end.z; k++)
                {
                    bool xBoundaryPml = (i < this->leftPmlBorder.x) || (i >= this->rightPmlBorder.x);
                    bool yBoundaryPml = (j < this->leftPmlBorder.y) || (j >= this->rightPmlBorder.y);
                    bool zBoundaryPml = (k < this->leftPmlBorder.z) || (k >= this->rightPmlBorder.z);
                    if (xBoundaryPml || yBoundaryPml || zBoundaryPml)
                    {
                        eIndex.push_back(Int3(i, j, k));
                    }
                }

        begin = this->fieldSolver->updateBAreaBegin;
        end = this->fieldSolver->updateBAreaEnd;

        for (int i = begin.x; i < end.x; i++)
            for (int j = begin.y; j < end.y; j++)
                for (int k = begin.z; k < end.z; k++)
                {
                    bool xBoundaryPml = (i < this->leftPmlBorder.x) || (i >= this->rightPmlBorder.x);
                    bool yBoundaryPml = (j < this->leftPmlBorder.y) || (j >= this->rightPmlBorder.y);
                    bool zBoundaryPml = (k < this->leftPmlBorder.z) || (k >= this->rightPmlBorder.z);
                    if (xBoundaryPml || yBoundaryPml || zBoundaryPml)
                    {
                        bIndex.push_back(Int3(i, j, k));
                    }
                }

        this->initializeSplitFieldsE(eIndex.size());
        this->initializeSplitFieldsB(bIndex.size());
    }

    template<GridTypes gridTypes>
    void Pml<gridTypes>::initializeSplitFieldsE(int sizeE) {
        exy.resize(sizeE, 0);
        exz.resize(sizeE, 0);
        eyx.resize(sizeE, 0);
        eyz.resize(sizeE, 0);
        ezx.resize(sizeE, 0);
        ezy.resize(sizeE, 0);
    }

    template<GridTypes gridTypes>
    void Pml<gridTypes>::initializeSplitFieldsB(int sizeB) {
        bxy.resize(sizeB, 0);
        bxz.resize(sizeB, 0);
        byx.resize(sizeB, 0);
        byz.resize(sizeB, 0);
        bzx.resize(sizeB, 0);
        bzy.resize(sizeB, 0);
    }

    template<GridTypes gridTypes>
    FP Pml<gridTypes>::computeSigma(FP coord, CoordinateEnum axis) const
    {
        Grid<FP, gridTypes>* grid = fieldSolver->grid;
        int d = (int)axis;

        FP globalLeftBorderCoord = grid->origin[d];
        FP globalRightBorderCoord = grid->origin[d] + grid->steps[d] * grid->numCells[d];
        FP leftPmlBorderCoord = grid->origin[d] + grid->steps[d] * leftPmlBorder[d];
        FP rightPmlBorderCoord = grid->origin[d] + grid->steps[d] * rightPmlBorder[d];

        FP distance;
        if (coord < leftPmlBorderCoord)
            distance = (leftPmlBorderCoord - coord) / (leftPmlBorderCoord - globalLeftBorderCoord);
        else if (coord >= rightPmlBorderCoord)
            distance = (coord - rightPmlBorderCoord) / (globalRightBorderCoord - rightPmlBorderCoord);
        else distance = 0;

        return this->maxSigma[d] * pow(distance, (FP)n);
    }


    template<GridTypes gridTypes>
    class PmlReal : public Pml<gridTypes>
    {
    public:

        PmlReal(RealFieldSolver<gridTypes>* fieldSolver, Int3 sizePml,
            FP nPmlParam = (FP)4.0, FP r0PmlParam = (FP)1e-8) :
            Pml<gridTypes>(static_cast<FieldSolver<gridTypes>*>(fieldSolver),
                sizePml, nPmlParam, r0PmlParam)
        {}

    };


    template<GridTypes gridTypes>
    class PmlSpectral : public Pml<gridTypes>
    {
    public:

        PmlSpectral(SpectralFieldSolver<gridTypes>* fieldSolver, Int3 sizePml,
            FP nPmlParam = (FP)4.0, FP r0PmlParam = (FP)1e-8) :
            Pml<gridTypes>(static_cast<FieldSolver<gridTypes>*>(fieldSolver),
                sizePml, nPmlParam, r0PmlParam)
        {}

        // first step of spectral pml: updating of split components
        virtual void updateExSplit() {};
        virtual void updateEySplit() {};
        virtual void updateEzSplit() {};
        virtual void updateBxSplit() {};
        virtual void updateBySplit() {};
        virtual void updateBzSplit() {};

        virtual void updateBSplit();
        virtual void updateESplit();

        // multiplication by exp(-sigma*dt)
        void multBySigmaE();
        void multBySigmaB();

        // summation of split components
        void sumFieldEx();
        void sumFieldEy();
        void sumFieldEz();
        void sumFieldBx();
        void sumFieldBy();
        void sumFieldBz();

        void sumFieldE();
        void sumFieldB();

        // second step of spectral pml: multiplication by exp(-sigma*dt) and summation of split components
        virtual void updateB();
        virtual void updateE();

    protected:

        void sumFieldComponent(ScalarField<FP>& field, const std::vector<Int3>& index,
            std::vector<FP>& splitField1, std::vector<FP>& splitField2);
    };

    template<GridTypes gridTypes>
    inline void PmlSpectral<gridTypes>::updateBSplit()
    {
        updateBxSplit();
        updateBySplit();
        updateBzSplit();
    }

    template<GridTypes gridTypes>
    inline void PmlSpectral<gridTypes>::updateESplit()
    {
        updateExSplit();
        updateEySplit();
        updateEzSplit();
    }

    template<GridTypes gridTypes>
    inline void PmlSpectral<gridTypes>::multBySigmaE()
    {
        Grid<FP, gridTypes>* grid = this->fieldSolver->grid;
        const int size = this->eIndex.size();
        const FP ñdt = constants::c * this->fieldSolver->dt;

        OMP_FOR()
            for (int idx = 0; idx < size; ++idx)
            {
                int i = eIndex[idx].x, j = eIndex[idx].y, k = eIndex[idx].z;

                FP sigmaX = this->computeSigma(grid->EyPosition(i, j, k).x, CoordinateEnum::x);
                FP sigmaY = this->computeSigma(grid->EzPosition(i, j, k).y, CoordinateEnum::y);
                FP sigmaZ = this->computeSigma(grid->ExPosition(i, j, k).z, CoordinateEnum::z);
                FP coeffX = exp(-sigmaX * ñdt), coeffY = exp(-sigmaY * ñdt), coeffZ = exp(-sigmaZ * ñdt);

                this->eyx[idx] *= coeffY;
                this->ezx[idx] *= coeffZ;
                this->ezy[idx] *= coeffZ;
                this->exy[idx] *= coeffX;
                this->exz[idx] *= coeffX;
                this->eyz[idx] *= coeffY;
            }
    }

    template<GridTypes gridTypes>
    inline void PmlSpectral<gridTypes>::multBySigmaB()
    {
        Grid<FP, gridTypes>* grid = this->fieldSolver->grid;
        const int size = this->bIndex.size();
        const FP ñdt = constants::c * this->fieldSolver->dt;

        OMP_FOR()
            for (int idx = 0; idx < size; ++idx)
            {
                int i = bIndex[idx].x, j = bIndex[idx].y, k = bIndex[idx].z;

                FP sigmaX = this->computeSigma(grid->ByPosition(i, j, k).x, CoordinateEnum::x);
                FP sigmaY = this->computeSigma(grid->BzPosition(i, j, k).y, CoordinateEnum::y);
                FP sigmaZ = this->computeSigma(grid->BxPosition(i, j, k).z, CoordinateEnum::z);
                FP coeffX = exp(-sigmaX * ñdt), coeffY = exp(-sigmaY * ñdt), coeffZ = exp(-sigmaZ * ñdt);

                this->byx[idx] *= coeffY;
                this->bzx[idx] *= coeffZ;
                this->bzy[idx] *= coeffZ;
                this->bxy[idx] *= coeffX;
                this->bxz[idx] *= coeffX;
                this->byz[idx] *= coeffY;
            }
    }

    template<GridTypes gridTypes>
    inline void PmlSpectral<gridTypes>::sumFieldEx()
    {
        Grid<FP, gridTypes>* grid = this->fieldSolver->grid;
        sumFieldComponent(grid->Ex, this->eIndex, this->eyx, this->ezx);
    }

    template<GridTypes gridTypes>
    inline void PmlSpectral<gridTypes>::sumFieldEy()
    {
        Grid<FP, gridTypes>* grid = this->fieldSolver->grid;
        sumFieldComponent(grid->Ey, this->eIndex, this->ezy, this->exy);
    }

    template<GridTypes gridTypes>
    inline void PmlSpectral<gridTypes>::sumFieldEz()
    {
        Grid<FP, gridTypes>* grid = this->fieldSolver->grid;
        sumFieldComponent(grid->Ez, this->eIndex, this->exz, this->eyz);
    }

    template<GridTypes gridTypes>
    inline void PmlSpectral<gridTypes>::sumFieldE()
    {
        this->sumFieldEx();
        this->sumFieldEy();
        this->sumFieldEz();
    }

    template<GridTypes gridTypes>
    inline void PmlSpectral<gridTypes>::sumFieldBx()
    {
        Grid<FP, gridTypes>* grid = this->fieldSolver->grid;
        sumFieldComponent(grid->Bx, this->bIndex, this->byx, this->bzx);
    }

    template<GridTypes gridTypes>
    inline void PmlSpectral<gridTypes>::sumFieldBy()
    {
        Grid<FP, gridTypes>* grid = this->fieldSolver->grid;
        sumFieldComponent(grid->By, this->bIndex, this->bzy, this->bxy);
    }

    template<GridTypes gridTypes>
    inline void PmlSpectral<gridTypes>::sumFieldBz()
    {
        Grid<FP, gridTypes>* grid = this->fieldSolver->grid;
        sumFieldComponent(grid->Bz, this->bIndex, this->bxz, this->byz);
    }

    template<GridTypes gridTypes>
    inline void PmlSpectral<gridTypes>::sumFieldB()
    {
        this->sumFieldBx();
        this->sumFieldBy();
        this->sumFieldBz();
    }

    template<GridTypes gridTypes>
    inline void PmlSpectral<gridTypes>::updateB()
    {
        multBySigmaB();
        sumFieldB();
    }

    template<GridTypes gridTypes>
    inline void PmlSpectral<gridTypes>::updateE()
    {
        multBySigmaE();
        sumFieldE();
    }

    template<GridTypes gridTypes>
    inline void PmlSpectral<gridTypes>::sumFieldComponent(ScalarField<FP>& field,
        const std::vector<Int3>& index,
        std::vector<FP>& splitField1, std::vector<FP>& splitField2)
    {
        const int size = index.size();

        OMP_FOR()
            for (int idx = 0; idx < size; ++idx)
            {
                int i = index[idx].x, j = index[idx].y, k = index[idx].z;

                field(i, j, k) = splitField1[idx] + splitField2[idx];
            }
    }
}
