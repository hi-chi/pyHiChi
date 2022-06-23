#pragma once
#include <omp.h>
#include "ScalarField.h"
#include "Grid.h"
#include "SpectralGrid.h"
#include "Enums.h"

#include "macros.h"

#ifdef __USE_FFT__
#include "fftw3.h"
#endif

namespace pfc
{
    namespace fourier_transform {

        inline Int3 getSizeOfComplexArray(Int3 sizeFP)
        {
            return Int3(sizeFP.x, sizeFP.y, sizeFP.z / 2 + 1);
        }

        enum Direction {
            RtoC, CtoR
        };
    }


    class ArrayFourierTransform3d {
#ifdef __USE_FFT__
        Int3 size, memSizeRealData;
        fftw_plan plans[2];  // RtoC/CtoR
        FP* realData;
        complexFP* complexData;
#endif

    public:

#ifdef __USE_FFT__
        ArrayFourierTransform3d()
        {
            plans[fourier_transform::Direction::RtoC] = 0;
            plans[fourier_transform::Direction::CtoR] = 0;
        }

        void initialize(FP* _realData, complexFP* _complexData,
            Int3 _size, Int3 _memSizeRealData)
        {
            size = _size;
            memSizeRealData = _memSizeRealData;
            realData = _realData;
            complexData = _complexData;
            createPlans();
        }

        ~ArrayFourierTransform3d() {
            destroyPlans();
        }

        void doDirectFourierTransform()
        {
            fftw_execute(plans[fourier_transform::Direction::RtoC]);
        }

        void doInverseFourierTransform()
        {
            fftw_execute(plans[fourier_transform::Direction::CtoR]);
            FP normCoeff = (FP)size.volume();
            int nx = memSizeRealData.x, ny = memSizeRealData.y, nz = memSizeRealData.z;
            OMP_FOR_COLLAPSE()
            for (int i = 0; i < size.x; i++)
                for (int j = 0; j < size.y; j++)
                    for (int k = 0; k < size.z; k++)
                        realData[k + (j + i * ny) * nz] /= normCoeff;
        }

#else
        ArrayFourierTransform3d() {}

        void initialize(FP* _realData, complexFP* _complexData,
            Int3 _size, Int3 _memSizeRealData) {}

        void doDirectFourierTransform() {}
        void doInverseFourierTransform() {}

        ~ArrayFourierTransform3d() {}
#endif

        void initialize(FP* _realData, complexFP* _complexData, Int3 _size) {
            initialize(_realData, _complexData, _size, _size);
        }

        ArrayFourierTransform3d(FP* _realData, complexFP* _complexData,
            Int3 _size, Int3 _memSizeRealData) {
            initialize(_realData, _complexData, _size, _memSizeRealData);
        }

        ArrayFourierTransform3d(FP* _realData, complexFP* _complexData, Int3 _size) {
            initialize(_realData, _complexData, _size);
        }

        void doFourierTransform(fourier_transform::Direction direction)
        {
            switch (direction) {
            case fourier_transform::Direction::RtoC:
                doDirectFourierTransform();
                break;
            case fourier_transform::Direction::CtoR:
                doInverseFourierTransform();
                break;
            default:
                break;
            }
        }

    private:

#ifdef __USE_FFT__
        void createPlans()
        {
            int Nx = size.x, Ny = size.y, Nz = size.z;

#ifdef __USE_OMP__
            fftw_plan_with_nthreads(omp_get_max_threads());
#endif
            plans[fourier_transform::Direction::RtoC] = fftw_plan_dft_r2c_3d(Nx, Ny, Nz,
                &(realData[0]), (fftw_complex*)&(complexData[0]), FFTW_ESTIMATE);
#ifdef __USE_OMP__
            fftw_plan_with_nthreads(omp_get_max_threads());
#endif
            plans[fourier_transform::Direction::CtoR] = fftw_plan_dft_c2r_3d(Nx, Ny, Nz,
                (fftw_complex*)&(complexData[0]), &(realData[0]), FFTW_ESTIMATE);
        }

        void destroyPlans()
        {
            if (plans[fourier_transform::Direction::RtoC] != 0)
                fftw_destroy_plan(plans[fourier_transform::Direction::RtoC]);
            if (plans[fourier_transform::Direction::CtoR] != 0)
                fftw_destroy_plan(plans[fourier_transform::Direction::CtoR]);
        }
#endif
    };


    class FourierTransformField : public ArrayFourierTransform3d {
    public:

        FourierTransformField() : ArrayFourierTransform3d() {}
        FourierTransformField(ScalarField<FP>* _realData,
            SpectralScalarField<FP, complexFP>* _complexData, Int3 _size) :
            ArrayFourierTransform3d(_realData->getData(), _complexData->getData(),
                _size, _realData->getSize()) {}

        void initialize(ScalarField<FP>* _realData,
            SpectralScalarField<FP, complexFP>* _complexData, Int3 _size) {
            ArrayFourierTransform3d::initialize(_realData->getData(),
                _complexData->getData(), _size, _realData->getSize());
        }
    };


    class FourierTransformGrid {

        FourierTransformField transform[3][3];  // field, coordinate

    public:

        FourierTransformGrid() {}
        
        template<GridTypes gridType>
        void initialize(Grid<FP, gridType>* gridFP, SpectralGrid<FP, complexFP>* gridCFP) {
            transform[(int)FieldEnum::E][(int)CoordinateEnum::x].initialize(&gridFP->Ex, &gridCFP->Ex, gridFP->numCells);
            transform[(int)FieldEnum::E][(int)CoordinateEnum::y].initialize(&gridFP->Ey, &gridCFP->Ey, gridFP->numCells);
            transform[(int)FieldEnum::E][(int)CoordinateEnum::z].initialize(&gridFP->Ez, &gridCFP->Ez, gridFP->numCells);

            transform[(int)FieldEnum::B][(int)CoordinateEnum::x].initialize(&gridFP->Bx, &gridCFP->Bx, gridFP->numCells);
            transform[(int)FieldEnum::B][(int)CoordinateEnum::y].initialize(&gridFP->By, &gridCFP->By, gridFP->numCells);
            transform[(int)FieldEnum::B][(int)CoordinateEnum::z].initialize(&gridFP->Bz, &gridCFP->Bz, gridFP->numCells);

            transform[(int)FieldEnum::J][(int)CoordinateEnum::x].initialize(&gridFP->Jx, &gridCFP->Jx, gridFP->numCells);
            transform[(int)FieldEnum::J][(int)CoordinateEnum::y].initialize(&gridFP->Jy, &gridCFP->Jy, gridFP->numCells);
            transform[(int)FieldEnum::J][(int)CoordinateEnum::z].initialize(&gridFP->Jz, &gridCFP->Jz, gridFP->numCells);
        }

        void doDirectFourierTransform(FieldEnum field, CoordinateEnum coord) {
            transform[(int)field][(int)coord].doDirectFourierTransform();
        }

        void doInverseFourierTransform(FieldEnum field, CoordinateEnum coord) {
            transform[(int)field][(int)coord].doInverseFourierTransform();
        }

        void doFourierTransform(FieldEnum field, CoordinateEnum coord,
            fourier_transform::Direction direction)
        {
            transform[(int)field][(int)coord].doFourierTransform(direction);
        }

    };

}