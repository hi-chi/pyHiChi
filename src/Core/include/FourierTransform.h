#pragma once
#include <omp.h>
#include "ScalarField.h"
#include "Grid.h"
#include "Enums.h"

#ifdef __USE_FFT__
#include "fftw3.h"
#endif

namespace pfc
{
    namespace fourier_transform {

        inline Int3 getSizeOfComplexArray(Int3 sizeOfFP)
        {
            return Int3(sizeOfFP.x, sizeOfFP.y, sizeOfFP.z / 2 + 1);
        }

        enum Direction {
            RtoC, CtoR
        };
    }


    class FourierTransformField {
#ifdef __USE_FFT__
        Int3 size;
        fftw_plan plans[2];  // RtoC/CtoR
        ScalarField<FP>* realField;
        ScalarField<complexFP>* complexField;
#endif

    public:

#ifdef __USE_FFT__
        FourierTransformField()
        {
            plans[fourier_transform::Direction::RtoC] = 0;
            plans[fourier_transform::Direction::CtoR] = 0;
        }

        void initialize(ScalarField<FP>* _realField, ScalarField<complexFP>* _complexField, Int3 _size)
        {
            size = _size;
            realField = _realField;
            complexField = _complexField;
            createPlans();
        }

        ~FourierTransformField() {
            destroyPlans();
        }

        void doDirectFourierTransform()
        {
            fftw_execute(plans[fourier_transform::Direction::RtoC]);
        }

        void doInverseFourierTransform()
        {
            fftw_execute(plans[fourier_transform::Direction::CtoR]);
            ScalarField<FP>& res = *realField;
#pragma omp parallel for
            for (int i = 0; i < size.x; i++)
                for (int j = 0; j < size.y; j++)
                    //#pragma omp simd
                    for (int k = 0; k < size.z; k++)
                        res(i, j, k) /= (FP)size.x*size.y*size.z;
        }

#else
        FourierTransformField() {}

        void initialize(ScalarField<FP>* _realField, ScalarField<complexFP>* _complexField, Int3 _size) {}

        void doDirectFourierTransform() {}
        void doInverseFourierTransform() {}

        ~FourierTransformField() {}
#endif

        FourierTransformField(ScalarField<FP>* _realField, ScalarField<complexFP>* _complexField, Int3 _size) {
            initialize(_realField, _complexField, _size);
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

            ScalarField<FP>& arrD = *(realField);
            ScalarField<complexFP>& arrC = *(complexField);

#ifdef __USE_OMP__
            fftw_plan_with_nthreads(omp_get_max_threads());
#endif
            plans[fourier_transform::Direction::RtoC] = fftw_plan_dft_r2c_3d(Nx, Ny, Nz,
                &(arrD(0, 0, 0)), (fftw_complex*)&(arrC(0, 0, 0)), FFTW_ESTIMATE);
#ifdef __USE_OMP__
            fftw_plan_with_nthreads(omp_get_max_threads());
#endif
            plans[fourier_transform::Direction::CtoR] = fftw_plan_dft_c2r_3d(Nx, Ny, Nz,
                (fftw_complex*)&(arrC(0, 0, 0)), &(arrD(0, 0, 0)), FFTW_ESTIMATE);
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


    class FourierTransformGrid {

        FourierTransformField transform[3][3];  // field, coordinate

    public:

        FourierTransformGrid() {}
        
        template<GridTypes gridType>
        void initialize(Grid<FP, gridType>* gridFP, Grid<complexFP, gridType>* gridCFP) {
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