#pragma once
#include <omp.h>
#include "ScalarField.h"
#include "Grid.h"

#ifdef __USE_FFT__
#include "fftw3.h"
#endif

namespace pfc
{
    enum Field {
        E, B, J
    };

    enum Coordinate {
        x, y, z
    };

    enum FourierTransformDirection {
        RtoC, CtoR
    };


    class FourierTransformGrid {

#ifdef __USE_FFT__
        Int3 size;
        fftw_plan plans[2][3][3];  // RtoC/CtoR, field, coordinate
        ScalarField<FP>* realFields[3][3];
        ScalarField<complexFP>* complexFields[3][3];
#endif

    public:

#ifdef __USE_FFT__
        FourierTransformGrid()
        {
            for (int f = 0; f < 3; f++)
                for (int d = 0; d < 3; d++) {
                    plans[FourierTransformDirection::RtoC][f][d] = 0;
                    plans[FourierTransformDirection::CtoR][f][d] = 0;
                }
        }

        template<GridTypes gridType>
        void initialize(Grid<FP, gridType>* gridFP, Grid<complexFP, gridType>* gridCFP)
        {
            size = gridFP->numCells;

            realFields[E][x] = &(gridFP->Ex), realFields[E][y] = &(gridFP->Ey), realFields[E][z] = &(gridFP->Ez);
            realFields[B][x] = &(gridFP->Bx), realFields[B][y] = &(gridFP->By), realFields[B][z] = &(gridFP->Bz);
            realFields[J][x] = &(gridFP->Jx), realFields[J][y] = &(gridFP->Jy), realFields[J][z] = &(gridFP->Jz);

            complexFields[E][x] = &(gridCFP->Ex), complexFields[E][y] = &(gridCFP->Ey), complexFields[E][z] = &(gridCFP->Ez);
            complexFields[B][x] = &(gridCFP->Bx), complexFields[B][y] = &(gridCFP->By), complexFields[B][z] = &(gridCFP->Bz);
            complexFields[J][x] = &(gridCFP->Jx), complexFields[J][y] = &(gridCFP->Jy), complexFields[J][z] = &(gridCFP->Jz);

            createPlans();
        }

        ~FourierTransformGrid() {
            destroyPlans();
        }

        void doDirectFourierTransform(Field field, Coordinate coord)
        {
            fftw_execute(plans[FourierTransformDirection::RtoC][field][coord]);
        }

        void doInverseFourierTransform(Field field, Coordinate coord)
        {
            fftw_execute(plans[FourierTransformDirection::CtoR][field][coord]);
            ScalarField<FP>& res = *realFields[field][coord];
#pragma omp parallel for
            for (int i = 0; i < size.x; i++)
                for (int j = 0; j < size.y; j++)
//#pragma omp simd
                    for (int k = 0; k < size.z; k++)
                        res(i, j, k) /= (FP)size.x*size.y*size.z;
        }

#else
        FourierTransformGrid() {}
        template<GridTypes gridType>
        FourierTransformGrid(Grid<FP, gridType>* gridFP, Grid<complexFP, gridType>* gridCFP) {}
        
        template<GridTypes gridType>
        void initialize(Grid<FP, gridType>* gridFP, Grid<complexFP, gridType>* gridCFP) {}

        void doDirectFourierTransform(Field field, Coordinate coord) {}
        void doInverseFourierTransform(Field field, Coordinate coord) {}
#endif

        void doFourierTransform(Field field, Coordinate coord,
            FourierTransformDirection direction)
        {
            switch (direction) {
            case RtoC:
                doDirectFourierTransform(field, coord);
                break;
            case CtoR:
                doInverseFourierTransform(field, coord);
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

            for (int f = 0; f < 3; f++)
                for (int d = 0; d < 3; d++) {
                    ScalarField<FP>& arrD = *(realFields[f][d]);
                    ScalarField<complexFP>& arrC = *(complexFields[f][d]);

#ifdef __USE_OMP__
                    fftw_plan_with_nthreads(omp_get_max_threads());
#endif
                    plans[FourierTransformDirection::RtoC][f][d] = fftw_plan_dft_r2c_3d(Nx, Ny, Nz,
                        &(arrD(0, 0, 0)), (fftw_complex*)&(arrC(0, 0, 0)), FFTW_ESTIMATE);
#ifdef __USE_OMP__
                    fftw_plan_with_nthreads(omp_get_max_threads());
#endif
                    plans[FourierTransformDirection::CtoR][f][d] = fftw_plan_dft_c2r_3d(Nx, Ny, Nz,
                        (fftw_complex*)&(arrC(0, 0, 0)), &(arrD(0, 0, 0)), FFTW_ESTIMATE);
                }
        }

        void destroyPlans()
        {
            for (int f = 0; f < 3; f++)
                for (int d = 0; d < 3; d++) {
                    if (plans[FourierTransformDirection::RtoC][f][d] != 0)
                        fftw_destroy_plan(plans[FourierTransformDirection::RtoC][f][d]);
                    if (plans[FourierTransformDirection::CtoR][f][d] != 0)
                        fftw_destroy_plan(plans[FourierTransformDirection::CtoR][f][d]);
                }
        }
#endif

    };


    class FourierTransform {

    public:

#ifdef __USE_FFT__
        
        static void doDirectFourierTransform(ScalarField<FP>& data, ScalarField<complexFP>& result)
        {
            fftw_plan plan = 0;
#ifdef __USE_OMP__
            fftw_plan_with_nthreads(omp_get_max_threads());
#endif
            plan = fftw_plan_dft_r2c_3d(data.getSize().x, data.getSize().y, data.getSize().z,
                (FP*)&(data(0, 0, 0)), (fftw_complex*)&(result(0, 0, 0)), FFTW_ESTIMATE);
            fftw_execute(plan);
            fftw_destroy_plan(plan);
        }

        static void doInverseFourierTransform(ScalarField<complexFP>& data, ScalarField<FP>& result)
        {
            fftw_plan plan = 0;
#ifdef __USE_OMP__
            fftw_plan_with_nthreads(omp_get_max_threads());
#endif
            plan = fftw_plan_dft_c2r_3d(result.getSize().x, result.getSize().y, result.getSize().z,
                (fftw_complex*)&(data(0, 0, 0)), (FP*)&(result(0, 0, 0)), FFTW_ESTIMATE);
            fftw_execute(plan);
            fftw_destroy_plan(plan);

#pragma omp parallel for
            for (int i = 0; i < result.getSize().x; i++)
                for (int j = 0; j < result.getSize().y; j++)
                    //#pragma omp simd
                    for (int k = 0; k < result.getSize().z; k++)
                        result(i, j, k) /= (FP)result.getSize().x*result.getSize().y*result.getSize().z;
        }
#else
        static void doDirectFourierTransform(ScalarField<FP>& data, ScalarField<complexFP>& result) {}
        static void doInverseFourierTransform(ScalarField<complexFP>& data, ScalarField<FP>& result) {}
#endif

        static void doFourierTransform(ScalarField<FP>& field1, ScalarField<complexFP>& field2,
            FourierTransformDirection direction)
        {
            switch (direction) {
            case RtoC:
                doDirectFourierTransform(field1, field2);
                break;
            case CtoR:
                doInverseFourierTransform(field2, field1);
                break;
            default:
                break;
            }
        }

        static Int3 getSizeOfComplex(const Int3& sizeOfFP)
        {
            return Int3(sizeOfFP.x, sizeOfFP.y, sizeOfFP.z / 2 + 1);
        }

    };
}