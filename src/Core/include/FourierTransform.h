#pragma once
#include <omp.h>
#include "ScalarField.h"

#ifdef __USE_FFT__
#include "fftw3.h"
#endif

namespace pfc
{

    enum FourierTransformDirection {
        RtoC = -1,  // FFTW_FORWARD, direct Fourier transform
        CtoR = +1  // FFTW_BACKWARD, inverse Fourier transform
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

            for (int i = 0; i < result.getSize().x; i++)
                for (int j = 0; j < result.getSize().y; j++)
                    for (int k = 0; k < result.getSize().z; k++)
                        result(i, j, k) /= (FP) result.getSize().x*result.getSize().y*result.getSize().z;
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

    private:

        class FFTWrapper {

        };

    };

}