#include "pyGrid.h"

#include "Constants.h"
#include "Psatd.h"
#include "Vectors.h"
#include "Enums.h"
#include "Mapping.h"
#include "FieldConfiguration.h"

#include <chrono>

using namespace pfc;


const FP FACTOR = 2.0;
const int NY = (int)(FACTOR * 256);
const int NZ = (int)(FACTOR * 256);
const int NX_BAND = int(56 * FACTOR);
const FP3 GRID_SIZE(NX_BAND, NY, NZ);

const FP WAVELENGTH = 1.0;
const FP PULSELENGTH = 2.0 * WAVELENGTH;
const FP PHASE = 0.0;
const FP R0 = 16 * WAVELENGTH;
const FP F_NUMBER = 0.3;
const FP EDGE_SMOOTHING_ANGLE = 0.1;
const FP TOTAL_POWER = constants::c;

const FP TIME_STEP = 8.0 / WAVELENGTH * constants::c;

const FP3 MIN_COORDS(-20.0, -20.0, -20.0);
const FP3 MAX_COORDS(20.0, 20.0, 20.0);

const FP D = 3.5 * PULSELENGTH;


int main(int argc, char **argv)
{

    TightFocusingField startConditions(F_NUMBER, R0, WAVELENGTH,
        PULSELENGTH, PHASE, TOTAL_POWER, EDGE_SMOOTHING_ANGLE);

    std::chrono::steady_clock::time_point startTimeInit = std::chrono::steady_clock::now();

    TightFocusingMapping mapping(R0, PULSELENGTH, D);

    FP xMin = mapping.getxMin();
    FP xMax = mapping.getxMax();

    FP3 gridMinCoords = FP3(xMin, MIN_COORDS.y, MIN_COORDS.z);
    FP3 gridMaxCoords = FP3(xMax, MAX_COORDS.y, MAX_COORDS.z);
    FP3 gridStep = (gridMaxCoords - gridMinCoords) / GRID_SIZE;

    pyPSATDGridMapping grid(GRID_SIZE, TIME_STEP, gridMinCoords, gridStep);
    grid.setMapping(&mapping);

    grid.setFieldConfiguration<TightFocusingField>(&startConditions);

    std::chrono::steady_clock::time_point endTimeInit = std::chrono::steady_clock::now();
    std::chrono::milliseconds timeInit =
        std::chrono::duration_cast<std::chrono::milliseconds>(endTimeInit - startTimeInit);

    std::cout << "Initialisation time is " << timeInit.count() << " ms" << std::endl;

    return 0;
}