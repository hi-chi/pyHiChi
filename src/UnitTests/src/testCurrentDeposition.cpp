#include "TestingUtility.h"

#include "CurrentDeposition.h"
#include "CurrentBoundaries.h"
#include "Constants.h"
#include "Fdtd.h"
#include "FieldGenerator.h"
#include "FieldValue.h"
#include "Particle.h"
#include "ParticleGenerator.h"
#include "ParticleBoundaries.h"
#include "ParticleSorting.h"
#include "Pusher.h"
#include <fstream>
#include <iostream>
#include <chrono>
#include <algorithm>

using namespace pfc;

TEST(FormfactorTest, TSCSumIsUnit) {

    FP stepx = 0.1, stepy = 0.2, stepz = 0.1;    
    Int3 numInternalCells = Int3(10, 7, 12);
    FP3 minCoords = FP3(0.0, 0.0, 0.0);
    FP3 steps = FP3(stepx, stepy, stepz);
    Int3 globalGridDims = numInternalCells;

    Int3 idxJx, idxJy, idxJz;
    FP3 internalCoordsX, internalCoordsY, internalCoordsZ;
    
    YeeGrid grid(numInternalCells, minCoords, steps, globalGridDims);
    FP x = (FP)0.3 * stepx;
    FP y = (FP)0.1 * stepy;
    FP z = (FP)0.4 * stepz;

    Particle3d::PositionType position(x, y, z);
    idxJx = grid.getClosestIndexJx(position);
    idxJy = grid.getClosestIndexJx(position);
    idxJz = grid.getClosestIndexJx(position);

    internalCoordsX = grid.getClosestInternalCoordsJx(position);
    internalCoordsY = grid.getClosestInternalCoordsJx(position);
    internalCoordsZ = grid.getClosestInternalCoordsJx(position);

    Int3 middleIdxX = truncate(internalCoordsX + FP3(0.5, 0.5, 0.5));
    Int3 middleIdxY = truncate(internalCoordsY + FP3(0.5, 0.5, 0.5));
    Int3 middleIdxZ = truncate(internalCoordsZ + FP3(0.5, 0.5, 0.5));
    FormFactorTSC formFactor;
    double sumX = 0, sumY = 0, sumZ = 0;

    formFactor(internalCoordsX - FP3(middleIdxX));
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                sumX += formFactor.c[0][i] * formFactor.c[1][j] * formFactor.c[2][k];
    }

    formFactor(internalCoordsY - FP3(middleIdxY));
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                sumY += formFactor.c[0][i] * formFactor.c[1][j] * formFactor.c[2][k];
    }

    formFactor(internalCoordsZ - FP3(middleIdxZ));
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                sumZ += formFactor.c[0][i] * formFactor.c[1][j] * formFactor.c[2][k];
    }
    FP maxError = (FP)1e-14;
    EXPECT_NEAR(sumX, 1.0, maxError);
    EXPECT_NEAR(sumY, 1.0, maxError);
    EXPECT_NEAR(sumZ, 1.0, maxError);
}

static FP testCurrentDepositionTestPlasmaOscillationTestParticleDensityFunc(FP x, FP y, FP z) {
    FP T0 = 0.01 * constants::electronMass * constants::c * constants::c;
    FP L = 1.0;
    FP dx = L / 64.0; //remember about nx=64
    FP D = T0 / (8 * constants::pi * constants::electronCharge * constants::electronCharge * 0.25 * dx * dx);
    return D * (1 + 0.05 * std::sin(2 * constants::pi * x / L));
}

static FP testCurrentDepositionTestPlasmaOscillationTestInitialTemperatureFunc(FP x, FP y, FP z) {
    FP T0 = 0.01 * constants::electronMass * constants::c * constants::c;
    return T0;
}

static FP3 testCurrentDepositionTestPlasmaOscillationTestInitialMomentumFunc(FP x, FP y, FP z) {
    return FP3(0.0, 0.0, 0.0);
}

TEST(CurrentDepositionTest, PlasmaOscillationTest) {
    int nx = 64, ny = 8, nz = 8;
    FP L = 1.0;
    FP dx = L / nx, dy = dx, dz = dx;
    int Nip = 256;
    int Np = nx / (std::sqrt(2) * constants::pi);  // number of periods
    int N = Nip * Np;
    FP T0 = 0.01 * constants::electronMass * constants::c * constants::c;
    FP L_Debay = dx * 0.5;
    FP D = T0 / (8 * constants::pi * constants::electronCharge * L_Debay * constants::electronCharge * L_Debay);
    FP wp = sqrt(4.0 * constants::pi * constants::electronCharge * constants::electronCharge * D
        / constants::electronMass);
    FP dt = 2 * (constants::pi / wp) / Nip;
    FP Nc = 30;
    Particle3d::WeightType w = D * dx * dx * dx / Nc;
    FP A = 0.05;
    FP3 p0 = FP3(0.0, 0.0, 0.0);

    BorisPusher scalarPusher;

    CurrentDepositionCIC<YeeGrid> currentDeposition(dt);

    Int3 numInternalCells = Int3(nx, ny, nz);
    FP3 minCoords = FP3(0.0, 0.0, 0.0);
    FP3 steps = FP3(dx, dy, dz);
    Int3 globalGridDims = numInternalCells;

    YeeGrid grid(numInternalCells, minCoords, steps, globalGridDims);

    for (int i = 0; i < grid.numCells.x; ++i)
        for (int j = 0; j < grid.numCells.y; ++j)
            for (int k = 0; k < grid.numCells.z; ++k) {
                grid.Ex(i, j, k) = -2 * L * D * A * constants::electronCharge 
                    * std::cos(2 * constants::pi * grid.getBaseCoords(i, j, k).x / L);
                grid.Ey(i, j, k) = 0.0;
                grid.Ez(i, j, k) = 0.0;
                grid.Bx(i, j, k) = 0.0;
                grid.By(i, j, k) = 0.0;
                grid.Bz(i, j, k) = 0.0;
            }

    ParticleArray3d particleArray;

    PeriodicalParticleBoundaryConditions particleSolver;
    FDTD fdtd(&grid, dt);
    fdtd.setPeriodicalBoundaryConditions();
    PeriodicalCurrentBoundaryConditionFdtd currentBoundary(&grid);

    ParticleGenerator particleGenerator;
    particleGenerator(&particleArray, &grid,
        testCurrentDepositionTestPlasmaOscillationTestParticleDensityFunc,
        testCurrentDepositionTestPlasmaOscillationTestInitialTemperatureFunc,
        testCurrentDepositionTestPlasmaOscillationTestInitialMomentumFunc,
        w, ParticleTypes::Electron
    );

    ParticleSorting<YeeGrid, ParticleArray3d> particleSorting(&grid, particleArray);

    std::ofstream fout("OscillationTestEx.txt");
    std::ofstream fout2("OscillationTestElectronDensity.txt");
    std::ofstream fout_energy("FieldEnergy.txt");
    

    double sort_time = 0.0;
    double current_time = 0.0;
    double fdtd_time = 0.0;
    double interpolate_time = 0.0;
    double pusher_time = 0.0;

    //-------------------------------------------initialize-field-energy-------------------------------------------------
    double energy = 0.0, energy_particles = 0.0;

    for (int i = 0; i < N; ++i) {
        //---------------------------------start-writing-Ex-density-energy-to-files--------------------------------------
         if (i % 16 == 0) {
            std::cout << "iter: " << i << "/"<< N << std::endl;
            //-----------------------------computing-current-energy------------------------------------------------------
            energy = 0.0; energy_particles = 0.0;
            Int3 numExternalLeftCells = grid.getNumExternalLeftCells();
            Int3 size = grid.numInternalCells + numExternalLeftCells;
            for (int i = numExternalLeftCells.x; i < size.x; i++)
                for (int j = numExternalLeftCells.y; j < size.y; j++)
                    for (int k = numExternalLeftCells.z; k < size.z; k++)
                        energy += (grid.Ex(i, j, k) * grid.Ex(i, j, k)
                            + grid.Ey(i, j, k) * grid.Ey(i, j, k)
                            + grid.Ez(i, j, k) * grid.Ez(i, j, k)
                            + grid.Bx(i, j, k) * grid.Bx(i, j, k)
                            + grid.By(i, j, k) * grid.By(i, j, k)
                            + grid.Bz(i, j, k) * grid.Bz(i, j, k));
            energy = energy * grid.steps.volume() * 1e-7 / ((FP)8 * constants::pi);
            fout_energy << (i) << " " << energy << std::endl;
            //-----------------------------Ex-depending-on-x------------------------------------------------------------
            for (int j = 0; j < grid.numInternalCells.x; ++j) {
                Int3 idx; FP3 internalCoords;
                fout << (i) << " " << minCoords.x + j * grid.steps.x << " " << grid.Ex(j
                    + grid.getNumExternalLeftCells().x, grid.getNumExternalLeftCells().y,
                    grid.getNumExternalLeftCells().z) << std::endl;
            }
            //----------------------------electron-density-depending-on-x-----------------------------------------------
            for (int k = 0; k < grid.numInternalCells.x; ++k) {
                int electronCount = 0;
                FP minCellCoords = minCoords.x + k * grid.steps.x;
                FP maxCellCoords = minCoords.x + (k + 1) * grid.steps.x;

                for (int j = 0; j < particleArray.size(); ++j) {
                    if ((particleArray[j].getPosition().x >= minCellCoords)
                        && (particleArray[j].getPosition().x <= maxCellCoords))
                        electronCount++;
                }
                // to sinchronize the output with Picador
                FP electronDensity = electronCount * w;  // density through plane
                fout2 << (i) << " " << minCellCoords << " " << electronDensity << std::endl;
            }
         }
        //---------------------------------end-writing-Ex-density-energy-to-files--------------------------------------

        auto interpolate_start = chrono::steady_clock::now();
        // interpolate
        std::vector<ValueField> fields;
        for (int i = 0; i < particleArray.size(); i++) {
            FP3 E, B;
            grid.getFieldsCIC(particleArray[i].getPosition(), E, B);
            fields.push_back(ValueField(E, B));
        }

        auto interpolate_end = chrono::steady_clock::now();
        double interpolate_local_time = (chrono::duration_cast<chrono::nanoseconds>(interpolate_end 
            - interpolate_start).count()) / 1e9;
        interpolate_time += interpolate_local_time;

        auto pusher_start = chrono::steady_clock::now();
        // pusher
        scalarPusher(&particleArray, fields, dt);
        auto pusher_end = chrono::steady_clock::now();
        double pusher_local_time = (chrono::duration_cast<chrono::nanoseconds>(pusher_end 
            - pusher_start).count()) / 1e9;
        pusher_time += pusher_local_time;

        // periodical particle position
        particleSolver.updateParticlePosition(&grid, &particleArray);
        // current deposition
        fdtd.updateDomainBorders();        
        auto sort_start = chrono::steady_clock::now();
        particleSorting.sortL(particleArray);

        auto sort_end = chrono::steady_clock::now();
        double sort_local_time = (chrono::duration_cast<chrono::nanoseconds>(sort_end
         - sort_start).count()) / 1e9;
        sort_time += sort_local_time;

        auto current_start = chrono::steady_clock::now();
        currentDeposition(&grid, &particleArray);
        auto current_end = chrono::steady_clock::now();
        double current_local_time = (chrono::duration_cast<chrono::nanoseconds>(current_end
         - current_start).count()) / 1e9;
        current_time += current_local_time;
        currentBoundary.updateCurrentBoundaries();

        auto fdtd_start = chrono::steady_clock::now();

        // fdtd
        fdtd.updateFields();

        auto fdtd_end = chrono::steady_clock::now();
        double fdtd_local_time = (chrono::duration_cast<chrono::nanoseconds>(fdtd_end - fdtd_start).count()) / 1e9;
        fdtd_time += fdtd_local_time;
    }

    std::cout << std::endl;
    std::cout << "AFTER: " << std::endl;
    for (int i = 0; i < particleArray.size(); ++i) {
        FP3 position = particleArray[i].getPosition();
        Int3 idx = grid.getBaseIndex(position);
        std::cout << idx << ", " << i << ", " << position << std::endl;                
    }
    fout.close();
    fout2.close();

    //--------------------------------computing-final-energy---------------------------------------------------------
    energy = 0.0;
    for (int i = grid.numInternalCells.x; i < grid.numCells.x; i++)
        for (int j = grid.numInternalCells.y; j < grid.numCells.y; j++)
            for (int k = grid.numInternalCells.z; k < grid.numCells.z; k++) {
                energy += grid.Ex(i, j, k) * grid.Ex(i, j, k);
                energy += grid.Ey(i, j, k) * grid.Ey(i, j, k);
                energy += grid.Ez(i, j, k) * grid.Ez(i, j, k);
                energy += grid.Bx(i, j, k) * grid.Bx(i, j, k);
                energy += grid.By(i, j, k) * grid.By(i, j, k);
                energy += grid.Bz(i, j, k) * grid.Bz(i, j, k);
            }
    energy = energy * dx * dy * dz * 1e-7 / ((FP)8 * constants::pi);
    std::cout << "final energy: " << energy << std::endl;
    fout_energy.close();

    std::cout << "interpolate_time= " << interpolate_time << std::endl;
    std::cout << "pusher_time= " << pusher_time << std::endl;
    std::cout << "sort time= " << sort_time << std::endl;
    std::cout << "current_time= " << current_time << std::endl;
    std::cout << "fdtd time= " << fdtd_time << std::endl;
    ASSERT_NO_THROW(true);
}
