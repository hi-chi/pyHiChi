#include "TestingUtility.h"

#include "Constants.h"
#include "Fdtd.h"
#include "Particle.h"
#include "ParticleGenerator.h"
#include "Species.h"

using namespace pfc;

static FP testParticleGeneratorTestParticlePositionIsInsideGridParticleDensityFunc(FP x, FP y, FP z) {
    return x;
}

static FP testParticleGeneratorTestParticlePositionIsInsideGridInitialTemperatureFunc(FP x, FP y, FP z) {
    return 0.0;
}

static FP3 testParticleGeneratorTestParticlePositionIsInsideGridInitialMomentumFunc(FP x, FP y, FP z) {
    return FP3(0.0, 0.0, 0.0);
}

TEST(ParticleGeneratorTest, ParticlePositionIsInsideGrid) {
    // create grid
    FP Ex = 0.0, Ey = 0.0, Ez = 0.0;
    FP Bx = 0.0, By = 0.0, Bz = 0.0;
    FP stepX = 1.0, stepY = 1.0, stepZ = 1.0;
    Int3 numInternalCells = Int3(10, 20, 5);
    FP3 minCoords = FP3(2.0, 3.0, 1.0);
    FP3 steps = FP3(stepX, stepY, stepZ);
    Int3 globalGridDims = numInternalCells;
    YeeGrid grid(numInternalCells, minCoords, steps, globalGridDims);
    FP3 maxCoords = minCoords + grid.numInternalCells * grid.steps;
    FP T0 = 1.0;

    for (int i = 0; i < grid.numCells.x; ++i)
        for (int j = 0; j < grid.numCells.y; ++j)
            for (int k = 0; k < grid.numCells.z; ++k) {
                grid.Ex(i, j, k) = Ex;
                grid.Ey(i, j, k) = Ey;
                grid.Ez(i, j, k) = Ez;
                grid.Bx(i, j, k) = Bx;
                grid.By(i, j, k) = By;
                grid.Bz(i, j, k) = Bz;
            }
    
    ParticleArray3d particleArray;

    ParticleGenerator particleGenerator;
    particleGenerator(&particleArray, &grid,
        testParticleGeneratorTestParticlePositionIsInsideGridParticleDensityFunc,
        testParticleGeneratorTestParticlePositionIsInsideGridInitialTemperatureFunc,
        testParticleGeneratorTestParticlePositionIsInsideGridInitialMomentumFunc
        );

    for (int i = 0; i < particleArray.size(); ++i)
        ASSERT_TRUE((particleArray[i].getPosition().x >= minCoords.x) && (particleArray[i].getPosition().x <= maxCoords.x)
                 && (particleArray[i].getPosition().y >= minCoords.y) && (particleArray[i].getPosition().y <= maxCoords.y)
                 && (particleArray[i].getPosition().z >= minCoords.z) && (particleArray[i].getPosition().z <= maxCoords.z));
}

