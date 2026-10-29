#include "TestingUtility.h"

#include "ParticleBoundaries.h"
#include "Constants.h"
#include "Fdtd.h"
#include "Particle.h"

using namespace pfc;

//TEST(ParticleSolverTest, PeriodicalParticleBoundariesIsCorrect) {
//    // create grid
//    FP Ex = 0.0, Ey = 0.0, Ez = 0.0;
//    FP Bx = 0.0, By = 0.0, Bz = 0.0;
//    FP stepX = 1.0, stepY = 1.0, stepZ = 1.0;
//    Int3 numInternalCells = Int3(10, 20, 5);
//    FP3 minCoords = FP3(2.0, 3.0, 1.0);
//    FP3 steps = FP3(stepX, stepY, stepZ);
//    Int3 globalGridDims = numInternalCells;
//    YeeGrid grid(numInternalCells, minCoords, steps, globalGridDims);
//
//    for (int i = 0; i < grid.numCells.x; ++i)
//        for (int j = 0; j < grid.numCells.y; ++j)
//            for (int k = 0; k < grid.numCells.z; ++k) {
//                grid.Ex(i, j, k) = Ex;
//                grid.Ey(i, j, k) = Ey;
//                grid.Ez(i, j, k) = Ez;
//                grid.Bx(i, j, k) = Bx;
//                grid.By(i, j, k) = By;
//                grid.Bz(i, j, k) = Bz;
//            }
//
//    //particles
//    FP x1 = 1.0, x2 = 3.0, x3 = 13.0;
//    FP y1 = 2.0, y2 = 4.0, y3 = 24.0;
//    FP z1 = 0.0, z2 = 2.0, z3 = 7.0;
//
//    Particle3d::PositionType position(x, y, z);
//    Particle3d::MomentumType momentum(FP3(1.0, 1.0, 1.0));
//
//    ParticleArray3d particleArray;
//    ParticleArray3d particleArrayCorrect;
//
//    particleArray.pushBack(Particle3d((FP3(x1, y1, z1)),momentum));
//    particleArray.pushBack(Particle3d((FP3(x1, y1, z2)), momentum));
//    particleArray.pushBack(Particle3d((FP3(x1, y1, z3)), momentum));
//    particleArray.pushBack(Particle3d((FP3(x1, y2, z1)), momentum));
//    particleArray.pushBack(Particle3d((FP3(x1, y2, z2)), momentum));
//    particleArray.pushBack(Particle3d((FP3(x1, y2, z3)), momentum));
//    particleArray.pushBack(Particle3d((FP3(x1, y3, z1)), momentum));
//    particleArray.pushBack(Particle3d((FP3(x1, y3, z2)), momentum));
//    particleArray.pushBack(Particle3d((FP3(x1, y3, z3)), momentum));
//
//    particleArray.pushBack(Particle3d((FP3(x2, y1, z1)), momentum));
//    particleArray.pushBack(Particle3d((FP3(x2, y1, z2)), momentum));
//    particleArray.pushBack(Particle3d((FP3(x2, y1, z3)), momentum));
//    particleArray.pushBack(Particle3d((FP3(x2, y2, z1)), momentum));
//    particleArray.pushBack(Particle3d((FP3(x2, y2, z2)), momentum));
//    particleArray.pushBack(Particle3d((FP3(x2, y2, z3)), momentum));
//    particleArray.pushBack(Particle3d((FP3(x2, y3, z1)), momentum));
//    particleArray.pushBack(Particle3d((FP3(x2, y3, z2)), momentum));
//    particleArray.pushBack(Particle3d((FP3(x2, y3, z3)), momentum));
//
//    particleArray.pushBack(Particle3d((FP3(x3, y1, z1)), momentum));
//    particleArray.pushBack(Particle3d((FP3(x3, y1, z2)), momentum));
//    particleArray.pushBack(Particle3d((FP3(x3, y1, z3)), momentum));
//    particleArray.pushBack(Particle3d((FP3(x3, y2, z1)), momentum));
//    particleArray.pushBack(Particle3d((FP3(x3, y2, z2)), momentum));
//    particleArray.pushBack(Particle3d((FP3(x3, y2, z3)), momentum));
//    particleArray.pushBack(Particle3d((FP3(x3, y3, z1)), momentum));
//    particleArray.pushBack(Particle3d((FP3(x3, y3, z2)), momentum));
//    particleArray.pushBack(Particle3d((FP3(x3, y3, z3)), momentum));
//
//    // correct position values
//    particleArrayCorrect.pushBack(Particle3d((FP3(11.0, 22.0, 5.0)), momentum));
//    particleArrayCorrect.pushBack(Particle3d((FP3(11.0, 22.0, z2)), momentum));
//    particleArrayCorrect.pushBack(Particle3d((FP3(11.0, 22.0, 2.0)), momentum));
//    particleArrayCorrect.pushBack(Particle3d((FP3(11.0, y2, 5.0)), momentum));
//    particleArrayCorrect.pushBack(Particle3d((FP3(11.0, y2, z2)), momentum));
//    particleArrayCorrect.pushBack(Particle3d((FP3(11.0, y2, 2.0)), momentum));
//    particleArrayCorrect.pushBack(Particle3d((FP3(11.0, 4.0, 5.0)), momentum));
//    particleArrayCorrect.pushBack(Particle3d((FP3(11.0, 4.0, z2)), momentum));
//    particleArrayCorrect.pushBack(Particle3d((FP3(11.0, 4.0, 2.0)), momentum));
//
//    particleArrayCorrect.pushBack(Particle3d((FP3(x2, 22.0, 5.0)), momentum));
//    particleArrayCorrect.pushBack(Particle3d((FP3(x2, 22.0, z2)), momentum));
//    particleArrayCorrect.pushBack(Particle3d((FP3(x2, 22.0, 2.0)), momentum));
//    particleArrayCorrect.pushBack(Particle3d((FP3(x2, y2, 5.0)), momentum));
//    particleArrayCorrect.pushBack(Particle3d((FP3(x2, y2, z2)), momentum));
//    particleArrayCorrect.pushBack(Particle3d((FP3(x2, y2, 2.0)), momentum));
//    particleArrayCorrect.pushBack(Particle3d((FP3(x2, 4.0, 5.0)), momentum));
//    particleArrayCorrect.pushBack(Particle3d((FP3(x2, 4.0, z2)), momentum));
//    particleArrayCorrect.pushBack(Particle3d((FP3(x2, 4.0, 2.0)), momentum));
//
//    particleArrayCorrect.pushBack(Particle3d((FP3(3.0, 22.0, 5.0)), momentum));
//    particleArrayCorrect.pushBack(Particle3d((FP3(3.0, 22.0, z2)), momentum));
//    particleArrayCorrect.pushBack(Particle3d((FP3(3.0, 22.0, 2.0)), momentum));
//    particleArrayCorrect.pushBack(Particle3d((FP3(3.0, y2, 5.0)), momentum));
//    particleArrayCorrect.pushBack(Particle3d((FP3(3.0, y2, z2)), momentum));
//    particleArrayCorrect.pushBack(Particle3d((FP3(3.0, y2, 2.0)), momentum));
//    particleArrayCorrect.pushBack(Particle3d((FP3(3.0, 4.0, 5.0)), momentum));
//    particleArrayCorrect.pushBack(Particle3d((FP3(3.0, 4.0, z2)), momentum));
//    particleArrayCorrect.pushBack(Particle3d((FP3(3.0, 4.0, 2.0)), momentum));
//
//    PeriodicalParticleBoundaryConditions particleSolver;
//    particleSolver.updateParticlePosition(&grid, &particleArray);
//
//    for (int i = 0; i < particleArray.size(); ++i)
//        ASSERT_EQ(particleArrayCorrect[i].getPosition(), particleArray[i].getPosition());
//}