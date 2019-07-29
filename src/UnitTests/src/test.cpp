#include "TestingUtility.h"
#include "Pml.h"
#include "Pstd.h"
#include "Fdtd.h"
#include <memory>

TEST(TEST_, test) {
    Int3 gridSize(32, 1, 1);
    Int3 minCoords = FP3(-1.0, 0.0, 0.0);
    Int3 maxCoords = FP3(1.0, 1.0, 1.0);
    FP3 steps((maxCoords.x - minCoords.x) / gridSize.x,
        (maxCoords.y - minCoords.y) / gridSize.y,
        (maxCoords.z - minCoords.z) / gridSize.z);
    Grid<FP, GridTypes::YeeGridType>* grid = new Grid<FP, GridTypes::YeeGridType>(gridSize, 0.5, minCoords, steps, gridSize);
    FDTD* pstd = new FDTD(grid);

    /*std::auto_ptr<Pml<GridTypes::YeeGridType>> pml;
    pml.reset(new PmlSpectral<GridTypes::YeeGridType>(&pstd, Int3(0, 0, 0)));
    pml.reset(new PmlPstd(&pstd, Int3(4, 0, 0)));*/

    //PmlSpectral<GridTypes::YeeGridType>* pml = new PmlSpectral<GridTypes::YeeGridType>(pstd, Int3(0, 0, 0));
    //delete pml;
    delete pstd;
    delete grid;
}