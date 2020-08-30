#pragma once
#include "GridTypes.h"
#include "Grid.h"
#include "Fdtd.h"

namespace pfc {


    struct Grid_Base
    {
        std::string name;
        virtual void *getPGrid() = 0;
    };

    class YeeGrid_Base : public Grid_Base
    {
        YeeGrid &grid;
    public:
        YeeGrid_Base(const Int3 & _numInternalCells, FP _dt, const FP3 & minCoords, const FP3 & _steps, const Int3 & globalGridDims) :
                     grid(YeeGrid(_numInternalCells, _dt, minCoords, _steps, globalGridDims)){ name = "YeeGrid"; }
        YeeGrid_Base(YeeGrid &grid): grid(grid) { name = "YeeGrid"; }
        YeeGrid& getGrid()
        {
            return grid;
        }
        void* getPGrid() override
        {
            return &grid;
        }
    };

    class SimpleGrid_Base : public Grid_Base
    {
        SimpleGrid& grid;
    public:
        SimpleGrid_Base(const Int3 & _numInternalCells, FP _dt,
            const FP3 & minCoords, const FP3 & _steps, const Int3 & globalGridDims) :
            grid(SimpleGrid(_numInternalCells, _dt, minCoords, _steps, globalGridDims))
        {
            name = "SimpleGrid";
        }
        SimpleGrid& getGrid()
        {
            return grid;
        }
        void* getPGrid() override
        {
            return &grid;
        }
    };

    class Module
    {
    public:
        std::string name;
        virtual void run() {};
    };

    class FDTD_Module : public Module
    {
        FDTD *fdtd;
        public:
            FDTD_Module(YeeGrid_Base* grid) : fdtd(new FDTD(&grid->getGrid())) { name = "FDTD_Module"; };
            FDTD_Module(FDTD* fdtd) : fdtd(fdtd) { name = "FDTD_Module"; };
            void run() override
            {
                fdtd->updateFields();
            }
    };

    class Simulation
    {
    private:
        Grid_Base *grid;
        std::vector<Module*> modules;
        const unsigned int numSteps;
        unsigned int currentStep = 0;
    public:
        Simulation(YeeGrid &g, unsigned int numSteps) : grid(new YeeGrid_Base(g)), numSteps(numSteps) {}
    
        void addModule(Module *module)
        {
            modules.push_back(module);
        }
        void addModule(FDTD *module)
        {
            FDTD_Module* m = new FDTD_Module(module);
            modules.push_back(m);
        }
    
        void run()
        {
            currentStep = 0;
            for (; currentStep < numSteps; ++currentStep)
            {
                for (int i = 0; i < modules.size(); i++)
                    modules[i]->run();
            }
        }
    
        Grid_Base* getGrid()
        {
            return grid;
        }
    };
}