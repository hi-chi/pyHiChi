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
        SimpleGrid_Base(SimpleGrid &grid) : grid(grid) { name = "SimpleGrid"; }
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
        FDTD &fdtd;
        public:
            FDTD_Module(FDTD& fdtd) : fdtd(fdtd) { name = "FDTD_Module"; };
            void run() override
            {
                fdtd.updateFields();
            }
    };

    class Simulation
    {
    private:
        Grid_Base &grid;
        std::vector<Module*> modules;
        const unsigned int numSteps;
        unsigned int currentStep = 0;
    public:
        Simulation(YeeGrid &g, unsigned int numSteps) : grid(YeeGrid_Base(g)), numSteps(numSteps) {}
    
        void addModule(FDTD &module)
        {
            FDTD_Module *tmp = new FDTD_Module(module);
            modules.push_back(tmp);
        }
    
        void run()
        {
            currentStep = 0; // need deleted
            for (; currentStep < numSteps; ++currentStep)
            {
                for (int i = 0; i < modules.size(); i++)
                    modules[i]->run();
            }
        }
    };
}