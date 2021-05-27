#include <sstream>
#include "Simulation.h"
#include "TestingUtility.h"

#include "Grid.h"
#include "Fdtd.h"
#include "FieldEntity.h"


TEST(SimulationTest, Simulation_save_load)
{
    std::shared_ptr<FieldEntity<YeeGrid, FDTD>> ptr1, ptr2;
    FieldEntity<YeeGrid, FDTD>* fieldEntity = new FieldEntity<YeeGrid, FDTD>(Int3(8, 9, 10), FP3(0.0, 0.0, 0.0), FP3(1.0, 1.0, 1.0), 2e-12);
    FieldEntity<YeeGrid, FDTD>* fieldTmp = new FieldEntity<YeeGrid, FDTD>(Int3(10, 7, 129), FP3(0.1, -0.2, 0.3), FP3(1.1, 1.2, 1.3), 4e-14);
    ptr1.reset(fieldEntity);
    ptr2.reset(fieldTmp);

    Simulation<YeeGrid, FDTD> simulation(ptr1);
    std::stringstream sstr;
    simulation.save(sstr);
    simulation = Simulation<YeeGrid, FDTD>(ptr2);
    simulation.load(sstr);
    ASSERT_EQ(fieldEntity->fieldSolver->dt, simulation.field->fieldSolver->dt);

}