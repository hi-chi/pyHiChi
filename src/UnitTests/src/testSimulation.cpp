#include "Simulation.h"
#include "TestingUtility.h"

TEST(SimulationTest, SimulationDefaultConstructor)
{
    ASSERT_NO_THROW(Simulation simulation());
}