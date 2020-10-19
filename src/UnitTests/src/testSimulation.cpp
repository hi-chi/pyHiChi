#include "TestingUtility.h"
#include "Fdtd.h"
#include "Simulation.h"
#include "DataManager.h"
using namespace pfc;

template <typename Real>
class SimulationTest : public ::testing::Test {
};
typedef ::testing::Types<float, double> types;
TYPED_TEST_CASE(SimulationTest, types);
TYPED_TEST(SimulationTest, Can_create_Simulation_with_Grid_and_FDTD)
{
    Int3 gridSize(11, 5, 6);
    FP3 minCoords = FP3(-1.0, 0.0, 0.0);
    FP3 maxCoords = FP3(1.0, 1.0, 1.0);
    FP3 steps((maxCoords.x - minCoords.x) / gridSize.x,
        (maxCoords.y - minCoords.y) / gridSize.y,
        (maxCoords.z - minCoords.z) / gridSize.z);
    FP timeStep = 1e-15;
    const unsigned int numSteps = 100;

    YeeGrid grid(gridSize, timeStep, minCoords, steps, gridSize);
    FDTD fdtd(&grid);
    Simulation simulation = Simulation(grid, numSteps);
    simulation.addModule(fdtd);
    ASSERT_NO_THROW(simulation.run());
}

template <typename Real>
class DataManagerTests : public ::testing::Test {};
typedef ::testing::Types<float, double> types;
TYPED_TEST_CASE(DataManagerTests, types);
TYPED_TEST(DataManagerTests, Test_Adios2)
{
    adios2::ADIOS adios;
    adios2::IO bpIO = adios.DeclareIO("BPFile_N2N");
    adios2::Engine engine = bpIO.Open("adios", adios2::Mode::Write);
    double w1 = 3.14;
    std::string s = "data";
    adios2::Variable<double> var1 = bpIO.DefineVariable<double>(s);
    ASSERT_NO_THROW(engine.Put(var1, w1));
    engine.Close();
    bpIO.RemoveAllVariables();

    adios2::Variable<double> var2 = bpIO.DefineVariable<double>(s);
    double r1;
    engine = bpIO.Open("adios", adios2::Mode::Read);
    ASSERT_NO_THROW(engine.Get(var2, r1));
    EXPECT_DOUBLE_EQ(w1, r1);
    engine.Close();
}
TYPED_TEST(DataManagerTests, Can_Put_Get_double)
{
    string var_name = "pi", path = "test_path";
    const double pi = 3.14;
    DataManager manager(path, IOType::Write);
    ASSERT_NO_THROW(manager.putVariable(var_name, pi));

    manager.setEngine(IOType::Read);
    double res;
    manager.getVariable(var_name, res);
    EXPECT_DOUBLE_EQ(res, pi);
}
TYPED_TEST(DataManagerTests, Can_Put_Get_Int3)
{
    string var_name = "int3Data", path = "test_path";
    DataManager manager(path, IOType::Write);
    const Int3 tmp(1, -2, 3);
    ASSERT_NO_THROW(manager.customPut(var_name, tmp));

    manager.setEngine(IOType::Read);
    Int3 res;
    manager.customGet(var_name, res);
    EXPECT_EQ(res, tmp);
}
TYPED_TEST(DataManagerTests, Can_Put_Get_ScalarField)
{
    string var_name = "field", path = "test_field";
    DataManager manager(path, IOType::Write);
    Int3 size(1, 2, 3);
    ScalarField<double> tmp(size);
    tmp(0, 0, 0) = 1.0;
    tmp(0, 1, 1) = 2.0;
    ASSERT_NO_THROW(manager.customPut(var_name, tmp));

    manager.setEngine(IOType::Read);
    ScalarField<double> res;
    manager.customGet(var_name, res);
    manager.endStep();
    EXPECT_EQ(res.getSize(), tmp.getSize());
    EXPECT_EQ(res(0, 0, 0), tmp(0, 0, 0));
    EXPECT_EQ(res(0, 1, 1), tmp(0, 1, 1));
}
