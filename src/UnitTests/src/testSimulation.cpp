#include <sstream>
#include "Simulation.h"
#include "TestingUtility.h"
#include "Grid.h"
#include "Fdtd.h"
#include "FieldEntity.h"

using namespace pfc;

template <class ParticleType>
class SimulationTest : public BaseParticleFixture<ParticleType> {
public:
    FieldEntity<YeeGrid, FDTD>* fieldEntity = new FieldEntity<YeeGrid, FDTD>(Int3(8, 9, 10), FP3(0.0, 0.0, 0.0), FP3(1.0, 1.0, 1.0), 2e-12);
    FieldEntity<YeeGrid, FDTD>* fieldTmp = new FieldEntity<YeeGrid, FDTD>(Int3(10, 7, 129), FP3(0.1, -0.2, 0.3), FP3(1.1, 1.2, 1.3), 4e-15);
    Ensemble<ParticleArray3d>* ensemble = new Ensemble<ParticleArray3d>();
    BorisPusher* pusher = new BorisPusher();

    std::shared_ptr<FieldEntity<YeeGrid, FDTD>> ptrField, ptrTmp;
    std::shared_ptr<Ensemble<ParticleArray3d>> ptrEnsemble;
    std::shared_ptr<BorisPusher> ptrPusher;

    virtual void SetUp()
    {
        ptrField.reset(fieldEntity);
        ptrTmp.reset(fieldTmp);
        ptrEnsemble.reset(ensemble);
        ptrPusher.reset(pusher);
        maxRelativeError = 1e-12;
    }
};
typedef ::testing::Types<Particle3d> types;
TYPED_TEST_CASE(SimulationTest, types);

TYPED_TEST(SimulationTest, Simulation_save_load)
{
    Simulation<YeeGrid, FDTD, ParticleArray3d> simulation(ptrField);
    std::stringstream sstr;
    simulation.save(sstr);
    //simulation = Simulation<YeeGrid, FDTD, ParticleArray3d>(ptrTmp);
    simulation.load(sstr);
    ASSERT_EQ(fieldEntity->fieldSolver->dt, simulation.field->fieldSolver->dt);
    ASSERT_EQ(fieldEntity->grid->numCells, simulation.field->grid->numCells);
}

TYPED_TEST(SimulationTest, Simulation_boris_pusher)
{
    for (int i = 0; i < 9; i++)
        ensemble->addParticle(this->randomParticle());
    ensemble->addParticle(this->randomParticle(Photon));
    ensemble->addParticle(this->randomParticle(Positron));
    ensemble->addParticle(this->randomParticle(Proton));
    Simulation<YeeGrid, FDTD, ParticleArray3d> simulation(ptrField, ptrEnsemble, ptrPusher);
    simulation.runIteration();
}

TYPED_TEST(SimulationTest, Simulation_boris_pusher_migration)
{
    for (int i = 0; i < 9; i++)
        ensemble->addParticle(this->randomParticle());
    ensemble->addParticle(this->randomParticle(Photon));
    ensemble->addParticle(this->randomParticle(Positron));
    ensemble->addParticle(this->randomParticle(Proton));
    Simulation<YeeGrid, FDTD, ParticleArray3d> simulation(ptrField, ptrEnsemble, ptrPusher);
    simulation.runIteration();   
    //ensemble->getMigrationParticles(fieldEntity->getGrid());
}