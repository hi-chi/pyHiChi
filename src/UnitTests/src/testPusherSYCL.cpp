#include "TestingUtility.h"

#include "Pusher.h"
#include "sycl/PusherSYCL.h"
#include "sycl/SpeciesSYCL.h"
#include "sycl/DeviceSYCL.h"

template <class SpeciesArrayType>
class PusherSYCLTest : public SpeciesTest<SpeciesArrayType> {
};


typedef ::testing::Types<
    sycl_pfc::Species<Three, Electron, ParticleRepresentation_AoS>,
    sycl_pfc::Species<Three, Positron, ParticleRepresentation_AoS>,
    sycl_pfc::Species<Three, Proton, ParticleRepresentation_AoS>,
    sycl_pfc::Species<Three, Electron, ParticleRepresentation_SoA>,
    sycl_pfc::Species<Three, Positron, ParticleRepresentation_SoA>,
    sycl_pfc::Species<Three, Proton, ParticleRepresentation_SoA>
> types;
TYPED_TEST_CASE(PusherSYCLTest, types);

TYPED_TEST(PusherSYCLTest, BorisPusherSYCL_CPU)
{
    typedef typename SpeciesTest<TypeParam>::SpeciesArray SpeciesArray;
    typedef typename SpeciesTest<TypeParam>::MomentumType MomentumType;

    SpeciesArray * speciesParticles = new SpeciesArray();
    SpeciesArray * speciesParticlesSYCL = new SpeciesArray();
    std::vector<ValueField> fields;
    sycl_pfc::sycl_vector<ValueField> * fieldsSYCL = new sycl_pfc::sycl_vector<ValueField>{ sycl_pfc::node.cpu_device };;
    FP Ex = 0.0, Ey = 0.0, Ez = 0.0;
    FP Bx = 1.0, By = 1.0, Bz = 1.0;
    int numParticles = 12;
    for (int i = 0; i < numParticles; i++)
    {
        speciesParticles->pushBack(this->randomParticle(speciesParticles->getType()));
        speciesParticlesSYCL->pushBack((*speciesParticles)[i]);
        fields.push_back(ValueField(Ex, Ey, Ez, Bx, By, Bz));
        fieldsSYCL->push_back(ValueField(Ex, Ey, Ez, Bx, By, Bz));
    }

    BorisPusher pusher;
    FP timeStep = 0.01;
    pusher(speciesParticles, fields, timeStep);


    sycl_pfc::BorisPusherSYCL pusherSYCL;
    pusherSYCL(sycl_pfc::Devices::CPU, &(speciesParticlesSYCL->particles), fieldsSYCL, timeStep);

    ASSERT_TRUE(this->eqParticleArrays(*speciesParticles, *speciesParticlesSYCL));
}

TYPED_TEST(PusherSYCLTest, BorisPusherSYCL_GPU)
{
    if (sizeof(FP) != 8)
    {
        typedef typename SpeciesTest<TypeParam>::SpeciesArray SpeciesArray;
        typedef typename SpeciesTest<TypeParam>::MomentumType MomentumType;

        SpeciesArray* speciesParticles = new SpeciesArray();
        SpeciesArray* speciesParticlesSYCL = new SpeciesArray();
        std::vector<ValueField> fields;
        sycl_pfc::sycl_vector<ValueField>* fieldsSYCL = new sycl_pfc::sycl_vector<ValueField>{ sycl_pfc::node.gpu_device };;
        FP Ex = 0.0, Ey = 0.0, Ez = 0.0;
        FP Bx = 1.0, By = 1.0, Bz = 1.0;
        int numParticles = 12;
        for (int i = 0; i < numParticles; i++)
        {
            speciesParticles->pushBack(this->randomParticle(speciesParticles->getType()));
            speciesParticlesSYCL->pushBack((*speciesParticles)[i]);
            fields.push_back(ValueField(Ex, Ey, Ez, Bx, By, Bz));
            fieldsSYCL->push_back(ValueField(Ex, Ey, Ez, Bx, By, Bz));
        }

        BorisPusher pusher;
        FP timeStep = 0.01;
        pusher(speciesParticles, fields, timeStep);


        sycl_pfc::BorisPusherSYCL pusherSYCL;
        pusherSYCL(sycl_pfc::Devices::GPU, &(speciesParticlesSYCL->particles), fieldsSYCL, timeStep);

        ASSERT_TRUE(this->eqParticleArrays(*speciesParticles, *speciesParticlesSYCL));
    }
}