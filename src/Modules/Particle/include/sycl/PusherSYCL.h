#pragma once
/*#include "../Pusher.h"
#include "../../../../Core/include/sycl/ParticleArray.h"

#include "../../../../Core/include/Particle.h"
#include "../../../../Core/include/ParticleArrayAccessor.h"
#include "../../../../Core/include/FieldValue.h"
*/
#include "Pusher.h"
#include "sycl/ParticleArraySYCL.h"

#include "Particle.h"
#include "ParticleArrayAccessor.h"
#include "FieldValue.h"

#include <sycl/CL/sycl.hpp>


#include <array>
#include <vector>

#include <iostream>

namespace pfc
{
    namespace sycl_pfc
    {
        class ParticlePusherSYCL
        {
        public:
            template<class T_Particle>
            inline void operator()(T_Particle* particle, ValueField& field, FP timeStep) {};

            template<class T_ParticleArray>
            inline void operator()(T_ParticleArray* particleArray, std::vector<ValueField>& fields, FP timeStep) { };
        };

        class BorisPusherSYCL : public ParticlePusherSYCL
        {
        public:
            template<class T_Particle>
            inline void operator()(T_Particle* particle, ValueField& field, FP timeStep)
            {
                FP3 e = field.getE();
                FP3 b = field.getB();
                FP eCoeff = timeStep * particle->getCharge() / (2 * particle->getMass() * Constants<FP>::lightVelocity());
                FP3 eMomentum = e * eCoeff;
                FP3 um = particle->getP() + eMomentum;
                FP3 t = b * eCoeff / sqrt((FP)1 + um.norm2());
                FP3 uprime = um + cross(um, t);
                FP3 s = t * (FP)2 / ((FP)1 + t.norm2());
                particle->setP(eMomentum + um + cross(uprime, s));
                particle->setPosition(particle->getPosition() + timeStep * particle->getVelocity());
            }

            template<class T_ParticleArray>
            void operator()(Devices device, T_ParticleArray* particleArray, sycl_vector<ValueField>* fields, FP timeStep)
            {
                typedef typename T_ParticleArray::ParticleProxyType ParticleProxyType;
                try
                {
                    
                    ParticleArrayAccessor3d particles_accessor(particleArray);

                    ValueField* ptr_fields = fields->data();
                    ParticleType* types = ParticleInfo::typesVector.data();

                    auto kernel =[&](sycl::handler& h) {
                        h.parallel_for(sycl::range<1>(particleArray->size()),
                            [=](sycl::id<1> i)
                        {
                            //using namespace sycl;
                            ParticleProxyType particle = particles_accessor[i];
                            FP3 e = ptr_fields[i].E;
                            FP3 b = ptr_fields[i].B;
                            FP eCoeff = timeStep * types[particle.getType()].charge
                                / ((FP)2 * types[particle.getType()].mass * Constants<FP>::lightVelocity());
                            FP3 eMomentum = e * eCoeff;
                            FP3 um = particle.getP() + eMomentum;
                            FP3 t = b * eCoeff / sycl::sqrt((FP)1 + um.norm2());
                            FP3 uprime = um + cross(um, t);
                            FP3 s = t * (FP)2 / ((FP)1 + t.norm2());
                            particle.setP_SYCL(eMomentum + um + cross(uprime, s));
                            particle.setPosition(particle.getPosition() + timeStep * particle.getVelocity());
                        });
                    };

                    if (device == Devices::CPU)
                        node.cpu_device.submit(kernel).wait_and_throw();
                    if (device == Devices::GPU)
                        node.gpu_device.submit(kernel).wait_and_throw();
                }
                catch (sycl::exception & e)
                {
                    std::cout << e.what() << std::endl;
                }
                catch (std::exception& e)
                {
                    std::cout << e.what() << std::endl;
                }
            } 
        };

    }
}