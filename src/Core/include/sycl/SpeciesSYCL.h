#pragma once
#include "Constants.h"
#include "Dimension.h"
#include "Particle.h"
#include "ParticleTypes.h"

#include "ParticleArraySYCL.h"

namespace pfc
{
    namespace sycl_pfc
    {
        template <Dimension dimension, ParticleTypes type, ParticleRepresentation storage>
        class Species
        {
        public:
            typedef typename ParticleArray<dimension, storage>::Type ParticleArray;

            typedef Particle<dimension>& ParticleRef;
            typedef const Particle<dimension>& ConstParticleRef;
            typedef ParticleProxy<dimension>& ParticleProxyRef;
            typedef const ParticleProxy<dimension>& ConstParticleProxyRef;

            typedef Particle<dimension> ParticleType;
            typedef ParticleProxy<dimension> ParticleProxyType;

            Species() :
                particles(type)
            {}

            inline int size() const { return particles.size(); }

            inline ParticleProxyType operator[](int idx)
            {
                return particles[idx];
            }

            inline ParticleProxyType back()
            {
                return particles.back();
            }

            inline void pushBack(ConstParticleRef particle)
            {
                if (particle.getType() == type) particles.pushBack(particle);
            }
            inline void popBack() { particles.pop_back(); }

            inline ParticleTypes getType()
            {
                return particles.getType();
            }
        
            ParticleArray particles;
        };
    }
}