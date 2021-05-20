#pragma once

#include "Dimension.h"
#include "Particle.h"
#include "ParticleArray.h"
#include "ParticleArrayAccessor.h"
#include "ParticleTypes.h"
#include "ParticleTraits.h"
#include "Vectors.h"
#include "VectorsProxy.h"

#include <map>
#include <vector>
#include <string>
#include <functional>

#include <sycl/CL/sycl.hpp>

#include "DeviceSYCL.h"

namespace pfc {

    namespace sycl_pfc
    {
        // Collection of particles with array-like semantics,
        // representation as array of structures
        template<Dimension dimension>
        class ParticleArrayAoS {
        public:

            typedef typename ParticleTraits<Particle<dimension>>::PositionType PositionType;
            typedef typename ParticleTraits<Particle<dimension>>::MomentumType MomentumType;
            typedef typename ParticleTraits<Particle<dimension>>::GammaType GammaType;
            typedef typename ParticleTraits<Particle<dimension>>::WeightType WeightType;
            typedef typename ParticleTraits<Particle<dimension>>::TypeIndexType TypeIndexType;
            typedef Particle<dimension> ParticleType;
            typedef Particle<dimension>& ParticleRef;
            typedef const Particle<dimension>& ConstParticleRef;
            static const int momentumDimension = VectorDimensionHelper<MomentumType>::dimension;
            typedef ParticleProxy<dimension> ParticleProxyType;


            typedef ParticleArrayAoS<dimension> typeArray;

            inline size_t size() const { return static_cast<size_t>(particles.size()); }

            ParticleArrayAoS(ParticleTypes type = Electron) : particles(node.default_device)
            {
                setType(type);
            }

            inline void setType(ParticleTypes type)
            {
                typeIndex = type;
            }

            inline ParticleTypes getType()
            {
                return static_cast<ParticleTypes>(typeIndex);
            }

            inline ParticleProxyType operator[](int idx)
            {
                return ParticleProxyType(particles.data()[idx]);
            }

            inline ParticleProxyType back()
            {
                return operator[](this->size() - 1);
            }

            inline void pushBack(ConstParticleRef particle)
            {
                if (particle.getType() == typeIndex)
                    particles.push_back(particle);
            }
            inline void popBack() { particles.pop_back(); }

            inline void deleteParticle(int idx)
            {
                if (idx < this->size())
                {
                    std::swap(particles[idx], particles[this->size() - 1]);
                    particles.pop_back();
                }
            }

            inline void clear()
            {
                particles.clear();
            }

        public:
            ParticleTypes typeIndex;
            sycl_vector<ParticleType> particles;
        };

        // Collection of particles with array-like semantics,
        // representation as structure of arrays
        template<Dimension dimension>
        class ParticleArraySoA {
        public:

            typedef typename ParticleTraits<Particle<dimension>>::PositionType PositionType;
            typedef typename ParticleTraits<Particle<dimension>>::MomentumType MomentumType;
            typedef typename ParticleTraits<Particle<dimension>>::GammaType GammaType;
            typedef typename ParticleTraits<Particle<dimension>>::WeightType WeightType;
            typedef typename ParticleTraits<Particle<dimension>>::TypeIndexType TypeIndexType;

            typedef typename ParticleTraits<Particle<dimension>>::PositionTypeProxy PositionTypeProxy;
            typedef typename ParticleTraits<Particle<dimension>>::MomentumTypeProxy MomentumTypeProxy;
            typedef typename ParticleTraits<Particle<dimension>>::GammaTypeProxy GammaTypeProxy;
            typedef typename ParticleTraits<Particle<dimension>>::WeightTypeProxy WeightTypeProxy;
            typedef typename ParticleTraits<Particle<dimension>>::TypeIndexTypeProxy TypeIndexTypeProxy;

            typedef Particle<dimension> ParticleType;
            typedef Particle<dimension>& ParticleRef;
            typedef const Particle<dimension>& ConstParticleRef;
            typedef ParticleProxy<dimension> ParticleProxyType;


            typedef ParticleArraySoA<dimension> typeArray;

            static const int positionDimension = VectorDimensionHelper<PositionType>::dimension;
            static const int momentumDimension = VectorDimensionHelper<MomentumType>::dimension;

            inline int size() const { return static_cast<int>(weights.size()); }

            ParticleArraySoA(ParticleTypes type = Electron);

            inline void setType(ParticleTypes type)
            {
                typeIndex = type;
            }

            inline ParticleTypes getType()
            {
                return static_cast<ParticleTypes>(typeIndex);
            }

            inline ParticleProxyType operator[](int idx);

            inline ParticleProxyType back()
            {
                return operator[](this->size() - 1);
            }

            inline void pushBack(ConstParticleRef particle)
            {
                if (particle.getType() == typeIndex)
                {
                    const PositionType position = particle.getPosition();
                    for (int d = 0; d < positionDimension; d++)
                        positions[d].push_back(position[d]);
                    const MomentumType p = particle.getP();
                    for (int d = 0; d < momentumDimension; d++)
                        ps[d].push_back(p[d]);
                    weights.push_back(particle.getWeight());
                    gammas.push_back(particle.getGamma());
                    typeIndexs.push_back(particle.getType());
                }

            }
            inline void popBack()
            {
                for (int d = 0; d < positionDimension; d++)
                    positions[d].pop_back();
                for (int d = 0; d < momentumDimension; d++)
                    ps[d].pop_back();
                weights.pop_back();
                gammas.pop_back();
                typeIndexs.pop_back();
            }
            
            inline void deleteParticle(int idx)
            {
                if (idx < this->size())
                {
                    int size = this->size();
                    for (int d = 0; d < positionDimension; d++)
                    {
                        std::swap(positions[d][idx], positions[d][size - 1]);
                        positions[d].pop_back();
                    }
                    for (int d = 0; d < momentumDimension; d++)
                    {
                        std::swap(ps[d][idx], ps[d][size - 1]);
                        ps[d].pop_back();
                    }

                    std::swap(weights[idx], weights[size - 1]);
                    weights.pop_back();
                    std::swap(gammas[idx], gammas[size - 1]);
                    gammas.pop_back();
                    std::swap(typeIndexs[idx], typeIndexs[size - 1]);
                    typeIndexs.pop_back();
                }
            }

            inline void clear()
            {
                for (int d = 0; d < positionDimension; d++)
                    positions[d].clear();
                for (int d = 0; d < momentumDimension; d++)
                    ps[d].clear();
                weights.clear();
                gammas.clear();
                typeIndexs.clear();
            }

        public:
            sycl_vector<typename ScalarType<PositionType>::Type> positions[positionDimension];
            sycl_vector<typename ScalarType<MomentumType>::Type> ps[momentumDimension];
            sycl_vector<WeightType> weights;
            sycl_vector<GammaType> gammas;
            sycl_vector<ParticleTypes> typeIndexs;
            ParticleTypes typeIndex;
        };

        template<>
        ParticleArraySoA<One>::ParticleArraySoA(ParticleTypes type) :
            positions{ sycl_vector<typename ScalarType<PositionType>::Type>(node.default_device) },
            ps{ sycl_vector<typename ScalarType<MomentumType>::Type>(node.default_device),
                sycl_vector<typename ScalarType<MomentumType>::Type>(node.default_device),
                sycl_vector<typename ScalarType<MomentumType>::Type>(node.default_device) },
            weights{ sycl_vector<WeightType>(node.default_device) },
            gammas{ sycl_vector<GammaType>(node.default_device) },
            typeIndexs{ sycl_vector<ParticleTypes>(node.default_device) }

        {
            setType(type);
        }

        template<>
        ParticleArraySoA<Two>::ParticleArraySoA(ParticleTypes type) :
            positions{ sycl_vector<typename ScalarType<PositionType>::Type>(node.default_device),
                       sycl_vector<typename ScalarType<PositionType>::Type>(node.default_device) },
            ps{ sycl_vector<typename ScalarType<MomentumType>::Type>(node.default_device),
                sycl_vector<typename ScalarType<MomentumType>::Type>(node.default_device),
                sycl_vector<typename ScalarType<MomentumType>::Type>(node.default_device) },
            weights{ sycl_vector<WeightType>(node.default_device) },
            gammas{ sycl_vector<GammaType>(node.default_device) },
            typeIndexs{ sycl_vector<ParticleTypes>(node.default_device) }

        {
            setType(type);
        }

        template<>
        ParticleArraySoA<Three>::ParticleArraySoA(ParticleTypes type) :
            positions{ sycl_vector<typename ScalarType<PositionType>::Type>(node.default_device),
                       sycl_vector<typename ScalarType<PositionType>::Type>(node.default_device),
                       sycl_vector<typename ScalarType<PositionType>::Type>(node.default_device) },
            ps{ sycl_vector<typename ScalarType<MomentumType>::Type>(node.default_device),
                sycl_vector<typename ScalarType<MomentumType>::Type>(node.default_device),
                sycl_vector<typename ScalarType<MomentumType>::Type>(node.default_device) },
            weights{ sycl_vector<WeightType>(node.default_device) },
            gammas{ sycl_vector<GammaType>(node.default_device) },
            typeIndexs{ sycl_vector<ParticleTypes>(node.default_device) }

        {
            setType(type);
        }

        template<>
        inline ParticleArraySoA<One>::ParticleProxyType ParticleArraySoA<One>::operator[](int idx)
        {
            PositionTypeProxy posProxy(positions[0][idx]);
            MomentumTypeProxy momProxy(ps[0][idx], ps[1][idx], ps[2][idx]);
            WeightTypeProxy weightRef = ref(weights[idx]);
            TypeIndexTypeProxy typeRef = ref(typeIndexs[idx]);
            GammaTypeProxy gammaRef = ref(gammas[idx]);

            return ParticleProxyType(posProxy, momProxy, weightRef, typeRef, gammaRef);
        }

        template<>
        inline ParticleArraySoA<Two>::ParticleProxyType ParticleArraySoA<Two>::operator[](int idx)
        {
            PositionTypeProxy posProxy(positions[0][idx], positions[1][idx]);
            MomentumTypeProxy momProxy(ps[0][idx], ps[1][idx], ps[2][idx]);
            WeightTypeProxy weightRef = ref(weights[idx]);
            TypeIndexTypeProxy typeRef = ref(typeIndexs[idx]);
            GammaTypeProxy gammaRef = ref(gammas[idx]);

            return ParticleProxyType(posProxy, momProxy, weightRef, typeRef, gammaRef);
        }

        template<>
        inline ParticleArraySoA<Three>::ParticleProxyType ParticleArraySoA<Three>::operator[](int idx)
        {
            PositionTypeProxy posProxy(positions[0][idx], positions[1][idx], positions[2][idx]);
            MomentumTypeProxy momProxy(ps[0][idx], ps[1][idx], ps[2][idx]);
            WeightTypeProxy weightRef = ref(weights[idx]);
            TypeIndexTypeProxy typeRef = ref(typeIndexs[idx]);
            GammaTypeProxy gammaRef = ref(gammas[idx]);

            return ParticleProxyType(posProxy, momProxy, weightRef, typeRef, gammaRef);
        }


        // Traits class to provide a Type corresponding to array of particles
        // according to the given representation
        template<Dimension dimension, ParticleRepresentation storage>
        struct ParticleArray {
        };

        template<Dimension dimension>
        struct ParticleArray<dimension, ParticleRepresentation_AoS> {
            typedef ParticleArrayAoS<dimension> Type;
        };

        template<Dimension dimension>
        struct ParticleArray<dimension, ParticleRepresentation_SoA> {
            typedef ParticleArraySoA<dimension> Type;
        };


        typedef typename ParticleArray<Three, ParticleRepresentation_SoA>::Type ParticleArray3d;
        typedef typename ParticleArray<Three, ParticleRepresentation_AoS>::Type ParticleArrayAoS3d;
    }  //namespace sycl_pfc

    template<>
    template<>
    inline ParticleArrayAccessor<One>::ParticleArrayAccessor(sycl_pfc::ParticleArrayAoS<One>* org)
    {
        positions[0] = (PositionElementType*)((char*)(org->particles.data()) + offsetof(ParticleType, position.x));
        positions_stride[0] = sizeof(ParticleType) / sizeof(PositionElementType);

        ps[0] = (MomentumElementType*)((char*)(org->particles.data()) + offsetof(ParticleType, p.x));
        ps[1] = (MomentumElementType*)((char*)(org->particles.data()) + offsetof(ParticleType, p.y));
        ps[2] = (MomentumElementType*)((char*)(org->particles.data()) + offsetof(ParticleType, p.z));
        ps_stride[0] = sizeof(ParticleType) / sizeof(MomentumElementType);
        ps_stride[1] = sizeof(ParticleType) / sizeof(MomentumElementType);
        ps_stride[2] = sizeof(ParticleType) / sizeof(MomentumElementType);

        weights = (WeightType*)((char*)(org->particles.data()) + offsetof(ParticleType, weight));
        gammas = (GammaType*)((char*)(org->particles.data()) + offsetof(ParticleType, gamma));
        typeIndex = (TypeIndexType*)((char*)(org->particles.data()) + offsetof(ParticleType, typeIndex));
        weights_stride = sizeof(ParticleType) / sizeof(WeightType);
        gammas_stride = sizeof(ParticleType) / sizeof(GammaType);
        typeIndex_stride = sizeof(ParticleType) / sizeof(TypeIndexType);
        _size = org->particles.size();

        if (sizeof(ParticleType) % sizeof(PositionElementType) ||
            sizeof(ParticleType) % sizeof(MomentumElementType) ||
            sizeof(ParticleType) % sizeof(WeightType) ||
            sizeof(ParticleType) % sizeof(GammaType) ||
            sizeof(ParticleType) % sizeof(TypeIndexType))
            std::cout << "Wrong create ParticleAccessor" << std::endl;
    };

    template<>
    template<>
    inline ParticleArrayAccessor<Two>::ParticleArrayAccessor(sycl_pfc::ParticleArrayAoS<Two>* org)
    {
        positions[0] = (PositionElementType*)((char*)(org->particles.data()) + offsetof(ParticleType, position.x));
        positions[1] = (PositionElementType*)((char*)(org->particles.data()) + offsetof(ParticleType, position.y));
        positions_stride[0] = sizeof(ParticleType) / sizeof(PositionElementType);
        positions_stride[1] = sizeof(ParticleType) / sizeof(PositionElementType);

        ps[0] = (MomentumElementType*)((char*)(org->particles.data()) + offsetof(ParticleType, p.x));
        ps[1] = (MomentumElementType*)((char*)(org->particles.data()) + offsetof(ParticleType, p.y));
        ps[2] = (MomentumElementType*)((char*)(org->particles.data()) + offsetof(ParticleType, p.z));
        ps_stride[0] = sizeof(ParticleType) / sizeof(MomentumElementType);
        ps_stride[1] = sizeof(ParticleType) / sizeof(MomentumElementType);
        ps_stride[2] = sizeof(ParticleType) / sizeof(MomentumElementType);

        weights = (WeightType*)((char*)(org->particles.data()) + offsetof(ParticleType, weight));
        gammas = (GammaType*)((char*)(org->particles.data()) + offsetof(ParticleType, gamma));
        typeIndex = (TypeIndexType*)((char*)(org->particles.data()) + offsetof(ParticleType, typeIndex));
        weights_stride = sizeof(ParticleType) / sizeof(WeightType);
        gammas_stride = sizeof(ParticleType) / sizeof(GammaType);
        typeIndex_stride = sizeof(ParticleType) / sizeof(TypeIndexType);
        _size = org->particles.size();

        if (sizeof(ParticleType) % sizeof(PositionElementType) ||
            sizeof(ParticleType) % sizeof(MomentumElementType) ||
            sizeof(ParticleType) % sizeof(WeightType) ||
            sizeof(ParticleType) % sizeof(GammaType) ||
            sizeof(ParticleType) % sizeof(TypeIndexType))
            std::cout << "Wrong create ParticleAccessor" << std::endl;
    };

    template<>
    template<>
    inline ParticleArrayAccessor<Three>::ParticleArrayAccessor(sycl_pfc::ParticleArrayAoS<Three>* org)
    {
        positions[0] = (PositionElementType*)((char*)(org->particles.data()) + offsetof(ParticleType, position.x));
        positions[1] = (PositionElementType*)((char*)(org->particles.data()) + offsetof(ParticleType, position.y));
        positions[2] = (PositionElementType*)((char*)(org->particles.data()) + offsetof(ParticleType, position.z));
        positions_stride[0] = sizeof(ParticleType) / sizeof(PositionElementType);
        positions_stride[1] = sizeof(ParticleType) / sizeof(PositionElementType);
        positions_stride[2] = sizeof(ParticleType) / sizeof(PositionElementType);

        ps[0] = (MomentumElementType*)((char*)(org->particles.data()) + offsetof(ParticleType, p.x));
        ps[1] = (MomentumElementType*)((char*)(org->particles.data()) + offsetof(ParticleType, p.y));
        ps[2] = (MomentumElementType*)((char*)(org->particles.data()) + offsetof(ParticleType, p.z));
        ps_stride[0] = sizeof(ParticleType) / sizeof(MomentumElementType);
        ps_stride[1] = sizeof(ParticleType) / sizeof(MomentumElementType);
        ps_stride[2] = sizeof(ParticleType) / sizeof(MomentumElementType);

        weights = (WeightType*)((char*)(org->particles.data()) + offsetof(ParticleType, weight));
        gammas = (GammaType*)((char*)(org->particles.data()) + offsetof(ParticleType, gamma));
        typeIndex = (TypeIndexType*)((char*)(org->particles.data()) + offsetof(ParticleType, typeIndex));
        weights_stride = sizeof(ParticleType) / sizeof(WeightType);
        gammas_stride = sizeof(ParticleType) / sizeof(GammaType);
        typeIndex_stride = sizeof(ParticleType) / sizeof(TypeIndexType);
        _size = org->particles.size();

        if (sizeof(ParticleType) % sizeof(PositionElementType) ||
            sizeof(ParticleType) % sizeof(MomentumElementType) ||
            sizeof(ParticleType) % sizeof(WeightType) ||
            sizeof(ParticleType) % sizeof(GammaType) ||
            sizeof(ParticleType) % sizeof(TypeIndexType))
            std::cout << "Wrong create ParticleAccessor" << std::endl;
    };

    template<>
    template<>
    inline ParticleArrayAccessor<One>::ParticleArrayAccessor(sycl_pfc::ParticleArraySoA<One>* org)
    {
        positions[0] = (PositionElementType*)(org->positions[0].data());
        positions_stride[0] = 1;

        ps[0] = (MomentumElementType*)(org->ps[0].data());
        ps[1] = (MomentumElementType*)(org->ps[1].data());
        ps[2] = (MomentumElementType*)(org->ps[2].data());
        ps_stride[0] = 1;
        ps_stride[1] = 1;
        ps_stride[2] = 1;

        weights = (WeightType*)(org->weights.data());
        gammas = (GammaType*)(org->gammas.data());
        typeIndex = (TypeIndexType*)(org->typeIndexs.data());
        weights_stride = 1;
        gammas_stride = 1;
        typeIndex_stride = 1;
        _size = org->weights.size();
    };

    template<>
    template<>
    inline ParticleArrayAccessor<Two>::ParticleArrayAccessor(sycl_pfc::ParticleArraySoA<Two>* org)
    {
        positions[0] = (PositionElementType*)(org->positions[0].data());
        positions[1] = (PositionElementType*)(org->positions[1].data());
        positions_stride[0] = 1;
        positions_stride[1] = 1;

        ps[0] = (MomentumElementType*)(org->ps[0].data());
        ps[1] = (MomentumElementType*)(org->ps[1].data());
        ps[2] = (MomentumElementType*)(org->ps[2].data());
        ps_stride[0] = 1;
        ps_stride[1] = 1;
        ps_stride[2] = 1;

        weights = (WeightType*)(org->weights.data());
        gammas = (GammaType*)(org->gammas.data());
        typeIndex = (TypeIndexType*)(org->typeIndexs.data());
        weights_stride = 1;
        gammas_stride = 1;
        typeIndex_stride = 1;
        _size = org->weights.size();
    };

    template<>
    template<>
    inline ParticleArrayAccessor<Three>::ParticleArrayAccessor(sycl_pfc::ParticleArraySoA<Three>* org)
    {
        positions[0] = (PositionElementType*)(org->positions[0].data());
        positions[1] = (PositionElementType*)(org->positions[1].data());
        positions[2] = (PositionElementType*)(org->positions[2].data());
        positions_stride[0] = 1;
        positions_stride[1] = 1;
        positions_stride[2] = 1;

        ps[0] = (MomentumElementType*)(org->ps[0].data());
        ps[1] = (MomentumElementType*)(org->ps[1].data());
        ps[2] = (MomentumElementType*)(org->ps[2].data());
        ps_stride[0] = 1;
        ps_stride[1] = 1;
        ps_stride[2] = 1;

        weights = (WeightType*)(org->weights.data());
        gammas = (GammaType*)(org->gammas.data());
        typeIndex = (TypeIndexType*)(org->typeIndexs.data());
        weights_stride = 1;
        gammas_stride = 1;
        typeIndex_stride = 1;
        _size = org->weights.size();
    };
} // namespace pfc