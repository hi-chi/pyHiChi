#pragma once

#include "Dimension.h"
#include "Particle.h"
#include "ParticleTypes.h"
#include "ParticleTraits.h"
#include "Vectors.h"
#include "VectorsProxy.h"

#include <map>
#include <vector>
#include <string>
#include <functional>

namespace pfc {

    template<Dimension dimension>
    class ParticleArrayAccessor {
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

        typedef typename ScalarType<PositionType>::Type PositionElementType;
        typedef typename ScalarType<MomentumType>::Type MomentumElementType;

        typedef Particle<dimension> ParticleType;
        typedef Particle<dimension>& ParticleRef;
        typedef const Particle<dimension>& ConstParticleRef;
        typedef ParticleProxy<dimension> ParticleProxyType;


        static const int positionDimension = VectorDimensionHelper<PositionType>::dimension;
        static const int momentumDimension = VectorDimensionHelper<MomentumType>::dimension;

        inline int size() const { return static_cast<int>(_size); }

        template<class pArray>
        inline ParticleArrayAccessor(pArray* org) {};

        inline ParticleProxyType operator[](size_t idx) const;

    public:
        PositionElementType* positions[positionDimension];
        MomentumElementType* ps[momentumDimension];
        WeightType* weights;
        GammaType* gammas;
        ParticleTypes* typeIndex;
        size_t positions_stride[positionDimension];
        size_t ps_stride[momentumDimension];
        size_t weights_stride;
        size_t gammas_stride;
        size_t typeIndex_stride;
        size_t _size;
    };

    template<>
    inline ParticleArrayAccessor<One>::ParticleProxyType ParticleArrayAccessor<One>::operator[](size_t idx) const
    {
        PositionTypeProxy posProxy(positions[0][idx * positions_stride[0]]);
        MomentumTypeProxy momProxy(ps[0][idx * ps_stride[0]],
            ps[1][idx * ps_stride[1]], ps[2][idx * ps_stride[2]]);
        WeightTypeProxy weightRef = ref(weights[idx * weights_stride]);
        TypeIndexTypeProxy typeRef = ref(typeIndex[idx * typeIndex_stride]);
        GammaTypeProxy gammaRef = ref(gammas[idx * gammas_stride]);

        return ParticleProxyType(posProxy, momProxy, weightRef, typeRef, gammaRef);
    }

    template<>
    inline ParticleArrayAccessor<Two>::ParticleProxyType ParticleArrayAccessor<Two>::operator[](size_t idx) const
    {
        PositionTypeProxy posProxy(positions[0][idx * positions_stride[0]],
            positions[1][idx * positions_stride[1]]);
        MomentumTypeProxy momProxy(ps[0][idx * ps_stride[0]],
            ps[1][idx * ps_stride[1]], ps[2][idx * ps_stride[2]]);
        WeightTypeProxy weightRef = ref(weights[idx * weights_stride]);
        TypeIndexTypeProxy typeRef = ref(typeIndex[idx * typeIndex_stride]);
        GammaTypeProxy gammaRef = ref(gammas[idx * gammas_stride]);

        return ParticleProxyType(posProxy, momProxy, weightRef, typeRef, gammaRef);
    }

    template<>
    inline ParticleArrayAccessor<Three>::ParticleProxyType ParticleArrayAccessor<Three>::operator[](size_t idx) const
    {
        PositionTypeProxy posProxy(positions[0][idx * positions_stride[0]],
            positions[1][idx * positions_stride[1]], positions[2][idx * positions_stride[2]]);
        MomentumTypeProxy momProxy(ps[0][idx * ps_stride[0]],
            ps[1][idx * ps_stride[1]], ps[2][idx * ps_stride[2]]);
        WeightTypeProxy weightRef = ref(weights[idx * weights_stride]);
        TypeIndexTypeProxy typeRef = ref(typeIndex[idx * typeIndex_stride]);
        GammaTypeProxy gammaRef = ref(gammas[idx * gammas_stride]);

        return ParticleProxyType(posProxy, momProxy, weightRef, typeRef, gammaRef);
    }

    typedef ParticleArrayAccessor<Three> ParticleArrayAccessor3d;

} // namespace pfc