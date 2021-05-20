#pragma once

#include "Dimension.h"
#include "Particle.h"
#include "ParticleArrayAccessor.h"
#include "ParticleTypes.h"
#include "ParticleTraits.h"
#include "Vectors.h"
#include "VectorsProxy.h"

#include <map>
#include <vector>
#include <string>
#include <functional>

#include <cstddef>

namespace pfc {

    template<typename pArray_t, typename ParticleType>
    class iteratorPArray : public std::iterator<std::random_access_iterator_tag, ParticleType, size_t>
    {
    private:
        pArray_t* pPArray;
        size_t index;
    public:

        iteratorPArray(pArray_t* _pArray) : pPArray(_pArray), index(0) {};
        iteratorPArray(pArray_t* _pArray, size_t _index) : pPArray(_pArray), index(_index) {};
        iteratorPArray(const iteratorPArray& it) : pPArray(it.pPArray), index(it.index) {};

        size_t getIdx() { return index; };

        ParticleType operator *() { return pPArray->operator[](index); }
        const iteratorPArray& operator ++() { ++index; return *this; }
        const iteratorPArray& operator --() { --index; return *this; }
        iteratorPArray operator ++(int)
        {
            iteratorPArray copy(*this);
            ++index;
            return copy;
        }
        iteratorPArray operator --(int)
        {
            iteratorPArray copy(*this);
            --index;
            return copy;
        }
        iteratorPArray& operator =(const iteratorPArray& other)
        {
            this->pParray = other.pPArray;
            this->index = other.index;
            return *this;
        }

        bool operator ==(const iteratorPArray& other) const
        {
            return index == other.index;
        }
        bool operator !=(const iteratorPArray& other) const
        {
            return index != other.index;
        }
        bool operator  <(const iteratorPArray& other) const
        {
            return index < other.index;
        }
        bool operator  >(const iteratorPArray& other) const
        {
            return index > other.index;
        }
        bool operator  <=(const iteratorPArray& other) const
        {
            return index <= other.index;
        }
        bool operator  >=(const iteratorPArray& other) const
        {
            return index >= other.index;
        }

        iteratorPArray& operator +(const size_t& add) const {
            iteratorPArray copy(*this);
            copy.index += add;
            return copy;
        }
        iteratorPArray& operator +=(const size_t& add) {
            index += add;
            return *this;
        }
        iteratorPArray& operator -(const size_t& add) const {
            iteratorPArray copy(*this);
            copy.index -= add;
            return copy;
        }
        iteratorPArray& operator -=(const long int& add) {
            index -= add;
            return *this;
        }

        ParticleType operator [](const size_t& n) const
        {
            return *pPArray[index + n];
        }
    };



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
        typedef iteratorPArray<typeArray, ParticleProxyType> iterator;

        inline size_t size() const { return static_cast<size_t>(particles.size()); }

        ParticleArrayAoS(ParticleTypes type = Electron)
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
            return ParticleProxyType(particles[idx]);
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

        inline void deleteParticle(iterator& idx)
        {
            deleteParticle(idx.getIdx());
            idx--;
        }

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

        inline iterator begin() { return iterator(this, 0); }
        inline iterator end() { return iterator(this, size()); }
        inline const iterator cbegin() { return begin(); }
        inline const iterator cend() { return end(); }


    public:
        ParticleTypes typeIndex;
        std::vector<ParticleType> particles;
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
        typedef iteratorPArray<typeArray, ParticleProxyType> iterator;

        static const int positionDimension = VectorDimensionHelper<PositionType>::dimension;
        static const int momentumDimension = VectorDimensionHelper<MomentumType>::dimension;

        inline int size() const { return static_cast<int>(weights.size()); }

        ParticleArraySoA(ParticleTypes type = Electron)
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
        }

        inline void deleteParticle(iterator& idx)
        {
            deleteParticle(idx.getIdx());
            idx--;
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
        }

        inline iterator begin() { return iterator(this, 0); }
        inline iterator end() { return iterator(this, size()); }
        inline const iterator cbegin() { return begin(); }
        inline const iterator cend() { return end(); }

    public:
        std::vector<typename ScalarType<PositionType>::Type> positions[positionDimension];
        std::vector<typename ScalarType<MomentumType>::Type> ps[momentumDimension];
        std::vector<WeightType> weights;
        std::vector<GammaType> gammas;
        ParticleTypes typeIndex;
    };

    template<>
    inline ParticleArraySoA<One>::ParticleProxyType ParticleArraySoA<One>::operator[](int idx)
    {
        PositionTypeProxy posProxy(positions[0][idx]);
        MomentumTypeProxy momProxy(ps[0][idx], ps[1][idx], ps[2][idx]);
        WeightTypeProxy weightRef = ref(weights[idx]);
        TypeIndexTypeProxy typeRef = ref(typeIndex);
        GammaTypeProxy gammaRef = ref(gammas[idx]);

        return ParticleProxyType(posProxy, momProxy, weightRef, typeRef, gammaRef);
    }

    template<>
    inline ParticleArraySoA<Two>::ParticleProxyType ParticleArraySoA<Two>::operator[](int idx)
    {
        PositionTypeProxy posProxy(positions[0][idx], positions[1][idx]);
        MomentumTypeProxy momProxy(ps[0][idx], ps[1][idx], ps[2][idx]);
        WeightTypeProxy weightRef = ref(weights[idx]);
        TypeIndexTypeProxy typeRef = ref(typeIndex);
        GammaTypeProxy gammaRef = ref(gammas[idx]);

        return ParticleProxyType(posProxy, momProxy, weightRef, typeRef, gammaRef);
    }

    template<>
    inline ParticleArraySoA<Three>::ParticleProxyType ParticleArraySoA<Three>::operator[](int idx)
    {
        PositionTypeProxy posProxy(positions[0][idx], positions[1][idx], positions[2][idx]);
        MomentumTypeProxy momProxy(ps[0][idx], ps[1][idx], ps[2][idx]);
        WeightTypeProxy weightRef = ref(weights[idx]);
        TypeIndexTypeProxy typeRef = ref(typeIndex);
        GammaTypeProxy gammaRef = ref(gammas[idx]);

        return ParticleProxyType(posProxy, momProxy, weightRef, typeRef, gammaRef);
    }





    template<>
    template<>
    inline ParticleArrayAccessor<One>::ParticleArrayAccessor(ParticleArrayAoS<One>* org)
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
    inline ParticleArrayAccessor<Two>::ParticleArrayAccessor(ParticleArrayAoS<Two>* org)
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
    inline ParticleArrayAccessor<Three>::ParticleArrayAccessor(ParticleArrayAoS<Three>* org)
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
    inline ParticleArrayAccessor<One>::ParticleArrayAccessor(ParticleArraySoA<One>* org)
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
        typeIndex = (TypeIndexType*)&(org->typeIndex);
        weights_stride = 1;
        gammas_stride = 1;
        typeIndex_stride = 0;
        _size = org->weights.size();
    };

    template<>
    template<>
    inline ParticleArrayAccessor<Two>::ParticleArrayAccessor(ParticleArraySoA<Two>* org)
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
        typeIndex = (TypeIndexType*)&(org->typeIndex);
        weights_stride = 1;
        gammas_stride = 1;
        typeIndex_stride = 0;
        _size = org->weights.size();
    };

    template<>
    template<>
    inline ParticleArrayAccessor<Three>::ParticleArrayAccessor(ParticleArraySoA<Three>* org)
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
        typeIndex = (TypeIndexType*)&(org->typeIndex);
        weights_stride = 1;
        gammas_stride = 1;
        typeIndex_stride = 0;
        _size = org->weights.size();
    };


    enum ParticleRepresentation { ParticleRepresentation_AoS, ParticleRepresentation_SoA };

    inline std::string toString(ParticleRepresentation particleRepresentation)
    {
        std::map<ParticleRepresentation, std::string> names;
        names[ParticleRepresentation_AoS] = "AoS";
        names[ParticleRepresentation_SoA] = "SoA";
        return names[particleRepresentation];
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

} // namespace pfc