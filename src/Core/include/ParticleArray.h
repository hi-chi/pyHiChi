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

    template<typename pArray_t, typename ParticleType>
    class iteratorPArray : public std::iterator<std::random_access_iterator_tag, ParticleType, size_t>
    {
    private:
        pArray_t* pPArray;
        size_t index;
    public:

        iteratorPArray(pArray_t* _pArray) : pPArray(_pArray), index(0) {};
        iteratorPArray(pArray_t* _pArray, size_t _index): pPArray(_pArray), index(_index) {};
        iteratorPArray(const iteratorPArray &it): pPArray(it.pPArray), index(it.index){};

        size_t getIdx() { return index; };

        ParticleType operator *() { return pPArray->operator[](index); }
        const iteratorPArray &operator ++() { ++index; return *this; }
        const iteratorPArray &operator --() { --index; return *this; }
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
        iteratorPArray & operator =(const iteratorPArray & other) 
        {
            this->pParray = other.pPArray;
            this->index = other.index;
            return *this; 
        }

        bool operator ==(const iteratorPArray &other) const
        { return index == other.index; }
        bool operator !=(const iteratorPArray &other) const
        { return index != other.index; }
        bool operator  <(const iteratorPArray &other) const
        { return index < other.index; }
        bool operator  >(const iteratorPArray &other) const
        { return index > other.index; }
        bool operator  <=(const iteratorPArray &other) const
        { return index <= other.index; }
        bool operator  >=(const iteratorPArray &other) const
        { return index >= other.index; }

        iteratorPArray & operator +(const size_t &add) const {
            iteratorPArray copy(*this);
            copy.index += add;
            return copy;
        }
        iteratorPArray & operator +=(const size_t &add) {
            index += add;
            return *this;
        }
        iteratorPArray & operator -(const size_t &add) const {
            iteratorPArray copy(*this);
            copy.index -= add;
            return copy;
        }
        iteratorPArray & operator -=(const long int &add) {
            index -= add;
            return *this;
        }

        ParticleType operator [](const size_t &n) const
        { return *pPArray[index + n]; }
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
        { return operator[](this->size() - 1); }

        inline void pushBack(ConstParticleRef particle) 
        { 
            if(particle.getType() == typeIndex)
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
        inline void save(std::ostream& os)
        {
            size_t tmp = size();
            os.write((char*)&tmp, sizeof(tmp));
            os.write((char*)&typeIndex, sizeof(typeIndex));
            os.write((char*)particles.data(), sizeof(Particle<dimension>)*tmp);
        }
        inline void load(std::istream& is)
        {
            size_t tmp = 0;
            is.read((char*)&tmp, sizeof(tmp));
            is.read((char*)&typeIndex, sizeof(typeIndex));
            particles.resize(tmp);
            is.read((char*)particles.data(), sizeof(Particle<dimension>) * size());
        }

        void cutMigratedParticles(std::vector<Particle3d> v[3][3][3], const FP3& minCoord, const FP3& maxCoord)
        {

        }

    private:
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

        inline void save(std::ostream& os)
        {
            size_t tmp = size();
            os.write((char*)&tmp, sizeof(tmp));
            for (int i = 0; i < positionDimension; i++)
                os.write((char*)positions[i].data(), sizeof(typename ScalarType<PositionType>::Type) * tmp);
            for (int i = 0; i < momentumDimension; i++)
                os.write((char*)ps[i].data(), sizeof(typename ScalarType<MomentumType>::Type) * tmp);
            os.write((char*)weights.data(), sizeof(WeightType) * tmp);
            os.write((char*)gammas.data(), sizeof(GammaType) * tmp);
            os.write((char*)&typeIndex, sizeof(typeIndex));
        }
        inline void load(std::istream& is)
        {
            size_t tmp = 0;
            is.read((char*)&tmp, sizeof(tmp));
            for (int i = 0; i < positionDimension; i++)
            {
                positions[i].resize(tmp);
                is.read((char*)positions[i].data(), sizeof(typename ScalarType<PositionType>::Type) * tmp);
            }
            for (int i = 0; i < momentumDimension; i++)
            {
                ps[i].resize(tmp);
                is.read((char*)ps[i].data(), sizeof(typename ScalarType<MomentumType>::Type) * tmp);
            }
            weights.resize(tmp);
            is.read((char*)weights.data(), sizeof(WeightType) * tmp);
            gammas.resize(tmp);
            is.read((char*)gammas.data(), sizeof(GammaType) * tmp);
            is.read((char*)&typeIndex, sizeof(typeIndex));
        }

        void checkParticlePos()
        {

        }

        void cutMigratedParticles(std::vector<Particle3d> v[3][3][3], const FP3& minCoord, const FP3& maxCoord)
        {
            for (int i = 0; i < this->size(); i++)
            {
                Int3 pos(1, 1, 1);
                for (int j = 0; j < positionDimension; j++)
                {
                    if (positions[j][i] > maxCoord[j]) pos[j] = 2;
                    else if (positions[j][i] < minCoord[j]) pos[j] = 0;
                }
                if (pos != Int3(1, 1, 1))
                {
                    v[pos.x][pos.y][pos.z].push_back(operator[](i));
                    deleteParticle(i);
                    i--;
                }
            }
        }

    private:
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