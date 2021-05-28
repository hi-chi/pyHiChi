#pragma once
#include "ParticleArray.h"
#include "ParticleTypes.h"

#include <map>

namespace pfc {

    // Collection of particles with array-like semantics,
    // representation as array of structures
    template<class pArray>
    class Ensemble {
    public:

        typedef typename pArray::PositionType PositionType;
        typedef typename pArray::MomentumType MomentumType;
        typedef typename pArray::GammaType GammaType;
        typedef typename pArray::WeightType WeightType;
        typedef typename pArray::TypeIndexType TypeIndexType;
        typedef typename pArray::ParticleType ParticleType;
        typedef typename pArray::ParticleRef ParticleRef;
        typedef typename pArray::ConstParticleRef ConstParticleRef;
        typedef typename pArray::ParticleProxyType ParticleProxyType;
        typedef pArray ParticleArray;

        inline int size() const 
        { 
            size_t size = 0;
            for (auto it = pArrays.begin(); it != pArrays.end(); it++)
            {
                size += it->second.size();
            }
            return static_cast<int>(size); 
        }

        Ensemble()
        {
            for (int t = 0; t < sizeParticleTypes; t++)
            {
                pArray newPArray(static_cast<ParticleTypes>(t));
                pArrays.insert(pair<string, pArray>(particleNames[t], newPArray));
            }
        }

        Ensemble(const Ensemble<pArray> & other)
        {
            pArrays = other.pArrays;
        }
        
        inline pArray& operator[](string name)
        {
            return pArrays.find(name)->second;
        }

        inline pArray& operator[](int ind)
        {
            return pArrays.find(particleNames[ind])->second;
        }
        
        inline void addParticle(ConstParticleRef particle)
        {
            pArrays[particleNames[particle.getType()]].pushBack(particle);
        }

        inline void clear() 
        {
            pArrays.clear();
        }

        inline void save(std::ostream& os)
        {
            size_t tmp = pArrays.size();
            os.write((char*)&tmp, sizeof(tmp));
            for (auto& el: pArrays)
            {
                size_t s_size = el.first.size();
                os.write((char*)&s_size, sizeof(s_size));
                os.write(el.first.data(), s_size);
                el.second.save(os);
            }
        }
        inline void load(std::istream& is)
        {
            size_t tmp = pArrays.size();
            is.read((char*)&tmp, sizeof(tmp));
            for (size_t i = 0; i < tmp; i++)
            {
                size_t s_size;
                is.read((char*)&s_size, sizeof(s_size));
                string s(s_size, '1');
                is.read(&s[0], s_size);
                pArray tmp_array;
                tmp_array.load(is);
                pArrays[s] = tmp_array;
            }
        }

    private:
        std::map<string, pArray> pArrays;
    };

    typedef Ensemble<ParticleArray3d> Ensemble3d;

    class NoParticleArray {};
    template<>
    class Ensemble<NoParticleArray>
    {
    public:
        inline void save(std::ostream& os) {}
        inline void load(std::istream& is) {}
    };
} // namespace pfc