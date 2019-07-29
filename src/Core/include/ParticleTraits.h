#pragma once

namespace pfc {


// This class defines types of particle data
// By default it relies on nested types in a particle class
// In case a particle class does not have these, ParticleTraits must be specialized
template<typename ParticleType>
struct ParticleTraits {
public:
    typedef typename ParticleType::PositionType PositionType;
    typedef typename ParticleType::MomentumType MomentumType;
    typedef typename ParticleType::GammaType GammaType;
    typedef typename ParticleType::WeightType WeightType;
    typedef typename ParticleType::TypeIndexType TypeIndexType;

	typedef typename ParticleType::PositionTypeProxy PositionTypeProxy;
	typedef typename ParticleType::MomentumTypeProxy MomentumTypeProxy;
	typedef typename ParticleType::GammaTypeProxy GammaTypeProxy;
	typedef typename ParticleType::WeightTypeProxy WeightTypeProxy;
	typedef typename ParticleType::TypeIndexTypeProxy TypeIndexTypeProxy;
};

} // namespace pfc
