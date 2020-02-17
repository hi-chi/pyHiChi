#pragma once
#include "Pusher.h"
#include "FieldValue.h"

#include <array>
#include <type_traits>
#include <vector>


template<class T_Handler>
struct IsVectorized : std::false_type {};

template<>
struct IsVectorized<VectorizedParticlePusher> : std::true_type {};
template<>
struct IsVectorized<VectorizedBorisPusher> : std::true_type {};

template< class T_Handler, class T_Particle>
void handleChunk(T_Handler& handler,T_Particle& particle, ValueField field, FP timeStep)
{
	handler(particle, field, timeStep);
}

template< class T_Handler, class T_Particle, size_t T_size>
typename std::enable_if<IsVectorized< T_Handler>::value>::type
handleChunk(T_Handler& handler, Chunk<T_Particle, T_size>& chunk, Chunk<ValueField, T_size> fields, FP timeStep)
{
	handler(chunk, fields, timeStep);
}

template< class T_Handler, class T_Particle, size_t T_size>
typename std::enable_if<!IsVectorized<T_Handler>::value>::type
handleChunk(T_Handler& handler, Chunk<T_Particle, T_size>& chunk, Chunk<ValueField, T_size> fields, FP timeStep)
{
	for (int i = 0; i < T_size; i++)
		handler(chunk[i], fields[i], timeStep);
}

class BorisPusher
{
public:
	template<class T_Particle>
	void operator()(T_Particle* particle, ValueField field, FP timeStep)
	{
		ScalarBorisPusher scalarPusher;
		handleChunk(scalarPusher, *particle, field, timeStep);
	}

	template<class T_ParticleArray>
	void operator()(T_ParticleArray* particleArray, std::vector<ValueField> fields, FP timeStep)
	{
		typedef typename T_ParticleArray::ParticleProxyType ParticleProxyType;
		ScalarBorisPusher scalarPusher;
		for (int i = 0; i < particleArray->size(); i++)
		{
			ParticleProxyType particle = (*particleArray)[i];
			ValueField field = fields[i];
			handleChunk(scalarPusher, particle, field, timeStep);
		}
	}
};

