#pragma once
#include "PicParticle.h"
#include <vector>

// not-thread-safe container for particles
class PicParticleContainer
{
	std::vector<PicParticle*> particles;

public:

	PicParticleContainer();
	~PicParticleContainer();

	void Push(const PicParticle& particle);
	void Pop(int index);  // swap and resize, changes vector

	// thread-safe method
	// deletes particle but does not change vector
	void DeleteContent(int index);

	PicParticle* operator[](int index);
	size_t size();

	PicParticleContainer& operator=(const PicParticleContainer&) = delete;
	PicParticleContainer(const PicParticleContainer&) = delete;

};

