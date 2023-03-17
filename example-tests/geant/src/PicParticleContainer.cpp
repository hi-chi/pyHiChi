#include "PicParticleContainer.h"
#include "MyException.h"


PicParticleContainer::PicParticleContainer() {}

PicParticleContainer::~PicParticleContainer()
{
	for (int i = 0; i < particles.size(); i++)
		if (particles[i]) DeleteContent(i);
}

void PicParticleContainer::Push(const PicParticle& particle)
{
	particles.push_back(new PicParticle(particle));
}

void PicParticleContainer::Pop(int index)
{
	if (index < 0 || index >= particles.size())
		throw MyException("Pop from PicParticleContainer is incorrect");
	DeleteContent(index);
	std::swap(particles[index], particles[particles.size() - 1]);
	particles.pop_back();
}

void PicParticleContainer::DeleteContent(int index)
{
	delete particles[index];
	particles[index] = nullptr;
}

PicParticle* PicParticleContainer::operator[](int index)
{
	return particles[index];
}

size_t PicParticleContainer::size()
{
	return particles.size();
}
