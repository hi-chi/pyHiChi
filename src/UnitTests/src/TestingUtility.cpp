#include "TestingUtility.h"

#include "Vectors.h"
#include "VectorsProxy.h"
#include "Particle.h"

#include <cmath>
#include <memory>

using namespace pfc;

namespace pfc {

	namespace ParticleInfo {
		std::vector<ParticleType> typesVector;
		const ParticleType* types;
		short numTypes;
	} // namespace ParticleInfo

} // namespace pica

void BaseFixture::SetUp() {
	srand(1);
	maxAbsoluteError = (FP)1e-4;
	maxRelativeError = (FP)1e-4;
}

void BaseFixture::TearDown() {
}

// Return whether two FP3s have coords differ by not larger than eps each.
bool BaseFixture::nearFP3(const FP3 & a, const FP3 & b, const FP eps) {
	return (fabs(a.x - b.x) <= eps) && (fabs(a.y - b.y) <= eps) &&
		(fabs(a.z - b.z) <= eps);
}

// Get uniformly distributed in [a, b) pseudo-random number.
FP BaseFixture::urand(FP a, FP b) const {
	return a + (b - a) * ((FP)rand()) / RAND_MAX;
}

int BaseFixture::urandInt(int a, int b) {
	return a + rand() % (b - a + 1);
}

// Get distributed in [a, b) pseudo-random vector.
FP3 BaseFixture::urandFP3(FP3 a, FP3 b) {
	FP3 result;
	result.x = urand(a.x, b.x);
	result.y = urand(a.y, b.y);
	result.z = urand(a.z, b.z);
	return result;
}

Int3 BaseFixture::urandInt3(Int3 a, Int3 b) {
	Int3 result;
	result.x = urandInt(a.x, b.x);
	result.y = urandInt(a.y, b.y);
	result.z = urandInt(a.z, b.z);
	return result;
}

// Get _n_ random vectors between _minValue_ and _maxValue_.
std::vector<FP3> BaseFixture::randomVectors(int n, FP3 & minValue, FP3 & maxValue) {
	std::vector<FP3> result(n);
	for (int i = 0; i < n; ++i)
		result[i] = urandFP3(minValue, maxValue);
	return result;
}