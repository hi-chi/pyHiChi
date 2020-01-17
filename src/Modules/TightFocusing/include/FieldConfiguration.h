#pragma once
#include "Vectors.h"

namespace pfc {

	class FieldConfiguration {

	public:

		virtual FP3 E(FP x, FP y, FP z) const { return FP3(); };
		virtual FP3 B(FP x, FP y, FP z) const { return FP3(); };
		virtual FP3 J(FP x, FP y, FP z) const { return FP3(); };
	};
}
