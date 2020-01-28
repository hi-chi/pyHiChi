#pragma once
#include "Vectors.h"

namespace pfc {

	enum FieldConfigurationType {
		NullFieldConfigurationType,
		TightFocusingFieldConfigurationType,
	};

	template<FieldConfigurationType fieldConfigurationType>
	class FieldConfiguration {
	};

	template<>
	class FieldConfiguration<FieldConfigurationType::NullFieldConfigurationType> {
	public:

		inline FP3 E(FP x, FP y, FP z) const { return FP3(); };
		inline FP3 B(FP x, FP y, FP z) const { return FP3(); };
	};

	typedef FieldConfiguration<FieldConfigurationType::NullFieldConfigurationType> NullField;
}
