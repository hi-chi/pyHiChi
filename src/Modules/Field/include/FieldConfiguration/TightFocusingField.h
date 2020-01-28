#pragma once
#include "Vectors.h"
#include "Constants.h"
#include "FieldConfiguration.h"
#include <cmath>
#include <functional>

namespace pfc {

	template<>
	class FieldConfiguration<FieldConfigurationType::TightFocusingFieldConfigurationType> {

	public:

		FP F_number;
		FP R0;
		FP wavelength;
		FP pulselength;
		FP phase;
		FP totalPower;
		FP openingAngle;
		FP edgeSmoothingAngle;
		FP3 polarisation;
		FP exclusionRadius;

		FieldConfiguration(FP F_number, FP R0, FP wavelength, FP pulselength, FP phase, FP totalPower,
			FP edgeSmoothingAngle, FP3 polarisation = FP3(0.0, 1.0, 0.0), FP exclusionRadius = 1e-5) :
			F_number(F_number), R0(R0), wavelength(wavelength), pulselength(pulselength), phase(phase),
			totalPower(totalPower), edgeSmoothingAngle(edgeSmoothingAngle), polarisation(polarisation),
			exclusionRadius(exclusionRadius), openingAngle(atan(1.0 / (2.0*F_number))) {}


		FP sign(FP x) const {
			return (x >= 0 ? 1.0 : -1.0);
		}

		FP block(FP x, FP xmin, FP xmax) const {
			return (sign(x - xmin) + sign(xmax - x))*0.5;
		}

		FP longitudinalFieldVariation(FP x_ct) const {
			return sin(2.0 * constants::pi*x_ct / wavelength + phase) *
				pow(cos(constants::pi*x_ct / pulselength), 2.0) *
				block(x_ct, -0.5*pulselength, 0.5*pulselength);
		}

		FP transverseShape(FP angle) const {
			FP tmp = block(angle, -1.0, openingAngle - edgeSmoothingAngle * 0.5);
			return (edgeSmoothingAngle == 0.0) ? tmp : (tmp +
				pow(cos(0.5*constants::pi*(angle - openingAngle + edgeSmoothingAngle * 0.5) / edgeSmoothingAngle), 2.0) *
				block(angle, openingAngle - edgeSmoothingAngle * 0.5, openingAngle + edgeSmoothingAngle * 0.5));
		}

		FP ampl() const {
			return sqrt(totalPower * 4.0 / (constants::c*(1.0 - cos(openingAngle))));
		}
		
		FP mask(FP x, FP y, FP z) const {
			FP R = sqrt(x * x + y * y + z * z);
			if (R > exclusionRadius) {
				FP angle = asin(sqrt(y * y + z * z) / R);
				return (ampl() / R)*longitudinalFieldVariation(R - R0) *transverseShape(angle)*(x < 0);
			}
			return 0.0;
		}


		FP3 E(FP x, FP y, FP z) const {
			FP3 r(x, y, z);
			FP3 s1 = cross(this->polarisation, r);
			s1.normalize();
			FP3 s0 = cross(r, s1);
			s0.normalize();
			return mask(x, y, z)*s0;
		};

		FP3 B(FP x, FP y, FP z) const {
			FP3 r(x, y, z);
			FP3 s1 = cross(this->polarisation, r);
			s1.normalize();
			return mask(x, y, z)*s1;
		};

	};

	typedef FieldConfiguration<FieldConfigurationType::TightFocusingFieldConfigurationType> TightFocusingField;

}