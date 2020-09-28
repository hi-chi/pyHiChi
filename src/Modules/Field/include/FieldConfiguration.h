#pragma once

#include "macros.h"

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

        forceinline FP3 getE(FP x, FP y, FP z) const { return FP3(); }
        forceinline FP3 getB(FP x, FP y, FP z) const { return FP3(); }

        // if E and B are computed in a similar way this function will work quicklier
        // than getE and getB separately
        forceinline void getEB(FP x, FP y, FP z, FP3* E, FP3* B) const {
            *E = FP3();
            *B = FP3();
        }
    };

    typedef FieldConfiguration<FieldConfigurationType::NullFieldConfigurationType> NullField;
    
    template<>
    class FieldConfiguration<FieldConfigurationType::TightFocusingFieldConfigurationType> {

    public:

        FP f_number;
        FP R0;
        FP wavelength;
        FP pulselength;
        FP totalPower;
        FP openingAngle;
        FP edgeSmoothingAngle;
        FP3 polarisation;
        FP exclusionRadius;
        FP ampl;

        FieldConfiguration(FP f_number, FP R0, FP wavelength, FP pulselength, FP totalPower,
            FP edgeSmoothingAngle, FP3 polarisation = FP3(0.0, 1.0, 0.0), FP exclusionRadius = 1e-5) :
            f_number(f_number), R0(R0), wavelength(wavelength), pulselength(pulselength),
            totalPower(totalPower), edgeSmoothingAngle(edgeSmoothingAngle), polarisation(polarisation),
            exclusionRadius(exclusionRadius), openingAngle(atan(1.0 / (2.0*f_number))),
            ampl(sqrt(totalPower * 4.0 / (1.0 - cos(openingAngle)) / constants::c))
        {}


        forceinline FP sign(FP x) const {
            return (x >= 0 ? 1.0 : -1.0);
        }

        forceinline FP block(FP x, FP xmin, FP xmax) const {
            return (sign(x - xmin) + sign(xmax - x))*0.5;
        }

        forceinline FP longitudinalFieldVariation(FP x_ct) const {
            FP cSin = sin(2.0 * constants::pi * x_ct / wavelength);
            FP cCos = cos(constants::pi * x_ct / pulselength);
            FP cBlock = block(x_ct, -0.5*pulselength, 0.5*pulselength);
            return cSin * cCos * cCos * cBlock;
        }

        forceinline FP transverseShape(FP angle) const {
            FP a = openingAngle - edgeSmoothingAngle * 0.5, b = openingAngle + edgeSmoothingAngle * 0.5;
            FP cMainBlock = block(angle, -1.0, a);
            if (edgeSmoothingAngle == 0.0)
                return cMainBlock;
            FP cEdgeSmAngleBlock = block(angle, a, b);
            FP cCos = cos(0.5 * constants::pi * (angle - a) / edgeSmoothingAngle);
            return cMainBlock + cCos * cCos * cEdgeSmAngleBlock;
        }
        
        forceinline FP mask(FP x, FP y, FP z) const {
            FP R = sqrt(x * x + y * y + z * z);
            if (R > exclusionRadius) {
                FP angle = asin(sqrt(y * y + z * z) / R);
                return (ampl / R) * longitudinalFieldVariation(R - R0) * transverseShape(angle) * (x < 0);
            }
            return 0.0;
        }


        forceinline FP3 getE(FP x, FP y, FP z) const {
            FP3 r(x, y, z);
            FP3 s1 = cross(this->polarisation, r);
            FP3 s0 = cross(r, s1);
            return s0.norm() != 0.0 ? (mask(x, y, z) / s0.norm()) * s0 : FP3(0.0, 0.0, 0.0);
        };

        forceinline FP3 getB(FP x, FP y, FP z) const {
            FP3 r(x, y, z);
            FP3 s1 = cross(this->polarisation, r);
            return s1.norm() != 0.0 ? (mask(x, y, z) / s1.norm()) * s1 : FP3(0.0, 0.0, 0.0);
        };

        forceinline void getEB(FP x, FP y, FP z, FP3* E, FP3* B) const {
            FP3 r(x, y, z);
            FP3 s1 = cross(this->polarisation, r);
            FP3 s0 = cross(r, s1);
            FP m = mask(x, y, z);
            *E = s0.norm() != 0.0 ? (m / s0.norm()) * s0 : FP3(0.0, 0.0, 0.0);
            *B = s1.norm() != 0.0 ? (m / s1.norm()) * s1 : FP3(0.0, 0.0, 0.0);
        };

    };

    typedef FieldConfiguration<FieldConfigurationType::TightFocusingFieldConfigurationType> TightFocusingField;
}
