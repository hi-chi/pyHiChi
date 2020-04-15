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

        FP3 getE(FP x, FP y, FP z) const { return FP3(); }
        FP3 getB(FP x, FP y, FP z) const { return FP3(); }

        // if E and B are computed in a similar way this function will work quicklier
        // than getE and getB separately
        void getEB(FP x, FP y, FP z, FP3* E, FP3* B) const {
            *E = FP3();
            *B = FP3();
        }
    };

    typedef FieldConfiguration<FieldConfigurationType::NullFieldConfigurationType> NullField;
    
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
        FP ampl;

        FieldConfiguration(FP F_number, FP R0, FP wavelength, FP pulselength, FP phase, FP totalPower,
            FP edgeSmoothingAngle, FP3 polarisation = FP3(0.0, 1.0, 0.0), FP exclusionRadius = 1e-5) :
            F_number(F_number), R0(R0), wavelength(wavelength), pulselength(pulselength), phase(phase),
            totalPower(totalPower), edgeSmoothingAngle(edgeSmoothingAngle), polarisation(polarisation),
            exclusionRadius(exclusionRadius), openingAngle(atan(1.0 / (2.0*F_number))),
            ampl(sqrt(totalPower * 4.0 / (constants::c*(1.0 - cos(openingAngle)))))
        {}


        FP sign(FP x) const {
            return (x >= 0 ? 1.0 : -1.0);
        }

        FP block(FP x, FP xmin, FP xmax) const {
            return (sign(x - xmin) + sign(xmax - x))*0.5;
        }

        FP longitudinalFieldVariation(FP x_ct) const {
            FP cSin = sin(2.0 * constants::pi * x_ct / wavelength + phase);
            FP cCos = cos(constants::pi * x_ct / pulselength);
            FP cBlock = block(x_ct, -0.5*pulselength, 0.5*pulselength);
            return cSin * cCos * cCos * cBlock;
        }

        FP transverseShape(FP angle) const {
            FP cMainBlock = block(angle, -1.0, openingAngle - edgeSmoothingAngle * 0.5);
            FP cEdgeSmAngleBlock = block(angle, openingAngle - edgeSmoothingAngle * 0.5, openingAngle + edgeSmoothingAngle * 0.5);
            FP cCos = cos(0.5*constants::pi*(angle - openingAngle + edgeSmoothingAngle * 0.5) / edgeSmoothingAngle);
            return (edgeSmoothingAngle == 0.0) ? cMainBlock : (cMainBlock + cCos * cCos * cEdgeSmAngleBlock);
        }
        
        FP mask(FP x, FP y, FP z) const {
            FP R = sqrt(x * x + y * y + z * z);
            if (R > exclusionRadius) {
                FP angle = asin(sqrt(y * y + z * z) / R);
                return (ampl / R) * longitudinalFieldVariation(R - R0) * transverseShape(angle) * (x < 0);
            }
            return 0.0;
        }


        FP3 getE(FP x, FP y, FP z) const {
            FP3 r(x, y, z);
            FP3 s1 = cross(this->polarisation, r);
            FP3 s0 = cross(r, s1);
            s0.normalize();
            return mask(x, y, z)*s0;
        };

        FP3 getB(FP x, FP y, FP z) const {
            FP3 r(x, y, z);
            FP3 s1 = cross(this->polarisation, r);
            s1.normalize();
            return mask(x, y, z)*s1;
        };

        void getEB(FP x, FP y, FP z, FP3* E, FP3* B) const {
            FP3 r(x, y, z);
            FP3 s1 = cross(this->polarisation, r);
            FP3 s0 = cross(r, s1);
            s0.normalize();
            s1.normalize();
            FP m = mask(x, y, z);
            *E = m*s0;
            *B = m*s1;
        };

    };

    typedef FieldConfiguration<FieldConfigurationType::TightFocusingFieldConfigurationType> TightFocusingField;
}
