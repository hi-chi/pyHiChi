#pragma once

namespace pfc {

// Mathematical and physical constants in CGS
namespace constants
{
    const double pi = 3.14159265358;
    const double c = 29979245800.0;
    const double lightVelocity = c;
    const double electronCharge = -4.80320427e-10;
    const double electronMass = 9.10938215e-28;
    const double protonMass = 1.672622964e-24;
    const double planck = 1.0545716818e-27;
    const double eV = 1.60217656535e-12;
    const double meV = 1e6 * eV;
} // namespace pfc::constants


// Provides same constants as in the constants namespace as static methods, but casted to the given type
template<typename Real>
struct Constants {
    static Real pi() { return static_cast<Real>(constants::pi); }
    static Real c() { return static_cast<Real>(constants::c); }
    static Real lightVelocity() { return static_cast<Real>(constants::lightVelocity); }
    static Real electronCharge() { return static_cast<Real>(constants::electronCharge); }
    static Real electronMass() { return static_cast<Real>(constants::electronMass); }
    static Real protonMass() { return static_cast<Real>(constants::protonMass); }
    static Real planck() { return static_cast<Real>(constants::planck); }
    static Real eV() { return static_cast<Real>(constants::eV); }
    static Real meV() { return static_cast<Real>(constants::meV); }
};

} // namespace pfc
