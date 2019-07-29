#pragma once

namespace pfc {

// Switch to use single or double precision as type FP, double by default.
#ifndef PFC_USE_SINGLE_PRECISION
    typedef double FP;
#else
    typedef float FP;
#endif

inline FP sqr(FP x) { return x * x; }

struct complexFP 
{
	FP real;
	FP imag;

    complexFP() {
        real = 0.0;
        imag = 0.0;
    }

	complexFP(FP _real, FP _imag) {
		real = _real;
		imag = _imag;
	}

	complexFP(FP _real) {
		real = _real;
		imag = 0.0;
	}

    static complexFP createInTrig(FP module, FP arg) {
        return complexFP(module * cos(arg), module * sin(arg));
    }

    friend int operator==(const complexFP& z1, const complexFP& z2) {
        return (z1.real == z2.real) && (z1.imag == z2.imag);
    }

    FP getModule() {
        return sqrt(real * real + imag * imag);
    }

    FP getArg() {
        return atan(imag / real);
    }

    complexFP getConj() {
        return complexFP(real, -imag);
    }

    complexFP operator-() {
        return complexFP(-real, -imag);
    }

    friend complexFP operator+(const complexFP& z1, const complexFP& z2) {
        return complexFP(z1.real + z2.real, z1.imag + z2.imag);
    }

    friend complexFP operator-(const complexFP& z1, const complexFP& z2) {
        return complexFP(z1.real - z2.real, z1.imag - z2.imag);
    }

    friend complexFP operator*(const complexFP& z1, const complexFP& z2) {
        return complexFP(z1.real * z2.real - z1.imag * z2.imag, z1.real * z2.imag + z2.real * z1.imag);
    }

    complexFP& operator+=(const complexFP& z) {
        real += z.real;
        imag += z.imag;
        return *this;
    }

    complexFP& operator-=(const complexFP& z) {
        real -= z.real;
        imag -= z.imag;
        return *this;
    }

    complexFP& operator*=(const complexFP& z) {
        *this = (*this)*z;
        return *this;
    }

    static inline complexFP i() {
        return complexFP(0, 1);
    }
};

typedef complexFP complex;

} // namespace pfc
