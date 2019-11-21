#pragma once
// This file contains a basic implementation of the Pseudo-Spectral Rotational Time-Dependent (PSRTD) 
// field solver, for details see Section 6.2. in A. Gonoskov, PhD thesis (2013): 
// http://www.diva-portal.org/smash/get/diva2:681092/FULLTEXT02.pdf 

#pragma once

#include "pyGrid.h"

#include "Constants.h"
#include "Dimension.h"
#include "Ensemble.h"
#include "Fdtd.h"
#include "FieldGenerator.h"
#include "FieldValue.h"
#include "Handler.h"
#include "Merging.h"
#include "Particle.h"
#include "ParticleArray.h"
#include "ParticleTypes.h"
#include "Pstd.h"
#include "Psatd.h"
#include "QED_AEG.h"
#include "Vectors.h"
#include "Thinning.h"

#include "fftw3.h"

#include <iostream>
#include <math.h>

#include "draft.h"

class complex3d // local implementation of a class for 3d arrays of complexFP supplied with fftw elemets
{
public:
	complexFP *data;
	int n;
	vector<int> dim;
	fftw_complex *fdata;
	fftw_plan fp, bp;
	complex3d(int d1, int d2, int d3, int plan = 0)
	{
		dim.push_back(d1);
		dim.push_back(d2);
		dim.push_back(d3);
		n = d1 * d2 * d3;
		fdata = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);
		data = (complexFP*)fdata;

		if (plan == 0) fp = fftw_plan_dft_3d(d3, d2, d1, fdata, fdata, FFTW_FORWARD, FFTW_ESTIMATE);
		if (plan == 0) bp = fftw_plan_dft_3d(d3, d2, d1, fdata, fdata, FFTW_BACKWARD, FFTW_ESTIMATE);

		if (plan == 1) fp = fftw_plan_dft_3d(d3, d2, d1, fdata, fdata, FFTW_FORWARD, FFTW_MEASURE);
		if (plan == 1) bp = fftw_plan_dft_3d(d3, d2, d1, fdata, fdata, FFTW_BACKWARD, FFTW_MEASURE);

		if (plan == 2) fp = fftw_plan_dft_3d(d3, d2, d1, fdata, fdata, FFTW_FORWARD, FFTW_PATIENT);
		if (plan == 2) bp = fftw_plan_dft_3d(d3, d2, d1, fdata, fdata, FFTW_BACKWARD, FFTW_PATIENT);
	}
	inline ~complex3d()
	{
		fftw_destroy_plan(fp);
		fftw_destroy_plan(bp);
		fftw_free(fdata);
	}
	complexFP& operator()(int i0, int i1, int i2)
	{
		return data[i0 + (i1 + i2 * dim[1]) * dim[0]];
	}
	void FFT(int direction)
	{
		if (direction == 1) fftw_execute(fp);
		if (direction == -1)
		{
			fftw_execute(bp);
			double inv_n = 1 / double(n);
#pragma omp parallel for 
			for (int i = 0; i < n; i++)
			{
				fdata[i][0] *= inv_n;
				fdata[i][1] *= inv_n;
			}
		}
	}
};

class PSRTD;

struct PSRTD_iterator : public fieldIterator
{
	PSRTD* _field;
	int x, y, z;
	PSRTD_iterator(PSRTD* _field_) : _field(_field_)
	{}
	void begin()
	{
		x = -1;
		y = 0;
		z = 0;
	}
	bool next(FP3& position);
	void setField(const FP3 E, const FP3 B);
};

class PSRTD : public fieldDraft
{
public:
	Int3 N; // size of the field matrix
	FP3 Min, Max; // physical limits of the grid
	complex3d Fx, Fy, Fz; // data arrays for the field F = E + i*B

	complex3d* qF[3];

	PSRTD(FP3 minCoords, FP3 maxCoords, Int3 gridSize, int plan = 0) : Min(minCoords), Max(maxCoords), N(gridSize),
		Fx(N.x, N.y, N.z, plan), Fy(N.x, N.y, N.z, plan), Fz(N.x, N.y, N.z, plan)
	{
		qF[0] = &Fx; qF[1] = &Fy; qF[2] = &Fz;
		std::cout << "PSRTD[" << Fx.dim[0] << ", " << Fx.dim[1] << ", " << Fx.dim[2] << "] allocated, fftw plans(" << plan << ") createed." << std::endl;
	}
	FP3 getPosition(int x, int y, int z)
	{
		FP3 coords;
		coords.x = Min.x + (Max.x - Min.x) * (x + 0.5) / double(N.x);
		coords.y = Min.y + (Max.y - Min.y) * (y + 0.5) / double(N.y);
		coords.z = Min.z + (Max.z - Min.z) * (z + 0.5) / double(N.z);
		return coords;
	}
	void setAnalyticalValues(int64_t funcEx, int64_t funcEy, int64_t funcEz, int64_t funcBx, int64_t funcBy, int64_t funcBz)
	{
		std::cout << "PSRTD: setting analytical values..." << std::endl;
		FP(*ex)(FP, FP, FP) = (FP(*)(FP, FP, FP))funcEx; // Does this take acceptibally short time?
		FP(*ey)(FP, FP, FP) = (FP(*)(FP, FP, FP))funcEy;
		FP(*ez)(FP, FP, FP) = (FP(*)(FP, FP, FP))funcEz;
		FP(*bx)(FP, FP, FP) = (FP(*)(FP, FP, FP))funcBx;
		FP(*by)(FP, FP, FP) = (FP(*)(FP, FP, FP))funcBy;
		FP(*bz)(FP, FP, FP) = (FP(*)(FP, FP, FP))funcBz;

		FP3 coords, E, B;
		for (int z = 0; z < N.z; z++)
		{
			for (int y = 0; y < N.y; y++)
				for (int x = 0; x < N.x; x++)
				{
					coords = getPosition(x, y, z);

					E.x = ex(coords.x, coords.y, coords.z);
					E.y = ey(coords.x, coords.y, coords.z);
					E.z = ez(coords.x, coords.y, coords.z);
					B.x = bx(coords.x, coords.y, coords.z);
					B.y = by(coords.x, coords.y, coords.z);
					B.z = bz(coords.x, coords.y, coords.z);

					Fx(x, y, z) = complexFP(E.x, B.x);
					Fy(x, y, z) = complexFP(E.y, B.y);
					Fz(x, y, z) = complexFP(E.z, B.z);
				}
			cout << "#";
		}
		cout << " completed." << endl;
	}
	void advance() // advances the field for the timeStep
	{
		//		Fx.FFT(1);
		//		Fy.FFT(1);
		//		Fz.FFT(1);
#pragma omp parallel for 
		for (int i = 0; i < 3; i++)qF[i]->FFT(1);

#pragma omp parallel for 
		for (int z = 0; z < N.z; z++)
			for (int y = 0; y < N.y; y++)
				for (int x = 0; x < N.x; x++)
				{
					FP P[3][3];
					complexFP Ft_[3];

					FP kx, ky, kz, K, cosa, sina;

					kx = (x <= N.x >> 1) ? 2 * pi*x / (Max.x - Min.x) : -2 * pi*(N.x - x) / (Max.x - Min.x);
					ky = (y <= N.y >> 1) ? 2 * pi*y / (Max.y - Min.y) : -2 * pi*(N.y - y) / (Max.y - Min.y);
					kz = (z <= N.z >> 1) ? 2 * pi*z / (Max.z - Min.z) : -2 * pi*(N.z - z) / (Max.z - Min.z);

					K = sqrt(kx * kx + ky * ky + kz * kz);

					cosa = cos(_timeStep*K*lightVelocity);
					sina = sin(_timeStep*K*lightVelocity);

					if ((x != 0) || (y != 0) || (z != 0)) { kx /= K; ky /= K; kz /= K; }

					P[0][0] = (1 - kx * kx)*cosa;	    P[0][1] = -kz * sina - kx * ky*cosa;	P[0][2] = ky * sina - kx * kz*cosa;
					P[1][0] = kz * sina - kx * ky*cosa; P[1][1] = (1 - ky * ky)*cosa;	    P[1][2] = -kx * sina - ky * kz*cosa;
					P[2][0] = -ky * sina - kx * kz*cosa;	P[2][1] = kx * sina - ky * kz*cosa; P[2][2] = (1 - kz * kz)*cosa;

					Ft_[0] = Fx(x, y, z);
					Ft_[1] = Fy(x, y, z);
					Ft_[2] = Fz(x, y, z);

					Fx(x, y, z) = P[0][0] * Ft_[0] + P[0][1] * Ft_[1] + P[0][2] * Ft_[2];
					Fy(x, y, z) = P[1][0] * Ft_[0] + P[1][1] * Ft_[1] + P[1][2] * Ft_[2];
					Fz(x, y, z) = P[2][0] * Ft_[0] + P[2][1] * Ft_[1] + P[2][2] * Ft_[2];
				}

#pragma omp parallel for 
		for (int i = 0; i < 3; i++)qF[i]->FFT(-1);

		//		Fx.FFT(-1);
		//		Fy.FFT(-1);
		//		Fz.FFT(-1);
	}
	FP3 getE(int x, int y, int z)
	{
		x = (x + N.x) % N.x;
		y = (y + N.y) % N.y;
		z = (z + N.z) % N.z;
		return FP3(Fx(x, y, z).real, Fy(x, y, z).real, Fz(x, y, z).real);
	}
	FP3 getB(int x, int y, int z)
	{
		x = (x + N.x) % N.x;
		y = (y + N.y) % N.y;
		z = (z + N.z) % N.z;
		return FP3(Fx(x, y, z).imag, Fy(x, y, z).imag, Fz(x, y, z).imag);
	}
	void get(FP3 const coords, FP3 &E, FP3 &B) // provides the field at position r; gives 0 outside ranges // this is just a simple, slow implementation
	{
		bool withinLimits = false;
		if ((coords.x >= Min.x) && (coords.x <= Max.x))
			if ((coords.y >= Min.y) && (coords.y <= Max.y))
				if ((coords.z >= Min.z) && (coords.z <= Max.z))
				{
					withinLimits = true;
					FP cx = N.x*(coords.x - Min.x) / (Max.x - Min.x) - 0.5; // can be precalculated
					int ix = int(cx + 1) - 1; cx -= ix;
					FP cy = N.y*(coords.y - Min.y) / (Max.y - Min.y) - 0.5; // can be precalculated
					int iy = int(cy + 1) - 1; cy -= iy;
					FP cz = N.z*(coords.z - Min.z) / (Max.z - Min.z) - 0.5; // can be precalculated
					int iz = int(cz + 1) - 1; cz -= iz;

					E = (1 - cx)*(1 - cy)*(1 - cz)*getE(ix, iy, iz) + (cx)*(1 - cy)*(1 - cz)*getE(ix + 1, iy, iz) +
						(1 - cx)*(cy)*(1 - cz)*getE(ix, iy + 1, iz) + (cx)*(cy)*(1 - cz)*getE(ix + 1, iy + 1, iz) +
						(1 - cx)*(1 - cy)*(cz)*getE(ix, iy, iz + 1) + (cx)*(1 - cy)*(cz)*getE(ix + 1, iy, iz + 1) +
						(1 - cx)*(cy)*(cz)*getE(ix, iy + 1, iz + 1) + (cx)*(cy)*(cz)*getE(ix + 1, iy + 1, iz + 1);

					B = (1 - cx)*(1 - cy)*(1 - cz)*getB(ix, iy, iz) + (cx)*(1 - cy)*(1 - cz)*getB(ix + 1, iy, iz) +
						(1 - cx)*(cy)*(1 - cz)*getB(ix, iy + 1, iz) + (cx)*(cy)*(1 - cz)*getB(ix + 1, iy + 1, iz) +
						(1 - cx)*(1 - cy)*(cz)*getB(ix, iy, iz + 1) + (cx)*(1 - cy)*(cz)*getB(ix + 1, iy, iz + 1) +
						(1 - cx)*(cy)*(cz)*getB(ix, iy + 1, iz + 1) + (cx)*(cy)*(cz)*getB(ix + 1, iy + 1, iz + 1);
				}
		if (!withinLimits)
		{
			E = FP3(0, 0, 0);
			B = FP3(0, 0, 0);
		}
	}
	pyField init(FP time, FP timeStep)
	{
		// here we can perform all necessary preparation and validations before issuing the 'field' for further use in HiChi
		return internalInit(time, timeStep);
	}
	fieldIterator* getIterator()
	{
		return new PSRTD_iterator(this);
	}
};

bool PSRTD_iterator::next(FP3& position)
{
	x++;
	if (x >= _field->N.x) { x = 0; y++; }
	if (y >= _field->N.y) { x = 0; y = 0; z++; }
	if (z >= _field->N.z) { x = 0; y = 0; z = 0; return false; }
	position = _field->getPosition(x, y, z);
	return true;
};
void PSRTD_iterator::setField(const FP3 E, const FP3 B)
{
	_field->Fx(x, y, z) = complexFP(E.x, B.x);
	_field->Fy(x, y, z) = complexFP(E.y, B.y);
	_field->Fz(x, y, z) = complexFP(E.z, B.z);
};


class PSRTD_test
{
public:
	PSRTD_test(int x, int y = 2)
	{
		std::cout << "hi from PSRTD" << x << ":" << y << std::endl;
	}
};