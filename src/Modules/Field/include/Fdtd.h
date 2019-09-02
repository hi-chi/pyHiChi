#pragma once
#include "Constants.h"
#include "FieldSolver.h"
#include "Grid.h"
#include "PmlFdtd.h"
#include "Vectors.h"

#include <algorithm>

namespace pfc {
	
	class FDTD : public RealFieldSolver<YeeGridType>
	{

	public:
		FDTD(YeeGrid* grid);

		void updateFields();

		void setPML(int sizePMLx, int sizePMLy, int sizePMLz);
		void setFieldGenerator(FieldGeneratorYee * _generator);

		void updateHalfB();
		void updateE();

        void setTimeStep(FP dt);

	private:

		void updateHalfB3D();
		void updateHalfB2D();
		void updateHalfB1D();
		void updateE3D();
		void updateE2D();
		void updateE1D();

		FP3 anisotropyCoeff;
		void setAnisotropy(const FP frequency, int axis);

	};

	inline FDTD::FDTD(YeeGrid* grid) :
		RealFieldSolver(grid)
	{
		updateDims();
		pml.reset(new Pml<GridTypes::YeeGridType>(this, Int3(0, 0, 0)));//pml.reset(new PmlFdtd(this));;
		generator.reset(new ReflectFieldGeneratorYee(this));
		updateInternalDims();
		anisotropyCoeff = FP3(1, 1, 1);
	}

	inline void FDTD::setPML(int sizePMLx, int sizePMLy, int sizePMLz)
	{
		pml.reset(new PmlFdtd(this, Int3(sizePMLx, sizePMLy, sizePMLz)));
		updateInternalDims();
	}

    inline void FDTD::setTimeStep(FP dt)
    {
        if (grid->setTimeStep(dt)) {
            if (pml->sizePML == Int3(0, 0, 0))
                pml.reset(new Pml<GridTypes::YeeGridType>(this, Int3(0, 0, 0)));
            else pml.reset(new PmlFdtd(this, pml->sizePML));
            generator.reset(new ReflectFieldGeneratorYee(this));
        }
    }

	inline void FDTD::setFieldGenerator(FieldGeneratorYee * _generator)
	{
		generator.reset(_generator);
	}

	inline void FDTD::setAnisotropy(FP frequency, int axis)
	{
		// We introduce artificial anisotropy, through one axis.
		// For this we upgrade Maxwell equations by coefficients,
		// which computes from major signal frequency. 
		// See more in Juntunen,Tsiboukis - Reduction of Numerical Dispersion in
		// FDTD Method Through Artificial Anisotropy.

		FP3 steps = grid->steps;
		FP WP = constants::pi * 2.0 * constants::c / frequency;
		FP R = WP / steps.norm();
		const FP q = 0.99;  // q - stability coefficient, 0 <= q <= 1
		FP Amax = constants::pi / (3 * R * asin(asin(constants::pi / (R * sqrt(3.0))) / sqrt(3.0)));
		FP Q = Amax - 1;
		FP c1 = 1 - Q / 2;
		int axis0 = axis;
		int axis1 = (axis + 1) % 3;
		int axis2 = (axis + 2) % 3;
		// equivalents of the variables
		// Z1 == Zy, Zz == Zz
		// Zy,Zz - designation from article
		FP Z1 = steps[axis0] / steps[axis1];
		FP Z2 = steps[axis0] / steps[axis2];
		// equivalents of the variables
		// CoeffA == K1, CoeffB == K2, a1 == a, a2 == b
		// K1, K2, a, b - designation from article
		FP CoeffA = constants::pi / (R * sqrt(1 + 1 / (Z1 * Z1) + 1 / (Z2 * Z2)));
		FP a1 = sin(CoeffA / c1) * sin(CoeffA / c1)
			/ (Z1 * Z1 * sin(CoeffA / (c1 * Z1)) * sin(CoeffA / (c1 * Z1)));
		FP a2 = sin(CoeffA / c1) * sin(CoeffA / c1)
			/ (Z2 * Z2 * sin(CoeffA / (c1 * Z2)) * sin(CoeffA / (c1 * Z2)));
		FP CoeffB = sqrt(1 + a1 * Z1 * Z1 + a2 * Z2 * Z2);
		anisotropyCoeff[axis0] = CoeffB / (CoeffA * q * sqrt(a1 * a2))
			* asin(q * sin(CoeffA / c1) / CoeffB);
		anisotropyCoeff[axis1] = a1 * anisotropyCoeff[axis0];
		anisotropyCoeff[axis2] = a2 * anisotropyCoeff[axis0];
	}

	inline void FDTD::updateFields()
	{
		updateHalfB();
		pml->updateB();
		generator->generateB();
		updateE();
		pml->updateE();
		generator->generateE();
		updateHalfB();
	}

	// Update grid values of magnetic field in FDTD.
	inline void FDTD::updateHalfB()
	{
		if (grid->dimensionality == 3)
			updateHalfB3D();
		else if (grid->dimensionality == 2)
			updateHalfB2D();
		else if (grid->dimensionality == 1)
			updateHalfB1D();
	}

	inline void FDTD::updateHalfB3D()
	{
		updateBAreaBegin = Int3(1, 1, 1);
		updateBAreaEnd = grid->numCells - Int3(1, 1, 1);
		for (int d = 0; d < 3; ++d)
		{
			internalBAreaBegin[d] = std::max(updateBAreaBegin[d], pml->leftDims[d]);
			internalBAreaEnd[d] = std::min(updateBAreaEnd[d],
				grid->numCells[d] - pml->rightDims[d]);
		}

		const FP cdt = constants::c * grid->dt * (FP)0.5;
		const FP coeffXY = cdt / (grid->steps.x * anisotropyCoeff.y);
		const FP coeffXZ = cdt / (grid->steps.x * anisotropyCoeff.z);
		const FP coeffYX = cdt / (grid->steps.y * anisotropyCoeff.x);
		const FP coeffYZ = cdt / (grid->steps.y * anisotropyCoeff.z);
		const FP coeffZX = cdt / (grid->steps.z * anisotropyCoeff.x);
		const FP coeffZY = cdt / (grid->steps.z * anisotropyCoeff.y);

		// In central area use b(i, j, k) += c * dt * -rot(e(i, j, k)), which is:
		// b.x(i, j, k) += c * dt * ((e.y(i, j, k) - e.y(i, j, k-1)) / eps_z * dz -
		//	 (e.z(i, j, k) - e.z(i, j-1, k)) / eps_y * dy),
		// b.y(i, j, k) += c * dt * ((e.z(i, j, k) - e.z(i-1, j, k)) / eps_x * dx -
		//	 (e.x(i, j, k) - e.x(i, j, k-1)) / eps_z * dz),
		// b.z(i, j, k) += c * dt * ((e.x(i, j, k) - e.x(i, j-1, k)) / eps_y * dy -
		//	 (e.y(i, j, k) - e.y(i-1, j, k)) / eps_x * dx),
		const Int3 begin = internalBAreaBegin;
		const Int3 end = internalBAreaEnd;
#pragma omp parallel for
		for (int i = begin.x; i < end.x; i++)
			for (int j = begin.y; j < end.y; j++)
			{
#pragma simd
				for (int k = begin.z; k < end.z; k++)
				{
					grid->Bx(i, j, k) += coeffZX * (grid->Ey(i, j, k) - grid->Ey(i, j, k - 1)) -
						coeffYX * (grid->Ez(i, j, k) - grid->Ez(i, j - 1, k));
					grid->By(i, j, k) += coeffXY * (grid->Ez(i, j, k) - grid->Ez(i - 1, j, k)) -
						coeffZY * (grid->Ex(i, j, k) - grid->Ex(i, j, k - 1));
					grid->Bz(i, j, k) += coeffYZ * (grid->Ex(i, j, k) - grid->Ex(i, j - 1, k)) -
						coeffXZ * (grid->Ey(i, j, k) - grid->Ey(i - 1, j, k));
				}
			}
	}

	inline void FDTD::updateHalfB2D()
	{
		updateBAreaBegin = Int3(1, 1, 0);
		updateBAreaEnd = grid->numCells - Int3(1, 1, 0);
		for (int d = 0; d < 2; ++d)
		{
			internalBAreaBegin[d] = std::max(updateBAreaBegin[d], pml->leftDims[d]);
			internalBAreaEnd[d] = std::min(updateBAreaEnd[d],
				grid->numCells[d] - pml->rightDims[d]);
		}

		const FP cdt = constants::c * grid->dt * (FP)0.5;
		const FP coeffXY = cdt / (grid->steps.x * anisotropyCoeff.y);
		const FP coeffXZ = cdt / (grid->steps.x * anisotropyCoeff.z);
		const FP coeffYX = cdt / (grid->steps.y * anisotropyCoeff.x);
		const FP coeffYZ = cdt / (grid->steps.y * anisotropyCoeff.z);

		// In central area use b(i, j, k) += c * dt * -rot(e(i, j, k)), which is:
		// b.x(i, j, k) += c * dt * ((e.y(i, j, k) - e.y(i, j, k-1)) / eps_z * dz -
		//	 (e.z(i, j, k) - e.z(i, j-1, k)) / eps_y * dy),
		// b.y(i, j, k) += c * dt * ((e.z(i, j, k) - e.z(i-1, j, k)) / eps_x * dx -
		//	 (e.x(i, j, k) - e.x(i, j, k-1)) / eps_z * dz),
		// b.z(i, j, k) += c * dt * ((e.x(i, j, k) - e.x(i, j-1, k)) / eps_y * dy -
		//	 (e.y(i, j, k) - e.y(i-1, j, k)) / eps_x * dx),
		const Int3 begin = internalBAreaBegin;
		const Int3 end = internalBAreaEnd;
#pragma omp parallel for
		for (int i = begin.x; i < end.x; i++) {
#pragma simd
			for (int j = begin.y; j < end.y; j++)
			{
				grid->Bx(i, j, 0) += -coeffYX * (grid->Ez(i, j, 0) - grid->Ez(i, j - 1, 0));
				grid->By(i, j, 0) += coeffXY * (grid->Ez(i, j, 0) - grid->Ez(i - 1, j, 0));
				grid->Bz(i, j, 0) += coeffYZ * (grid->Ex(i, j, 0) - grid->Ex(i, j - 1, 0)) -
					coeffXZ * (grid->Ey(i, j, 0) - grid->Ey(i - 1, j, 0));
			}
		}
	}

	inline void FDTD::updateHalfB1D()
	{
		updateBAreaBegin = Int3(1, 0, 0);
		updateBAreaEnd = grid->numCells - Int3(1, 0, 0);
		for (int d = 0; d < 1; ++d)
		{
			internalBAreaBegin[d] = std::max(updateBAreaBegin[d], pml->leftDims[d]);
			internalBAreaEnd[d] = std::min(updateBAreaEnd[d],
				grid->numCells[d] - pml->rightDims[d]);
		}

		const FP cdt = constants::c * grid->dt * (FP)0.5;
		const FP coeffXY = cdt / (grid->steps.x * anisotropyCoeff.y);
		const FP coeffXZ = cdt / (grid->steps.x * anisotropyCoeff.z);

		// In central area use b(i, j, k) += c * dt * -rot(e(i, j, k)), which is:
		// b.x(i, j, k) += c * dt * ((e.y(i, j, k) - e.y(i, j, k-1)) / eps_z * dz -
		//	 (e.z(i, j, k) - e.z(i, j-1, k)) / eps_y * dy),
		// b.y(i, j, k) += c * dt * ((e.z(i, j, k) - e.z(i-1, j, k)) / eps_x * dx -
		//	 (e.x(i, j, k) - e.x(i, j, k-1)) / eps_z * dz),
		// b.z(i, j, k) += c * dt * ((e.x(i, j, k) - e.x(i, j-1, k)) / eps_y * dy -
		//	 (e.y(i, j, k) - e.y(i-1, j, k)) / eps_x * dx),
		const Int3 begin = internalBAreaBegin;
		const Int3 end = internalBAreaEnd;
#pragma omp parallel for
		for (int i = begin.x; i < end.x; i++) {
			grid->By(i, 0, 0) += coeffXY * (grid->Ez(i, 0, 0) - grid->Ez(i - 1, 0, 0));
			grid->Bz(i, 0, 0) += -coeffXZ * (grid->Ey(i, 0, 0) - grid->Ey(i - 1, 0, 0));
		}
	}

	// Update grid values of electric field in FDTD.
	inline void FDTD::updateE()
	{
		if (grid->dimensionality == 3)
			updateE3D();
		else if (grid->dimensionality == 2)
			updateE2D();
		else if (grid->dimensionality == 1)
			updateE1D();
	}

	inline void FDTD::updateE3D()
	{
		updateEAreaBegin = Int3(0, 0, 0);
		updateEAreaEnd = grid->numCells - Int3(1, 1, 1);
		for (int d = 0; d < 3; ++d)
		{
			internalEAreaBegin[d] = std::max(updateEAreaBegin[d], pml->leftDims[d]);
			internalEAreaEnd[d] = std::min(updateEAreaEnd[d],
				grid->numCells[d] - pml->rightDims[d]);
		}

		const FP coeffCurrent = -(FP)4 * constants::pi * grid->dt;
		const FP cdt = constants::c * grid->dt;
		const FP coeffXY = cdt / (grid->steps.x * anisotropyCoeff.y);
		const FP coeffXZ = cdt / (grid->steps.x * anisotropyCoeff.z);
		const FP coeffYX = cdt / (grid->steps.y * anisotropyCoeff.x);
		const FP coeffYZ = cdt / (grid->steps.y * anisotropyCoeff.z);
		const FP coeffZX = cdt / (grid->steps.z * anisotropyCoeff.x);
		const FP coeffZY = cdt / (grid->steps.z * anisotropyCoeff.y);

		// In internal area use:
		// e.x(i, j, k) += dt * -4pi * j.x(i, j, k) + c * dt * ((b.z(i, j+1, k) -
		//	 b.z(i, j, k)) / eps_y * dy - (b.y(i, j, k+1) - b.y(i, j, k)) / eps_z * dz),
		// e.y(i, j, k) += dt * -4pi * j.y(i, j, k) + c * dt * ((b.x(i, j, k+1) -
		//	 b.x(i, j, k)) / eps_z * dz - (b.z(i+1, j, k) - b.z(i, j, k)) / eps_x * dx),
		// e.z(i, j, k) += dt * -4pi * j.z(i, j, k) + c * dt * ((b.y(i+1, j, k) -
		//	 b.y(i, j, k)) / eps_x * dx - (b.x(i, j+1, k) - b.x(i, j, k)) / eps_y * dy),
		const Int3 begin = internalEAreaBegin;
		const Int3 end = internalEAreaEnd;
#pragma omp parallel for
		for (int i = begin.x; i < end.x; i++)
			for (int j = begin.y; j < end.y; j++)
			{
#pragma simd
				for (int k = begin.z; k < end.z; k++)
				{
					grid->Ex(i, j, k) += coeffCurrent * grid->Jx(i, j, k) +
						coeffYX * (grid->Bz(i, j + 1, k) - grid->Bz(i, j, k)) -
						coeffZX * (grid->By(i, j, k + 1) - grid->By(i, j, k));
					grid->Ey(i, j, k) += coeffCurrent * grid->Jy(i, j, k) +
						coeffZY * (grid->Bx(i, j, k + 1) - grid->Bx(i, j, k)) -
						coeffXY * (grid->Bz(i + 1, j, k) - grid->Bz(i, j, k));
					grid->Ez(i, j, k) += coeffCurrent * grid->Jz(i, j, k) +
						coeffXZ * (grid->By(i + 1, j, k) - grid->By(i, j, k)) -
						coeffYZ * (grid->Bx(i, j + 1, k) - grid->Bx(i, j, k));
				}
			}

		// Process edge values
		if (updateEAreaEnd.x == grid->numCells.x - 1)
		{
			int i = updateEAreaEnd.x;
#pragma omp parallel for
			for (int j = begin.y; j < end.y; j++)
				for (int k = begin.z; k < end.z; k++)
					grid->Ex(i, j, k) += coeffCurrent * grid->Jx(i, j, k) +
					coeffYX * (grid->Bz(i, j + 1, k) - grid->Bz(i, j, k)) -
					coeffZX * (grid->By(i, j, k + 1) - grid->By(i, j, k));
		}
		if (updateEAreaEnd.y == grid->numCells.y - 1)
		{
			int j = updateEAreaEnd.y;
#pragma omp parallel for
			for (int i = begin.x; i < end.x; i++)
				for (int k = begin.z; k < end.z; k++)
					grid->Ey(i, j, k) += coeffCurrent * grid->Jy(i, j, k) +
					coeffZY * (grid->Bx(i, j, k + 1) - grid->Bx(i, j, k)) -
					coeffXY * (grid->Bz(i + 1, j, k) - grid->Bz(i, j, k));
		}
		if (updateEAreaEnd.z == grid->numCells.z - 1)
		{
			int k = updateEAreaEnd.z;
#pragma omp parallel for
			for (int i = begin.x; i < end.x; i++)
				for (int j = begin.y; j < end.y; j++)
					grid->Ez(i, j, k) += coeffCurrent * grid->Jz(i, j, k) +
					coeffXZ * (grid->By(i + 1, j, k) - grid->By(i, j, k)) -
					coeffYZ * (grid->Bx(i, j + 1, k) - grid->Bx(i, j, k));
		}
	}

	inline void FDTD::updateE2D()
	{
		updateEAreaBegin = Int3(0, 0, 0);
		updateEAreaEnd = grid->numCells - Int3(1, 1, 0);
		for (int d = 0; d < 2; ++d)
		{
			internalEAreaBegin[d] = std::max(updateEAreaBegin[d], pml->leftDims[d]);
			internalEAreaEnd[d] = std::min(updateEAreaEnd[d],
				grid->numCells[d] - pml->rightDims[d]);
		}

		const FP coeffCurrent = -(FP)4 * constants::pi * grid->dt;
		const FP cdt = constants::c * grid->dt;
		const FP coeffXY = cdt / (grid->steps.x * anisotropyCoeff.y);
		const FP coeffXZ = cdt / (grid->steps.x * anisotropyCoeff.z);
		const FP coeffYX = cdt / (grid->steps.y * anisotropyCoeff.x);
		const FP coeffYZ = cdt / (grid->steps.y * anisotropyCoeff.z);

		// In internal area use:
		// e.x(i, j, k) += dt * -4pi * j.x(i, j, k) + c * dt * ((b.z(i, j+1, k) -
		//	 b.z(i, j, k)) / eps_y * dy - (b.y(i, j, k+1) - b.y(i, j, k)) / eps_z * dz),
		// e.y(i, j, k) += dt * -4pi * j.y(i, j, k) + c * dt * ((b.x(i, j, k+1) -
		//	 b.x(i, j, k)) / eps_z * dz - (b.z(i+1, j, k) - b.z(i, j, k)) / eps_x * dx),
		// e.z(i, j, k) += dt * -4pi * j.z(i, j, k) + c * dt * ((b.y(i+1, j, k) -
		//	 b.y(i, j, k)) / eps_x * dx - (b.x(i, j+1, k) - b.x(i, j, k)) / eps_y * dy),
		const Int3 begin = internalEAreaBegin;
		const Int3 end = internalEAreaEnd;
#pragma omp parallel for
		for (int i = begin.x; i < end.x; i++) {
#pragma simd
			for (int j = begin.y; j < end.y; j++) {
				grid->Ex(i, j, 0) += coeffCurrent * grid->Jx(i, j, 0) +
					coeffYX * (grid->Bz(i, j + 1, 0) - grid->Bz(i, j, 0));
				grid->Ey(i, j, 0) += coeffCurrent * grid->Jy(i, j, 0) -
					coeffXY * (grid->Bz(i + 1, j, 0) - grid->Bz(i, j, 0));
				grid->Ez(i, j, 0) += coeffCurrent * grid->Jz(i, j, 0) +
					coeffXZ * (grid->By(i + 1, j, 0) - grid->By(i, j, 0)) -
					coeffYZ * (grid->Bx(i, j + 1, 0) - grid->Bx(i, j, 0));
			}
		}

		// Process edge values
		if (updateEAreaEnd.x == grid->numCells.x - 1)
		{
			int i = updateEAreaEnd.x;
#pragma omp parallel for
			for (int j = begin.y; j < end.y; j++)
				grid->Ex(i, j, 0) += coeffCurrent * grid->Jx(i, j, 0) +
				coeffYX * (grid->Bz(i, j + 1, 0) - grid->Bz(i, j, 0));
		}
		if (updateEAreaEnd.y == grid->numCells.y - 1)
		{
			int j = updateEAreaEnd.y;
#pragma omp parallel for
			for (int i = begin.x; i < end.x; i++)
				grid->Ey(i, j, 0) += coeffCurrent * grid->Jy(i, j, 0) -
				coeffXY * (grid->Bz(i + 1, j, 0) - grid->Bz(i, j, 0));
		}
	}

	inline void FDTD::updateE1D()
	{
		updateEAreaBegin = Int3(0, 0, 0);
		updateEAreaEnd = grid->numCells - Int3(1, 0, 0);
		for (int d = 0; d < 1; ++d)
		{
			internalEAreaBegin[d] = std::max(updateEAreaBegin[d], pml->leftDims[d]);
			internalEAreaEnd[d] = std::min(updateEAreaEnd[d],
				grid->numCells[d] - pml->rightDims[d]);
		}

		const FP coeffCurrent = -(FP)4 * constants::pi * grid->dt;
		const FP cdt = constants::c * grid->dt;
		const FP coeffXY = cdt / (grid->steps.x * anisotropyCoeff.y);
		const FP coeffXZ = cdt / (grid->steps.x * anisotropyCoeff.z);

		// In internal area use:
		// e.x(i, j, k) += dt * -4pi * j.x(i, j, k) + c * dt * ((b.z(i, j+1, k) -
		//	 b.z(i, j, k)) / eps_y * dy - (b.y(i, j, k+1) - b.y(i, j, k)) / eps_z * dz),
		// e.y(i, j, k) += dt * -4pi * j.y(i, j, k) + c * dt * ((b.x(i, j, k+1) -
		//	 b.x(i, j, k)) / eps_z * dz - (b.z(i+1, j, k) - b.z(i, j, k)) / eps_x * dx),
		// e.z(i, j, k) += dt * -4pi * j.z(i, j, k) + c * dt * ((b.y(i+1, j, k) -
		//	 b.y(i, j, k)) / eps_x * dx - (b.x(i, j+1, k) - b.x(i, j, k)) / eps_y * dy),
		const Int3 begin = internalEAreaBegin;
		const Int3 end = internalEAreaEnd;
#pragma omp parallel for
		for (int i = begin.x; i < end.x; i++) {
			grid->Ex(i, 0, 0) += coeffCurrent * grid->Jx(i, 0, 0);
			grid->Ey(i, 0, 0) += coeffCurrent * grid->Jy(i, 0, 0) -
				coeffXY * (grid->Bz(i + 1, 0, 0) - grid->Bz(i, 0, 0));
			grid->Ez(i, 0, 0) += coeffCurrent * grid->Jz(i, 0, 0) +
				coeffXZ * (grid->By(i + 1, 0, 0) - grid->By(i, 0, 0));
		}
	}

}