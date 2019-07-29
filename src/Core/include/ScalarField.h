#pragma once
#include <vector>

#include "FormFactor.h"
#include "Vectors.h"
#include "VectorsProxy.h"

namespace pfc {
	/* Class for storing 3d scalar field on a regular grid.
	Provides index-wise access, interpolation and deposition.*/
	template <typename Data>
	class ScalarField
	{
	public:

		ScalarField(const Int3& size);
		ScalarField(const ScalarField& field);
		ScalarField& operator =(const ScalarField& field);

		std::vector<Data>& toVector() {
			return elements;
		}

		Int3 getSize() const
		{
			return size;
		}

		/* Read-only access by scalar indexes */
		Data operator()(int i, int j, int k) const
		{
			return raw[k + (j + i * size.y) * size.z];
		}

		/* Read-write access by scalar indexes */
		Data& operator()(int i, int j, int k)
		{
			return raw[k + (j + i * size.y) * size.z];
		}

		/* Read-only access by vector index */
		Data operator()(const Int3& index) const
		{
			return raw[index.z + (index.y + index.x * size.y) * size.z];
		}

		/* Read-write access by vector index */
		Data& operator()(const Int3& index)
		{
			return raw[index.z + (index.y + index.x * size.y) * size.z];
		}

		/* Interpolation: with given base index and coefficients */
		FP interpolateCIC(const Int3& baseIdx, const FP3& coeffs) const;
		FP interpolateTSC(const Int3& baseIdx, const FP3& coeffs) const;
		FP interpolateSecondOrder(const Int3& baseIdx, const FP3& coeffs) const;
		FP interpolateFourthOrder(const Int3& baseIdx, const FP3& coeffs) const;
		FP interpolatePCS(const Int3& baseIdx, const FP3& coeffs) const;

		/* Make all elements zero */
		inline void zeroize();

	private:

		FP interpolateThreePoints(const Int3& baseIdx, FP c[3][3]) const;

		std::vector<Data> elements; // storage
		Data* raw; // raw pointer to elements vector
		Int3 size; // size of each dimension
		Int3 dimensionCoeffInt; // 0 for fake dimensions, 1 otherwise
		FP3 dimensionCoeffFP; // 0 for fake dimensions, 1 otherwise
	};

	template <class Data>
	inline ScalarField<Data>::ScalarField(const Int3& _size)
	{
		size = _size;
		elements.resize(size.volume());
		raw = elements.data();
		for (int d = 0; d < 3; d++) {
			dimensionCoeffInt[d] = (size[d] > 1) ? 1 : 0;
			dimensionCoeffFP[d] = (FP)dimensionCoeffInt[d];
		}
		zeroize();
	}

	template <class Data>
	inline ScalarField<Data>::ScalarField(const ScalarField& field)
	{
		size = field.size;
		elements = field.elements;
		raw = elements.data();
		dimensionCoeffInt = field.dimensionCoeffInt;
		dimensionCoeffFP = field.dimensionCoeffFP;
	}

	template <class Data>
	inline ScalarField<Data>& ScalarField<Data>::operator =(const ScalarField& field)
	{
		size = field.size;
		elements = field.elements;
		raw = elements.data();
		dimensionCoeffInt = field.dimensionCoeffInt;
		dimensionCoeffFP = field.dimensionCoeffFP;
		return *this;
	}

	template <>
	inline FP ScalarField<FP>::interpolateCIC(const Int3& baseIdx, const FP3& coeffs) const
	{
		FP3 c = coeffs * dimensionCoeffFP;
		FP3 invC = FP3(1, 1, 1) - c;
		Int3 base = baseIdx * dimensionCoeffInt;
		Int3 next = base + dimensionCoeffInt;
		return invC.x * (invC.y * (invC.z * (*this)(base.x, base.y, base.z) + c.z * (*this)(base.x, base.y, next.z)) +
							c.y * (invC.z * (*this)(base.x, next.y, base.z) + c.z * (*this)(base.x, next.y, next.z))) +
				  c.x * (invC.y * (invC.z * (*this)(next.x, base.y, base.z) + c.z * (*this)(next.x, base.y, next.z)) +
							c.y * (invC.z * (*this)(next.x, next.y, base.z) + c.z * (*this)(next.x, next.y, next.z)));
	}
	
	template <class Data>
	inline FP ScalarField<Data>::interpolateTSC(const Int3& baseIdx, const FP3& coeffs) const
	{
		FP c[3][3];
		for (int i = 0; i < 3; i++)
			c[0][i] = formfactorTSC(FP(i - 1) - coeffs.x);
		for (int j = 0; j < 3; j++)
			c[1][j] = formfactorTSC(FP(j - 1) - coeffs.y);
		for (int k = 0; k < 3; k++)
			c[2][k] = formfactorTSC(FP(k - 1) - coeffs.z);
		return interpolateThreePoints(baseIdx, c);
	}
	
	template <class Data>
	inline FP ScalarField<Data>::interpolateSecondOrder(const Int3& baseIdx, const FP3& coeffs) const
	{
		FP c[3][3];
		c[0][0] = (FP)0.5 * (coeffs.x * (coeffs.x - (FP)1));
		c[0][1] = (FP)1 - coeffs.x * coeffs.x;
		c[0][2] = (FP)0.5 * (coeffs.x * (coeffs.x + (FP)1));
		c[1][0] = (FP)0.5 * (coeffs.y * (coeffs.y - (FP)1));
		c[1][1] = (FP)1 - coeffs.y * coeffs.y;
		c[1][2] = (FP)0.5 * (coeffs.y * (coeffs.y + (FP)1));
		c[2][0] = (FP)0.5 * (coeffs.z * (coeffs.z - (FP)1));
		c[2][1] = (FP)1 - coeffs.z * coeffs.z;
		c[2][2] = (FP)0.5 * (coeffs.z * (coeffs.z + (FP)1));
		return interpolateThreePoints(baseIdx, c);
	}
	
	template <>
	inline FP ScalarField<FP>::interpolateThreePoints(const Int3& baseIdx, FP c[3][3]) const
	{
		for (int d = 0; d < 3; d++)
			if (!dimensionCoeffInt[d]) {
				c[d][0] = 0;
				c[d][1] = 1.0;
				c[d][2] = 0;
			}
		FP result = 0;
		Int3 minIndex = Int3(-1, -1, -1) * dimensionCoeffInt;
		Int3 maxIndex = Int3(1, 1, 1) * dimensionCoeffInt;
		Int3 base = baseIdx * dimensionCoeffInt;
		for (int ii = minIndex.x; ii <= maxIndex.x; ii++)
			for (int jj = minIndex.y; jj <= maxIndex.y; jj++)
				for (int kk = minIndex.z; kk <= maxIndex.z; kk++)
					result += c[0][ii + 1] * c[1][jj + 1] * c[2][kk + 1] * (*this)(base.x + ii, base.y + jj, base.z + kk);
		return result;
	}

	template <>
	inline FP ScalarField<FP>::interpolateFourthOrder(const Int3& baseIdx, const FP3& coeffs) const
	{
		Int3 base = baseIdx * dimensionCoeffInt;
		const Int3 minAllowedIdx = Int3(2, 2, 2) * dimensionCoeffInt;
		const Int3 maxAllowedIdx = (size - Int3(4, 4, 4)) * dimensionCoeffInt;
		if ((base >= minAllowedIdx) && (base <= maxAllowedIdx) == false)
			return interpolateTSC(baseIdx, coeffs);
		FP c[3][5];
		formfactorFourthOrder(coeffs.x, c[0]);
		formfactorFourthOrder(coeffs.y, c[1]);
		formfactorFourthOrder(coeffs.z, c[2]);
		for (int d = 0; d < 3; d++)
			if (!dimensionCoeffInt[d]) {
				c[d][0] = 0;
				c[d][1] = 0;
				c[d][2] = 1.0;
				c[d][3] = 0;
				c[d][4] = 0;
			}
		FP result = 0;
		Int3 minIndex = Int3(-2, -2, -2) * dimensionCoeffInt;
		Int3 maxIndex = Int3(2, 2, 2) * dimensionCoeffInt;
		for (int ii = minIndex.x; ii <= maxIndex.x; ii++)
			for (int jj = minIndex.y; jj <= maxIndex.y; jj++)
				for (int kk = minIndex.z; kk <= maxIndex.z; kk++)
					result += c[0][ii + 2] * c[1][jj + 2] * c[2][kk + 2] * (*this)(base.x + ii, base.y + jj, base.z + kk);
		return result;
	}

	template <>
	inline FP ScalarField<FP>::interpolatePCS(const Int3& baseIdx, const FP3& coeffs) const
	{
		FP c[3][4];
		for (int i = 0; i < 4; i++)
			c[0][i] = formfactorPCS(FP(i - 1) - coeffs.x);
		for (int j = 0; j < 4; j++)
			c[1][j] = formfactorPCS(FP(j - 1) - coeffs.y);
		for (int k = 0; k < 4; k++)
			c[2][k] = formfactorPCS(FP(k - 1) - coeffs.z);
		for (int d = 0; d < 3; d++)
			if (!dimensionCoeffInt[d]) {
				c[d][0] = 1.0;
				c[d][1] = 0;
				c[d][2] = 0;
				c[d][3] = 0;
			}
		FP result = 0;
		Int3 base = (baseIdx - Int3(1, 1, 1)) * dimensionCoeffInt;
		Int3 minIndex = Int3(0, 0, 0);
		Int3 maxIndex = Int3(3, 3, 3) * dimensionCoeffInt;
		for (int ii = minIndex.x; ii <= maxIndex.x; ii++)
			for (int jj = minIndex.y; jj <= maxIndex.y; jj++)
				for (int kk = minIndex.z; kk <= maxIndex.z; kk++)
					result += c[0][ii] * c[1][jj] * c[2][kk] * (*this)(base.x + ii, base.y + jj, base.z + kk);
		return result;
	}

	template <>
	inline void ScalarField<FP>::zeroize()
	{
#pragma omp parallel for
		for (int idx = 0; idx < (int)elements.size(); idx++)
			elements[idx] = 0;
	}

	template <>
	inline FP ScalarField<complex>::interpolateCIC(const Int3& baseIdx, const FP3& coeffs) const
	{
		FP3 c = coeffs * dimensionCoeffFP;
		FP3 invC = FP3(1, 1, 1) - c;
		Int3 base = baseIdx * dimensionCoeffInt;
		Int3 next = base + dimensionCoeffInt;
		return invC.x * (invC.y * (invC.z * (*this)(base.x, base.y, base.z).real + c.z * (*this)(base.x, base.y, next.z).real) +
			c.y * (invC.z * (*this)(base.x, next.y, base.z).real + c.z * (*this)(base.x, next.y, next.z).real)) +
			c.x * (invC.y * (invC.z * (*this)(next.x, base.y, base.z).real + c.z * (*this)(next.x, base.y, next.z).real) +
				c.y * (invC.z * (*this)(next.x, next.y, base.z).real + c.z * (*this)(next.x, next.y, next.z).real));
	}
	
	template <>
	inline FP ScalarField<complex>::interpolateThreePoints(const Int3& baseIdx, FP c[3][3]) const
	{
		for (int d = 0; d < 3; d++)
			if (!dimensionCoeffInt[d]) {
				c[d][0] = 0;
				c[d][1] = 1.0;
				c[d][2] = 0;
			}
		FP result = 0;
		Int3 minIndex = Int3(-1, -1, -1) * dimensionCoeffInt;
		Int3 maxIndex = Int3(1, 1, 1) * dimensionCoeffInt;
		Int3 base = baseIdx * dimensionCoeffInt;
		for (int ii = minIndex.x; ii <= maxIndex.x; ii++)
			for (int jj = minIndex.y; jj <= maxIndex.y; jj++)
				for (int kk = minIndex.z; kk <= maxIndex.z; kk++)
					result += c[0][ii + 1] * c[1][jj + 1] * c[2][kk + 1] * (*this)(base.x + ii, base.y + jj, base.z + kk).real;
		return result;
	}

	template <>
	inline FP ScalarField<complex>::interpolateFourthOrder(const Int3& baseIdx, const FP3& coeffs) const
	{
		Int3 base = baseIdx * dimensionCoeffInt;
		const Int3 minAllowedIdx = Int3(2, 2, 2) * dimensionCoeffInt;
		const Int3 maxAllowedIdx = (size - Int3(4, 4, 4)) * dimensionCoeffInt;
		if ((base >= minAllowedIdx) && (base <= maxAllowedIdx) == false)
			return interpolateTSC(baseIdx, coeffs);
		FP c[3][5];
		formfactorFourthOrder(coeffs.x, c[0]);
		formfactorFourthOrder(coeffs.y, c[1]);
		formfactorFourthOrder(coeffs.z, c[2]);
		for (int d = 0; d < 3; d++)
			if (!dimensionCoeffInt[d]) {
				c[d][0] = 0;
				c[d][1] = 0;
				c[d][2] = 1.0;
				c[d][3] = 0;
				c[d][4] = 0;
			}
		FP result = 0;
		Int3 minIndex = Int3(-2, -2, -2) * dimensionCoeffInt;
		Int3 maxIndex = Int3(2, 2, 2) * dimensionCoeffInt;
		for (int ii = minIndex.x; ii <= maxIndex.x; ii++)
			for (int jj = minIndex.y; jj <= maxIndex.y; jj++)
				for (int kk = minIndex.z; kk <= maxIndex.z; kk++)
					result += c[0][ii + 2] * c[1][jj + 2] * c[2][kk + 2] * (*this)(base.x + ii, base.y + jj, base.z + kk).real;
		return result;
	}

	template <>
	inline FP ScalarField<complex>::interpolatePCS(const Int3& baseIdx, const FP3& coeffs) const
	{
		FP c[3][4];
		for (int i = 0; i < 4; i++)
			c[0][i] = formfactorPCS(FP(i - 1) - coeffs.x);
		for (int j = 0; j < 4; j++)
			c[1][j] = formfactorPCS(FP(j - 1) - coeffs.y);
		for (int k = 0; k < 4; k++)
			c[2][k] = formfactorPCS(FP(k - 1) - coeffs.z);
		for (int d = 0; d < 3; d++)
			if (!dimensionCoeffInt[d]) {
				c[d][0] = 1.0;
				c[d][1] = 0;
				c[d][2] = 0;
				c[d][3] = 0;
			}
		FP result = 0;
		Int3 base = (baseIdx - Int3(1, 1, 1)) * dimensionCoeffInt;
		Int3 minIndex = Int3(0, 0, 0);
		Int3 maxIndex = Int3(3, 3, 3) * dimensionCoeffInt;
		for (int ii = minIndex.x; ii <= maxIndex.x; ii++)
			for (int jj = minIndex.y; jj <= maxIndex.y; jj++)
				for (int kk = minIndex.z; kk <= maxIndex.z; kk++)
					result += c[0][ii] * c[1][jj] * c[2][kk] * (*this)(base.x + ii, base.y + jj, base.z + kk).real;
		return result;
	}

	template <>
	inline void ScalarField<complex>::zeroize()
	{
#pragma omp parallel for
		for (int idx = 0; idx < (int)elements.size(); idx++)
		{
			elements[idx].real = 0;
			elements[idx].imag = 0;
		}
	}
}