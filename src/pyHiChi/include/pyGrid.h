#pragma once
#include "Grid.h"
#include "Mapping.h"

#include "pybind11/pybind11.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace pfc
{
    
	enum pyGridTypes {
		pyGridType,
		pyGridMappingType
	};

    template<class TypeGrid, class TDerived>
	class pyGridAttributes : public TypeGrid
	{
	public:

		pyGridAttributes(const Int3 & _numInternalCells, FP _dt,
			const FP3 & minCoords, const FP3 & _steps) :
			TypeGrid(Int3(_numInternalCells), _dt, minCoords, _steps, _numInternalCells)
		{
			static_assert(std::is_member_function_pointer<decltype(&TDerived::convertCoords)>::value,
				"Wrong instance of pyGrid");
			fEt[0] = 0; fEt[1] = 0; fEt[2] = 0;
			fBt[0] = 0; fBt[1] = 0; fBt[2] = 0;
			isAnalytical = false;
		}

		void setTime(FP time) { globalT = time; }

		void setAnalytical(int64_t _fEx, int64_t _fEy, int64_t _fEz, int64_t _fBx, int64_t _fBy, int64_t _fBz)
		{
			fEt[0] = _fEx; fEt[1] = _fEy; fEt[2] = _fEz;
			fBt[0] = _fBx; fBt[1] = _fBy; fBt[2] = _fBz;
			isAnalytical = true;
		}

		FP3 getE(const FP3& coords) const
		{
			FP3 result;
			if (isAnalytical)
			{
				FP(*fx)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))fEt[0];
				FP(*fy)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))fEt[1];
				FP(*fz)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))fEt[2];
				result[0] = fx(coords.x, coords.y, coords.z, globalT + this->timeShiftE);
				result[1] = fy(coords.x, coords.y, coords.z, globalT + this->timeShiftE);
				result[2] = fz(coords.x, coords.y, coords.z, globalT + this->timeShiftE);
			}
			else {
				result = TypeGrid::getE(coords);
			}
			return result;
		}

		FP3 getB(const FP3& coords) const
		{
			FP3 result;
			if (isAnalytical)
			{
				FP(*fx)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))fBt[0];
				FP(*fy)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))fBt[1];
				FP(*fz)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))fBt[2];
				result[0] = fx(coords.x, coords.y, coords.z, globalT + this->timeShiftB);
				result[1] = fy(coords.x, coords.y, coords.z, globalT + this->timeShiftB);
				result[2] = fz(coords.x, coords.y, coords.z, globalT + this->timeShiftB);
			}
			else {
				result = TypeGrid::getB(coords);
			}
			return result;
		}


		void analyticalUpdateFields(FP t)
		{
			if (isAnalytical)
			{
				setExyzt(fEt[0], fEt[1], fEt[2], t);
				setBxyzt(fBt[0], fBt[1], fBt[2], t);
			}
		}

		template <class FieldConficurationType>
		void setFieldConfiguration(const FieldConficurationType* fieldConf) {
			TDerived* derived = static_cast<TDerived*>(this);
#pragma omp parallel for
			for (int i = 0; i < this->numCells.x; i++)
				for (int j = 0; j < this->numCells.y; j++)
					for (int k = 0; k < this->numCells.z; k++)
					{
						FP3 cEx, cEy, cEz, cBx, cBy, cBz, cJx, cJy, cJz;

						cEx = derived->convertCoords(this->ExPosition(i, j, k));
						cEy = derived->convertCoords(this->EyPosition(i, j, k));
						cEz = derived->convertCoords(this->EzPosition(i, j, k));
						this->Ex(i, j, k) = fieldConf->E(cEx.x, cEx.y, cEx.z).x;
						this->Ey(i, j, k) = fieldConf->E(cEy.x, cEy.y, cEy.z).y;
						this->Ez(i, j, k) = fieldConf->E(cEz.x, cEz.y, cEz.z).z;

						cBx = derived->convertCoords(this->ExPosition(i, j, k));
						cBy = derived->convertCoords(this->EyPosition(i, j, k));
						cBz = derived->convertCoords(this->EzPosition(i, j, k));
						this->Bx(i, j, k) = fieldConf->B(cBx.x, cBx.y, cBx.z).x;
						this->By(i, j, k) = fieldConf->B(cBy.x, cBy.y, cBy.z).y;
						this->Bz(i, j, k) = fieldConf->B(cBz.x, cBz.y, cBz.z).z;
					}
		}

		void pySetExyz(py::function fEx, py::function fEy, py::function fEz)
		{
			TDerived* derived = static_cast<TDerived*>(this);
			for (int i = 0; i < this->numCells.x; i++)
				for (int j = 0; j < this->numCells.y; j++)
					for (int k = 0; k < this->numCells.z; k++)
					{
						FP3 cEx, cEy, cEz;
						cEx = derived->convertCoords(this->ExPosition(i, j, k));
						cEy = derived->convertCoords(this->EyPosition(i, j, k));
						cEz = derived->convertCoords(this->EzPosition(i, j, k));
						this->Ex(i, j, k) = fEx("x"_a = cEx.x, "y"_a = cEx.y, "z"_a = cEx.z).template cast<FP>();
						this->Ey(i, j, k) = fEy("x"_a = cEy.x, "y"_a = cEy.y, "z"_a = cEy.z).template cast<FP>();
						this->Ez(i, j, k) = fEz("x"_a = cEz.x, "y"_a = cEz.y, "z"_a = cEz.z).template cast<FP>();
					}
		}

		void pySetE(py::function fE)
		{
			TDerived* derived = static_cast<TDerived*>(this);
			for (int i = 0; i < this->numCells.x; i++)
				for (int j = 0; j < this->numCells.y; j++)
					for (int k = 0; k < this->numCells.z; k++)
					{
						FP3 cEx, cEy, cEz;
						cEx = derived->convertCoords(this->ExPosition(i, j, k));
						cEy = derived->convertCoords(this->EyPosition(i, j, k));
						cEz = derived->convertCoords(this->EzPosition(i, j, k));
						this->Ex(i, j, k) = fE("x"_a = cEx.x, "y"_a = cEx.y, "z"_a = cEx.z).template cast<FP3>().x;
						this->Ey(i, j, k) = fE("x"_a = cEy.x, "y"_a = cEy.y, "z"_a = cEy.z).template cast<FP3>().y;
						this->Ez(i, j, k) = fE("x"_a = cEz.x, "y"_a = cEz.y, "z"_a = cEz.z).template cast<FP3>().z;
					}
		}

		void setExyz(int64_t _fEx, int64_t _fEy, int64_t _fEz)
		{
			TDerived* derived = static_cast<TDerived*>(this);
			FP(*fEx)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fEx;
			FP(*fEy)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fEy;
			FP(*fEz)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fEz;
#pragma omp parallel for
			for (int i = 0; i < this->numCells.x; i++)
				for (int j = 0; j < this->numCells.y; j++)
					for (int k = 0; k < this->numCells.z; k++)
					{
						FP3 cEx, cEy, cEz;
						cEx = derived->convertCoords(this->ExPosition(i, j, k));
						cEy = derived->convertCoords(this->EyPosition(i, j, k));
						cEz = derived->convertCoords(this->EzPosition(i, j, k));
						this->Ex(i, j, k) = fEx(cEx.x, cEx.y, cEx.z);
						this->Ey(i, j, k) = fEy(cEy.x, cEy.y, cEy.z);
						this->Ez(i, j, k) = fEz(cEz.x, cEz.y, cEz.z);
					}
		}

		void setExyzt(int64_t _fEx, int64_t _fEy, int64_t _fEz, FP t)
		{
			TDerived* derived = static_cast<TDerived*>(this);
			FP(*fEx)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fEx;
			FP(*fEy)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fEy;
			FP(*fEz)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fEz;
#pragma omp parallel for
			for (int i = 0; i < this->numCells.x; i++)
				for (int j = 0; j < this->numCells.y; j++)
					for (int k = 0; k < this->numCells.z; k++)
					{
						FP3 cEx, cEy, cEz;
						cEx = derived->convertCoords(this->ExPosition(i, j, k));
						cEy = derived->convertCoords(this->EyPosition(i, j, k));
						cEz = derived->convertCoords(this->EzPosition(i, j, k));
						this->Ex(i, j, k) = fEx(cEx.x, cEx.y, cEx.z, t + this->timeShiftE);
						this->Ey(i, j, k) = fEy(cEy.x, cEy.y, cEy.z, t + this->timeShiftE);
						this->Ez(i, j, k) = fEz(cEz.x, cEz.y, cEz.z, t + this->timeShiftE);
					}
		}

		void setE(int64_t _fE)
		{
			TDerived* derived = static_cast<TDerived*>(this);
			FP3(*fE)(FP, FP, FP) = (FP3(*)(FP, FP, FP))_fE;
#pragma omp parallel for
			for (int i = 0; i < this->numCells.x; i++)
				for (int j = 0; j < this->numCells.y; j++)
					for (int k = 0; k < this->numCells.z; k++)
					{
						FP3 cEx, cEy, cEz;
						cEx = derived->convertCoords(this->ExPosition(i, j, k));
						cEy = derived->convertCoords(this->EyPosition(i, j, k));
						cEz = derived->convertCoords(this->EzPosition(i, j, k));
						this->Ex(i, j, k) = fE(cEx.x, cEx.y, cEx.z).x;
						this->Ey(i, j, k) = fE(cEy.x, cEy.y, cEy.z).y;
						this->Ez(i, j, k) = fE(cEz.x, cEz.y, cEz.z).z;
					}
		}

		void pySetBxyz(py::function fBx, py::function fBy, py::function fBz)
		{
			TDerived* derived = static_cast<TDerived*>(this);
			for (int i = 0; i < this->numCells.x; i++)
				for (int j = 0; j < this->numCells.y; j++)
					for (int k = 0; k < this->numCells.z; k++)
					{
						FP3 cBx, cBy, cBz;
						cBx = derived->convertCoords(this->BxPosition(i, j, k));
						cBy = derived->convertCoords(this->ByPosition(i, j, k));
						cBz = derived->convertCoords(this->BzPosition(i, j, k));
						this->Bx(i, j, k) = fBx("x"_a = cBx.x, "y"_a = cBx.y, "z"_a = cBx.z).template cast<FP>();
						this->By(i, j, k) = fBy("x"_a = cBy.x, "y"_a = cBy.y, "z"_a = cBy.z).template cast<FP>();
						this->Bz(i, j, k) = fBz("x"_a = cBz.x, "y"_a = cBz.y, "z"_a = cBz.z).template cast<FP>();
					}
		}

		void pySetB(py::function fB)
		{
			TDerived* derived = static_cast<TDerived*>(this);
			for (int i = 0; i < this->numCells.x; i++)
				for (int j = 0; j < this->numCells.y; j++)
					for (int k = 0; k < this->numCells.z; k++)
					{
						FP3 cBx, cBy, cBz;
						cBx = derived->convertCoords(this->BxPosition(i, j, k));
						cBy = derived->convertCoords(this->ByPosition(i, j, k));
						cBz = derived->convertCoords(this->BzPosition(i, j, k));
						this->Bx(i, j, k) = fB("x"_a = cBx.x, "y"_a = cBx.y, "z"_a = cBx.z).template cast<FP3>().x;
						this->By(i, j, k) = fB("x"_a = cBy.x, "y"_a = cBy.y, "z"_a = cBy.z).template cast<FP3>().y;
						this->Bz(i, j, k) = fB("x"_a = cBz.x, "y"_a = cBz.y, "z"_a = cBz.z).template cast<FP3>().z;
					}
		}

		void setBxyz(int64_t _fBx, int64_t _fBy, int64_t _fBz)
		{
			TDerived* derived = static_cast<TDerived*>(this);
			FP(*fBx)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fBx;
			FP(*fBy)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fBy;
			FP(*fBz)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fBz;
#pragma omp parallel for
			for (int i = 0; i < this->numCells.x; i++)
				for (int j = 0; j < this->numCells.y; j++)
					for (int k = 0; k < this->numCells.z; k++)
					{
						FP3 cBx, cBy, cBz;
						cBx = derived->convertCoords(this->BxPosition(i, j, k));
						cBy = derived->convertCoords(this->ByPosition(i, j, k));
						cBz = derived->convertCoords(this->BzPosition(i, j, k));
						this->Bx(i, j, k) = fBx(cBx.x, cBx.y, cBx.z);
						this->By(i, j, k) = fBy(cBy.x, cBy.y, cBy.z);
						this->Bz(i, j, k) = fBz(cBz.x, cBz.y, cBz.z);
					}
		}

		void setBxyzt(int64_t _fBx, int64_t _fBy, int64_t _fBz, FP t)
		{
			TDerived* derived = static_cast<TDerived*>(this);
			FP(*fBx)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fBx;
			FP(*fBy)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fBy;
			FP(*fBz)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fBz;
#pragma omp parallel for
			for (int i = 0; i < this->numCells.x; i++)
				for (int j = 0; j < this->numCells.y; j++)
					for (int k = 0; k < this->numCells.z; k++)
					{
						FP3 cBx, cBy, cBz;
						cBx = derived->convertCoords(this->BxPosition(i, j, k));
						cBy = derived->convertCoords(this->ByPosition(i, j, k));
						cBz = derived->convertCoords(this->BzPosition(i, j, k));
						this->Bx(i, j, k) = fBx(cBx.x, cBx.y, cBx.z, t + this->timeShiftB);
						this->By(i, j, k) = fBy(cBy.x, cBy.y, cBy.z, t + this->timeShiftB);
						this->Bz(i, j, k) = fBz(cBz.x, cBz.y, cBz.z, t + this->timeShiftB);
					}
		}

		void setB(int64_t _fB)
		{
			TDerived* derived = static_cast<TDerived*>(this);
			FP3(*fB)(FP, FP, FP) = (FP3(*)(FP, FP, FP))_fB;
#pragma omp parallel for
			for (int i = 0; i < this->numCells.x; i++)
				for (int j = 0; j < this->numCells.y; j++)
					for (int k = 0; k < this->numCells.z; k++)
					{
						FP3 cBx, cBy, cBz;
						cBx = derived->convertCoords(this->BxPosition(i, j, k));
						cBy = derived->convertCoords(this->ByPosition(i, j, k));
						cBz = derived->convertCoords(this->BzPosition(i, j, k));
						this->Bx(i, j, k) = fB(cBx.x, cBx.y, cBx.z).x;
						this->By(i, j, k) = fB(cBy.x, cBy.y, cBy.z).y;
						this->Bz(i, j, k) = fB(cBz.x, cBz.y, cBz.z).z;
					}
		}

		void pySetJxyz(py::function fJx, py::function fJy, py::function fJz)
		{
			TDerived* derived = static_cast<TDerived*>(this);
			for (int i = 0; i < this->numCells.x; i++)
				for (int j = 0; j < this->numCells.y; j++)
					for (int k = 0; k < this->numCells.z; k++)
					{
						FP3 cJx, cJy, cJz;
						cJx = derived->convertCoords(this->JxPosition(i, j, k));
						cJy = derived->convertCoords(this->JyPosition(i, j, k));
						cJz = derived->convertCoords(this->JzPosition(i, j, k));
						this->Jx(i, j, k) = fJx("x"_a = cJx.x, "y"_a = cJx.y, "z"_a = cJx.z).template cast<FP>();
						this->Jy(i, j, k) = fJy("x"_a = cJy.x, "y"_a = cJy.y, "z"_a = cJy.z).template cast<FP>();
						this->Jz(i, j, k) = fJz("x"_a = cJz.x, "y"_a = cJz.y, "z"_a = cJz.z).template cast<FP>();
					}
		}

		void pySetJ(py::function fJ)
		{
			TDerived* derived = static_cast<TDerived*>(this);
			for (int i = 0; i < this->numCells.x; i++)
				for (int j = 0; j < this->numCells.y; j++)
					for (int k = 0; k < this->numCells.z; k++)
					{
						FP3 cJx, cJy, cJz;
						cJx = derived->convertCoords(this->JxPosition(i, j, k));
						cJy = derived->convertCoords(this->JyPosition(i, j, k));
						cJz = derived->convertCoords(this->JzPosition(i, j, k));
						this->Jx(i, j, k) = fJ("x"_a = cJx.x, "y"_a = cJx.y, "z"_a = cJx.z).template cast<FP3>().x;
						this->Jy(i, j, k) = fJ("x"_a = cJy.x, "y"_a = cJy.y, "z"_a = cJy.z).template cast<FP3>().y;
						this->Jz(i, j, k) = fJ("x"_a = cJz.x, "y"_a = cJz.y, "z"_a = cJz.z).template cast<FP3>().z;
					}
		}

		void setJxyz(int64_t _fJx, int64_t _fJy, int64_t _fJz)
		{
			TDerived* derived = static_cast<TDerived*>(this);
			FP(*fJx)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fJx;
			FP(*fJy)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fJy;
			FP(*fJz)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fJz;
#pragma omp parallel for
			for (int i = 0; i < this->numCells.x; i++)
				for (int j = 0; j < this->numCells.y; j++)
					for (int k = 0; k < this->numCells.z; k++)
					{
						FP3 cJx, cJy, cJz;
						cJx = derived->convertCoords(this->JxPosition(i, j, k));
						cJy = derived->convertCoords(this->JyPosition(i, j, k));
						cJz = derived->convertCoords(this->JzPosition(i, j, k));
						this->Jx(i, j, k) = fJx(cJx.x, cJx.y, cJx.z);
						this->Jy(i, j, k) = fJy(cJy.x, cJy.y, cJy.z);
						this->Jz(i, j, k) = fJz(cJz.x, cJz.y, cJz.z);
					}
		}

		void setJxyzt(int64_t _fJx, int64_t _fJy, int64_t _fJz, FP t)
		{
			TDerived* derived = static_cast<TDerived*>(this);
			FP(*fJx)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fJx;
			FP(*fJy)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fJy;
			FP(*fJz)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fJz;
#pragma omp parallel for
			for (int i = 0; i < this->numCells.x; i++)
				for (int j = 0; j < this->numCells.y; j++)
					for (int k = 0; k < this->numCells.z; k++)
					{
						FP3 cJx, cJy, cJz;
						cJx = derived->convertCoords(this->JxPosition(i, j, k));
						cJy = derived->convertCoords(this->JyPosition(i, j, k));
						cJz = derived->convertCoords(this->JzPosition(i, j, k));
						this->Jx(i, j, k) = fJx(cJx.x, cJx.y, cJx.z, t + this->timeShiftJ);
						this->Jy(i, j, k) = fJy(cJy.x, cJy.y, cJy.z, t + this->timeShiftJ);
						this->Jz(i, j, k) = fJz(cJz.x, cJz.y, cJz.z, t + this->timeShiftJ);
					}
		}

		void setJ(int64_t _fJ)
		{
			TDerived* derived = static_cast<TDerived*>(this);
			FP3(*fJ)(FP, FP, FP) = (FP3(*)(FP, FP, FP))_fJ;
#pragma omp parallel for
			for (int i = 0; i < this->numCells.x; i++)
				for (int j = 0; j < this->numCells.y; j++)
					for (int k = 0; k < this->numCells.z; k++)
					{
						FP3 cJx, cJy, cJz;
						cJx = derived->convertCoords(this->JxPosition(i, j, k));
						cJy = derived->convertCoords(this->JyPosition(i, j, k));
						cJz = derived->convertCoords(this->JzPosition(i, j, k));
						this->Jx(i, j, k) = fJ(cJx.x, cJx.y, cJx.z).x;
						this->Jy(i, j, k) = fJ(cJy.x, cJy.y, cJy.z).y;
						this->Jz(i, j, k) = fJ(cJz.x, cJz.y, cJz.z).z;
					}
		}

	private:

		int64_t fEt[3], fBt[3];
		FP globalT;
		bool isAnalytical;

	};


	template<class TypeGrid>
	class pyGrid : public pyGridAttributes<TypeGrid, pyGrid<TypeGrid>>
	{
	public:

		pyGrid(const Int3 & _numInternalCells, FP _dt,
			const FP3 & minCoords, const FP3 & _steps) :
			pyGridAttributes<TypeGrid, pyGrid<TypeGrid>>(_numInternalCells, _dt, minCoords, _steps) {}

		inline FP3 convertCoords(const FP3& coords) const {
			return coords;
		}
	};


	typedef pyGrid<YeeGrid> pyYeeGrid;
	typedef pyGrid<PSTDGrid> pyPSTDGrid;
	typedef pyGrid<PSATDGrid> pyPSATDGrid;


	template<class TypeGrid>
    class pyGridMapping : public pyGridAttributes<TypeGrid, pyGridMapping<TypeGrid>>
    {
    public:

        pyGridMapping(const Int3 & _numInternalCells, FP _dt,
            const FP3 & minCoords, const FP3 & _steps) :
			pyGridAttributes<TypeGrid, pyGridMapping<TypeGrid>>(_numInternalCells, _dt, minCoords, _steps) {}

        FP3 getE(const FP3& coords) const {
            bool status = false;
            FP3 inverseCoords = getInverseCoords(coords, &status);
            if (!status)
                return FP3(0, 0, 0);
            return pyGridAttributes::getE(inverseCoords);
        }

        FP3 getB(const FP3& coords) const {
            bool status = false;
            FP3 inverseCoords = getInverseCoords(coords, &status);
            if (!status)
                return FP3(0, 0, 0);
            return pyGridAttributes::getB(inverseCoords);
        }

		void setMapping(Mapping* mapping) {
			mappings.push_back(mapping);
		}

		inline FP3 convertCoords(const FP3& coords) const {
			bool status = true;
			return getDirectCoords(coords, &status);
		}

    private:

        std::vector<Mapping*> mappings;
		
		inline FP3 getDirectCoords(const FP3& coords, bool* status) const {
			FP3 coords_ = coords;
			bool status_ = true;
			*status = true;
			for (size_t i = 0; i < mappings.size(); i++) {
				coords_ = mappings[i]->getDirectCoords(coords_, &status_);
				*status = (*status) && status_;
			}
			return coords_;
			//return mappings[0]->getDirectCoords(coords, status);
		}

		inline FP3 getInverseCoords(const FP3& coords, bool* status) const {
			FP3 coords_ = coords;
			bool status_ = true;
			*status = true;
			for (size_t i = mappings.size(); i >= 1; i--) {
				coords_ = mappings[i-1]->getInverseCoords(coords_, &status_);
				*status = (*status) && status_;
			}
			return coords_;
			//return mappings[0]->getInverseCoords(coords, status);
		}

    };

	typedef pyGridMapping<YeeGrid> pyYeeGridMapping;
	typedef pyGridMapping<PSTDGrid> pyPSTDGridMapping;
	typedef pyGridMapping<PSATDGrid> pyPSATDGridMapping;
}
