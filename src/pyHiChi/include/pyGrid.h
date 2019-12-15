#pragma once
#include "Grid.h"
#include "Mapping.h"

#include "pybind11/pybind11.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace pfc
{
    template<class TypeGrid>
	class pyGrid : public TypeGrid
	{
	public:

		pyGrid(const Int3 & _numInternalCells, FP _dt,
			const FP3 & minCoords, const FP3 & _steps) :
			TypeGrid(Int3(_numInternalCells), _dt, minCoords, _steps, _numInternalCells)
		{
			fEt[0] = 0; fEt[1] = 0; fEt[2] = 0;
			fBt[0] = 0; fBt[1] = 0; fBt[2] = 0;
			isAnalytical = false;
            mapping = 0;
		}

		void setTime(FP time) { globalT = time; }
        void setMapping(Mapping* mapping) { this->mapping = mapping; }

		void setAnalytical(int64_t _fEx, int64_t _fEy, int64_t _fEz, int64_t _fBx, int64_t _fBy, int64_t _fBz)
		{
			fEt[0] = _fEx; fEt[1] = _fEy; fEt[2] = _fEz;
			fBt[0] = _fBx; fBt[1] = _fBy; fBt[2] = _fBz;
			isAnalytical = true;
            mapping = 0;
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
            else if (mapping) {
                result = TypeGrid::getE(mapping->getInverseCoords(coords));
            }
			else
			{
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
            else if (mapping) {
                result = TypeGrid::getB(mapping->getInverseCoords(coords));
            }
			else
			{
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

        void pySetExyz(py::function fEx, py::function fEy, py::function fEz)
        {
            for (int i = 0; i < this->numCells.x; i++)
            for (int j = 0; j < this->numCells.y; j++)
            for (int k = 0; k < this->numCells.z; k++)
            {
                FP3 cEx, cEy, cEz;
                cEx = mapping ? mapping->getDirectCoords(this->ExPosition(i, j, k)) : this->ExPosition(i, j, k);
                cEy = mapping ? mapping->getDirectCoords(this->EyPosition(i, j, k)) : this->EyPosition(i, j, k);
                cEz = mapping ? mapping->getDirectCoords(this->EzPosition(i, j, k)) : this->EzPosition(i, j, k);
                this->Ex(i, j, k) = fEx("x"_a = cEx.x, "y"_a = cEx.y, "z"_a = cEx.z).template cast<FP>();
                this->Ey(i, j, k) = fEy("x"_a = cEy.x, "y"_a = cEy.y, "z"_a = cEy.z).template cast<FP>();
                this->Ez(i, j, k) = fEz("x"_a = cEz.x, "y"_a = cEz.y, "z"_a = cEz.z).template cast<FP>();
            }
        }

		void pySetE(py::function fE)
		{
			for (int i = 0; i < this->numCells.x; i++)
			for (int j = 0; j < this->numCells.y; j++)
			for (int k = 0; k < this->numCells.z; k++)
			{
				FP3 cEx, cEy, cEz;
                cEx = mapping ? mapping->getDirectCoords(this->ExPosition(i, j, k)) : this->ExPosition(i, j, k);
                cEy = mapping ? mapping->getDirectCoords(this->EyPosition(i, j, k)) : this->EyPosition(i, j, k);
                cEz = mapping ? mapping->getDirectCoords(this->EzPosition(i, j, k)) : this->EzPosition(i, j, k);
                this->Ex(i, j, k) = fE("x"_a = cEx.x, "y"_a = cEx.y, "z"_a = cEx.z).template cast<FP3>().x;
				this->Ey(i, j, k) = fE("x"_a = cEy.x, "y"_a = cEy.y, "z"_a = cEy.z).template cast<FP3>().y;
				this->Ez(i, j, k) = fE("x"_a = cEz.x, "y"_a = cEz.y, "z"_a = cEz.z).template cast<FP3>().z;
			}
		}

		void setExyz(int64_t _fEx, int64_t _fEy, int64_t _fEz)
		{
			FP(*fEx)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fEx;
			FP(*fEy)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fEy;
			FP(*fEz)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fEz;
#pragma omp parallel for
			for (int i = 0; i < this->numCells.x; i++)
			for (int j = 0; j < this->numCells.y; j++)
			for (int k = 0; k < this->numCells.z; k++)
			{
				FP3 cEx, cEy, cEz;
                cEx = mapping ? mapping->getDirectCoords(this->ExPosition(i, j, k)) : this->ExPosition(i, j, k);
                cEy = mapping ? mapping->getDirectCoords(this->EyPosition(i, j, k)) : this->EyPosition(i, j, k);
                cEz = mapping ? mapping->getDirectCoords(this->EzPosition(i, j, k)) : this->EzPosition(i, j, k);
                this->Ex(i, j, k) = fEx(cEx.x, cEx.y, cEx.z);
				this->Ey(i, j, k) = fEy(cEy.x, cEy.y, cEy.z);
				this->Ez(i, j, k) = fEz(cEz.x, cEz.y, cEz.z);
			}
		}

		void setExyzt(int64_t _fEx, int64_t _fEy, int64_t _fEz, FP t)
		{
			
			FP(*fEx)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fEx;
			FP(*fEy)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fEy;
			FP(*fEz)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fEz;
#pragma omp parallel for
			for (int i = 0; i < this->numCells.x; i++)
			for (int j = 0; j < this->numCells.y; j++)
			for (int k = 0; k < this->numCells.z; k++)
			{
				FP3 cEx, cEy, cEz;
                cEx = mapping ? mapping->getDirectCoords(this->ExPosition(i, j, k)) : this->ExPosition(i, j, k);
                cEy = mapping ? mapping->getDirectCoords(this->EyPosition(i, j, k)) : this->EyPosition(i, j, k);
                cEz = mapping ? mapping->getDirectCoords(this->EzPosition(i, j, k)) : this->EzPosition(i, j, k);
                this->Ex(i, j, k) = fEx(cEx.x, cEx.y, cEx.z, t + this->timeShiftE);
				this->Ey(i, j, k) = fEy(cEy.x, cEy.y, cEy.z, t + this->timeShiftE);
				this->Ez(i, j, k) = fEz(cEz.x, cEz.y, cEz.z, t + this->timeShiftE);
			}
		}

		void setE(int64_t _fE)
		{
			FP3(*fE)(FP, FP, FP) = (FP3(*)(FP, FP, FP))_fE;
#pragma omp parallel for
			for (int i = 0; i < this->numCells.x; i++)
			for (int j = 0; j < this->numCells.y; j++)
			for (int k = 0; k < this->numCells.z; k++)
			{
				FP3 cEx, cEy, cEz;
                cEx = mapping ? mapping->getDirectCoords(this->ExPosition(i, j, k)) : this->ExPosition(i, j, k);
                cEy = mapping ? mapping->getDirectCoords(this->EyPosition(i, j, k)) : this->EyPosition(i, j, k);
                cEz = mapping ? mapping->getDirectCoords(this->EzPosition(i, j, k)) : this->EzPosition(i, j, k);
                this->Ex(i, j, k) = fE(cEx.x, cEx.y, cEx.z).x;
				this->Ey(i, j, k) = fE(cEy.x, cEy.y, cEy.z).y;
				this->Ez(i, j, k) = fE(cEz.x, cEz.y, cEz.z).z;
			}
		}

		void pySetBxyz(py::function fBx, py::function fBy, py::function fBz)
		{
			for (int i = 0; i < this->numCells.x; i++)
			for (int j = 0; j < this->numCells.y; j++)
			for (int k = 0; k < this->numCells.z; k++)
			{
				FP3 cBx, cBy, cBz;
                cBx = mapping ? mapping->getDirectCoords(this->BxPosition(i, j, k)) : this->BxPosition(i, j, k);
                cBy = mapping ? mapping->getDirectCoords(this->ByPosition(i, j, k)) : this->ByPosition(i, j, k);
                cBz = mapping ? mapping->getDirectCoords(this->BzPosition(i, j, k)) : this->BzPosition(i, j, k);
                this->Bx(i, j, k) = fBx("x"_a = cBx.x, "y"_a = cBx.y, "z"_a = cBx.z).template cast<FP>();
				this->By(i, j, k) = fBy("x"_a = cBy.x, "y"_a = cBy.y, "z"_a = cBy.z).template cast<FP>();
				this->Bz(i, j, k) = fBz("x"_a = cBz.x, "y"_a = cBz.y, "z"_a = cBz.z).template cast<FP>();
			}
		}

		void pySetB(py::function fB)
		{
			for (int i = 0; i < this->numCells.x; i++)
			for (int j = 0; j < this->numCells.y; j++)
			for (int k = 0; k < this->numCells.z; k++)
			{
				FP3 cBx, cBy, cBz;
                cBx = mapping ? mapping->getDirectCoords(this->BxPosition(i, j, k)) : this->BxPosition(i, j, k);
                cBy = mapping ? mapping->getDirectCoords(this->ByPosition(i, j, k)) : this->ByPosition(i, j, k);
                cBz = mapping ? mapping->getDirectCoords(this->BzPosition(i, j, k)) : this->BzPosition(i, j, k);
                this->Bx(i, j, k) = fB("x"_a = cBx.x, "y"_a = cBx.y, "z"_a = cBx.z).template cast<FP3>().x;
				this->By(i, j, k) = fB("x"_a = cBy.x, "y"_a = cBy.y, "z"_a = cBy.z).template cast<FP3>().y;
				this->Bz(i, j, k) = fB("x"_a = cBz.x, "y"_a = cBz.y, "z"_a = cBz.z).template cast<FP3>().z;
			}
		}

		void setBxyz(int64_t _fBx, int64_t _fBy, int64_t _fBz)
		{
			FP(*fBx)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fBx;
			FP(*fBy)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fBy;
			FP(*fBz)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fBz;
#pragma omp parallel for
			for (int i = 0; i < this->numCells.x; i++)
			for (int j = 0; j < this->numCells.y; j++)
			for (int k = 0; k < this->numCells.z; k++)
			{
				FP3 cBx, cBy, cBz;
                cBx = mapping ? mapping->getDirectCoords(this->BxPosition(i, j, k)) : this->BxPosition(i, j, k);
                cBy = mapping ? mapping->getDirectCoords(this->ByPosition(i, j, k)) : this->ByPosition(i, j, k);
                cBz = mapping ? mapping->getDirectCoords(this->BzPosition(i, j, k)) : this->BzPosition(i, j, k);
                this->Bx(i, j, k) = fBx(cBx.x, cBx.y, cBx.z);
				this->By(i, j, k) = fBy(cBy.x, cBy.y, cBy.z);
				this->Bz(i, j, k) = fBz(cBz.x, cBz.y, cBz.z);
			}
		}

		void setBxyzt(int64_t _fBx, int64_t _fBy, int64_t _fBz, FP t)
		{
			FP(*fBx)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fBx;
			FP(*fBy)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fBy;
			FP(*fBz)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fBz;
#pragma omp parallel for
			for (int i = 0; i < this->numCells.x; i++)
			for (int j = 0; j < this->numCells.y; j++)
		    for (int k = 0; k < this->numCells.z; k++)
		    {
		    	FP3 cBx, cBy, cBz;
                cBx = mapping ? mapping->getDirectCoords(this->BxPosition(i, j, k)) : this->BxPosition(i, j, k);
                cBy = mapping ? mapping->getDirectCoords(this->ByPosition(i, j, k)) : this->ByPosition(i, j, k);
                cBz = mapping ? mapping->getDirectCoords(this->BzPosition(i, j, k)) : this->BzPosition(i, j, k);
                this->Bx(i, j, k) = fBx(cBx.x, cBx.y, cBx.z, t + this->timeShiftB);
		    	this->By(i, j, k) = fBy(cBy.x, cBy.y, cBy.z, t + this->timeShiftB);
		    	this->Bz(i, j, k) = fBz(cBz.x, cBz.y, cBz.z, t + this->timeShiftB);
		    }
		}

		void setB(int64_t _fB)
		{
			FP3(*fB)(FP, FP, FP) = (FP3(*)(FP, FP, FP))_fB;
#pragma omp parallel for
			for (int i = 0; i < this->numCells.x; i++)
			for (int j = 0; j < this->numCells.y; j++)
			for (int k = 0; k < this->numCells.z; k++)
			{
				FP3 cBx, cBy, cBz;
                cBx = mapping ? mapping->getDirectCoords(this->BxPosition(i, j, k)) : this->BxPosition(i, j, k);
                cBy = mapping ? mapping->getDirectCoords(this->ByPosition(i, j, k)) : this->ByPosition(i, j, k);
                cBz = mapping ? mapping->getDirectCoords(this->BzPosition(i, j, k)) : this->BzPosition(i, j, k);
                this->Bx(i, j, k) = fB(cBx.x, cBx.y, cBx.z).x;
				this->By(i, j, k) = fB(cBy.x, cBy.y, cBy.z).y;
				this->Bz(i, j, k) = fB(cBz.x, cBz.y, cBz.z).z;
			}
		}

		void pySetJxyz(py::function fJx, py::function fJy, py::function fJz)
		{
			for (int i = 0; i < this->numCells.x; i++)
			for (int j = 0; j < this->numCells.y; j++)
			for (int k = 0; k < this->numCells.z; k++)
			{
				FP3 cJx, cJy, cJz;
                cJx = mapping ? mapping->getDirectCoords(this->JxPosition(i, j, k)) : this->JxPosition(i, j, k);
                cJy = mapping ? mapping->getDirectCoords(this->JyPosition(i, j, k)) : this->JyPosition(i, j, k);
                cJz = mapping ? mapping->getDirectCoords(this->JzPosition(i, j, k)) : this->JzPosition(i, j, k);
                this->Jx(i, j, k) = fJx("x"_a = cJx.x, "y"_a = cJx.y, "z"_a = cJx.z).template cast<FP>();
				this->Jy(i, j, k) = fJy("x"_a = cJy.x, "y"_a = cJy.y, "z"_a = cJy.z).template cast<FP>();
				this->Jz(i, j, k) = fJz("x"_a = cJz.x, "y"_a = cJz.y, "z"_a = cJz.z).template cast<FP>();
			}
		}

		void pySetJ(py::function fJ)
		{
			for (int i = 0; i < this->numCells.x; i++)
			for (int j = 0; j < this->numCells.y; j++)
			for (int k = 0; k < this->numCells.z; k++)
			{
				FP3 cJx, cJy, cJz;
                cJx = mapping ? mapping->getDirectCoords(this->JxPosition(i, j, k)) : this->JxPosition(i, j, k);
                cJy = mapping ? mapping->getDirectCoords(this->JyPosition(i, j, k)) : this->JyPosition(i, j, k);
                cJz = mapping ? mapping->getDirectCoords(this->JzPosition(i, j, k)) : this->JzPosition(i, j, k);
                this->Jx(i, j, k) = fJ("x"_a = cJx.x, "y"_a = cJx.y, "z"_a = cJx.z).template cast<FP3>().x;
				this->Jy(i, j, k) = fJ("x"_a = cJy.x, "y"_a = cJy.y, "z"_a = cJy.z).template cast<FP3>().y;
				this->Jz(i, j, k) = fJ("x"_a = cJz.x, "y"_a = cJz.y, "z"_a = cJz.z).template cast<FP3>().z;
			}
		}

		void setJxyz(int64_t _fJx, int64_t _fJy, int64_t _fJz)
		{
			FP(*fJx)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fJx;
			FP(*fJy)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fJy;
			FP(*fJz)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fJz;
#pragma omp parallel for
			for (int i = 0; i < this->numCells.x; i++)
			for (int j = 0; j < this->numCells.y; j++)
			for (int k = 0; k < this->numCells.z; k++)
			{
				FP3 cJx, cJy, cJz;
                cJx = mapping ? mapping->getDirectCoords(this->JxPosition(i, j, k)) : this->JxPosition(i, j, k);
                cJy = mapping ? mapping->getDirectCoords(this->JyPosition(i, j, k)) : this->JyPosition(i, j, k);
                cJz = mapping ? mapping->getDirectCoords(this->JzPosition(i, j, k)) : this->JzPosition(i, j, k);
                this->Jx(i, j, k) = fJx(cJx.x, cJx.y, cJx.z);
				this->Jy(i, j, k) = fJy(cJy.x, cJy.y, cJy.z);
				this->Jz(i, j, k) = fJz(cJz.x, cJz.y, cJz.z);
			}
		}

		void setJxyzt(int64_t _fJx, int64_t _fJy, int64_t _fJz, FP t)
		{
			FP(*fJx)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fJx;
			FP(*fJy)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fJy;
			FP(*fJz)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fJz;
#pragma omp parallel for
			for (int i = 0; i < this->numCells.x; i++)
			for (int j = 0; j < this->numCells.y; j++)
		    for (int k = 0; k < this->numCells.z; k++)
		    {
		    	FP3 cJx, cJy, cJz;
                         cJx = mapping ? mapping->getDirectCoords(this->JxPosition(i, j, k)) : this->JxPosition(i, j, k);
                         cJy = mapping ? mapping->getDirectCoords(this->JyPosition(i, j, k)) : this->JyPosition(i, j, k);
                         cJz = mapping ? mapping->getDirectCoords(this->JzPosition(i, j, k)) : this->JzPosition(i, j, k);
                         this->Jx(i, j, k) = fJx(cJx.x, cJx.y, cJx.z, t + this->timeShiftJ);
		    	this->Jy(i, j, k) = fJy(cJy.x, cJy.y, cJy.z, t + this->timeShiftJ);
		    	this->Jz(i, j, k) = fJz(cJz.x, cJz.y, cJz.z, t + this->timeShiftJ);
		    }
		}

		void setJ(int64_t _fJ)
		{
			FP3(*fJ)(FP, FP, FP) = (FP3(*)(FP, FP, FP))_fJ;
#pragma omp parallel for
			for (int i = 0; i < this->numCells.x; i++)
			for (int j = 0; j < this->numCells.y; j++)
			for (int k = 0; k < this->numCells.z; k++)
			{
				FP3 cJx, cJy, cJz;
                cJx = mapping ? mapping->getDirectCoords(this->JxPosition(i, j, k)) : this->JxPosition(i, j, k);
                cJy = mapping ? mapping->getDirectCoords(this->JyPosition(i, j, k)) : this->JyPosition(i, j, k);
                cJz = mapping ? mapping->getDirectCoords(this->JzPosition(i, j, k)) : this->JzPosition(i, j, k);
                this->Jx(i, j, k) = fJ(cJx.x, cJx.y, cJx.z).x;
				this->Jy(i, j, k) = fJ(cJy.x, cJy.y, cJy.z).y;
				this->Jz(i, j, k) = fJ(cJz.x, cJz.y, cJz.z).z;
			}
		}
	private:
		int64_t fEt[3],  fBt[3];
		FP globalT;
		bool isAnalytical;
        Mapping* mapping = 0;

	};

    typedef pyGrid<YeeGrid> pyYeeGrid;
    typedef pyGrid<PSTDGrid> pyPSTDGrid;
    typedef pyGrid<PSATDGrid> pyPSATDGrid;
}
