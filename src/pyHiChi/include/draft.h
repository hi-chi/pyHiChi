//This is a temporary(!!!) immplementation (draft) for the structure of classes that handle fields.

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

#include <iostream>
#include <cmath>
#include <windows.h>

using namespace std;

struct fieldIterator
{
	virtual void begin() = 0;
	virtual bool next(FP3& position) = 0;
	virtual void setField(const FP3 E, const FP3 B) = 0;
};
struct pyFieldIterator //  wrapper for Python
{
	fieldIterator* iterator;
	pyFieldIterator(fieldIterator* _iterator) : iterator(_iterator)
	{}
	void begin() { iterator->begin(); }
	bool next(FP3& position) { return iterator->next(position); }
	void setField(const FP3 E, const FP3 B) { iterator->setField(E, B); }
	inline pyFieldIterator& operator= (pyFieldIterator& v)
	{
		iterator = v.iterator;
		return *this;
	}
	void selectPositions(int64_t condition);
	void thinOut(int interval);
};

class fieldBase // abstract base class for fields (can be either transformedField or newField) 
{
public:
	virtual void get(FP3 const coords, FP3 &E, FP3 &B) = 0;
	virtual void advance() = 0;
	virtual FP getTime() = 0;
	virtual FP* getTimeReferece() = 0; // for movingwindow, should be possible to remove by managing type transformation in python
	virtual fieldIterator* getIterator() = 0;
	virtual void resetTime(FP time) = 0;
};

class pyField // class 'field' for Python, a wrapper for fieldBase (this is to implement virtual functions in Python)
{
	fieldBase* _fieldBase;
public:
	double advanceTimeSpent; // counter for the time spent on advance()
	double getAdvanceTimeSpent() { return advanceTimeSpent; }
	pyField()
	{
		advanceTimeSpent = 0;
	}
	pyField(fieldBase& __fieldBase) : _fieldBase(&__fieldBase)
	{
		advanceTimeSpent = 0;
	}
	virtual void get(FP3 const coords, FP3 &E, FP3 &B)
	{
		_fieldBase->get(coords, E, B);
	}
	virtual void advance()
	{
		LARGE_INTEGER frequency, t1, t2;
		QueryPerformanceFrequency(&frequency);
		QueryPerformanceCounter(&t1);

		_fieldBase->advance();

		QueryPerformanceCounter(&t2);
		advanceTimeSpent += (t2.QuadPart - t1.QuadPart) / double(frequency.QuadPart);
	}
	virtual FP getTime()
	{
		return _fieldBase->getTime();
	}
	inline pyField& operator= (pyField& v)
	{
		_fieldBase = v._fieldBase;
		return *this;
	}
	FP* getTimeReferece()
	{
		return _fieldBase->getTimeReferece();
	}
	pyFieldIterator getIterator()
	{
		return pyFieldIterator(_fieldBase->getIterator());
	}
	void resetTime(FP time) { _fieldBase->resetTime(time); }
	void interpolateFrom(pyField fromField);
	friend class transformedField;
	friend class newField;
};

class fieldDraft // class for the draft of a field, to be further used to initiate a 'field' via init(time, timeStep) in Python
{
	bool initiated; // this is to restrict to only one call of init()
public:
	fieldDraft() : initiated(false)
	{}
	FP _time, _timeStep; // !!! either move time and timeStep to the fieldDraft or prohibit a reassignment
	virtual void get(FP3 const coords, FP3 &E, FP3 &B) = 0;
	virtual void advance() = 0;
	virtual fieldIterator* getIterator() = 0;
	virtual pyField init(FP time, FP timeStep) = 0; // can be used to check that the drafted field is well-defined and can be initialized for further use (calling advance() and get()).
													// to amend: we should not give time to fieldDarf childs, since time is controlled here.
	pyField internalInit(FP time, FP timeStep); // internal function to init
	void resetTime(FP time) { _time = time; }
};

class transform_ // base class for transforms
{
public:
	virtual pyField operator()(pyField& field) = 0; // It would be better to define it here but I do not know how to make the function of parent class available for its childs in Python.
	virtual void directCoords(FP3& vect) = 0;
	virtual bool inverseCoords(FP3& vect) = 0; // returns false if converted position is outside the ranges
	virtual void directField(FP3& vect) = 0;
	virtual void inverseField(FP3& vect) = 0;
	//we can define other types of transforms for pseudoVectors, spinors, tensors, pseudoTensors and so on.  
};

struct transformedFieldIterator : public fieldIterator
{
	fieldIterator* _iterator;
	reference_wrapper<transform_> _transform;
	transformedFieldIterator(transform_& _transform_, fieldIterator* _iterator_) : _transform(_transform_), _iterator(_iterator_)
	{}
	void begin() { _iterator->begin(); }
	bool next(FP3& position)
	{
		bool _next = _iterator->next(position);
		_transform.get().directCoords(position);
		return _next;
	}
	void setField(const FP3 E, const FP3 B)
	{
		FP3 Et(E), Bt(B);
		_transform.get().inverseField(Et);
		_transform.get().inverseField(Bt); //can be optimized by calling only requested entity
		_iterator->setField(Et, Bt);
	}
};

class transformedField : public fieldBase // class for transformed field (wrapped with a given transform)
{
public:
	fieldBase* _fieldBase;
	reference_wrapper<transform_> _transform;
	transformedField(transform_& _transform_, pyField& field) : _transform(_transform_), _fieldBase(field._fieldBase)
	{}
	virtual void get(FP3 const coords, FP3 &E, FP3 &B)
	{
		FP3 _coords(coords);
		bool insideRanges = _transform.get().inverseCoords(_coords);
		if (insideRanges)
		{
			_fieldBase->get(_coords, E, B);
			_transform.get().directField(E);
			_transform.get().directField(B);
		}
		else
		{
			E = FP3(0, 0, 0);
			B = FP3(0, 0, 0);
		}
	}
	void advance()
	{
		_fieldBase->advance();
	}
	FP getTime()
	{
		return _fieldBase->getTime();
	}
	FP* getTimeReferece()
	{
		return _fieldBase->getTimeReferece();
	}
	fieldIterator* getIterator()
	{
		return new transformedFieldIterator(_transform.get(), _fieldBase->getIterator());
	}
	void resetTime(FP time) { _fieldBase->resetTime(time); }
};

class newField : public fieldBase // newly initiated field, with time controller (see advance())
{
	fieldDraft* _fieldDraft;
public:
	newField(fieldDraft& field) : _fieldDraft(&field)
	{}
	void get(FP3 const coords, FP3 &E, FP3 &B)
	{
		_fieldDraft->get(coords, E, B);
	}
	void advance()
	{
		_fieldDraft->advance();
		_fieldDraft->_time += _fieldDraft->_timeStep;
	}
	FP getTime()
	{
		return _fieldDraft->_time;
	}
	FP* getTimeReferece()
	{
		return &(_fieldDraft->_time);
	}
	fieldIterator* getIterator()
	{
		return _fieldDraft->getIterator();
	}
	void resetTime(FP time) { _fieldDraft->resetTime(time); }
};

pyField fieldDraft::internalInit(FP time, FP timeStep)
{
	if (initiated)
	{
		std::cout << "fieldDraft: Error: second call of init() is prohibited (init() binds the draft to the only field instance)." << endl;
		exit(0);
	}
	initiated = true;
	_time = time;
	_timeStep = timeStep;
	newField* _newField = new newField(*this);
	return pyField(*_newField);
};


class analyticalField : public fieldDraft
{
	int64_t _funcEx, _funcEy, _funcEz, _funcBx, _funcBy, _funcBz; // analytical functions for all field components
public:
	analyticalField(int64_t funcEx, int64_t funcEy, int64_t funcEz, int64_t funcBx, int64_t funcBy, int64_t funcBz)
	{
		_funcEx = funcEx;
		_funcEy = funcEy;
		_funcEz = funcEz;
		_funcBx = funcBx;
		_funcBy = funcBy;
		_funcBz = funcBz;
	}
	void get(FP3 const coords, FP3 &E, FP3 &B)
	{
		FP(*ex)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_funcEx; // Does this take acceptibally short time?
		FP(*ey)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_funcEy;
		FP(*ez)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_funcEz;
		FP(*bx)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_funcBx;
		FP(*by)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_funcBy;
		FP(*bz)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_funcBz;
		E.x = ex(coords.x, coords.y, coords.z, _time);
		E.y = ey(coords.x, coords.y, coords.z, _time);
		E.z = ez(coords.x, coords.y, coords.z, _time);
		B.x = bx(coords.x, coords.y, coords.z, _time);
		B.y = by(coords.x, coords.y, coords.z, _time);
		B.z = bz(coords.x, coords.y, coords.z, _time);
	}
	void advance()
	{
		std::cout << "analyticalField: advance: do nothing" << std::endl;
	}
	pyField init(FP time, FP timeStep)
	{
		// here we can perform all necessary preparation and validations before issuing the 'field' for further use in HiChi
		return internalInit(time, timeStep);
	}
	fieldIterator* getIterator()
	{
		cout << "analyticalField: Error: calling getIterator() is inappropriate." << endl;
		exit(0);
		return nullptr;
	}
};

class rotateZ : public transform_ // rotation about z axis
{
	FP cos_angle, sin_angle;
public:
	rotateZ(FP angle)
	{
		cos_angle = cos(angle);
		sin_angle = sin(angle);
	}
	pyField operator()(pyField& field)
	{
		transformedField* v = new transformedField(*this, field);
		return pyField(*v);
	}
	void directCoords(FP3& vect)
	{
		FP3 out;
		out.x = vect.x * cos_angle - vect.y * sin_angle;
		out.y = vect.x * sin_angle + vect.y * cos_angle;
		out.z = vect.z;
		vect = out;
	}
	bool inverseCoords(FP3& vect)
	{
		FP3 out;
		out.x = vect.x * cos_angle + vect.y * sin_angle;
		out.y = -vect.x * sin_angle + vect.y * cos_angle;
		out.z = vect.z;
		vect = out;
		return true;
	}
	void directField(FP3& vect)
	{
		directCoords(vect);
	}
	void inverseField(FP3& vect)
	{
		inverseCoords(vect);
	}
	//we can define other types of transforms for pseudoVectors, spinors, tensors, pseudoTensors and so on.  
};

class movingWindowX : public transform_ // mimicks window moving in positive x direction, requires field solver with periodic boundaries along this axis. 
{
	FP _xMin, _xMax, invLx;
	FP* _time; // here
public:
	/*movingWindowX(FP xMin, FP xMax, fieldDraft& fieldToSyncWith) : _xMin(xMin), _xMax(xMax)
	{
		_time = &(fieldToSyncWith._time);
		invLx = 1 / (_xMax - _xMin);
	}*/
	movingWindowX(FP xMin, FP xMax, pyField& fieldToSyncWith) : _xMin(xMin), _xMax(xMax)
	{
		_time = fieldToSyncWith.getTimeReferece();
		invLx = 1 / (_xMax - _xMin);
	}
	pyField operator()(pyField& field)
	{
		transformedField* v = new transformedField(*this, field);
		return pyField(*v);
	}
	void directCoords(FP3& vect)
	{
		FP d = (constants::c*(*_time) + _xMax - vect.x) / (_xMax - _xMin);
		vect.x += (_xMax - _xMin)*(int(d) - (d < 0));
	}
	bool inverseCoords(FP3& vect)
	{
		if ((vect.x - constants::c*(*_time) < _xMin) || (vect.x - constants::c*(*_time) > _xMax)) return false;
		double tmp;
		vect.x = _xMin + modf((vect.x - _xMin)*invLx, &tmp)*(_xMax - _xMin);
		return true;
	}
	void directField(FP3& vect)
	{
		// do nothing
	}
	void inverseField(FP3& vect)
	{
		// do nothing
	}
};
class movingWindowCX : public transform_ // mimicks window moving in positive x direction, requires field solver with periodic boundaries along this axis. 
{
	FP _xMin, _xMax, _xD, invLx;
	FP* _time; // here
public:
	/*movingWindowX(FP xMin, FP xMax, fieldDraft& fieldToSyncWith) : _xMin(xMin), _xMax(xMax)
	{
		_time = &(fieldToSyncWith._time);
		invLx = 1 / (_xMax - _xMin);
	}*/
	movingWindowCX(FP xMin, FP xMax, FP xD, pyField& fieldToSyncWith) : _xMin(xMin), _xMax(xMax), _xD(xD)
	{
		_time = fieldToSyncWith.getTimeReferece();
		invLx = 1 / (_xMax - _xMin);
	}
	pyField operator()(pyField& field)
	{
		transformedField* v = new transformedField(*this, field);
		return pyField(*v);
	}
	void directCoords(FP3& vect) // will need to rethink and optimize this
	{
		FP t = *_time;
		bool None = false;
		bool C = false;
		int i = 0;

		for (; !C; i++)
		{
			FP3 vv = vect - FP3(-i * (_xMax - _xMin), 0, 0);
			FP r = vv.norm();
			if (_xMax + constants::c*t < 0)
				C = ((r >= -_xMax - constants::c*t) && (r < -_xD - constants::c*t) && (vv.x <= 0));

			if ((_xMax + constants::c*t >= 0) && (_xD + constants::c*t <= 0))
			{
				FP t1 = t - (-_xMax) / constants::c;
				if (vv.x <= 0) C = (r < (_xMax - _xD));
				else C = (r < t1*constants::c);
			}
			if (_xD + constants::c*t > 0)
			{
				FP t2 = t - (-_xD) / constants::c;
				C = ((r >= constants::c*t2) && (r < constants::c*t2 + (_xMax - _xD)) && (vv.x >= 0));
			}
			if (i == 32) { None = true; C = true; } //Will need to fix this later.
		}
		vect.x += (_xMax - _xMin)*(i - 1);
		if (None) vect.x += 1e+10;//Will need to fix this later.
	}
	bool inverseCoords(FP3& vect)
	{
		FP t = *_time;
		FP r = vect.norm();
		if (_xMax + constants::c*t < 0)
			if ((r >= -_xD - constants::c*t) || (r < -_xMax - constants::c*t) || (vect.x > 0)) return false;

		if ((_xMax + constants::c*t >= 0) && (_xD + constants::c*t <= 0))
		{
			if (vect.x < 0) if (r > (_xMax - _xD)) return false;
			if (vect.x >= 0) if (r > constants::c*(t - (-_xMax) / constants::c)) return false;
		}

		if (_xD + constants::c*t > 0)
		{
			FP t1 = t - (-_xD) / constants::c;
			if ((r < t1*constants::c) || (r > t1*constants::c + (_xMax - _xD)) || (vect.x < 0)) return false;
		}


		//if()
		//{
		//	FP t1 = t - -_xMax/constants::c
		//	//if ((vect.x < _xD + constants::c*t) || (vect.x > _xMax - constants::c*t)) return false;
		//	if ((r >= -_xD - constants::c*t) || (r < -_xMax - constants::c*t) || (vect.x > 0)) return false;
		//}

		double tmp;
		vect.x = _xMin + modf((vect.x - _xMin)*invLx, &tmp)*(_xMax - _xMin);
		return true;
	}
	void directField(FP3& vect)
	{
		// do nothing
	}
	void inverseField(FP3& vect)
	{
		// do nothing
	}
};

//----------fieldIterator options

struct fieldIteratorSelectPositions : public fieldIterator
{
	fieldIterator* iterator;
	int64_t condition;
	fieldIteratorSelectPositions(int64_t _condition, fieldIterator* _iterator) : iterator(_iterator)
	{
		condition = _condition;
	}
	void begin() { iterator->begin(); }
	bool next(FP3& position)
	{
		FP(*cond)(FP, FP, FP) = (FP(*)(FP, FP, FP))condition; // Does this take acceptibally short time?
		while (iterator->next(position))
		{
			if (cond(position.x, position.y, position.z) >= 0) return true;
		}
		return false;
	}
	void setField(const FP3 E, const FP3 B) { iterator->setField(E, B); }
};

void pyFieldIterator::selectPositions(int64_t condition)
{
	iterator = new fieldIteratorSelectPositions(condition, this->iterator);
}

struct fieldIteratorThinOut : public fieldIterator
{
	fieldIterator* iterator;
	int interval;
	fieldIteratorThinOut(int _interval, fieldIterator* _iterator) : iterator(_iterator), interval(_interval)
	{}
	void begin() { iterator->begin(); }
	bool next(FP3& position)
	{
		int k = 0;
		while (iterator->next(position))
		{
			k++;
			if (k == interval) return true;
		}
		return false;
	}
	void setField(const FP3 E, const FP3 B) { iterator->setField(E, B); }
};

void pyFieldIterator::thinOut(int interval)
{
	iterator = new fieldIteratorThinOut(interval, this->iterator);
}

void pyField::interpolateFrom(pyField fromField)
{
	fieldIterator* it = _fieldBase->getIterator();
	FP3 position, E, B;
	it->begin();
	int counter = 0;
	while (it->next(position))
	{
		fromField.get(position, E, B);
		it->setField(E, B);
		counter++;
		if (counter % (1024 * 1024) == 0) cout << "#";
	}
}

void arTest()
{
	cout << "arTest()" << endl;
}

struct draft
{
	FP _version;
public:
	draft()
	{
		_version = 0.01;
	}
	FP version()
	{
		return _version;
	}
};
