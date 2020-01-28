#pragma once
#include "Vectors.h"
#include "Mapping.h"
#include "Constants.h"
#include <functional>

namespace pfc {

    class TightFocusingMapping : public PeriodicalXMapping {

    public:

		TightFocusingMapping(FP R0, FP L, FP D, FP cutAngle = 0.5*constants::pi) :
			PeriodicalXMapping(-R0 - D + 0.5*L, -R0 + 0.5*L),
            xL(-R0 - 0.5*L), cutAngle(cutAngle), ifCut(true) {}

		void setIfCut(bool ifCut = true) {
			this->ifCut = ifCut;
		}

        FP3 getDirectCoords(const FP3& coords, bool* status = 0) {
            FP ct = constants::c*time;
            FP r = coords.norm();

            FP3 directCoords = coords;
            if (ifCut) setFailStatus(status);
			else setOkStatus(status);

            if (coords.x < xMin + ct || coords.x >= xMax + ct)
                return coords;

            int nPeriods = 0;
            int shiftSign = 0;

            if (xMax + ct < 0) {
                nPeriods = int(((xMax + ct)*cos(cutAngle) - (xMin + ct)) / D) + 1;  // целая часть сверху
                shiftSign = 1;
            }
            else if (xMin + ct > 0) {
                nPeriods = int(((xMax + ct) - (xMin + ct)*cos(cutAngle)) / D) + 1;  // целая часть сверху
                shiftSign = -1;
            }
            else nPeriods = 1;

            for (int i = 0; i < nPeriods; i++) {
                FP shift = i * D * shiftSign;

                FP3 coordsShift(coords.x + shift, coords.y, coords.z);
                FP rShift = coordsShift.norm();

                if (ifInArea(coordsShift)) {
                    directCoords = coordsShift;
                    setOkStatus(status);
                    break;
                }
            }

            return directCoords;
        }

        FP3 getInverseCoords(const FP3& coords, bool* status = 0) {
            
			setOkStatus(status);

            if (!ifInArea(coords) && ifCut) {
				setFailStatus(status);
            }

            return PeriodicalXMapping::getInverseCoords(coords);

        }

        FP getxMin() const { return xMin; }
        FP getxMax() const { return xMax; }

        bool ifInArea(const FP3& coords) {
            FP ct = constants::c*time;
            FP r = coords.norm();
            FP angle = atan(abs(sqrt(coords.y*coords.y + coords.z*coords.z)/coords.x));

            if (xMax + ct < 0) {
                if ((r >= -xL - ct) || (r < -xMax - ct) || (coords.x > 0))
                    return false;
                if (angle > cutAngle)
                    return false;
            }
            else if ((xMax + ct >= 0) && (xL + ct <= 0))
            {
                if (coords.x < 0 && r > xMax - xL)
                    return false;
                if (coords.x >= 0 && r > xMax + ct)
                    return false;
            }
            else if (xL + ct > 0)
            {
                if ((r <= xL + ct) || (r > xMax + ct) || (coords.x < 0))
                    return false;
                if (angle > cutAngle)
                    return false;
            }

            return true;
        }

        FP cutAngle, xL;
		bool ifCut = true;

    };

}