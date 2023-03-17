
#ifndef ExG4DetectorConstruction01_h
#define ExG4DetectorConstruction01_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4UserLimits.hh"
#include "G4Material.hh"

#include "PicParticle.h"

class G4VPhysicalVolume;
class G4LogicalVolume;


//  Класс геометрии установки, объявление материалов и детекторов
class MyGeometry : public G4VUserDetectorConstruction
{
public:

    MyGeometry();
    virtual ~MyGeometry();

    virtual G4VPhysicalVolume* Construct();

    virtual void ConstructSDandField();
};
#endif
