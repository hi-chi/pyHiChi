#include "MyGeometry.h"

#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4UserLimits.hh"

#include "MySensitiveDetector.h"

MyGeometry::MyGeometry() : G4VUserDetectorConstruction() {
}

MyGeometry::~MyGeometry() {
}

G4VPhysicalVolume* MyGeometry::Construct()
{
    G4NistManager* nist = G4NistManager::Instance();
    G4Material* worldMat = nist->FindOrBuildMaterial("G4_Galactic");
    G4Material* detMat = nist->FindOrBuildMaterial("G4_Galactic");

    G4Box* solidWorld = new G4Box("World", 0.5 * m, 0.5 * m, 0.5 * m);

    G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, worldMat, "World");

    G4VPhysicalVolume* physWorld = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.),
        logicWorld, "World", 0, false, 0, true);

    G4Isotope* isotopeRad = new G4Isotope("Au197", 79, 197);

    G4Element* elementRad = new G4Element("Gold", "Au", 1);
    elementRad->AddIsotope(isotopeRad, 100. * perCent);

    G4Material* materialRad = new G4Material("Gold", 19.32 * g / cm3, 1);
    materialRad->AddElement(elementRad, 100. * perCent);

    G4Sphere* solidRadiator = new G4Sphere("Radiator", 0.2 * m, 0.23 * m, 0, 360, 0, 360);
    //G4Box *solidRadiator = new G4Box("solidRadiator", 0.2*m, 0.2*m, 0.1*m);
    G4LogicalVolume* logicRadiator = new G4LogicalVolume(solidRadiator, materialRad, "Radiator");
    //fScoringVolume = logicRadiator;
    G4VPhysicalVolume* physRadiator = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicRadiator,
        "Radiator", logicWorld, false, 0, true);

    G4Sphere* solidDetector = new G4Sphere("Detector", 0.48 * m, 0.49 * m, 0, 360, 0, 360);
    G4LogicalVolume* logicDetector = new G4LogicalVolume(solidDetector, detMat, "Detector");
    G4VPhysicalVolume* physDetector = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicDetector,
        "Detector", logicWorld, false, 0, true);

    return physWorld;
}

void MyGeometry::ConstructSDandField()
{
    MySensitiveDetector* sensDet = new MySensitiveDetector("Detector");
    G4SDManager::GetSDMpointer()->AddNewDetector(sensDet);
    SetSensitiveDetector("Detector", sensDet, true);
}
