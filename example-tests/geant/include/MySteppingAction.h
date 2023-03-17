#pragma once

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class MySteppingAction : public G4UserSteppingAction {
public:
    MySteppingAction();
    virtual ~MySteppingAction();

    // גלוסעמ processHits
    virtual void UserSteppingAction(const G4Step*);

private:

    void Save(G4int trackID, G4int parentTrackID, G4int picID, G4String ptype,
        const G4ThreeVector& posPart, G4double timePart, const G4ThreeVector& dirPart,
        const G4ThreeVector& momPart, G4double velocity, G4double weight) const;
};

