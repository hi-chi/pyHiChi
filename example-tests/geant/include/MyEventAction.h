#pragma once
#include "G4UserEventAction.hh"

class MyEventAction : public G4UserEventAction {
public:
    MyEventAction() {}

    virtual void BeginOfEventAction(const G4Event* anEvent);
    virtual void EndOfEventAction(const G4Event* anEvent);
};
