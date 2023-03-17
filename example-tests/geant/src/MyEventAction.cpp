#include "MyEventAction.h"
#include "G4Event.hh"

void MyEventAction::BeginOfEventAction(const G4Event* anEvent) {
    //G4cout << "Start event: " << anEvent->GetEventID() << G4endl;
}

void MyEventAction::EndOfEventAction(const G4Event* anEvent) {
    //G4cout << "End event: " << anEvent->GetEventID() << G4endl;
}
