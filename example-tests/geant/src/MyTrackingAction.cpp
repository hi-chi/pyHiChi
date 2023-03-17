#include "MyTrackingAction.h"

#include "G4ios.hh"
#include "G4Track.hh"
#include "G4TrackingManager.hh"
#include "G4PrimaryParticle.hh"

#include "MyUserPrimaryParticleInformation.h"
#include "MyUserTrackInformation.h"


MyTrackingAction::MyTrackingAction() {}

MyTrackingAction::~MyTrackingAction() {}

void MyTrackingAction::PreUserTrackingAction(const G4Track* aTrack) {
    // преобразование "начальной" (MyUserPrimaryParticleInformation) информации о треке в "текущую" (MyUserTrackInformation)
    if (aTrack->GetUserInformation() == nullptr && aTrack->GetDynamicParticle()->GetPrimaryParticle()) {
        MyUserPrimaryParticleInformation* startInfo = dynamic_cast<MyUserPrimaryParticleInformation*>(
            aTrack->GetDynamicParticle()->GetPrimaryParticle()->GetUserInformation());
        if (startInfo) {
            aTrack->SetUserInformation(new MyUserTrackInformation(*startInfo));
        }
    }
}

void MyTrackingAction::PostUserTrackingAction(const G4Track* aTrack) {
    // передача информации от родительской частицы в дочерние
    G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
    if (secondaries)
    {
        MyUserTrackInformation* info = dynamic_cast<MyUserTrackInformation*>(aTrack->GetUserInformation());
        if (info) {
            for (size_t i = 0; i < secondaries->size(); i++) {
                (*secondaries)[i]->SetUserInformation(new MyUserTrackInformation(*info));
            }
        }
    }
}
