#include "MySteppingAction.h"

#include "G4ThreeVector.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4PrimaryParticle.hh"

#include "MyAnalysisManager.h"
#include "MyUserTrackInformation.h"


MySteppingAction::MySteppingAction()
{
}

MySteppingAction::~MySteppingAction()
{
}

void MySteppingAction::UserSteppingAction(const G4Step* aStep)
{
	G4Track* track = aStep->GetTrack();

	if (track->GetVolume()->GetName() == "Detector") {
		G4String ptype = track->GetDefinition()->GetParticleName();

		G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
		G4ThreeVector posPart = preStepPoint->GetPosition();
		G4ThreeVector momPart = preStepPoint->GetMomentum();
		G4ThreeVector dirPart = preStepPoint->GetMomentumDirection();
		G4double timePart = preStepPoint->GetGlobalTime();
		G4double vPart = preStepPoint->GetVelocity();

		MyUserTrackInformation* userInfo = dynamic_cast<MyUserTrackInformation*>(aStep->GetTrack()->GetUserInformation());
		G4double weight = 0;
		G4int picID = -1;
		if (userInfo) {
			weight = userInfo->GetInfo().GetWeigth();
			picID = userInfo->GetInfo().GetPicID();
		}

		G4int trackID = track->GetTrackID();
		G4int parentID = track->GetParentID();

		Save(trackID, parentID, picID, ptype, posPart, timePart, dirPart, momPart, vPart, weight);

		track->SetTrackStatus(G4TrackStatus::fStopAndKill);
	}
}

void MySteppingAction::Save(G4int trackID, G4int parentTrackID, G4int picID, G4String ptype,
	const G4ThreeVector& posPart, G4double timePart, const G4ThreeVector& dirPart,
	const G4ThreeVector& momPart, G4double velocity, G4double weight) const
{
	MyAnalysisManager* man = MyAnalysisManager::Instance();
	man->FillNtupleIColumn(0, trackID);
	man->FillNtupleIColumn(1, parentTrackID);
	man->FillNtupleIColumn(2, picID);
	man->FillNtupleSColumn(3, ptype);
	man->FillNtupleDColumn(4, posPart[0]);
	man->FillNtupleDColumn(5, posPart[1]);
	man->FillNtupleDColumn(6, posPart[2]);
	man->FillNtupleDColumn(7, timePart);
	man->FillNtupleDColumn(8, dirPart[0]);
	man->FillNtupleDColumn(9, dirPart[1]);
	man->FillNtupleDColumn(10, dirPart[2]);
	man->FillNtupleDColumn(11, momPart[0]);
	man->FillNtupleDColumn(12, momPart[1]);
	man->FillNtupleDColumn(13, momPart[2]);
	man->FillNtupleDColumn(14, velocity);
	man->FillNtupleDColumn(15, weight);
	man->AddNtupleRow();
}
