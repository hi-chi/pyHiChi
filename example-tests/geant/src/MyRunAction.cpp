#include "MyRunAction.h"
#include "MyAnalysisManager.h"

MyRunAction::MyRunAction()
{}

MyRunAction::~MyRunAction()
{}

void MyRunAction::BeginOfRunAction(const G4Run*)
{
	MyAnalysisManager*man = MyAnalysisManager::Instance();
	man->SetNtupleMerging(true);
	man->OpenFile(outFileName);

	man->CreateNtuple("Hits","Hits");
	man->CreateNtupleIColumn("fTrackID");
	man->CreateNtupleIColumn("fParentID");
	man->CreateNtupleIColumn("fPicID");
	man->CreateNtupleSColumn("fType");
	man->CreateNtupleDColumn("fX");
	man->CreateNtupleDColumn("fY");
	man->CreateNtupleDColumn("fZ");
	man->CreateNtupleDColumn("fTime");
	man->CreateNtupleDColumn("fDirX");
	man->CreateNtupleDColumn("fDirY");
	man->CreateNtupleDColumn("fDirZ");
	man->CreateNtupleDColumn("fPx");
	man->CreateNtupleDColumn("fPy");
	man->CreateNtupleDColumn("fPz");
	man->CreateNtupleDColumn("fV");
	man->CreateNtupleDColumn("fW");
	man->FinishNtuple();
}

void MyRunAction::EndOfRunAction(const G4Run*)
{
	MyAnalysisManager* man = MyAnalysisManager::Instance();

	man->Write();
	man->CloseFile();
}
