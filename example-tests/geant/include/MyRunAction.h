#ifndef RUN_HH
#define RUN_HH

#include "G4UserRunAction.hh"
#include "G4String.hh"


class MyRunAction : public G4UserRunAction
{
	const G4String outFileName = "output";
public:
	MyRunAction();
	~MyRunAction();
	
	virtual void BeginOfRunAction(const G4Run*);
	virtual void EndOfRunAction(const G4Run*);
};

#endif
