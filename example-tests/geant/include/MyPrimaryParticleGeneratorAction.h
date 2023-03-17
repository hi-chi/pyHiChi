
#ifndef ExG4PrimaryGeneratorAction_h
#define ExG4PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

#include "MyParticleGun.h"
//#include "G4ParticleGun.hh"
//using MyParticleGun = G4ParticleGun;

#include "PicParticleContainer.h"

class G4ParticleGun;
class G4Event;
class G4Box;

/// Класс определения источника первичных частиц
class MyPrimaryParticleGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    MyPrimaryParticleGeneratorAction(PicParticleContainer* particles);
    virtual ~MyPrimaryParticleGeneratorAction();

    // Метод из базового класса, задает параметры источника начальных частиц
    virtual void GeneratePrimaries(G4Event*);

private:
    MyParticleGun* CreatePrimary(const PicParticle* picParticle);

    PicParticleContainer* particles;
    std::vector<MyParticleGun*> fParticleGuns; //указатель на источники частиц
};
#endif
