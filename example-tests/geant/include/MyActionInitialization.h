
#ifndef ExG4ActionInitialization01_h
#define ExG4ActionInitialization01_h 1

#include "G4VUserActionInitialization.hh"

#include "PicParticleContainer.h"

#include "MyPrimaryParticleGeneratorAction.h"

// Обязательный класс, который должен быть объявлен в проекте Geant4
// Имя класса может быть другим, и он должен наследоваться от
// класса G4VUserActionInitialization
class MyActionInitialization : public G4VUserActionInitialization
{
public:
    MyActionInitialization(PicParticleContainer* particles);
    virtual ~MyActionInitialization();
    virtual void Build() const;
    virtual void BuildForMaster() const;

private:
    PicParticleContainer* particles;

};

#endif