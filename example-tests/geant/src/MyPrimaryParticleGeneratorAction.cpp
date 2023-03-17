#include "MyPrimaryParticleGeneratorAction.h"

#include "MyUserPrimaryParticleInformation.h"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

using namespace pfc;

// Класс, в котором описывается положение, тип, энергия, направление вылета
// и распределение начальных частиц
MyPrimaryParticleGeneratorAction::MyPrimaryParticleGeneratorAction(
    PicParticleContainer* particles)
    : G4VUserPrimaryGeneratorAction(), fParticleGuns(particles->size(), 0),
    particles(particles)
{}

// Деструктор, удаляем созданный в конструкторе экземпляр класса источника G4ParticleGun
MyPrimaryParticleGeneratorAction::~MyPrimaryParticleGeneratorAction()
{}

// Эта функция вызывается в начале каждого события
// каждое событие отвечает за генерацию и обработку соответствующей частицы
void MyPrimaryParticleGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    G4int index = anEvent->GetEventID();
    const PicParticle* p = (*particles)[index];
    // создаем источник
    fParticleGuns[index] = CreatePrimary(p);
    // регистрируем источник
    fParticleGuns[index]->GeneratePrimaryVertex(anEvent);
    // прикрепляем к излучаемой частице некоторую информацию об исходной частице
    anEvent->GetPrimaryVertex()->GetPrimary()->SetUserInformation(
        new MyUserPrimaryParticleInformation(p->factor, p->id));
}

MyParticleGun* MyPrimaryParticleGeneratorAction::CreatePrimary(const PicParticle* picParticle)
{
    MyParticleGun* fParticleGun = new MyParticleGun();

    // Создаем частицу
    G4ParticleDefinition* particleDef = picParticle->GetParticleDefinition();
    // Устанавливаем полученную частицу в качестве испускаемого типа начальных частиц в источнике
    fParticleGun->SetParticleDefinition(particleDef);
    // Устанавливаем направление движение частицы по (x,y,z)
    FP3 direction = picParticle->GetDirection();
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(
        direction.x, direction.y, direction.z
    ));
    // Установка начальной энергии испускаемых частиц
    FP energy = picParticle->GetEnergy();
    fParticleGun->SetParticleEnergy(energy * erg);
    // Задаем координаты, где будет источник начальных частиц
    G4double x0 = picParticle->position.x * cm;
    G4double y0 = picParticle->position.y * cm;
    G4double z0 = picParticle->position.z * cm;
    // Устанавливаем позицию источника начальных частиц
    fParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
    // Устанавливаем время испускания частицы
    fParticleGun->SetParticleTime(picParticle->emissionTime * second);

    return fParticleGun;
}
