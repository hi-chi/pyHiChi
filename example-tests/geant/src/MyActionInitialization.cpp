#include "MyActionInitialization.h"
#include "PicParticleContainer.h"

#include "MyPrimaryParticleGeneratorAction.h"
#include "MyTrackingAction.h"
#include "MyRunAction.h"
#include "MySteppingAction.h"
#include "MyEventAction.h"


MyActionInitialization::MyActionInitialization(PicParticleContainer* particles)
 : G4VUserActionInitialization(), particles(particles)
{}


MyActionInitialization::~MyActionInitialization()
{}

void MyActionInitialization::Build() const
{
    // создаем частицы
    SetUserAction(new MyPrimaryParticleGeneratorAction(particles));
    // регистрация потока
    SetUserAction(new MyRunAction());
    // Создание обработчика треков
    SetUserAction(new MyTrackingAction());
    // вывод в файл пошагово
    SetUserAction(new MySteppingAction());
    // вызывается в начале и конце обработки каждой частицы Picador
    SetUserAction(new MyEventAction());
}

void MyActionInitialization::BuildForMaster() const {
     // Просто создаем и регистрируем поток
     SetUserAction(new MyRunAction());
}

