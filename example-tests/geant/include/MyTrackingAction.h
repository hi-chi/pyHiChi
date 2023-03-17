#pragma once

#include "G4UserTrackingAction.hh"

// класс, отвечающий за обратотку трека каждой частицы в начале и в конце
class MyTrackingAction : public G4UserTrackingAction {
public:

    MyTrackingAction();
    virtual ~MyTrackingAction();

    // копирует исходную информацию о частице в текущий трек
    virtual void PreUserTrackingAction(const G4Track*);
    // копирует информацию о частице в дочерние частицы
    virtual void PostUserTrackingAction(const G4Track*);
};