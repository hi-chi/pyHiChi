#pragma once

#include "G4VUserTrackInformation.hh"

#include "MyUserPrimaryParticleInformation.h"


// класс, отвечающий за хранение дополнительной пользовательской информации о треке
class MyUserTrackInformation : public G4VUserTrackInformation {
    // храним исходную информацию, которая была получена вместе с испускаемой генератором частицей
    MyUserPrimaryParticleInformation info;

public:

    MyUserTrackInformation();
    MyUserTrackInformation(const MyUserPrimaryParticleInformation& info);

    void Print() const override;

    MyUserPrimaryParticleInformation GetInfo() const;
};