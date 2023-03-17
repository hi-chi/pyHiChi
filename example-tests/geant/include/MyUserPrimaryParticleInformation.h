#pragma once

#include "G4VUserPrimaryParticleInformation.hh"

// класс, хранящий пользовательскую информацию об испускаемой генератором частице
class MyUserPrimaryParticleInformation : public G4VUserPrimaryParticleInformation {
    double weight = 0; // factor
    int picID = -1;

public:

    MyUserPrimaryParticleInformation();
    MyUserPrimaryParticleInformation(double weigth, int picID);

    virtual ~MyUserPrimaryParticleInformation();

    virtual void Print() const;

    double GetWeigth() const { return weight; }
    double GetPicID() const { return picID; }
};