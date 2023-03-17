#include "MyUserPrimaryParticleInformation.h"

#include "G4ios.hh"


MyUserPrimaryParticleInformation::MyUserPrimaryParticleInformation() {}

MyUserPrimaryParticleInformation::MyUserPrimaryParticleInformation(double weight, int picID) :
    weight(weight), picID(picID)
{}

MyUserPrimaryParticleInformation::~MyUserPrimaryParticleInformation() {}

void MyUserPrimaryParticleInformation::Print() const {
    // ������� ������ �������
    G4cout << "  Weight: " << this->weight;
};