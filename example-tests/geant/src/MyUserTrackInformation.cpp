#include "MyUserTrackInformation.h"

MyUserTrackInformation::MyUserTrackInformation() {}
MyUserTrackInformation::MyUserTrackInformation(const MyUserPrimaryParticleInformation& info) : info(info) {}

void MyUserTrackInformation::Print() const {
    info.Print();
}

MyUserPrimaryParticleInformation MyUserTrackInformation::GetInfo() const {
    return info;
}