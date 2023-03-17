#pragma once

#include "G4VUserTrackInformation.hh"

#include "MyUserPrimaryParticleInformation.h"


// �����, ���������� �� �������� �������������� ���������������� ���������� � �����
class MyUserTrackInformation : public G4VUserTrackInformation {
    // ������ �������� ����������, ������� ���� �������� ������ � ����������� ����������� ��������
    MyUserPrimaryParticleInformation info;

public:

    MyUserTrackInformation();
    MyUserTrackInformation(const MyUserPrimaryParticleInformation& info);

    void Print() const override;

    MyUserPrimaryParticleInformation GetInfo() const;
};