#pragma once

#include "G4VSensitiveDetector.hh"

class G4Step;
class G4HCofThisEvent;

// класс для обработки треков в процессе моделирования
class MySensitiveDetector : public G4VSensitiveDetector
{
public:
    MySensitiveDetector(const G4String& name);
    ~MySensitiveDetector() override;

    G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override {
        return true;
    }
};

