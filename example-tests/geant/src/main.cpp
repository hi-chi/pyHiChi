#include "ArgParser.h"

#include "PicParticleHandler.h"
#include "PicParticleReader.h"
#include "PicParticleContainer.h"

#include <iostream>
#include <tuple>
#include <string>
#include <sstream>

#include "G4UImanager.hh"

#include "MyGeometry.h"
#include "MyActionInitialization.h"
#include "Shielding.hh"

#include "G4VisManager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif


G4RunManager* createAndInitRunManager(PicParticleContainer* particles, int nThreads) {

#ifdef G4MULTITHREADED
    G4cout << "multithreaded mode: " << nThreads << " threads" << G4endl;

    G4MTRunManager* runManager = new G4MTRunManager;
    runManager->SetNumberOfThreads(nThreads);
#else
    G4cout << "single-threaded mode" << G4endl;

    G4RunManager* runManager = new G4RunManager;
#endif

    // инициализация Geant
    runManager->SetUserInitialization(new MyGeometry());
    runManager->SetUserInitialization(new Shielding());
    runManager->SetUserInitialization(new MyActionInitialization(particles));
    runManager->Initialize();

    return runManager;
}

void runGeantGUI(PicParticleContainer* particles,
    const std::vector<std::string>& macFiles, int nThreads, int argc, char** argv)
{
    G4RunManager* runManager = createAndInitRunManager(particles, nThreads);

    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();

    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    std::for_each(macFiles.begin(), macFiles.end(), [&](const std::string& fileName) {
        UImanager->ApplyCommand("/control/execute " + fileName);
        });

    runManager->BeamOn((int)particles->size());

    ui->SessionStart();
    delete ui;
    delete visManager;

    delete runManager;
}

void runGeantCMD(PicParticleContainer* particles,
    const std::vector<std::string>& macFiles, int nThreads)
{
    G4RunManager* runManager = createAndInitRunManager(particles, nThreads);

    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    std::for_each(macFiles.begin(), macFiles.end(), [&](const std::string& fileName) {
        UImanager->ApplyCommand("/control/execute " + fileName);
        });

    runManager->BeamOn((int)particles->size());

    delete runManager;
}


int main(int argc, char** argv) {
    
    try {

        ArgParser parser;
        ArgParser::ParsedArgs params = parser.ParseArgs(argc, argv);

        if (params.status != ArgParser::Status::SUCCESS)
            return 0;

        // чтение частиц из файла
        PicParticleContainer particles;
        PicParticleReader::ReadFromFile(params.picParticleFiles, particles);

        // в main используются единицы измерения Picador (СГС)
        // здесь должна быть запущена фильтрация частиц из Picador, которые не попадают в область Geant
        // но геометрия такая, что все частицы попадут
        // поэтому следующие строчки запускать не обязательно
        PicParticleHandler particleHandler;
        particleHandler.FilterAndMoveParticles(particles, pfc::FP3(-40, -40, -40), pfc::FP3(40, 40, 40), 1.0);

        // запуск оболочки Geant
        // устанавливаем зерно для генератора
        CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
        CLHEP::HepRandom::setTheSeed(0);

        switch (params.mode) {
        case ArgParser::Mode::GUI:
            runGeantGUI(&particles, params.macFiles, params.nThreads, argc, argv);
            break;
        case ArgParser::Mode::CMD:
            runGeantCMD(&particles, params.macFiles, params.nThreads);
            break;
        default: break;
        }

    }
    catch (std::exception& e) {
        std::cout << "ERROR: An exception is catched in the main() function" << std::endl;
        std::cout << e.what() << std::endl;
    }
    catch (...) {
        std::cout << "ERROR: An undefined exception is catched in the main() function" << std::endl;
    }

    return 0;
}