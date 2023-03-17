#include "ArgParser.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <G4Threading.hh>


bool ArgParser::IsFileExist(std::string fileName) {
    std::ifstream file(fileName);
    if (!file.is_open()) return false;
    file.close();
    return true;
}

void ArgParser::ShowHelp() {
    std::cout << "TestGeant [-h] [-v files_vis.mac or -r files_run.mac] [-p files_of_picador_output.txt] [-nt number_of_threads]" << std::endl;
    std::cout << "-h, --help: show help" << std::endl;
    std::cout << "-v files_vis.mac: start GUI mode and run commands from files_vis.mac" << std::endl;
    std::cout << "-r files_run.mac: start command line mode and run commands from files_run.mac" << std::endl;
    std::cout << "-p files_output.txt: list of Hi-Chi particles" << std::endl;
    std::cout << "-nt number_of_threads: number of threads in a multithreaded mode" << std::endl;
    std::cout << "Example 1: \"TestGeant -r run.mac -p hichi_output1.txt hichi_output2.txt -nt 4\"" << std::endl;
    std::cout << "Example 2: \"TestGeant -v init_vis.mac run.mac -p hichi_output.txt" << std::endl;
    std::cout << "Note: in the single-threaded mode the -nt parameter is ignored" << std::endl;
}

void ArgParser::ReadFileLine(int argc, char** argv, int& i, std::vector<std::string>& res) {
    while (i + 1 < argc && argv[i + 1][0] != '-') {
        i++;
        if (!IsFileExist(std::string(argv[i])))
            throw MyExceptionFileNotExist(std::string(argv[i]));
        res.push_back(std::string(argv[i]));
    }
}

void ArgParser::SetDefaultParams(ArgParser::ParsedArgs& params) {
#ifdef G4MULTITHREADED
    params.nThreads = G4Threading::G4GetNumberOfCores();
#else
    params.nThreads = 1;
#endif
}

ArgParser::ParsedArgs ArgParser::ParseArgs(int argc, char** argv)
{
    ArgParser::ParsedArgs params;
    SetDefaultParams(params);

    try {
        int i = 1;
        while (i < argc) {
            if (std::string(argv[i]) == "-h" || std::string(argv[i]) == "--help") {
                ShowHelp();
                params.status = ArgParser::Status::HELP;
                return params;
            }
            else if (std::string(argv[i]) == "-r" && params.mode != ArgParser::Mode::GUI) {
                params.mode = ArgParser::Mode::CMD;
                ReadFileLine(argc, argv, i, params.macFiles);
            }
            else if (std::string(argv[i]) == "-v") {
                params.mode = ArgParser::Mode::GUI;
                ReadFileLine(argc, argv, i, params.macFiles);
            }
            else if (std::string(argv[i]) == "-p") {
                ReadFileLine(argc, argv, i, params.picParticleFiles);
            }
            else if (std::string(argv[i]) == "-nt") {
                if (++i >= argc) throw MyExceptionParameterNotSet("number of threads");
                params.nThreads = std::atoi(argv[i]);
                if (params.nThreads < 1) throw MyExceptionParameterInvalid("number of threads", "1 or more");
            }
            i++;
        }
        if (params.mode == ArgParser::Mode::NONE)
            throw MyExceptionParameterNotSet("\"-v\" or \"-r\"");
        if (params.mode == ArgParser::Mode::CMD && params.macFiles.size() < 1)
            throw MyExceptionParameterInvalid("number of files_run", "1 or more");
        if (params.picParticleFiles.size() < 1)
            throw MyExceptionParameterInvalid("number of files_of_picador_output", "1 or more");

        params.status = ArgParser::Status::SUCCESS;
    }
    catch (std::exception& e) {
        std::cout << "ERROR: " << e.what() << std::endl;
        params.status = ArgParser::Status::FAILED;
    }

    return params;
}
