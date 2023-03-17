#pragma once
#include <exception>
#include <string>
#include <vector>

#include "MyException.h"


class ArgParser
{
public:

    enum class Mode {
        NONE,
        GUI,
        CMD
    };

    enum class Status {
        HELP,
        FAILED,
        SUCCESS
    };

    struct ParsedArgs {
        Status status = Status::FAILED;
        Mode mode = Mode::NONE;
        std::vector<std::string> macFiles;
        std::vector<std::string> picParticleFiles;
        int nThreads;
    };

    void ShowHelp();
    ParsedArgs ParseArgs(int argc, char** argv);

protected:

    bool IsFileExist(std::string fileName);
    void ReadFileLine(int argc, char** argv, int& i, std::vector<std::string>& res);
    void SetDefaultParams(ArgParser::ParsedArgs& params);

};

