#include "MyException.h"
#include <sstream>

MyExceptionFileNotExist::MyExceptionFileNotExist(std::string fileName) {
    std::stringstream ss;
    ss << "file " << std::string(fileName) << " does not exist" << std::endl;
    this->whatString = ss.str();
}

MyExceptionParameterNotSet::MyExceptionParameterNotSet(std::string parameterName) {
    std::stringstream ss;
    ss << "expected " << parameterName << ". Use \"--help\" to obtain more information" << std::endl;
    this->whatString = ss.str();
}

MyExceptionParameterInvalid::MyExceptionParameterInvalid(std::string parameterName, std::string expValue) {
    std::stringstream ss;
    ss << parameterName << " should be " << expValue << ". Use \"--help\" to obtain more information" << std::endl;
    this->whatString = ss.str();
}
