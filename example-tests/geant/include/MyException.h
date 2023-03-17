#pragma once
#include <string>


class MyException : public std::exception {
protected:
    std::string whatString;
public:
    MyException() {}
    MyException(std::string whatString) : whatString(whatString) {}

    const char* what() const noexcept override {
        return whatString.data();
    }
};


class MyExceptionFileNotExist : public MyException {
public:
    MyExceptionFileNotExist(std::string fileName);
};

class MyExceptionParameterNotSet : public MyException {
public:
    MyExceptionParameterNotSet(std::string parameterName);
};

class MyExceptionParameterInvalid : public MyException {
public:
    MyExceptionParameterInvalid(std::string parameterName, std::string expValue);
};
