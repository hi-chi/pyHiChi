#pragma once
#include <string>
#include "Grid.h"

#include <adios2.h>
namespace pfc {

#define put_var(x) PutVariable((#x), (x));
#define put_arr(x) PutArray((#x), (x));
#define get_arr(x) GetArray((#x), (x));
#define get_arr(x) GetArray((#x), (x));
class DataManager
{
    std::string path; //откуда брать данные и куда их писать
    int iteration = 0;
    enum Mode { read, write, append };
    Mode mode;
    adios2::IO bpIO;
public:
    //DataManager() {}
    DataManager(std::string path, Mode mode = Mode::write, int it = 0): path(path), mode(mode), iteration(it)
    {
        adios2::ADIOS adios;
        bpIO = adios.DeclareIO(path);
		mode = Mode::read;
    }
    void BeginIteration() //create new file
    {
		PutVariable(path, path);
    }
    void EndIteration() //close file
    {

    }
    void BeginStep() //logical step 
    {

    }
    void EndStep() //flush buffer
    {

    }

    template<typename T>
    void PutVariable(std::string name, const T var)
    {

    }
    template<typename T>
    void PutArray(std::string name, const T* arr, size_t size)
    {

    }
    //void PutGrid(int *arr, size_t size, std::string name)
    //{
    //
    //}

    template<typename T>
    void GetVariable(std::string name, T& var)
    {

    }
    template<typename T>
    void GetArray(std::string name, T* arr, size_t size)
    {

    }

    template<typename T>
    void PutGrid(const T &grid)
    {

    }
    template<typename T>
    void getGrid(T &grid)
    {

    }
};
} // namespace pfc