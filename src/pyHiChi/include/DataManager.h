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
    adios2::Mode mode = adios2::Mode::Undefined;
    adios2::IO ioFactory = adios.DeclareIO("hiChiIO");; //
    adios2::ADIOS adios;

    adios2::Engine engine;
    int iteration = 0;
public:
    DataManager(){}
    DataManager(std::string path, int it = 0): path(path), iteration(it) {}
    void setEngine(adios2::Mode new_mode)
    {
        if (engine && mode != new_mode)
            engine.Close();
        engine = ioFactory.Open(path, mode);
        ioFactory.SetEngine("BP4");
    }
    void beginIteration() //create new file
    {
        //PutVariable(path, path);
    }
    void endIteration() //close file
    {

    }
    void beginStep() //logical step 
    {

    }
    void endStep() //flush buffer
    {

    }

    template<typename T>
    void PutVariable(const std::string name, const T& val)
    {
        adios2::Variable<T> var = ioFactory.DefineVariable<T>(name);
        engine.Put(var, val);
    }
    template<typename T>
    void PutArray(const std::string name, const T* arr, size_t size)
    {

    }
    //void PutGrid(int *arr, size_t size, std::string name)
    //{
    //
    //}

    template<typename T>
    void GetVariable(const std::string name, T& var)
    {

    }
    template<typename T>
    void GetArray(const std::string name, T* arr, size_t size)
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