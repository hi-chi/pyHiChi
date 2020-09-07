#pragma once
#include <string>
#include "Grid.h"

#include <adios2.h>
namespace pfc {

#define put_var(x) PutVariable((#x), (x));
#define put_arr(x) PutArray((#x), (x));
#define get_arr(x) GetArray((#x), (x));
#define get_arr(x) GetArray((#x), (x));
enum IOType { Undefined, Write, Read };
class DataManager
{
    std::string path; //откуда брать данные и куда их писать
    IOType mode;

    adios2::ADIOS adios;
    adios2::IO ioFactory = adios.DeclareIO("hiChiIO");
    adios2::Engine engine;
    int iteration = 0;
public:
    DataManager(){}
    DataManager(std::string path, int it = 0): path(path), iteration(it) {}
    DataManager(std::string path, IOType mode, int it = 0): path(path), mode(mode), iteration(it)
    {
        setEngine(mode);
    }
    void setEngine(IOType new_mode)
    {
        if (engine && mode != new_mode)
            engine.Close();
        if (mode == IOType::Read || mode == IOType::Write)
        {
            engine = ioFactory.Open(path, (adios2::Mode)mode);
            ioFactory.SetEngine("BP4");
        }
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
    void putVariable(const std::string name, const T& val)
    {
        adios2::Variable<T> var = ioFactory.DefineVariable<T>(name);
        engine.Put(var, val);
    }
    template<typename T>
    void putArray(const std::string name, const T* arr, size_t size)
    {
        adios2::Variable<T> var = bpIO.DefineVariable<T>(name, {},
                                  {}, { size }, adios2::ConstantDims);
        engine.Put(var, arr);
    }
    template<typename T>
    void getVariable(const std::string name, T& val)
    {
        adios2::Variable<T> var = bpIO.DefineVariable<T>(name);
        engine.Get(var, val);
    }
    template<typename T>
    void getArray(const std::string name, T* arr, size_t size)
    {
        adios2::Variable<T> var = bpIO.DefineVariable<T>(name, {},
                                  {}, { size }, adios2::ConstantDims);
        engine.Get(var, arr);
    }
    
    template<typename Data, GridTypes gridType>
    void customPut(const Grid<Data, gridType> &grid)
    {
        Put(grid.dt);
    }
};
} // namespace pfc