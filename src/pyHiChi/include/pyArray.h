#pragma once
#include <vector>
#include <FP.h>

namespace pfc
{
    class pyArray1d {
    public:

        pyArray1d(size_t nx): data(nx, 0) {}

        FP* getData() const {
            return data.data();
        }

        size_t getSize() const {
            return data.size();
        }

        // read-write accessors
        FP& get(size_t index) {
            return data[index];
        }

    private:

        std::vector<FP> data;
    };
}