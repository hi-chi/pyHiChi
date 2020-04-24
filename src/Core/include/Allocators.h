#pragma once
#include <vector>

namespace pfc {

    // NUMA allocator for ScalarField
    template <class Data>
    class NUMA_Allocator {
    public:

        using value_type = Data;
        using propagate_on_container_move_assignment = true_type;
        using is_always_equal = true_type;

        constexpr NUMA_Allocator() noexcept {}
        constexpr NUMA_Allocator(const NUMA_Allocator&) noexcept = default;
        template<class OtherData>
        constexpr NUMA_Allocator(const NUMA_Allocator<OtherData>&) noexcept {}

        value_type * allocate(const size_t num)
        {
            const size_t size_vt = sizeof(value_type);
            const size_t len = num * size_vt;
            char * p = reinterpret_cast<char*>(std::malloc(len));
#pragma omp parallel for simd
            for (int i = 0; i < len; i += size_vt)
#pragma unroll(16)
                for (int j = 0; j < size_vt; j++)
                    p[i + j] = 0;
            return reinterpret_cast<value_type*>(p);
        }

        void deallocate(value_type * const p, const size_t num)
        {
            std::free(p);
        }

        friend int operator==(const NUMA_Allocator& a1, const NUMA_Allocator& a2) {
            return true;
        }

        friend int operator!=(const NUMA_Allocator& a1, const NUMA_Allocator& a2) {
            return false;
        }

    };

}