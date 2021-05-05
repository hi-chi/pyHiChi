#pragma once
#include <cstring>
#include <vector>
#include <omp.h>

namespace pfc {

    // NUMA allocator for ScalarField
    // Based upon ideas by Georg Hager and Gerhard Wellein
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
            const size_t num_threads = omp_get_max_threads();
            if (num_threads > num) {
                char * p = reinterpret_cast<char*>(std::malloc(len));
                std::memset(p, 0, len);
                return reinterpret_cast<value_type*>(p);
            }
            const size_t block_size = (num / num_threads) * size_vt;
            const size_t block_size_rem = len - block_size * (num_threads - 1);
            char * p = reinterpret_cast<char*>(std::malloc(len));
#pragma omp parallel for
            for (int thr = 0; thr < num_threads; thr++) {
                const size_t cur_block_size = thr == num_threads - 1 ? block_size_rem : block_size;
                std::memset(p + thr * block_size, 0, cur_block_size);
            }
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
