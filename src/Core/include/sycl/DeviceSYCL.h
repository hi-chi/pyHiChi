#pragma once

#include <sycl/CL/sycl.hpp>

namespace pfc {

    namespace sycl_pfc {

        enum Devices { CPU = 0, GPU = 1 };

        template <typename T>
        using sycl_vector = std::vector<T, sycl::usm_allocator<T, sycl::usm::alloc::shared>>;

        class AllDevices
        {
        public:
            AllDevices() : default_device{ sycl::default_selector{} },
                cpu_device{ sycl::cpu_selector{}, sycl::async_handler{} }                
            {
                try {
                    gpu_device = sycl::queue { sycl::gpu_selector{}, sycl::async_handler{} };
                } catch (sycl::exception const& e) {
                    std::cout << "Cannot select a GPU\n";
                    std::cout << "Using a CPU device\n";
                    gpu_device = sycl::queue { sycl::cpu_selector{}, sycl::async_handler{} };
                }
                
                std::cout << "Host Device: " << default_device.get_device().get_info<sycl::info::device::name>() << std::endl;
                std::cout << "CPU Device: " << cpu_device.get_device().get_info<sycl::info::device::name>() << std::endl;
                std::cout << "GPU Device: " << gpu_device.get_device().get_info<sycl::info::device::name>() << std::endl;
                std::cout << std::fflush;
            };
            sycl::queue default_device;
            sycl::queue cpu_device;
            sycl::queue gpu_device;
        };

        extern AllDevices node;
    }
}