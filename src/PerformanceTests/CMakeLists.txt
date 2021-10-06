include_directories(
    include
    ../../3rdparty/google/benchmark/include
    ../Core/include
    ../Modules/Particle/include
    ../Modules/Field/include
    ${FFT_INCLUDES})

add_executable(ptests
    src/ptestFDTD.cpp
    src/ptestFDTDSYCL.cpp
    src/ptestPusher.cpp
    src/ptestPusherSYCL.cpp
    src/Main.cpp)

if (APPLE)
    target_link_libraries(ptests
        -Wl,-all_load
        benchmark::benchmark
        ${FFT_LIBS}
        -Wl,-noall_load)
elseif (UNIX)
    target_link_libraries(ptests
        -Wl,-whole-archive
        benchmark::benchmark
        -Wl,--start-group ${FFT_LIBS} -Wl,--end-group
        -Wl,-no-whole-archive)
else()
    target_link_libraries(ptests
        benchmark::benchmark
        ${FFT_LIBS})
endif()

set_target_properties(benchmark PROPERTIES FOLDER PerformanceTests)
set_target_properties(ptests PROPERTIES FOLDER PerformanceTests)