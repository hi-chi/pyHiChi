
include_directories(
    include
    ../Core/include
    ../Modules/Particle/include
    ../Modules/Field/include
    ../pyHiChi/include
    ${FFT_INCLUDES})

add_executable(application
    src/main.cpp)

if (APPLE)
    target_link_libraries(application
        -Wl,-all_load
        ${FFT_LIBS}
        -Wl,-noall_load)
elseif (UNIX)
    target_link_libraries(application
        -Wl,-whole-archive
        -Wl,--start-group ${FFT_LIBS} -Wl,--end-group
        -Wl,-no-whole-archive)
else()
    target_link_libraries(application
        ${FFT_LIBS})
endif()

set_target_properties(application PROPERTIES FOLDER Application)