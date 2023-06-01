include_directories(
    include
    ../../3rdparty/google/googletest/googletest/include
    ../Core/include
    ../Modules/Particle/include
    ../Modules/Field/include
    ${FFT_INCLUDES})

add_executable(tests
    src/testConstants.cpp
    src/testDimension.cpp
    src/testEnsemble.cpp
    src/testFieldSolvers.cpp
    src/testFourierTransform.cpp
    src/testFP.cpp
    src/testGrid.cpp
    src/testMerging.cpp
    src/testParticle.cpp
    src/testParticleArray.cpp
    src/testParticleProxy.cpp
    src/testPML.cpp
    src/testPusherAndHandler.cpp
    src/testSaveLoadParticle.cpp
    src/testSaveLoadGrid.cpp
    src/testScalarField.cpp
    src/testSpecies.cpp
    src/testThinning.cpp
    src/testVectors.cpp
    src/testVectorsProxy.cpp
    src/TestingUtility.cpp
    src/Main.cpp)

if (APPLE)
    target_link_libraries(tests
        -Wl,-all_load
        gtest
        ${FFT_LIBS}
        -Wl,-noall_load)
elseif (UNIX)
    target_link_libraries(tests
        -Wl,-whole-archive
        gtest
        -Wl,--start-group ${FFT_LIBS} -Wl,--end-group
        -Wl,-no-whole-archive)
else()
    target_link_libraries(tests
        gtest
        ${FFT_LIBS})
endif()

set_target_properties(gtest PROPERTIES FOLDER UnitTests)
set_target_properties(tests PROPERTIES FOLDER UnitTests)