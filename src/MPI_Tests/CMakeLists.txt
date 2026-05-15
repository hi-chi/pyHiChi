include_directories(
    include
    ../../3rdparty/google/googletest/googletest/include
    ../Core/include
    ../Modules/Particle/include
    ../Modules/Field/include
    ../MPI/include
    ${FFT_INCLUDES})


file(GLOB_RECURSE MPI_SOURCES "../MPI/src/*.cpp")

add_executable(mpi_tests
    src/testMPI_Topology.cpp
    src/testMPI_FieldExchanger.cpp
    src/testMPI_GridSlicer.cpp
    src/testMPI_FieldSolver.cpp
    src/testMPI_BoundaryConditions.cpp

    src/TestingUtility.cpp

    ${MPI_SOURCES}

    src/MPI_Main.cpp
)

target_link_libraries(mpi_tests
    gtest
    MPI::MPI_CXX    
    ${FFT_LIBS}
)

set_target_properties(gtest PROPERTIES FOLDER UnitTests)
set_target_properties(mpi_tests PROPERTIES FOLDER UnitTests)