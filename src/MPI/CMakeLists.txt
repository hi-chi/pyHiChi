include_directories(
    include
    ../Core/include
    ../Modules/Particle/include
    ../Modules/Field/include
    ${FFT_INCLUDES})

add_executable(mpi
    src/MPI_Utilities.cpp
    src/MPI_Topology.cpp
    src/MPI_FieldUtils.cpp
    src/MPI_FieldExchanger.cpp
    src/MPI_GridSlicer.cpp
    Main.cpp)

target_link_libraries(mpi
    MPI::MPI_CXX)

set_target_properties(mpi PROPERTIES FOLDER MPI)