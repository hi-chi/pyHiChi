# valarray in GCC 9.5.0 is incompatible with strict C++ compilers (icpx, clang)
# applying patch from 3rdparty/patches/libstdcxx
set(CHECK_SRC "${CMAKE_BINARY_DIR}/CMakeFiles/check_valarray.cpp")
if(NOT EXISTS "${CHECK_SRC}")
	file(WRITE "${CHECK_SRC}" "#include <valarray>\nint main(){return 0;}\n")
endif()

try_compile(CHECK_RES_OK
	${CMAKE_BINARY_DIR}/CMakeFiles/check_valarray
	"${CHECK_SRC}"
	CXX_STANDARD 17
	OUTPUT_VARIABLE CHECK_LOG
)

if(NOT CHECK_RES_OK)
	message(STATUS "libstdc++ <valarray> noexcept bug detected, applying patch")
	include_directories(BEFORE SYSTEM "${CMAKE_SOURCE_DIR}/3rdparty/patches/libstdcxx")
endif()
