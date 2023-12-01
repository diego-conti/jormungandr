cmake_minimum_required(VERSION 3.10)
project(jormungandr VERSION 0.3)
include(FetchContent)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
add_compile_options(-g -O3 -Wctor-dtor-privacy -Wreorder -Wold-style-cast -Wsign-promo -Wchar-subscripts -Winit-self -Wmissing-braces -Wparentheses -Wreturn-type -Wswitch -Wtrigraphs -Wextra -Wno-sign-compare -Wno-narrowing -Wno-attributes)


FetchContent_Declare(demonblast
	GIT_REPOSITORY https://github.com/diego-conti/DEMONbLAST
	GIT_TAG c2aaa9a7337b8e27c80d9eceb8ddc6406b54afaf
)
FetchContent_Populate(demonblast)

set(DEMONBLAST_SOURCES partitions.cpp tree.cpp labeled_tree.cpp weightbasis.cpp log.cpp niceliegroup.cpp liegroupsfromdiagram.cpp gauss.cpp niceeinsteinliegroup.cpp ricci.cpp filter.cpp permutations.cpp weightmatrix.cpp implicitmetric.cpp antidiagonal.cpp adinvariantobstruction.cpp parsetree.cpp automorphisms.cpp partitionprocessor.cpp diagramprocessor.cpp nicediagramsinpartition.cpp)
list(TRANSFORM DEMONBLAST_SOURCES PREPEND ${demonblast_SOURCE_DIR}/src/)

set(SOURCES src/jormungandr.cpp ${DEMONBLAST_SOURCES})

add_executable(jormungandr ${SOURCES})

target_link_libraries(jormungandr PUBLIC ginac wedge cocoa gmp boost_program_options)
target_link_directories(jormungandr PUBLIC ${GINAC_DIR}/ginac/.libs $ENV{WEDGE_PATH}/lib)
target_include_directories(jormungandr PUBLIC "${PROJECT_BINARY_DIR}" ${PROJECT_SOURCE_DIR}/src ${demonblast_SOURCE_DIR}/src/ $ENV{WEDGE_PATH}/include)


