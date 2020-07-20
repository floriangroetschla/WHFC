cmake_minimum_required(VERSION 3.10)
project(WHFC)
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -pthread -Wall -Wpedantic -Wunused -lstdc++fs")
set(CMAKE_CXX_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -g -D_GLIBCXX_DEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -mtune=native -march=native -g")
set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -lm")
#set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -fno-omit-frame-pointer")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(TBB REQUIRED)
include_directories(SYSTEM ${TBB_INCLUDE_DIRS})

#add_executable(TESTS run_tests.cpp)

add_executable(PaToH test_patoh.cpp)
target_link_libraries(PaToH ${CMAKE_SOURCE_DIR}/extern/libpatoh.a)
target_link_libraries(PaToH ${TBB_LIBRARIES})

#add_executable(TestCSR test_csr.cpp)
#target_link_libraries(TestCSR ${CMAKE_SOURCE_DIR}/extern/libpatoh.a)

add_executable(RecursiveBisection test_rb.cpp)
target_link_libraries(RecursiveBisection ${CMAKE_SOURCE_DIR}/extern/libpatoh.a)
target_link_libraries(RecursiveBisection ${TBB_LIBRARIES})

add_executable(KWayRefinement test_kway.cpp)
target_link_libraries(KWayRefinement ${CMAKE_SOURCE_DIR}/extern/libpatoh.a)
target_link_libraries(KWayRefinement ${TBB_LIBRARIES})

add_executable(KWayRefinementParallel test_kway_parallel.cpp)
target_link_libraries(KWayRefinementParallel ${CMAKE_SOURCE_DIR}/extern/libpatoh.a)
target_link_libraries(KWayRefinementParallel ${TBB_LIBRARIES})

add_executable(SnapshotRunner snapshot_runner.cpp)
target_link_libraries(SnapshotRunner ${TBB_LIBRARIES})

