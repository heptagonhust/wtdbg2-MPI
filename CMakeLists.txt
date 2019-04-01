cmake_minimum_required(VERSION 2.8)
SET(CMAKE_C_COMPILER icc )
SET(CMAKE_CXX_COMPILER icc)

set(CMAKE_C_FLAGS "-Wall -Wno-unused-but-set-variable -mpopcnt -march=native -Wno-pointer-arith -Wno-unused-function -Wno-write-strings -Wno-format")
set(CMAKE_CXX_FLAGS ${CMAKE_C_FLAGS})

add_library(ksw ksw.cpp kbm_opt.cpp kbm_core.cpp debrain.h)
target_link_libraries(ksw pthread m rt z glog)

add_executable(kbm kbm.cpp)
add_executable(wtdbg wtdbg.cpp)
add_executable(wtdbg-cns wtdbg-cns.cpp)
add_executable(wtpoa-cns wtpoa-cns.cpp)

target_link_libraries(kbm ksw)
target_link_libraries(wtdbg ksw)
target_link_libraries(wtdbg-cns ksw)
target_link_libraries(wtpoa-cns ksw)




add_executable(pgzf pgzf.cpp)
target_link_libraries(pgzf pthread m rt z)
