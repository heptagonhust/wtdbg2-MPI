cmake_minimum_required(VERSION 3.12)

set(CMAKE_C_FLAGS "-Wall -Wno-unused-but-set-variable -mpopcnt -march=native -Wno-pointer-arith -Wno-unused-function")
set(CMAKE_CXX_FLAGS ${CMAKE_C_FLAGS})

add_library(ksw ksw.c kbm_opt.cpp kbm_core.c)
target_link_libraries(ksw pthread m rt z)

add_executable(kbm kbm.c)
add_executable(wtdbg wtdbg.c)
add_executable(wtdbg-cns wtdbg-cns.c)
add_executable(wtpoa-cns wtpoa-cns.c)

target_link_libraries(kbm ksw)
target_link_libraries(wtdbg ksw)
target_link_libraries(wtdbg-cns ksw)
target_link_libraries(wtpoa-cns ksw)




add_executable(pgzf pgzf.c)
target_link_libraries(pgzf pthread m rt z)