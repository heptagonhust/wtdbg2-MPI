cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_FLAGS "-Wall -Wno-unused-but-set-variable -mpopcnt -march=native")

include_directories(include)
add_library(ksw ksw.cpp)
target_link_libraries(ksw pthread m rt z)


# add_executable(kbm kbm.cpp)
# add_executable(wtdbg wtdbg.cpp)
#add_executable(wtdbg-cns wtdbg-cns.cpp)
#add_executable(wtpoa-cns wtpoa-cns.cpp)

#target_link_libraries(kbm ksw)
# target_link_libraries(wtdbg ksw)
#target_link_libraries(wtdbg-cns ksw)
#target_link_libraries(wtpoa-cns ksw)




#add_executable(pgzf pgzf.cpp)
#target_link_libraries(pgzf pthread m rt z)
