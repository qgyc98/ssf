# 设置所需的 CMake 最低版本
cmake_minimum_required(VERSION 3.10)

# 设置项目名称
project(SIFEL)

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# 包含头文件目录
include_directories(${CMAKE_SOURCE_DIR}/SIFEL/GEFEL/CMLFILE)

# 包含子目录
add_subdirectory(GEFEL)
add_subdirectory(HPARAL)
# 移除 MAN 目录
add_subdirectory(MAN)
add_subdirectory(MEFEL)
add_subdirectory(METR)
add_subdirectory(PARGEF)
add_subdirectory(PARMEF)
add_subdirectory(PARMETR)
add_subdirectory(PARTRF)
add_subdirectory(PREP)
add_subdirectory(TRFEL)