# 设置目标输出目录
set(OUTPUTPATH ${CMAKE_BINARY_DIR}/MAN)

# 包含头文件目录
include_directories(${CMAKE_SOURCE_DIR}/MAN)

# 添加源文件
file(GLOB_RECURSE SOURCES "${CMAKE_SOURCE_DIR}/MAN/*.cpp")

# 如果没有找到源文件，则不创建库
if(NOT SOURCES)
    message(STATUS "No source files found for target: man, skipping library creation")
else()
    # 添加库
    add_library(man STATIC ${SOURCES})
endif()