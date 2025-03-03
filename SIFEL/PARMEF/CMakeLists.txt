# 设置目标输出目录
set(OUTPUTPATH ${CMAKE_BINARY_DIR}/PARMEF)

# 包含头文件目录
include_directories(${CMAKE_SOURCE_DIR}/PARMEF)

# 添加源文件
file(GLOB_RECURSE SOURCES "${CMAKE_SOURCE_DIR}/PARMEF/*.cpp")

# 添加库
add_library(parmef STATIC ${SOURCES})

# 链接库（如果有其他库需要链接，在这里添加）
# target_link_libraries(parmef <library_name>)