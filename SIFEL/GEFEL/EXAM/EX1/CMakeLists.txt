# 设置目标输出目录
set(OUTPUTPATH ${CMAKE_BINARY_DIR}/EX1)

# 包含头文件目录
include_directories(${CMAKE_SOURCE_DIR}/GEFEL/EXAM/EX1)

# 添加源文件
file(GLOB_RECURSE SOURCES "${CMAKE_SOURCE_DIR}/GEFEL/EXAM/EX1/*.cpp")

# 添加可执行文件
add_executable(ex1 ${SOURCES})

# 链接库（如果有其他库需要链接，在这里添加）
# target_link_libraries(ex1 <library_name>)