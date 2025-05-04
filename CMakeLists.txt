cmake_minimum_required(VERSION 3.10)
project(XIXILBJ_hw1 C)  # 项目名称需与目标一致

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")

include_directories(inc)  # 头文件目录

# 指定源文件路径
add_executable(XIXILBJ_hw1   # 目标名称与 project() 一致
    src/algebra.c
    src/main.c               # 明确指定 main.c 的路径
)

target_link_libraries(XIXILBJ_hw1 m)