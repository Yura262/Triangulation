cmake_minimum_required(VERSION 3.14)
project(triangulation LANGUAGES CXX)


add_executable(main main.cpp)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find pybind11
set(PYTHON_EXECUTABLE "C:/Users/yura/AppData/Local/Programs/Python/Python311/python.exe")
set(PYTHON_INCLUDE_DIR "C:/Users/yura/AppData/Local/Programs/Python/Python311/include")
set(PYTHON_LIBRARY "C:/Users/yura/AppData/Local/Programs/Python/Python311/libs/python311.lib") 
set(pybind11_DIR "C:/Users/yura/AppData/Local/Programs/Python/Python311/Lib/site-packages/pybind11/share/cmake/pybind11") 

find_package(pybind11 CONFIG REQUIRED)
find_package(Python 3 COMPONENTS Interpreter REQUIRED)

# Build Python module
pybind11_add_module(triangulation bindings.cpp )