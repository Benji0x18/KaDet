cmake_minimum_required(VERSION 3.10)

project(Kadet)

set(CMAKE_CXX_STANDARD 17)
set(source_dir "${PROJECT_SOURCE_DIR}/src/")
set(include_dir "${PROJECT_SOURCE_DIR}/include/")

find_package(Threads REQUIRED)
include_directories("/usr/include" "/usr/local/include")
link_directories("/usr/lib/cantera" "usr/local/lib")

add_library(kadet SHARED ${source_dir}/Kadet.cpp)
target_include_directories(kadet PRIVATE ${include_dir})
target_link_libraries(kadet cantera)

add_executable(main "${source_dir}/main.cpp")
target_include_directories(main PRIVATE ${include_dir})
target_link_libraries(main kadet)
target_link_libraries(main cantera)



