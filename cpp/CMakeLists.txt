cmake_minimum_required(VERSION 3.13)
project(herring_model)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_BUILD_TYPE Release)

find_package(Boost COMPONENTS program_options)
IF(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

find_package(HDF5 COMPONENTS CXX)
IF(HDF5_FOUND)
    include_directories(${HDF5_INCLUDE_DIRS})
endif()

include_directories(./include)

add_executable(herring_model programs/herring_model.cpp src/def.cpp src/utils.cpp src/state_type.cpp src/observer_type.cpp src/solver_type.cpp include/herring_func.hpp)
target_link_libraries(herring_model ${Boost_LIBRARIES} ${HDF5_LIBRARIES})
