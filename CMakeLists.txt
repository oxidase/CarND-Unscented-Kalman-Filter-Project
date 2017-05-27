project(UnscentedKF)

cmake_minimum_required (VERSION 3.5)

add_definitions(-std=c++11)

set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")

set(sources src/ukf.cpp src/main.cpp src/tools.cpp)

set(LIBUV_LIBRARY "uv")
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

include_directories(/usr/local/include)
include_directories(/usr/local/opt/openssl/include)
link_directories(/usr/local/lib)
link_directories(/usr/local/opt/openssl/lib)
link_directories(/usr/local/Cellar/libuv/1.11.0/lib)

elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

  find_library(LIBUV_LIBRARY NAMES uv)
  if(NOT LIBUV_LIBRARY)
    message(FATAL_ERROR "uv library is not found")
  endif()
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")


add_executable(UnscentedKF ${sources})

target_link_libraries(UnscentedKF z ssl ${LIBUV_LIBRARY} uWS)
