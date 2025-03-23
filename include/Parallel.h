#pragma once

#include <omp.h>
#include <thread>

#if defined(__cpp_lib_hardware_interference_size)
// default cacheline size from runtime
constexpr std::size_t CL_SIZE = std::hardware_constructive_interference_size;
constexpr std::size_t CL_NUM_PADDING = (std::hardware_constructive_interference_size / 8) - 1;
#else
// most common cacheline size otherwise
constexpr std::size_t CL_SIZE = 64;
constexpr std::size_t CL_NUM_PADDING = 7;
#endif
