#pragma once

#include "logger.h"
#include <chrono>

namespace sfem::common
{
    class Timer
    {
    public:
        Timer(const std::string &_process_name);
        ~Timer();

    private:
        std::chrono::time_point<std::chrono::high_resolution_clock> _start;
        std::chrono::time_point<std::chrono::high_resolution_clock> _stop;
        std::string _function_name;
    };
}