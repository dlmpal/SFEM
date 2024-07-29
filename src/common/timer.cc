#include "timer.h"

namespace sfem::common
{
    //=============================================================================
    Timer::Timer(const std::string &_process_name)
        : _function_name(_process_name)
    {
        _start = std::chrono::high_resolution_clock::now();
    }
    //=============================================================================
    Timer::~Timer()
    {
        _stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(_stop - _start);

        std::string message = _function_name + " completed in: " +
                              std::to_string(duration.count()) + " milliseconds\n";
        Logger::instance().log_message(message, Logger::Level::all);
    }
}