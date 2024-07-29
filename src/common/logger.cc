#include "logger.h"
#include <mpi.h>
#include <iostream>

namespace sfem
{
    //=============================================================================
    Logger::Logger(const std::string &application_name,
                   const std::string &path,
                   Level level)
        : level_(level)
    {
        MPI_Comm_rank(SFEM_COMM_WORLD, &proc_rank_);
        MPI_Comm_size(SFEM_COMM_WORLD, &n_procs_);

        if (!path.empty())
        {
            log_file_ = std::ofstream(path + "_proc_" + std::to_string(proc_rank_));
        }

        std::string message = "Application: " + application_name + "\n";
        message += "Number of processes launched: " + std::to_string(n_procs_) + "\n";
        log_message(message);
    }
    //=============================================================================
    Logger &Logger::instance(const std::string &application_name,
                             const std::string &path,
                             Level level)
    {
        static Logger logger(application_name, path, level);
        return logger;
    }
    //=============================================================================
    int Logger::proc_rank() const
    {
        return proc_rank_;
    }
    //=============================================================================
    int Logger::n_procs() const
    {
        return n_procs_;
    }
    //=============================================================================
    void Logger::log_message(const std::string &message,
                             Level level,
                             const std::string &file,
                             int line)
    {
        if (level < level_)
        {
            return;
        }

        std::string _message;
        _message += "#=============================================================================#\n";
        _message += "SFEM-PROCESS-" + std::to_string(proc_rank_) + "-" + level_to_string(level) + "\n";
        if (level > Level::info)
        {
            _message += "File: " + file + "\nLine: " + std::to_string(line) + "\n";
        }
        _message += message;
        _message += "#=============================================================================#\n";

        std::cout << _message;

        if (log_file_.is_open())
        {
            log_file_ << _message;
        }

        if (level == Level::error)
        {
            abort();
        }
    }
    //=============================================================================
    void Logger::info(const std::string &message)
    {
        log_message(message);
    }
    //=============================================================================
    void Logger::warn(const std::string &message, const std::string &file, int line)
    {
        log_message(message, Level::warn, file, line);
    }
    //=============================================================================
    void Logger::error(const std::string &message, const std::string &file, int line)
    {
        log_message(message, Level::error, file, line);
    }
    //=============================================================================
    void Logger::abort() const
    {
        MPI_Abort(SFEM_COMM_WORLD, (int)Level::error);
    }
    //=============================================================================
    std::string Logger::level_to_string(Level level) const
    {
        switch (level)
        {
        case Level::all:
            return std::string("INFO");
        case Level::info:
            return std::string("INFO");
        case Level::warn:
            return std::string("WARNING");
        default:
            return std::string("ERROR");
        }
    }
}