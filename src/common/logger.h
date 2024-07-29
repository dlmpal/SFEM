#pragma once

#include "config.h"
#include <string>
#include <fstream>

namespace sfem
{
    class Logger
    {
    public:
        enum class Level : int
        {
            all = 0,
            info = 1,
            warn = 2,
            error = 3
        };

        Logger(const Logger &) = delete;
        Logger &operator=(const Logger &) = delete;

        /// @brief Get the (static) Logger instance
        /// @note All input parameters only matter for initialization
        /// @param appliaction_name The name of the SFEM application
        /// @param path The path for the log file
        static Logger &instance(const std::string &appliaction_name = "",
                                const std::string &path = "",
                                Level level = Level::info);

        ///@brief Get the MPI process rank
        int proc_rank() const;

        /// @brief Get the number of MPI processes
        int n_procs() const;

        /// @brief Log a message
        /// @param message The messaged to be logged
        /// @param level The log level
        /// @param file The file from which the message originates (e.g __FILE__)
        /// @param line The line in the file (e.g __LINE__)
        void log_message(const std::string &message,
                         Level level = Level::info,
                         const std::string &file = "",
                         int line = -1);

        /// @brief Log an info message
        /// @note See log_message
        void info(const std::string &message);

        /// @brief Log a warning
        /// @note See log_message
        void warn(const std::string &message,
                  const std::string &file = "",
                  int line = -1);

        /// @brief Log an error
        /// @note Also aborts the program
        /// @note See log_message
        void error(const std::string &message,
                   const std::string &file = "",
                   int line = -1);

        /// @brief Abort the application
        void abort() const;

        /// @brief Get the corresponding string for a Level
        std::string level_to_string(Level level) const;

    private:
        /// @brief Constructor
        /// @note See GetInstance
        Logger(const std::string &application_name = "SFEM_PROGRAM",
               const std::string &path = "",
               Level level = Level::info);

        Level level_;
        int proc_rank_;
        int n_procs_;
        std::ofstream log_file_;
    };
}