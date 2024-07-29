#pragma once

#include "logger.h"

namespace sfem
{
    /// @brief Should be called at the start of each sfem application.
    void initialize(int *argc,
                    char ***argv,
                    const std::string &application_name = "SFEM_PROGRAM",
                    const std::string &log_path = "",
                    Logger::Level level = Logger::Level::info);

    /// @brief Should be called at the end of each sfem application
    void finalize();
}