#include "init.h"

#ifdef SFEM_HAS_PETSC
#include <petsc.h>
#endif

#ifdef SFEM_HAS_SLEPC
#include <slepc.h>
#endif

namespace sfem
{
    //=============================================================================
    void initialize(int *argc,
                    char ***argv,
                    const std::string &application_name,
                    const std::string &log_path,
                    Logger::Level level)
    {
        if (argc != nullptr)
        {
#ifdef SFEM_HAS_PETSC
            PetscInitialize(argc, argv, nullptr, nullptr);
#endif
#ifdef SFEM_HAS_SLEPC
            SlepcInitialize(argc, argv, nullptr, nullptr);
#endif
        }
        Logger::instance(application_name, log_path, level);
    }
    //=============================================================================
    void finalize()
    {
#ifdef SFEM_HAS_PETSC
        PetscBool petsc_init;
        PetscInitialized(&petsc_init);
        if (petsc_init)
        {
            PetscFinalize();
        }
#endif
    }
}