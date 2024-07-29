"""
Installation script for SFEM
python3 install.py --remove_previous_build --mpi_dir=/usr/lib/x86_64-linux-gnu/openmpi/ 
--mpi_cxx_compiler=/usr/bin/mpic++ --petsc_dir=$PETSC_DIR 
--petsc_arch=$PETSC_ARCH --metis_dir=/home/user/local 
--slepc_dir=$SLEPC_DIR --install_dir=$SFEM_DIR --with_apps --with_pysfem --with_stubs
"""

import os
import argparse
from subprocess import run

# ==============================================================================
parser = argparse.ArgumentParser(
    description='Build and install the SFEM finite-element package')
parser.add_argument('--install_dir',
                    type=str,
                    default="install",
                    help='The path of the installation directory for SFEM')
parser.add_argument('--build_dir',
                    type=str,
                    default="build",
                    help='The path of the build directory for SFEM')
parser.add_argument('--mpi_dir',
                    type=str,
                    required=True,
                    help='The path of the root directory for MPI')
parser.add_argument('--mpi_cxx_compiler',
                    type=str,
                    required=True,
                    help='The path of the MPI C++ compiler')
parser.add_argument('--petsc_dir',
                    type=str,
                    help='The path of the root directory for PETSc')
parser.add_argument('--petsc_arch',
                    type=str,
                    help='The PETSC-ARCH directory')
parser.add_argument('--slepc_dir',
                    type=str,
                    required=False,
                    help='The path of the root directory for SLEPc')
parser.add_argument('--metis_dir',
                    type=str,
                    required=False,
                    help='The path of the root directory for METIS')
parser.add_argument("--build_type",
                    type=str,
                    default="RELEASE",
                    required=False,
                    help="The build type for cmake. Either RELEASE or DEBUG.")
parser.add_argument('--with_apps',
                    action="store_true",
                    help="Build the applications")
parser.add_argument('--with_pysfem',
                    action="store_true",
                    help="Build the Python bindings for SFEM")
parser.add_argument('--with_stubs',
                    action="store_true",
                    help="Whether to generate stubs (.pyi) for the Python bindings")
parser.add_argument('--remove_previous_build',
                    action="store_true",
                    help="Remove the previous build, assuming the path has not changed")
args = parser.parse_args()
# ==============================================================================
########################## DO NOT EDIT BELOW ###################################
# ==============================================================================
# Options passed to cmake
CMAKE_INSTALL_OPTIONS = {
    "CMAKE_INSTALL_PREFIX": args.install_dir,
    "CMAKE_BUILD_TYPE": args.build_type,
    "CMAKE_CXX_COMPILER": args.mpi_cxx_compiler,
    "MPI_DIR": args.mpi_dir,
    "PETSC_DIR":  args.petsc_dir,
    "PETSC_ARCH": args.petsc_arch,
    "SLEPC_DIR": args.slepc_dir,
    "METIS_DIR": args.metis_dir,
    "WITH_APPS": "On" if args.with_apps is True else "Off",
    "WITH_PYSFEM": "On" if args.with_pysfem is True else "Off"
}
# ==============================================================================
# Paths of the build/config/install directories
BUILD_DIR = args.build_dir
CMAKE_DIR = "cmake"
INSTALL_DIR = CMAKE_INSTALL_OPTIONS["CMAKE_INSTALL_PREFIX"]
# ==============================================================================
# Number of job for make to run
N_JOBS_MAKE = os.cpu_count()
# ==============================================================================
# Remove folders from previous build/config/install
if args.remove_previous_build:
    remove_prev_build = f"rm -rf {BUILD_DIR}"
    run(remove_prev_build, shell=True)

    remove_prev_config = f"rm -rf cmake"
    run(remove_prev_config, shell=True)

remove_prev_install = f"rm -rf {INSTALL_DIR}"
run(remove_prev_install, shell=True)
# ==============================================================================
# Create the Config.cmake.in
cmake_config_filepah = os.path.join(CMAKE_DIR, "Config.cmake.in")
os.makedirs(os.path.dirname(cmake_config_filepah), exist_ok=True)
with open(cmake_config_filepah, "w") as file:
    file.writelines(
        "@PACKAGE_INIT@\ninclude(${CMAKE_CURRENT_LIST_DIR}/sfemTargets.cmake)\ncheck_required_components(sfem)")
# ==============================================================================
# Add the install options, and run cmake
run_cmake = f"cmake -S . -B {BUILD_DIR} "
for k, v in CMAKE_INSTALL_OPTIONS.items():
    if v is not None:
        run_cmake += f"-D{k}={v} "
run(run_cmake, shell=True)
# ==============================================================================
# Build
run_make = f"cd {BUILD_DIR}; make -j {N_JOBS_MAKE}"
run(run_make, shell=True)
# ==============================================================================
# Installation
run_make_install = f"cd {BUILD_DIR}; make install"
run(run_make_install, shell=True)
# ==============================================================================
# Create the stub .pyi files for the Python bindings
if args.with_pysfem:
    pysfem_dir = os.path.join(args.install_dir, "lib", "pysfem")
    with open(os.path.join(pysfem_dir, "__init__.py"), "w") as file:
        file.write("from .pysfem import *")

    if args.with_stubs:
        create_stubs = f"cd {pysfem_dir}; python3 -m nanobind.stubgen -m pysfem --recursive"
        run(create_stubs, shell=True)
        move_root_pyi = f"cd {pysfem_dir}; mv pysfem.pyi __init__.pyi"
        run(move_root_pyi, shell=True)
