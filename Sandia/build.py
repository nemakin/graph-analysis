#!/usr/bin/env python3
"""
Build script for sandia_spla.cpp
Compiles the Sandia Delta-Stepping SSSP algorithm using SPLA library
"""

import os
import sys
import subprocess
import shutil
from pathlib import Path

# Color codes for terminal output
class Colors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def log_info(msg):
    print(f"{Colors.OKBLUE}[INFO]{Colors.ENDC} {msg}")

def log_success(msg):
    print(f"{Colors.OKGREEN}[SUCCESS]{Colors.ENDC} {msg}")

def log_warning(msg):
    print(f"{Colors.WARNING}[WARNING]{Colors.ENDC} {msg}")

def log_error(msg):
    print(f"{Colors.FAIL}[ERROR]{Colors.ENDC} {msg}")

def run_command(cmd, cwd=None, description=None):
    """Run a shell command and handle errors"""
    if description:
        log_info(f"{description}...")
    try:
        result = subprocess.run(cmd, shell=True, cwd=cwd, check=True)
        return result.returncode == 0
    except subprocess.CalledProcessError as e:
        log_error(f"Command failed: {cmd}")
        return False

def main():
    # Setup paths
    script_dir = Path(__file__).parent.resolve()
    project_root = script_dir.parent
    spla_dir = project_root / "spla"
    build_dir = script_dir / "build"
    output_dir = build_dir / "bin"
    cpp_file = script_dir / "sandia_spla.cpp"

    print(f"{Colors.BOLD}{Colors.WARNING}=== Sandia SPLA Compilation Script ==={Colors.ENDC}")
    log_info(f"Script directory: {script_dir}")
    log_info(f"Project root: {project_root}")
    log_info(f"SPLA directory: {spla_dir}")
    log_info(f"Build directory: {build_dir}")

    # Validation
    if not spla_dir.exists():
        log_error(f"SPLA library not found at {spla_dir}")
        return 1

    if not cpp_file.exists():
        log_error(f"Source file not found at {cpp_file}")
        return 1

    # Create build directory
    build_dir.mkdir(exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Build using CMake
    log_info("Configuring build with CMake...")
    
    cmake_cmd = f"""
    cmake -B {build_dir} \
          -S {script_dir} \
          -DCMAKE_BUILD_TYPE=Release \
          -DSPLA_BUILD_TESTS=OFF \
          -DSPLA_BUILD_EXAMPLES=OFF \
          -DSPLA_BUILD_OPENCL=ON
    """
    
    if not run_command(cmake_cmd, description="CMake configuration"):
        log_error("CMake configuration failed")
        return 1

    log_info("Building project...")
    
    # Get number of processors
    try:
        import multiprocessing
        num_procs = multiprocessing.cpu_count()
    except:
        num_procs = 1

    build_cmd = f"cmake --build {build_dir} -j {num_procs}"
    
    if not run_command(build_cmd, description="Building with CMake"):
        log_error("Build failed")
        return 1

    # Check if executable was created
    exe_path = output_dir / "sandia_spla"
    if exe_path.exists():
        log_success(f"Compilation successful!")
        log_success(f"Output: {exe_path}")
        return 0
    else:
        log_error(f"Executable not found at {exe_path}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
