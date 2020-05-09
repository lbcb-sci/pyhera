# Scaffolding with ScaRa
ScaRa is our scaffolding tool that uses come concepts from the HERA scaffolder paper (https://www.biorxiv.org/content/early/2018/06/13/345983). This repository contains both Python and C++ implementations.

## Installation

  1. Clone the repository, and include all submodules.
  
    git clone --recursive https://github.com/kkrizanovic/ScaRa.git
  
  __Note:__ if you omitted `<--recursive>` from `<git clone>`, run `<git submodule update --init --recursive>` before proceeding.
  
  After cloning, Python scripts are ready to go.
  
  2. Building C++ version
  
    mkdir build
    cd build
    cmake ..
    make

### Dependencies
Python scripts require PYthon2.7. C++ version requires CMake 3.5.

## Running Python version
ScaRa is run using `<bridger.py>` script from the root folder of the repository. Running the script without arguments will print out a help message.

## Running C++ version
C++ executable should be in the build folder. Running it without arguments will print out a help message.

## OLD ScaRa 
Old ScaRa README (containing only the Python version) can be found at [OLD ScaRa](https://github.com/lbcb-sci/ScaRa/blob/master/SCARA_OLD.md)
