# Scaffolding with ScaRa
This repository contains OLD python scripts and binaries comprising our ScaRa tools for genome scaffolding.

ScaRa is based on Ezra scaffolding tool (http://complex.zesoi.fer.hr/data/pdf/Ivan_Krpelnik_diplomski.pdf) and uses come concepts from the HERA scaffolder paper (https://www.biorxiv.org/content/early/2018/06/13/345983). It consists of two phases, extension phase and bridging phase, each can be run multiple time. Extension phase uses MSA (Multiple Sequence Alignment) and POA (Partial OrderAlignment) graphs, while bridging phase constructs paths betwewen contigs using reads that overlap and then tries to determine the best path.

ScaRa uses Minimap2 to generate everlaps as needed, and can use Racon to correct the contigs at the start, et the end od after each iteration of the scaffolding process.

## Installation

  1. Clone the repository, and include all submodules.
  
    git clone --recursive https://github.com/kkrizanovic/ScaRa.git
  
  __Note:__ if you omitted `<--recursive>` from `<git clone>`, run `<git submodule update --init --recursive>` before proceeding.
  
  2. Build Ezra
  
    cd Ezra
    mkdir build
    cd build
    cmake ..
    make

  3. Build Minimap2
  
    cd Minimap2
    make

  4. Build Racon
  
    cd racon
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make

Pythons scripts, such as Bridger, ScaRa and Samscripts tool do not need to be installed.

### Dependencies
Python scripts require PYthon2.7. Ezra and Racon require CMake 3.5.

## Running ScaRa
ScaRa is run using `<scara.py>` script from the root folder of the repository.

A basic command for running ScaRa:

    python scara.py contigs.fasta reads.fastq --results results_folder

Available options will be printed out when running the script without arguments.

The scaffolding will run  in several iterations for more complete and correct scaffolding. The script will generate a separate folder for each iteration of scffolding. The final result will be in the root results folder named `<scara_scaffolds_final.fasta>`.

__Racon__
ScaRa script can also be instructed run Racon on the initial contigs, on the final scaffolds and after each iteration of scaffolding. This can be specified using options:
  - `<--racon-start>` - run racon on the initial contigs, before scaffolding
  - `<--racon-end>`   - run racon on the final scaffold
  - `<--racon>`       - run racon after each iteration

__Minimap2__
ScaRa will run Minimap2 to calcualte overlaps as needed for the scaffolding and also for Racon.
