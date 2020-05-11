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
ScaRa is run using `<bridger.py>` script from the root folder of the repository. Running the script without arguments will print out a help message. The script should be run wtith the first parameter set to `scaffold` in the following way:

    bridger.py scaffold <contigs FASTA> <reads FASTA> <reads-contigs overlaps PAF/SAM> <reads-reads overlaps PAF/SAM options
    options:"
    -o (--output) <file> : output file to which the report will be written

## Running C++ version
C++ executable should be in the build folder. Running it without arguments will print out following help message.

    How to run ScaRa :
    ScaRa -f (--folder) <Input folder>
        or
    ScaRa -r <Reads file> -c <Contigs file> -o <Reads to Contigs Overlaps file> -s <Reads to Reds Overlaps file>
        The program will perform one iteration of the algorithm
        and output contigs to the standard output!.
        <Input foder> must contain the following files:
        - reads.fastq - reads in FASTQ/FASTA format
        - readsToContigs.paf - overlaps between reads and contigs
        - readsToReads.paf - overlaps between reads
        - contigs.fasta - contigs in FASTA format
    If Reads file, Contigs file or Overlaps are specified,
    they will not be looked for in the Input folder!
    Options:
    -f (--folder)     specify input folder for ScaRa
    -r (--reads)      specify reads file for ScaRa
    -c (--contigs)    specify contigs file ScaRa
    -o (--overlapsRC)   specify contig-read overlaps file for ScaRa
    -s (--overlapsRR)   specify read self overlaps file for ScaRa
    -m (--multithreading)   use multithreading
    -v (--version)    print program version
    -h (--help)       print this help message


## OLD ScaRa 
Old ScaRa README (containing only the Python version) can be found at [OLD ScaRa](https://github.com/lbcb-sci/ScaRa/blob/master/SCARA_OLD.md)
