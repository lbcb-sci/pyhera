# Scaffolding with pyhera
This repository contains python scripts and binaries comprising our tools for genome scaffolding.

The  tools are currently under development and any feedback on them is very welcome and will be greatly appreciated.

## Tools
### PyHera
PyHera is a python implementation of HERA scaffolder. HERA uses overlaps between contigs and reads and also self-overalps between reads to generate a graph. It then generates numerous paths between contigs, groups them according to their length and choses the best representative path.
The paper on the HERA is available at bioRchiv: https://www.biorxiv.org/content/early/2018/06/13/345983.

### Ezra
Ezra is a result of a graduation thesis work done by Ivan Krpelnik on University of Zareb Faculty of Electrical Engineering and Computing. Ezra uses MSA (Multiple Sequence Alignment) to extend contigs on each end and is also able to connect two contigs using long reads that overlap with both of them. To make the best use of the extension step, Ezra should be run in several iterations.
This thesis is available here: http://complex.zesoi.fer.hr/data/pdf/Ivan_Krpelnik_diplomski.pdf.

Ezra is included as a submodule, with main repository at: https://gitlab.com/Krpa/ezra.

### Scaffolding script
Scaffolding script combines PYHera and Ezra, according to a given scaffolding plan, to iteratively perform the scaffolding, using output of the previous iteration as input for the next one. The script will use Minimap2 to produce overlaps needed for the scaffolding.

## Installation

  1. Clone the repository, and include all submodules.
  
    git clone --recursive https://github.com/kkrizanovic/pyhera.git
  
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

  Pythons scripts, such as PyHera, Scaffolder script and samscripts tool do not need to be installed.

### Dependencies
Python scripts require PYthon2.7. Ezra requires CMake 3.5.

## Running the scripts
### Running Ezra
After building, Ezra executable will be created in `<Ezra/build>` folder. Ezra receives a single folder as an argument and outputs scaffolds to the standard output. The folder must contain the following files:
- reads.fastq - file with read sequences, must be in FASTQ format
- contigs.fasta - file with contig sequences
- readsToContigs.paf - file with overlaps between reads and contigs, when calculating overlaps, contigs must be used as reference

Typical Ezra run command:

    ezra folder > output.fasta

### Running PyHera
PyHera is run by running `<pyhera.py>` script, from the root folder of the repository. PyHera can be run with several options and various arguments, which are printed to the screen when running the script without arguments.

A typical command to run PyHera for scaffolding:
    
    python pyhera.py scaffold contigs.fasta reads.fastq readsToContigs.paf readsToReads.paf -o output.fasta

Aside from overlaps between reads and contigs, PyHera also requires overlaps between reads and reads (self-overlaps).

### Running Scaffolding script
Scaffolding script combines Ezra and PyHera and is run using `<scaffold.py>` script from the root folder of the repository.

A typical command to run the scaffolding script:

    python scaffold.py contigs.fasta reads.fastq --results results_folder --plan P1E3

The script will combine Ezra and PyHera according to a given scaffolding plan. A scaffolding plan is a series of letter-number pairs, where letter defines the scaffolder and the number defines the number of iterations for that scaffolder. The scaffolding script will run Minimap2 to calcualte needed overlaps and will generate a separate folder for each iteration of scffolding. The final result will be in the folder relating to the last iteration.

If the scffolding plan is not specified, the default plan will be used. The default plan is E3P1 - run Ezra three times and then PyHera once. Ezra usually gives better results after a few iterations so we suggest running Ezra at least three times.

## Testing
Scaffolding script was tested on several Salmonella Enterica datasets. [Minimap2](https://github.com/lh3/minimap2) was used to generate overalps. Initial assembly was done using [Rala](https://github.com/rvaser/rala). Obtained contigs were then scaffolded using scaffolding script and several different scaffolding plans. Finally scaffolds were compared to the reference using [gepard tool](http://cube.univie.ac.at/gepard).

In general, all plans resulted in the improvement of the assembly, but the final quality heavily depended on the plan.

### Example
On one dataset Rala produced highly fragmented assembly containing 18 contigs. Gepard image is given below:
![Rala contigs](https://github.com/lbcb-sci/pyhera/blob/master/images/NCTC10384_rala_dotplot.jpeg|width=400)

![Scaffolding with plan P1E3](https://github.com/lbcb-sci/pyhera/blob/master/images/NCTC10384_P1E3_dotplot.jpeg)

![Scaffolding with plan E3P1](https://github.com/lbcb-sci/pyhera/blob/master/images/NCTC10384_E3P1_dotplot.jpeg)


