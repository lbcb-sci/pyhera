# Scaffolding with ScaRa
This repository contains python scripts and binaries comprising our tools for genome scaffolding.

The tools are currently under development and any feedback on them is very welcome and will be greatly appreciated.

## Tools
### ScaRa
ScaRa is a python implementation of our scaffolding algorithm which uses some concepts from the HERA scaffolder paper (https://www.biorxiv.org/content/early/2018/06/13/345983). It uses overlaps between contigs and reads and also self-overalps between reads to generate a graph. It then generates numerous paths between contigs, groups them according to their length and choses the best representative path.

### Ezra
Ezra is a result of a graduation thesis work done by Ivan Krpelnik on University of Zareb Faculty of Electrical Engineering and Computing. Ezra uses MSA (Multiple Sequence Alignment) to extend contigs on each end and is also able to connect two contigs using long reads that overlap with both of them. To make the best use of the extension step, Ezra should be run in several iterations.
This thesis is available here: http://complex.zesoi.fer.hr/data/pdf/Ivan_Krpelnik_diplomski.pdf.

Ezra is included as a submodule, with main repository at: https://gitlab.com/Krpa/ezra.

### Scaffolding script
Scaffolding script combines ScaRa and Ezra, according to a given scaffolding plan, to iteratively perform the scaffolding, using output of the previous iteration as input for the next one. The script will use Minimap2 to produce overlaps needed for the scaffolding.

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

  4. Nuild Racon
  
    cd racon
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make

  Pythons scripts, such as ScaRa, Scaffolder script and samscripts tool do not need to be installed.

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

### Running ScaRa
ScaRa is run by running `<scara.py>` script, from the root folder of the repository. ScaRa can be run with several options and various arguments, which are printed to the screen when running the script without arguments.

A typical command to run ScaRa for scaffolding:
    
    python scara.py scaffold contigs.fasta reads.fastq readsToContigs.paf readsToReads.paf -o output.fasta

Aside from overlaps between reads and contigs, ScaRa also requires overlaps between reads and reads (self-overlaps).

### Running Scaffolding script
Scaffolding script combines Ezra and ScaRa and is run using `<scaffold.py>` script from the root folder of the repository.

A basic command to run the scaffolding script:

    python scaffold.py contigs.fasta reads.fastq --results results_folder

Available options will be printed out when running the script without arguments.

The scaffolding script will combine Ezra and ScaRa tools and run them in several iterations for more complete and correct scaffolding. The script will generate a separate folder for each iteration of scffolding. The final result will be in the folder relating to the last iteration. By default, the script will run Ezra 3 times and then ScarRa 1 time. We found that this combination worked reasonably well in our testing.

__Racon__
The scaffolding script can also be instructed run racon on the initial contigs, on the final scaffolds and after each iteration of scaffolding. This can be specified using options:
  - `<--racon-start>` - run racon on the initial contigs, before scaffolding
  - `<--racon-end>`   - run racon on the final scaffold
  - `<--racon>`       - run racon after each iteration

__Minimap2__
The scaffolding script will run Minimap2 to calcualte overlaps as needed by Ezra, Scara and possibly Racon.

## Testing
Scaffolding script was tested on two Salmonella Enterica datasets (NCTC10384 and NCTC12417). [Minimap2](https://github.com/lh3/minimap2) was used to generate overalps. Initial assembly was done using [Rala](https://github.com/rvaser/rala). Initial assembly was first corrected using [Racon](https://github.com/isovic/racon). Obtained contigs were then scaffolded using our scaffolding script. After each scaffolding iteration, Racon was also run to correct the results. Finally scaffolds were compared to the reference using using N50 measure and visually using [gepard tool](http://cube.univie.ac.at/gepard).

In both cases, scaffolding resulted in the improvement of the assembly, but the final quality is heavily dependant on the initial assembly.

### Example

|Dataset|No. sequences|Total length|N50|
|---|---|---|---|
|Reference| 3|5133713|5133713|
|NCTC10384 before ScaRa| 18|4770612|391888|
|NCTC10384 after ScaRa| 3|4775953|4478216|
|NCTC12417 before ScaRa| 5|5041995|3311937|
|NCTC12417 after ScaRa| 4|5031200|3481081|


<p float="left">
  <img src=NCTC10384_rala_racon_dotplot.jpeg width="400" />
  <img src=NCTC10384_scara_dotplot.jpeg width="400" /> 
</p>


<p float="left">
  <img src=NCTC12417_rala_dotplot.jpeg width="400" />
  <img src=NCTC12417_scara_dotplot.jpeg width="400" /> 
</p>


