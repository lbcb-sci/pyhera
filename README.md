# Scaffolding with ScaRa
This repository contains python scripts and binaries comprising our ScaRa tools for genome scaffolding.

The tool is still under development and any feedback on it is very welcome and will be greatly appreciated. We also do not recommend running it on larger genome since it is run on a single thread, is not optimised and might run for a long time. We are currently working on parallelizing all major processes.

ScaRa is based on Ezra scaffolding tool (http://complex.zesoi.fer.hr/data/pdf/Ivan_Krpelnik_diplomski.pdf) and uses come concepts from the HERA scaffolder paper (https://www.biorxiv.org/content/early/2018/06/13/345983). It consists of two phases, extension phase abd bridging phase, which can be run multiple time. Extension phase uses MSA (Multiple Sequence Alignment) and POA (Partial OrderAlignment) graphs, while bridging phase constructs paths betwewen contigs using reads that overlap and then tries to determine the best path.

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

## Testing
Scaffolding script was tested on two Salmonella Enterica datasets (NCTC10384 and NCTC12417). [Minimap2](https://github.com/lh3/minimap2) was used to generate overalps. Initial assembly was done using [Rala](https://github.com/rvaser/rala). Initial assembly was first corrected using [Racon](https://github.com/isovic/racon). Obtained contigs were then scaffolded using our scaffolding script. After each scaffolding iteration, Racon was also run to correct the results. Finally scaffolds were compared to the reference using using N50 measure and visually using [gepard tool](http://cube.univie.ac.at/gepard).

In both cases, scaffolding resulted in the improvement of the assembly, but the final quality is heavily dependant on the initial assembly. For the first dataset, the assembly is improved significantly, and very close to the reference, both in the number of contigs and their length. For the second dataset, scaffolding only connected two contigs, increasing the assembly quality only slightly.

Table with the results:

|Dataset|No. sequences|Total length|N50|
|---|---|---|---|
|Reference| 3|5133713|5133713|
|NCTC10384 before ScaRa| 18|4770612|391888|
|NCTC10384 after ScaRa| 3|4775953|4478216|
|NCTC12417 before ScaRa| 5|5041995|3311937|
|NCTC12417 after ScaRa| 4|5031200|3481081|

Gepard dotplots for dataset NCTC10384 against the reference (before the scaffolding - left, after the scaffolding - right):

<p float="left">
  <img src=images/NCTC10384_rala_racon_dotplot.jpeg width="400" />
  <img src=images/NCTC10384_scara_dotplot.jpeg width="400" /> 
</p>

Gepard dotplots for dataset NCTC12417 against the reference (before the scaffolding - left, after the scaffolding - right):

<p float="left">
  <img src=images/NCTC12417_rala_dotplot.jpeg width="400" />
  <img src=images/NCTC12417_scara_dotplot.jpeg width="400" /> 
</p>
