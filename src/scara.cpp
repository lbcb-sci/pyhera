#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <map>
#include <fstream>
#include <algorithm>
#include <set>
#include <iterator>
#include <unordered_set>

#include <unistd.h>
#include <ios>
#include <ctime>

#include <sys/time.h>
#include <sys/resource.h>

#include <getopt.h>

#include "scara.h"
#include "SBridger.h"
#include "globals.h"

using namespace std;
using namespace scara;

namespace scara {

int multithreading;

uint32_t MinMCPaths, HardNodeLimit, SoftNodeLimit;
uint32_t NumDFSNodes, MaxMCIterations;
uint32_t MinPathsinGroup;
float SImin, OHmax;

bool test_short_length, test_contained_reads, test_low_quality;
bool print_output;
double PathGroupHalfSize;

std::string logFile;

DebugLevel globalDebugLevel;
}

// For tracking memory usage withing the program
// Taken from:
// https://www.tutorialspoint.com/how-to-get-memory-usage-at-runtime-using-cplusplus
void mem_usage(double& vm_usage, double& resident_set) {
   vm_usage = 0.0;
   resident_set = 0.0;
   ifstream stat_stream("/proc/self/stat",ios_base::in); //get info from proc directory
   //create some variables to get info
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;
   unsigned long vsize;
   long rss;
   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
   >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
   >> utime >> stime >> cutime >> cstime >> priority >> nice
   >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest
   stat_stream.close();
   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // for x86-64 is configured to use 2MB pages
   vm_usage = vsize / 1024.0;
   resident_set = rss * page_size_kb;
}

void print_mem_usage(std::string location) {
   double vm, rss;
   mem_usage(vm, rss);
   cerr << "\nMemory usage at: " << location;
   cerr << "\nVirtual Memory: " << vm << "\nResident set size: " << rss << endl;
}

// Setting global parameters, variables are defined in globals.h
void setGlobalParameters() {

  // Global variable that determines whether to use multithreading or not
  // Not using multithreading unless specified in the parameters
  scara::multithreading = 0;

  // A minimum number of paths generated by Monte Carlo method
  scara::MinMCPaths = 40;
  // A maximum number of nodes that a path can contain, longer paths will not be generated
  scara::HardNodeLimit = 1000;
  // A number of nodes in a path that will cause a warning to be displayed
  scara::SoftNodeLimit = 100;

  // A number of nodes added to the stack in each step of DFS graph traversal
  scara::NumDFSNodes = 5;
  // Maximum number of iterations using Monte Carlo approach
  scara::MaxMCIterations = 100000;

  // A minimum number of paths in a group that will generate a scaffold
  scara::MinPathsinGroup = 3;

  // Minimum sequence identity for filtering overlaps
  scara::SImin = 0.60;
  // Maximum allowed overhang percentage for filtering overlaps
  scara::OHmax = 0.25;

  scara::test_short_length = true;
  scara::test_contained_reads = true;
  scara::test_low_quality = true;

  // Determines if more verbose messages are printed
  // TODO: replace this with DEBUG_LEVEL with multiple levels
  scara::print_output = true;

  // A path is placed in a group if its length falls within pathGtoupHalfSize of groups representative length
  scara::PathGroupHalfSize = 5000;

  // Setting default Debugg level
  scara::globalDebugLevel = DL_INFO;

  // Setting Log file
  scara::logFile = "scaraLog.txt";
}

DebugLevel debugLevelFromString(std::string str) {
	if (str == "0") return DL_NONE;
	else if (str == "1") return DL_INFO;
	else if (str == "2") return DL_VERBOSE;
	else if (str == "3") return DL_DEBUG;
	else return DL_NONE;
}

void print_version_message_and_exit() {
  const char* versionmessage = "\nScaRa version 1.3!\n";
  std::cerr << versionmessage;
  exit(0);
}

void print_help_message_and_exit() {
  const char* helpmessage = "\nHow to run ScaRa :"
    "\nScaRa -f (--folder) <Input folder>"
    "\n    or"
    "\nScaRa -r <Reads file> -c <Contigs file> -o <Reads to Contigs Overlaps file> -s <Reads to Reds Overlaps file>"
    "\n    The program will perform one iteration of the algorithm"
    "\n    and output contigs to the standard output!."
    "\n    <Input foder> must contain the following files:"
    "\n    - reads.fastq - reads in FASTQ/FASTA format"
    "\n    - readsToContigs.paf - overlaps between reads and contigs"
    "\n    - readsToReads.paf - overlaps between reads"
    "\n    - contigs.fasta - contigs in FASTA format"
    "\nIf Reads file, Contigs file or Overlaps are specified,"
	  "\nthey will not be looked for in the Input folder!"
    "\nGeneral options:"
    "\n-f (--folder)     specify input folder for ScaRa"
    "\n-r (--reads)      specify reads file for ScaRa"
    "\n-c (--contigs)    specify contigs file ScaRa"
    "\n-o (--overlapsRC)   specify contig-read overlaps file for ScaRa"
    "\n-s (--overlapsRR)   specify read self overlaps file for ScaRa"
    "\n-m (--multithreading)   use multithreading"
    "\n-D (--debug_level) [level] set a debugg level which determines "
    "\n 					the amount of output the program generates to stderr"
    "\n 					level can be set to values 0 - 3, with 0 being the least"
    "\n           and 3 being the most verbose"
    "\n________________________________________________________________________"
    "\nAlgorithm parameter options:"
    "\npMinMCPaths - minimum number of paths that will try to be generated using"
    "\n              Monte Carlo method (defult 40)."
    "\npMAXMCIterations - maximum nbumber of iterations during Monte Carlo (defualt 100000)"
    "\npHardNodeLimit - a maximum number of nodes in a path (defult 1000)"
    "\npNumDFSNodes - a number of nodes that will be placed on the stack during"
    "\n               DSF search of the graph (defult 5)"
    "\npMinPathsInGroup - a minimum number of paths in a path group (default 3)"
    "\npPathGroupHalfSize - a size of a bucket in which paths are grouped"
    "\n                     according to their length (default 5000)"
    "\n________________________________________________________________________"
    "\n-v (--version)    print program version"
    "\n-h (--help)       print this help message\n";

  std::cerr << helpmessage;
  exit(0);
}

int main(int argc, char **argv)
{

  // KK: Defining basic program nOptions
  int pMinMCPaths, pMaxMCIterations, pHardNodeLimit, pNumDFSNodes, pMinPathsInGroup, pPathGroupHalfSize;
  const char* const short_opts = "hvr:c:o:s:f:mD:";
  const option long_opts[] = {
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'v'},
    {"folder", required_argument, NULL, 'f'},
    {"reads", required_argument, NULL, 'r'},
    {"contigs", required_argument, NULL, 'c'},
    {"overlapsRC", required_argument, NULL, 'o'},
    {"overlapsRR", required_argument, NULL, 's'},
    {"multithreading", no_argument, NULL, 'm'},
    {"debug_level", required_argument, NULL, 'D'},
    {"pMinMCPaths", required_argument, &pMinMCPaths, 0},
    {"pMAXMCIterations", required_argument, &pMaxMCIterations, 0},
    {"pHardNodeLimit", required_argument, &pHardNodeLimit, 0},
    {"pNumDFSNodes", required_argument, &pNumDFSNodes, 0},
    {"pMinPathsInGroup", required_argument, &pMinPathsInGroup, 0},
    {"pPathGroupHalfSize", required_argument, &pPathGroupHalfSize, 0},
    {NULL, no_argument, NULL, 0}
  };

  int readsSet, contigsSet, RCOverlapsSet, RROverlapsSet;
  readsSet = contigsSet = RCOverlapsSet = RROverlapsSet = 0;
  string strData, strReadsFasta, strContigsFasta, strR2COvlPaf, strR2ROvlPaf;

  setGlobalParameters();

  int opt;
  while ((opt = getopt_long(argc, argv, short_opts, long_opts, NULL)) != -1)
    switch(opt) {
    case 'h':
      print_help_message_and_exit();
      break;
    case 'v':
      print_version_message_and_exit();
      break;
    case 'f':
      readsSet = contigsSet = RCOverlapsSet = RROverlapsSet = 1;
      strData = optarg;
      strReadsFasta = strData + "/reads.fastq";
      strR2COvlPaf = strData + "/readsToContigs.paf";
      strR2ROvlPaf = strData + "/readsToReads.paf";
      strContigsFasta = strData + "/contigs.fasta";
      break;
    case 'r':
      readsSet = 1;
      strReadsFasta = optarg;
      break;
    case 'c':
      contigsSet = 1;
      strContigsFasta = optarg;
      break;
    case 'o':
      RCOverlapsSet = 1;
      strR2COvlPaf = optarg;
      break;
    case 's':
      RROverlapsSet = 1;
      strR2ROvlPaf = optarg;
      break;
    case 'm':
      scara::multithreading = 1;
      break;
    case 'D':
      scara::globalDebugLevel = debugLevelFromString(optarg);
      break;
    case 0:
      if (pMinMCPaths) scara::MinMCPaths = stoi(optarg);
      if (pMaxMCIterations) scara::MaxMCIterations = stoi(optarg);
      if (pHardNodeLimit) scara::HardNodeLimit = stoi(optarg);
      if (pNumDFSNodes) scara::NumDFSNodes = stoi(optarg);
      if (pMinPathsInGroup) scara::MinPathsinGroup = stoi(optarg);
      if (pPathGroupHalfSize) scara::PathGroupHalfSize = stoi(optarg);
      break;
    default:
      print_help_message_and_exit();
    }


  if (argc < 2 || readsSet == 0 || contigsSet == 0 || RCOverlapsSet == 0 || RROverlapsSet == 0) {
     std::cerr << "\nNot all arguments specified!";
     print_help_message_and_exit();
  }

  std::time_t start_time = std::time(nullptr);
  std::cerr << "\nSCARA: Starting ScaRa: " << std::asctime(std::localtime(&start_time)) << start_time << " seconds since the Epoch";
  if (scara::globalDebugLevel == DL_DEBUG) print_mem_usage("Start");

  std::time_t current_time = std::time(nullptr);
  std::cerr << "\nSCARA: Finished loading: " << std::asctime(std::localtime(&current_time)) << (current_time - start_time) << " seconds since the start\n";

  scara::SBridger sbridger(strReadsFasta, strContigsFasta, strR2COvlPaf, strR2ROvlPaf);
  if (scara::globalDebugLevel >= DL_VERBOSE) {
    print_mem_usage("After initialization");
    sbridger.printState();
  }

  current_time = std::time(nullptr);
  std::cerr << "\nSCARA: Finished initializing bridger: " << std::asctime(std::localtime(&current_time)) << (current_time - start_time) << " seconds since the start\n";
  std::cerr << "\nSCARA: Generating graph:";

  sbridger.generateGraph();
  sbridger.print();

  current_time = std::time(nullptr);
  std::cerr << "\nSCARA: Finished generating graph: " << std::asctime(std::localtime(&current_time)) << (current_time - start_time) << " seconds since the start\n";
  std::cerr << "\nSCARA: Cleaning up the graph:";

  sbridger.cleanupGraph();
  sbridger.print();
  if (scara::globalDebugLevel >= DL_VERBOSE) {
    print_mem_usage("After graph generation");
    sbridger.printState();
  }

  current_time = std::time(nullptr);
  std::cerr << "\nSCARA: Finished cleaning up the graph: " << std::asctime(std::localtime(&current_time)) << (current_time - start_time) << " seconds since the start\n";
  std::cerr << "\nSCARA: Generating paths:";

  int numPaths = sbridger.generatePaths();

  if (scara::globalDebugLevel >= DL_VERBOSE) {
    print_mem_usage("After path generation");
    sbridger.printState();
  }
  current_time = std::time(nullptr);
  std::cerr << "\nSCARA: Finished generating paths: " << std::asctime(std::localtime(&current_time)) << (current_time - start_time) << " seconds since the start\n";
  std::cerr << "\nSCARA: Paths generated: " << numPaths;

  std::cerr << "\nSCARA: printing paths:";
  sbridger.printPaths();

  std::cerr << "\nSCARA: Grouping and processing paths:";

  int numGroups = sbridger.groupAndProcessPaths();

  if (scara::globalDebugLevel >= DL_VERBOSE) {
    print_mem_usage("After path grouping");
    sbridger.printState();
  }
  current_time = std::time(nullptr);
  std::cerr << "\nSCARA: Finished grouping and processing paths: " << std::asctime(std::localtime(&current_time)) << (current_time - start_time) << " seconds since the start\n";
  std::cerr << "\nSCARA: Final number of path groups: " << numGroups;
  std::cerr << "\nSCARA: Generating sequences:";

  int numSeq = sbridger.generateSequences();

  current_time = std::time(nullptr);
  std::cerr << "\nSCARA: Finished generating sequences: " << std::asctime(std::localtime(&current_time)) << (current_time - start_time) << " seconds since the start\n";
  std::cerr << "\nSCARA: Final number of sequences generated: " << numSeq;

  if (scara::globalDebugLevel >= DL_VERBOSE) {
    print_mem_usage("After sequence generation");
    sbridger.printState();
  }
  current_time = std::time(nullptr);
  std::cerr << "\nSCARA: Finished ScaRa run: " << std::asctime(std::localtime(&current_time)) << (current_time - start_time) << " seconds since the start\n";

  return 0;
}
