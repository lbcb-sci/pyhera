#pragma once 

namespace scara {

enum DebugLevel {
    DL_NONE = 0,
    DL_INFO = 1 << 0,
    DL_VERBOSE = 1 << 1,
    DL_DEBUG = 1<< 2,
  };

extern DebugLevel globalDebugLevel;

extern int multithreading;

extern int MinMCPaths, HardNodeLimit, SoftNodeLimit;
extern int numDFSNodes, maxMCIterations;
extern float SImin, OHmax;

extern bool test_short_length, test_contained_reads, test_low_quality;

extern bool print_output;

extern double pathGroupHalfSize;

}
