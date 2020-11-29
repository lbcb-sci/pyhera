#include <iostream>
#include <string>

#include <getopt.h>

#include "Types.h"
#include "Loader.h"


void print_version_message_and_exit() {
  const char* versionmessage = "\nLoad_HiC version 0.1!\n";
  std::cerr << versionmessage;
  exit(0);
}

void print_help_message_and_exit() {
  const char* helpmessage = "\nHow to run Load_HiC :"
    "\nLoad_HiC -1 (--Hi-C1) <First Hi-C file> -2 (--Hi-C2) <Second Hi-C file> [-r (--reference) <Reference file>]"
    "\n    or"
    "\nLoad_HiC -c (--Hi-C) <Combined Hi-C file> [-r (--reference) <Reference file>]"
    "\n"
    "\n-v (--version)    print program version"
    "\n-h (--help)       print this help message\n";

  std::cerr << helpmessage;
  exit(0);
}

void analyzeCombinedHiC(scara::MapIdToSeq mIdToReadHiC, scara::MapIdToSeq mIdToReadReference, int refSet) {
  std::cout << "\n\nAnalyzing combined Hi-C reads ...";
  std::cout << "\nHi-C reads file size: " << mIdToReadHiC.size();
  if (refSet) {
    std::cout << "\nReference file size: " << mIdToReadReference.size();
  }
  else {
  	std::cout << "\nReference NOT SET!";
  }
}

void analyzeSeparateHiC(scara::MapIdToSeq mIdToReadHiC1, scara::MapIdToSeq mIdToReadHiC2, scara::MapIdToSeq mIdToReadReference, int refSet) {
  std::cout << "\n\nAnalyzing separate Hi-C reads ...";
  std::cout << "\nFirst Hi-C reads file size: " << mIdToReadHiC1.size();
  std::cout << "\nSecond Hi-C reads file size: " << mIdToReadHiC2.size();
  if (refSet) {
    std::cout << "\nReference file size: " << mIdToReadReference.size();
  }
  else {
  	std::cout << "\nReference NOT SET!";
  }
}

int main(int argc, char **argv)
{

  // KK: Defining basic program options
  const char* const short_opts = "hvr:c:1:2:";
  const option long_opts[] = {
    {"help", no_argument, NULL, 'h'},                   // option_index = 0
    {"version", no_argument, NULL, 'v'},                // option_index = 1
    {"Hi-C", required_argument, NULL, 'c'},             // option_index = 2
    {"Hi-C1", required_argument, NULL, '1'},            // option_index = 3
    {"Hi-C2", required_argument, NULL, '2'},            // option_index = 4
    {"reference", required_argument, NULL, 'r'},        // option_index = 5
    {NULL, no_argument, NULL, 0}
  };

  int HiCSet, HiC1Set, HiC2Set, refSet;
  HiCSet = HiC1Set = HiC2Set = refSet = 0;
  std::string strHiCFQ, strHiC1FQ, strHiC2FQ, strRefFA;

  int opt, option_index;
  while ((opt = getopt_long(argc, argv, short_opts, long_opts, &option_index)) != -1)
    switch(opt) {
    case 'h':
      print_help_message_and_exit();
      break;
    case 'v':
      print_version_message_and_exit();
      break;
    case 'r':
      refSet = 1;
      strRefFA = optarg;
      break;
    case 'c':
      HiCSet = 1;
      strHiCFQ = optarg;
      break;
    case '1':
      HiC1Set = 1;
      strHiC1FQ = optarg;
      break;
    case '2':
      HiC2Set = 1;
      strHiC2FQ = optarg;
      break;
    default:
      print_help_message_and_exit();
    }


  if (HiCSet == 0 && HiC1Set == 0 && HiC2Set == 0) {
     std::cerr << "\nHi-C file or files not specified!\n";
     print_help_message_and_exit();
  }

  if (HiC1Set + HiC2Set == 1) {
     std::cerr << "\nBoth Hi-C files must be specified!\n";
     print_help_message_and_exit();
  }

  scara::MapIdToSeq mIdToReadHiC;
  scara::MapIdToSeq mIdToReadHiC1;
  scara::MapIdToSeq mIdToReadHiC2;
  scara::MapIdToSeq mIdToReadReference;

  if (refSet) {
    parseProcessFasta(strRefFA, mIdToReadReference); 
  }

  if (HiCSet) {
    parseProcessFastq(strHiCFQ, mIdToReadHiC);
    analyzeCombinedHiC(mIdToReadHiC, mIdToReadReference, refSet);
  }
  else if (HiC1Set && HiC2Set) {
    parseProcessFastq(strHiC1FQ, mIdToReadHiC1);
    parseProcessFastq(strHiC2FQ, mIdToReadHiC2);
    analyzeSeparateHiC(mIdToReadHiC1, mIdToReadHiC2, mIdToReadReference, refSet);
  }
  else {
  	 std::cerr << "\nHi-C file or files not specified!\n";
     print_help_message_and_exit();
  }

  return 0;
}
