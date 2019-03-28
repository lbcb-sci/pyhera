#pragma once

#include <string>

#include "globals.h"

namespace scara {

  enum SequenceStrand {
    SS_INVALID = 0,
    SS_FORWARD = 1 << 0,
    SS_REVERSE = 1 << 1
  };

  enum OverlapType {
    ET_INVALID = 0,
    ET_PREFIX = 1 << 0,
    ET_SUFFIX = 1 << 1,
    ET_CONTAINED = 1 << 2
  };

  extern SequenceStrand reverseStrand(SequenceStrand strand);

  extern OverlapType GetOtherExtension(OverlapType oType);

  extern std::string StringFromExtensionType(OverlapType oType);

  class Overlap {
  public:
    OverlapType ext_oType;
    std::string ext_strName;
    std::string ext_strTarget;
    uint32_t ext_ulQBegin;
    uint32_t ext_ulQEnd;
    uint32_t ext_ulQLen;
    uint32_t ext_ulTBegin;
    uint32_t ext_ulTEnd;
    uint32_t ext_ulTLen;
    bool ext_bOrientation; // true '+'
    uint32_t ext_ulTargetDelta;
    uint32_t ext_ulOverlapLength;
    uint32_t ext_ulExtLen;

    uint32_t paf_matching_bases;
    uint32_t paf_overlap_length;
    uint32_t paf_mapping_quality;

    Overlap(
      const char* q_name, uint32_t q_name_length,
      uint32_t q_length,
      uint32_t q_begin,
      uint32_t q_end,
      char orientation,
      const char* t_name, uint32_t t_name_length,
      uint32_t t_length,
      uint32_t t_begin,
      uint32_t t_end,
      uint32_t matching_bases,
      uint32_t overlap_length,
      uint32_t mapping_quality
    );

    void GetExtensionPosition(uint32_t& ulBegin, uint32_t& ulEnd);

    bool Test(void);
  };


  // KK: copied from Ezra, probably not needed
  /*
  struct ExtBridge {
    uint32_t eb_ulSufBegin; // where it begins on suffix (prev contig)
    uint32_t eb_ulPreEnd; // where it ends on prefix (next contig)
    uint32_t eb_ulReadBegin;
    uint32_t eb_ulReadEnd;
    std::string eb_strRead; // read id
    std::string eb_strPrevContig;
    std::string eb_strNextContig;
    bool eb_bOrientation;
    uint32_t eb_ulCoverage;

    ExtBridge(const Extension& extSuf, const Extension& extPre);
  };
  */

}
