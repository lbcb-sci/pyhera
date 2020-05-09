#include "Overlap.h"
#include <string>

namespace scara {


	SequenceStrand reverseStrand(SequenceStrand strand) {
		if (strand == SS_FORWARD) return SS_REVERSE;
		else if (strand == SS_REVERSE) return SS_FORWARD;
		else return strand;
	}

  std::string StringFromOverlapType(OverlapType eType) {
    if (eType == ET_PREFIX) return "EXT_PREFIX";
    if (eType == ET_SUFFIX) return "EXT_SUFFIX";
    return "EXT_INVALID";
  }

  Overlap::Overlap(
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
  ) : ext_oType(ET_INVALID), ext_strName(q_name, q_name_length), ext_strTarget(t_name, t_name_length),
        ext_ulQBegin(q_begin), ext_ulQEnd(q_end), ext_ulQLen(q_length), 
        ext_ulTBegin(t_begin), ext_ulTEnd(t_end), ext_ulTLen(t_length), 
        ext_bOrientation(orientation == '+'), ext_ulOverlapLength(overlap_length),
        paf_matching_bases(matching_bases), paf_overlap_length(overlap_length), paf_mapping_quality(mapping_quality)
  {
    const uint32_t ulTOverlap = t_end - t_begin;
    const uint32_t ulQOverlap = q_end - q_begin;
    const uint32_t ulTohR = t_length - t_end;
    const uint32_t ulTohL = t_begin;
    const uint32_t ulQohR = ext_bOrientation ? q_length - q_end : q_begin;
    const uint32_t ulQohL = ext_bOrientation ? q_begin : q_length - q_end;

    const uint32_t ulThreshold = std::min(ulTOverlap, ulQOverlap) * 0.125;
    if (ulQohR > ulQohL && ulQohR > ulTohR && std::max(ulTohR, ulQohL) < ulThreshold) {
      ext_oType = ET_SUFFIX;
      ext_ulTargetDelta = t_length - t_begin;
      ext_ulExtLen = ulQohR;
    } else if (ulQohL > ulQohR && ulQohL > ulTohL && std::max(ulTohL, ulQohR) < ulThreshold) {
      ext_oType = ET_PREFIX;
      ext_ulTargetDelta = t_end;
      ext_ulExtLen = ulQohL;
    } else if(ulQOverlap > ulQohL + ulQohR && ulQohR < ulTohR && ulQohL < ulTohL) {
      ext_oType = ET_CONTAINED;
    }
  }


   void Overlap::GetExtensionPosition(uint32_t& ulBegin, uint32_t& ulEnd) {
    if ((ext_oType == ET_PREFIX && ext_bOrientation) || (ext_oType == ET_SUFFIX && !ext_bOrientation)) {
      ulBegin = 0; ulEnd = ext_ulQEnd;
    } else {
      ulBegin = ext_ulQBegin; ulEnd = ext_ulQLen;
    }
  }

  // Check is an Overlap is suitable to be used as an Edge
  bool Overlap::Test() {
  	return true;
  }

  /*
  ExtBridge::ExtBridge(const Extension& extSuf, const Extension& extPre) :
    eb_ulSufBegin(extSuf.ext_ulTBegin), eb_ulPreEnd(extPre.ext_ulTEnd), eb_strRead(extSuf.ext_strName),
    eb_strPrevContig(extSuf.ext_strTarget), eb_strNextContig(extPre.ext_strTarget),
    eb_bOrientation(extSuf.ext_bOrientation), eb_ulCoverage(0) {
    if (extSuf.ext_bOrientation) {
      eb_ulReadBegin = extSuf.ext_ulQBegin;
      eb_ulReadEnd = extPre.ext_ulQEnd;
    } else {
      eb_ulReadBegin = extPre.ext_ulQBegin;
      eb_ulReadEnd = extSuf.ext_ulQEnd;
    }
  }
  */

}
