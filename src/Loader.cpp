#include "Loader.h"
#include "bioparser/bioparser.hpp"
#include "Sequence.h"
#include "Overlap.h"

#include <iostream>

#include <unordered_set>

namespace scara {
  using namespace std;

  void parseProcessFastq(const string& strFastq, MapIdToSeq& mIdToSeq) {
    vector<unique_ptr<Sequence>> aReads;
    auto fastqParser = bioparser::createParser<bioparser::FastqParser, scara::Sequence>(strFastq);
    fastqParser->parse_objects(aReads, -1);
    for (uint32_t i = 0; i < aReads.size(); i++) {
      // Loosing quals, to save space
      // KK: Since I'm a newb and don know how exactly this works, I'm extending length by 1, for \0 in case its needed
      unique_ptr<Sequence> newSeq_ptr = make_unique<Sequence>(aReads[i]->seq_strName.c_str(), (uint32_t)(aReads[i]->seq_strName.length()),
                                                              aReads[i]->seq_strData.c_str(), (uint32_t)(aReads[i]->seq_strData.length()));

      mIdToSeq.emplace(newSeq_ptr->seq_strName, std::move(newSeq_ptr));
    }
  }

  void parseProcessFasta(const string& strFasta, MapIdToSeq& mIdToSeq) {
    vector<unique_ptr<Sequence>> aReads;
    auto fastaParser = bioparser::createParser<bioparser::FastaParser, scara::Sequence>(strFasta);
    fastaParser->parse_objects(aReads, -1);
    for (uint32_t i = 0; i < aReads.size(); i++) {
      mIdToSeq.emplace(aReads[i]->seq_strName, std::move(aReads[i]));
    }
  }

  void parseProcessPaf(const string& strPaf, MapIdToOvl& mIdToOvl) {
    vector<unique_ptr<scara::Overlap>> aExtensions;
    auto pafParser = bioparser::createParser<bioparser::PafParser, Overlap>(strPaf);
    pafParser->parse_objects(aExtensions, -1);

    // preprocessing - mark contained reads
    unordered_set<string> usContained;
    for (const auto& ext : aExtensions) {
      if (ext->ext_oType == ET_CONTAINED) {
        usContained.emplace(ext->ext_strName);
      }
    }

    for (uint32_t i = 0; i < aExtensions.size(); i++) {
      if (aExtensions[i]->ext_oType != ET_INVALID && usContained.find(aExtensions[i]->ext_strName) == usContained.end()) {
        mIdToOvl[{aExtensions[i]->ext_strTarget, aExtensions[i]->ext_oType}].emplace_back(std::move(aExtensions[i]));
      }
    }
  }

  extern void parseProcessPaf(const std::string& strPaf, VecOvl& vOvl) {
    vector<unique_ptr<scara::Overlap>> aExtensions;
    auto pafParser = bioparser::createParser<bioparser::PafParser, Overlap>(strPaf);
    pafParser->parse_objects(aExtensions, -1);

    for (uint32_t i = 0; i < aExtensions.size(); i++) {
      vOvl.emplace_back(std::move(aExtensions[i]));
    }
  }
}
