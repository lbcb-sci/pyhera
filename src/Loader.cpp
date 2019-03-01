#include "Loader.h"
#include "bioparser/bioparser.hpp"
#include "Sequence.h"
#include "Overlap.h"

#include <unordered_set>

namespace scara {
  using namespace std;

  void parseProcessFastq(const string& strFastq, MapIdToSeq& mIdToSeq) {
    vector<unique_ptr<Sequence>> aReads;
    auto fastqParser = bioparser::createParser<bioparser::FastqParser, scara::Sequence>(strFastq);
    fastqParser->parse_objects(aReads, -1);
    for (uint32_t i = 0; i < aReads.size(); i++) {
      mIdToSeq.emplace(aReads[i]->seq_strName, std::move(aReads[i]));
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
}
