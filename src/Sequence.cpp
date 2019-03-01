#include "Sequence.h"
#include <string>

namespace scara {

  char _bioBaseComplement(char c) {
    if (c == 'A') return 'T';
    if (c == 'T') return 'A';
    if (c == 'G') return 'C';
    if (c == 'C') return 'G';

    return c;
  }

  std::string _bioReverseComplement(const std::string& strData, uint32_t ulStart, uint32_t ulEnd)
  {
    std::string strRc(ulEnd - ulStart, 'x');
    for(uint32_t ulPos = ulStart; ulPos < ulEnd; ulPos++) {
      strRc[ulEnd - ulPos - 1] = _bioBaseComplement(strData[ulPos]);
    }

    return strRc;
  }

  std::string _bioReverseComplement(const std::string& strData)
  {
    return _bioReverseComplement(strData, 0, strData.size());
  }

  Sequence::Sequence(const char* name, uint32_t name_length,
    const char* sequence, uint32_t sequence_length
  ) : seq_strName(name, name_length), seq_strData(sequence, sequence_length)
  {
  }

  Sequence::Sequence(const char* name, uint32_t name_length,
    const char* sequence, uint32_t sequence_length,
    const char* quality, uint32_t quality_length
  ) : seq_strName(name, name_length), seq_strData(sequence, sequence_length) , seq_strQuality(quality, quality_length)
  {
  }

}
