#pragma once

#include <string>

namespace scara {

  extern char _bioBaseComplement(char c);

  extern std::string _bioReverseComplement(const std::string& strData, uint32_t ulStart, uint32_t ulEnd);

  extern std::string _bioReverseComplement(const std::string& strData);

  class Sequence {
  public:
    std::string seq_strName;
    std::string seq_strData;
    std::string seq_strQuality;

  public:
    Sequence(const char* name, uint32_t name_length,
      const char* sequence, uint32_t sequence_length
    );

    Sequence(const char* name, uint32_t name_length,
      const char* sequence, uint32_t sequence_length,
      const char* quality, uint32_t quality_length
    );
  };

}
