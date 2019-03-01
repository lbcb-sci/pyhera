#include "Types.h"

namespace scara {

	extern void parseProcessFastq(const std::string& strFastq, MapIdToSeq& mIdToSeq);
	extern void parseProcessFasta(const std::string& strFasta, MapIdToSeq& mIdToSeq);
	extern void parseProcessPaf(const std::string& strPaf, MapIdToOvl& mIdToOvl);

}
