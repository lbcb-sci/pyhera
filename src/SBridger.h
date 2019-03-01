#include "Types.h"
#include "Graph.h"
#include <string>

namespace scara {

	using namespace std;

	class SBridger {
	private:
		MapIdToOvl mIdToOvlR2C;
	  	MapIdToOvl mIdToOvlR2R;
	  	MapIdToSeq mIdToContig;
	  	MapIdToSeq mIdToRead;

	  	// Statistical information
	  	uint32_t numANodes;
	  	uint32_t numRNodes;

	  	uint32_t numAREdges_all;
	  	uint32_t numAREdges_usable;
	  	uint32_t numAREdges_contained;
	  	uint32_t numAREdges_short;
	  	uint32_t numAREdges_lowqual;
	  	uint32_t numAREdges_zero;

	  	uint32_t numRREdges_all;
	  	uint32_t numRREdges_usable;
	  	uint32_t numRREdges_contained;
	  	uint32_t numRREdges_short;
	  	uint32_t numRREdges_lowqual;
	  	uint32_t numRREdges_zero;

	  	uint32_t isolatedANodes = 0;
	  	uint32_t isolatedRNodes = 0;

	  	// Paths generated through the graph
	  	std::vector<shared_ptr<Path>> vPaths;

	  	// Path info vector, for faster analysis
	  	std::vector<shared_ptr<PathInfo>> vPathInfos;


	public:
		int bGraphCreated;
		MapIdToNode mAnchorNodes;
		MapIdToNode mReadNodes;
		std::vector<shared_ptr<Edge>> vEdges;

		SBridger(const string& strReadsFasta, const string& strContigsFasta, const string& strR2Cpaf, const string& strR2Rpaf);

	  	void Initialize(const string& strReadsFasta, const string& strContigsFasta, const string& strR2Cpaf, const string& strR2CRaf);

	  	void printData(void);

	  	void printGraph(void);

	  	void print(void);

	  	void generateGraph(void);

	  	void cleanupGraph(void);

	  	int generatePaths(void);

	  	int groupAndProcessPaths(void);

	  	int generateSequences(void);

	  	void Execute(void);
	};

}
