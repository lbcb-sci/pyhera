#include "Types.h"
#include "Graph.h"
#include "Overlap.h"
#include "Loader.h"
#include <string>
#include <fstream>

namespace scara {

	using namespace std;

	class SBridger {
	private:
		VecOvl vOvlR2C;
	  	VecOvl vOvlR2R;
	  	MapIdToSeq mIdToContig;
	  	MapIdToSeq mIdToRead;

	  	// Statistical information
	  	uint32_t numANodes;
	  	uint32_t numRNodes;

	  	uint32_t numEdges_all;
	  	uint32_t numEdges_usable;
	  	uint32_t numEdges_contained;
	  	uint32_t numEdges_short;
	  	uint32_t numEdges_lowqual;
	  	uint32_t numEdges_zero;

	  	uint32_t isolatedANodes = 0;
	  	uint32_t isolatedRNodes = 0;

	  	// Paths generated through the graph
	  	std::vector<shared_ptr<Path>> vPaths;

	  	// Path info vector, for faster analysis
	  	std::vector<shared_ptr<PathInfo>> vPathInfos;

	  	// Path group vector, for faster analysis
	  	std::vector<shared_ptr<PathGroup>> vPathGroups;

	  	// Final scaffolds
	  	std::vector<shared_ptr<std::vector<shared_ptr<PathInfo>>>> scaffolds;


	public:
		int bGraphCreated;
		MapIdToNode mAnchorNodes;
		MapIdToNode mReadNodes;
		// std::vector<shared_ptr<Edge>> vEdges;

		SBridger(const string& strReadsFasta, const string& strContigsFasta, const string& strR2Cpaf, const string& strR2Rpaf);

	  	void Initialize(const string& strReadsFasta, const string& strContigsFasta, const string& strR2Cpaf, const string& strR2CRaf);

	  	void printData(void);

	  	void printGraph(void);

	  	void printPaths(void);

	  	void print(void);

	  	void generateGraph(void);

	  	void cleanupGraph(void);

	  	int generatePaths(void);

	  	int groupAndProcessPaths(void);

	  	int generateSequences(void);

	  	void Execute(void);

	  	void printState();

	private:
		shared_ptr<PathInfo> getBestPath_AvgSI();

		void printOvlToStream(VecOvl &vOvl, ofstream& outStream);
		void printNodeToStream(MapIdToNode &map, ofstream& outStream);

		bool scaffoldsEqual(shared_ptr<std::vector<shared_ptr<PathGroup>>> scaff1, shared_ptr<std::vector<shared_ptr<PathGroup>>> scaff2);
		// int compareScaffolds(shared_ptr<std::vector<shared_ptr<PathGroup>>> scaff1, shared_ptr<std::vector<shared_ptr<PathGroup>>> scaff2);
	};

}
