#include "SBridger.h"
#include "globals.h"
#include <vector>
// #include <string>
#include <iostream>
#include <set>

namespace scara {

  using namespace std;

  extern int multithreading;

  SBridger::SBridger(const string& strReadsFasta, const string& strContigsFasta, const string& strR2Cpaf, const string& strR2Rpaf) {
    Initialize(strReadsFasta, strContigsFasta, strR2Cpaf, strR2Rpaf);
    bGraphCreated = 0;
  }

  void SBridger::Initialize(const string& strReadsFasta, const string& strContigsFasta, const string& strR2Cpaf, const string& strR2CRaf) {
      parseProcessFastq(strReadsFasta, mIdToRead);
      parseProcessFasta(strContigsFasta, mIdToContig);
      parseProcessPaf(strR2Cpaf, mIdToOvlR2C);
      parseProcessPaf(strR2CRaf, mIdToOvlR2R);
  }

	void SBridger::printData(void){
	  std::cerr << "\nLoaded data\n";
      std::cerr << "Number of contigs:" << mIdToContig.size() << '\n';
      std::cerr << "Number of reads:" << mIdToRead.size() << '\n';
      std::cerr << "Number of read-to-contig overlaps:" << mIdToOvlR2C.size() << '\n';
      std::cerr << "Number of read-to-read overlaps:" << mIdToOvlR2R.size() << '\n';
	}


	void SBridger::printGraph(void){
		std::cerr << "\nConstructed graph:\n";
		if (bGraphCreated == 0) {
			std::cerr << "Graph is not yet constructed!\n";
		}
		else {
	        std::cerr << "Number of anchor nodes:" << mAnchorNodes.size() << '\n';
	        std::cerr << "Number of isolated anchor nodes:" << isolatedANodes << '\n';
	        std::cerr << "Number of read nodes:" << mReadNodes.size() << '\n';
	        std::cerr << "Number of isolated read nodes:" << isolatedRNodes << '\n';
	        std::cerr << "Number of edges:" << vEdges.size() << '\n';

	        std::cerr << "\nAnchor node - read node edges:" << '\n';
			std::cerr << "Usable: " << numAREdges_usable << '\n';
			std::cerr << "Contained: " << numAREdges_contained << '\n';
			std::cerr << "Short: " << numAREdges_short << '\n';
			std::cerr << "Low quality: " << numAREdges_lowqual << '\n';
			std::cerr << "Zero extension: " << numAREdges_zero << '\n';

			std::cerr << "\nRead node - read node edges:" << '\n';
			std::cerr << "Usable: " << numRREdges_usable << '\n';
			std::cerr << "Contained: " << numRREdges_contained << '\n';
			std::cerr << "Short: " << numRREdges_short << '\n';
			std::cerr << "Low quality: " << numRREdges_lowqual << '\n';
			std::cerr << "Zero extension: " << numRREdges_zero << '\n';
		}
	}

  void SBridger::print(void) {
      std::cerr << "\nSBridger:\n";
      printData();
      printGraph();
  }

  void SBridger::generateGraph(void) {

  	 numANodes = numRNodes = 0;

	 numAREdges_all = numAREdges_usable = numAREdges_contained = numAREdges_short = numAREdges_lowqual = numAREdges_zero = 0;
	 numRREdges_all = numRREdges_usable = numRREdges_contained = numRREdges_short = numRREdges_lowqual = numRREdges_zero = 0;

	// 1. Generate anchor nodes
	for (auto const& it : mIdToContig) {
		auto node_ptr = make_shared<Node>(it.second, NT_ANCHOR);
		mAnchorNodes.emplace(it.first, node_ptr);
	}
	numANodes = mAnchorNodes.size();

	// 2. Generate read nodes
	for (auto const& it : mIdToRead) {
		auto node_ptr = make_shared<Node>(it.second, NT_READ);
		mReadNodes.emplace(it.first, node_ptr);
	}
	numRNodes = mReadNodes.size();

	// 3. Generate edges, function Overlap::Test() is used for filtering
	for (auto const& it1 : mIdToOvlR2C) {
		for (auto const& it2 : it1.second) {
			if (it2->Test()) {
				auto edge_ptr = make_shared<Edge>(it2, mAnchorNodes, mReadNodes);
				int test_val = edge_ptr->test();
				switch (test_val) {
					case (-1):
						numAREdges_contained += 1;
						break;
					case (-2):
						numAREdges_short += 1;
						break;
					case (-3):
						numAREdges_lowqual += 1;
						break;
					case (-4):
						numAREdges_zero += 1;
						break;
					default:
						break;
				}

				if (test_val > 0) {
					numAREdges_usable += 1;
					vEdges.emplace_back(edge_ptr);

					// Creating the second Edge in the opposite direction
					auto edge_ptr2 = make_shared<Edge>(it2, mAnchorNodes, mReadNodes);
					edge_ptr2->reverseNodes();
					vEdges.emplace_back(edge_ptr2);
				}
			}
		}
	}
	
	for (auto const& it1 : mIdToOvlR2R) {
		for (auto const& it2 : it1.second) {
			if (it2->Test()) {
				auto edge_ptr = make_shared<Edge>(it2, mAnchorNodes, mReadNodes);
				int test_val = edge_ptr->test();
				switch (test_val) {
					case (-1):
						numRREdges_contained += 1;
						break;
					case (-2):
						numRREdges_short += 1;
						break;
					case (-3):
						numRREdges_lowqual += 1;
						break;
					case (-4):
						numRREdges_zero += 1;
						break;
					default:
						break;
				}
				if (test_val > 0) {
					numRREdges_usable += 1;
					vEdges.emplace_back(edge_ptr);

					// Creating the second Edge in the opposite direction
					auto edge_ptr2 = make_shared<Edge>(it2, mAnchorNodes, mReadNodes);
					edge_ptr2->reverseNodes();
					vEdges.emplace_back(edge_ptr2);
				}
			}
		}
	}

	// 3.1 Connect Edges to Nodes
	for (auto const& edge : vEdges) {
		std::shared_ptr<Node> startNode = edge->startNode;
		std::shared_ptr<Node> endNode = edge->endNode;
		startNode->vOutEdges.emplace_back(edge);
		// KK: Creating two Edges per overlap, in both directions
		// The fist one already exists, creating the second one now
		// endNode->vOutEdges.emplace_back(edge);
	}

	// 4. Filter nodes
	// Remove isolated and contained read nodes
	// TODO:
	// Currently only calculating isolated Anchor and Read Nodes
	for (auto const& itANode : mAnchorNodes) {
		std::string aNodeName = itANode.first;
		std::shared_ptr<Node> aNode = itANode.second;
		if (aNode->vOutEdges.size() == 0) isolatedANodes += 1;
	}
	for (auto const& itRNode : mReadNodes) {
		std::string rNodeName = itRNode.first;
		std::shared_ptr<Node> rNode = itRNode.second;
		if (rNode->vOutEdges.size() == 0) isolatedRNodes += 1;
	}

	bGraphCreated = 1;
  }


  void SBridger::cleanupGraph(void) {
  	// TODO: this is currently a placeholder
  }

  int SBridger::generatePaths(void) {
  	int numPaths_maxOvl = scara::generatePathsDeterministic(vPaths, mAnchorNodes, PGT_MAXOS);
  	if (scara::print_output)
  		std::cerr << "\nSCARA: Generating paths using maximum overlap score. Number of paths generated: " << numPaths_maxOvl;
    int numPaths_maxExt = scara::generatePathsDeterministic(vPaths, mAnchorNodes, PGT_MAXES);
    if (scara::print_output)
    	std::cerr << "\nSCARA: Generating paths using maximum extension score. Number of paths generated: " << numPaths_maxExt;
    int minMCPaths = numPaths_maxExt + numPaths_maxOvl;
    if (minMCPaths < scara::MinMCPaths) minMCPaths = scara::MinMCPaths;
    int numPaths_MC = scara::generatePaths_MC(vPaths, mAnchorNodes, minMCPaths);
    if (scara::print_output)
  		std::cerr << "\nSCARA: Generating paths using Monte Carlo approach. Number of paths generated: " << numPaths_MC;


    return vPaths.size();
  }

  int SBridger::groupAndProcessPaths(void) {
  	int numGroups = 0;
  	
  	/* KK: Was only for testing */
  	// Printing paths before processing
  	if (scara::print_output) {
  		std::cerr << "\n\nSCARA: paths before processing:";
  		for (auto const& path_ptr : vPaths) {
  			shared_ptr<PathInfo> pathinfo_ptr = make_shared<PathInfo>(path_ptr);
  			std::cerr << "\nPATHINFO: SNODE(" << pathinfo_ptr->startNodeName << "), ";
  			std::cerr << "ENODE(" << pathinfo_ptr->endNodeName << "), ";
 			std::cerr << "DIRECTION(" << Direction2String(pathinfo_ptr->pathDir) << "), ";
  			std::cerr << "NODES(" << pathinfo_ptr->numNodes << "), ";
  			std::cerr << "BASES(" << pathinfo_ptr->length << "), ";
  			std::cerr << "BASES2(" << pathinfo_ptr->length2 << "), ";
  			std::cerr << "AVG SI(" << pathinfo_ptr->avgSI << "), ";
  			std::cerr << "CONSISTENT(" << checkPath(pathinfo_ptr->path_ptr) << ")";
  		}
  	}
  	/**/

  	// Path extending to the left are reversed so that all paths extend to the right
  	// Simulaneously path are grouped into buckets of set size
  	if (scara::print_output) {
  		std::cerr << "\n\nSCARA: paths after processing:";
  	}
  	for (auto const& path_ptr : vPaths) {
  		shared_ptr<Edge> firstEdge = path_ptr->edges[0];
  		Direction dir = D_LEFT;
        if (firstEdge->QES2 > firstEdge->QES1) dir = D_RIGHT;
        shared_ptr<PathInfo> pathinfo_ptr;
        if (dir == D_RIGHT) pathinfo_ptr = make_shared<PathInfo>(path_ptr);
        else {
        	shared_ptr<Path> revPath = path_ptr->reversedPath();
        	pathinfo_ptr = make_shared<PathInfo>(revPath);
        }

  		vPathInfos.emplace_back(pathinfo_ptr);
  		if (scara::print_output) {
  			std::cerr << "\nPATHINFO: SNODE(" << pathinfo_ptr->startNodeName << "), ";
  			std::cerr << "ENODE(" << pathinfo_ptr->endNodeName << "), ";
 			std::cerr << "DIRECTION(" << Direction2String(pathinfo_ptr->pathDir) << "), ";
  			std::cerr << "NODES(" << pathinfo_ptr->numNodes << "), ";
  			std::cerr << "BASES(" << pathinfo_ptr->length << "), ";
  			std::cerr << "BASES2(" << pathinfo_ptr->length2 << "), ";
  			std::cerr << "AVG SI(" << pathinfo_ptr->avgSI << "), ";
  			std::cerr << "CONSISTENT(" << checkPath(pathinfo_ptr->path_ptr) << ")";
  		}

  		// Grouping the path
  		bool grouped = false;
  		for (auto const& pgroup_ptr : vPathGroups) {
  			if (pgroup_ptr->addPathInfo(pathinfo_ptr)) {
  				grouped = true;
  				break;
  			}
  		}

  		if (!grouped) {
  			numGroups += 1;
  			shared_ptr<PathGroup> pgroup_ptr = make_shared<PathGroup>(pathinfo_ptr);
  			vPathGroups.emplace_back(pgroup_ptr);
  		}
  	}

  	if (scara::print_output) {
  		std::cerr << "\n\nSCARA: Groups before processing:";
	  	for (auto const& pgroup_ptr : vPathGroups) {
	  		std::cerr << "\nPATHGROUP: SNODE(" << pgroup_ptr->startNodeName << "), ";
  			std::cerr << "ENODE(" << pgroup_ptr->endNodeName << "), ";
  			std::cerr << "BASES(" << pgroup_ptr->length << "), ";
  			std::cerr << "PATHS(" << pgroup_ptr->numPaths << "), ";
	  	}
	}

	// For each node that acts as a starting node preserve only the best group
	// Currently this is a group with the largest number of paths
	// 1. Construct a Map with StartNode name as key and a vector of corresponding groups as value
	std::map<std::string, shared_ptr<vector<shared_ptr<PathGroup>>>> mGroups;
	for (auto const& pgroup_ptr : vPathGroups) {
		if (mGroups.find(pgroup_ptr->startNodeName) == mGroups.end()) {
			shared_ptr<vector<shared_ptr<PathGroup>>> val = make_shared<vector<shared_ptr<PathGroup>>>();
			val->emplace_back(pgroup_ptr);
			mGroups[pgroup_ptr->startNodeName] = val;
		}
		else {
			shared_ptr<vector<shared_ptr<PathGroup>>> val = mGroups[pgroup_ptr->startNodeName];
			val->emplace_back(pgroup_ptr);
		}
	}

	// 2. For each startNode in the map use only the best PathGroup
	std::map<std::string, shared_ptr<PathGroup>> vFilteredGroups;
	for (auto const& it : mGroups) {
		shared_ptr<vector<shared_ptr<PathGroup>>> val = it.second;
		shared_ptr<PathGroup> bestGroup = val->front();
		double bestSize = bestGroup->numPaths;
		for (uint32_t i=1; i<val->size(); i++) {
			shared_ptr<PathGroup> pgroup_ptr = (*val)[i];
			if (pgroup_ptr->numPaths > bestSize) {
				bestGroup = pgroup_ptr;
				bestSize = pgroup_ptr->numPaths;
			}
		}
		vFilteredGroups.emplace(bestGroup->startNodeName, bestGroup);
	}

	// Join groups that contain the same anchoring node, i.e groups 1->2 and 2->3, should be joined
	// to connect all three anchoring nodes in a single scaffold
	// First find all nodes that are starting nodes for a path, but are not ending nodes of any path
	// Those nodes will start a scaffold
	std::set<std::string> startNodes;
	std::set<std::string> endNodes;
	std::vector<shared_ptr<std::vector<shared_ptr<PathGroup>>>> scaffolds_temp;
	for (auto const& it : vFilteredGroups) {
		startNodes.insert(it.second->startNodeName);
		endNodes.insert(it.second->endNodeName);
	}

	// Only look at start nodes which are not in the end nodes set!
	for (auto const& startNode : startNodes) {
		if (endNodes.find(startNode) == endNodes.end()) {
			// Create a scaffold (a series of connected paths)
			auto newVec = make_shared<std::vector<shared_ptr<PathGroup>>>();
			// We can assume that the path with the considered startNode exists!
			std::string newStartNode = startNode;
			do {
				auto pgroup_ptr = vFilteredGroups[newStartNode];
				newVec->emplace_back(pgroup_ptr);
				newStartNode = pgroup_ptr->endNodeName;
			} while (vFilteredGroups.find(newStartNode) != vFilteredGroups.end());
			scaffolds_temp.emplace_back(newVec);
		}
	}

	if (scara::print_output) {
  		std::cerr << "\n\nSCARA: Final scaffolds before sequence generation:";
  		int i= 0;
	  	for (auto const& vec_ptr : scaffolds_temp) {
	  		i++;
	  		std::cerr << "\nSCAFFOLD " << i << ": ";
	  		for (auto const& pgroup_ptr : (*vec_ptr)) {
	  			std::cerr << pgroup_ptr->startNodeName << " --> " << pgroup_ptr->endNodeName << "(" << pgroup_ptr->numPaths << "), ";
	  		}
	  	}
	}

	// Processing scaffolds by chosing a best path within PathGroup
	// scaffolds_temp contains vectors of PathGroups
	// scaffolds contains vectors of PathInfos, pointing to the best path for each group
	for (auto const&  vec_ptr: scaffolds_temp) {
		auto newVec = make_shared<std::vector<shared_ptr<PathInfo>>>();
		for (auto const& pgroup_ptr : (*vec_ptr)) {
			auto best_pinfo_ptr = pgroup_ptr->vPathInfos[0];
	  		auto best_avgSI = best_pinfo_ptr->avgSI;
	  		for (auto const& pinfo_ptr : pgroup_ptr->vPathInfos) {		// KK: looking at the first element again, lazy to write it better
				if (pinfo_ptr->avgSI > best_avgSI) {
					best_avgSI = pinfo_ptr->avgSI;
					best_pinfo_ptr = pinfo_ptr;
				}
			}
			newVec->emplace_back(best_pinfo_ptr);
		}
		scaffolds.emplace_back(newVec);
	}

  	return scaffolds.size();
  }


  /* Each scaffold is represented as a series of paths, following one another, with the end node of a previous path
   * being the start node of the next one. Start and end nodes represent contigs, and they should be used completely
   * for generating a scaffold sequence, with the holes filled up using reads.
   * 1. Determine scaffold size and allocate enough space for a compelte sequence.
   * 2. Itterate over paths and generate sequence for each path
   *	Take care to use each contig only once (since they are present twice, except for the first and last one)
   */
  int SBridger::generateSequences(void) {
  	int i = 0;
  	for (auto const&  vec_ptr: scaffolds) {
  		// bool fistContigUsed = false;
  		i++;
  		// Generate header and calculate scaffold length
  		std::string header = ">Scaffold_" + to_string(i);
  		cerr << "SCARA: Generating sequence and header for scaffold " << i << endl;
  		cerr << "SCARA: sacffold edges: ";

  		uint32_t slength = 0;
		uint32_t lastNodeLength;
		Path scaffoldPath;
  		for (auto const& pinfo_ptr : (*vec_ptr)) {
  			header += ' ' + pinfo_ptr->startNodeName;
  			slength += pinfo_ptr->length;
  			// Remove the length of the endNode (as not to be added twice)
  			auto lastEdge = pinfo_ptr->path_ptr->edges.back();
  			if (!lastEdge->reversed) lastNodeLength = lastEdge->ovl_ptr->ext_ulTLen;
  			else lastNodeLength = lastEdge->ovl_ptr->ext_ulQLen;
  			slength -= lastNodeLength;

  			// Construct a joint path from the edges of all paths in the scaffold
  			cerr << pinfo_ptr->path_ptr->edges.size() << " (" << pinfo_ptr->length << ") ";
  			for (auto const& edge_ptr: pinfo_ptr->path_ptr->edges) scaffoldPath.appendEdge(edge_ptr);  
  		}
  		header += ' ' + vec_ptr->back()->endNodeName;		// Add the last endNode
  		slength += lastNodeLength;		// For the last path, add the endNode length
  		// uint32_t adjusted_slength = (uint32_t)(slength*1.05);		// Add 5%, just ot be sure we didnt miscalculate something

  		// Output the scaffold header to cout!
  		cerr << "SCARA generated header " << header << endl;
  		cout << header << endl;

  		cerr << "SCARA generating sequence of length " << slength << " from " << scaffoldPath.edges.size() << " nodes!" << endl;

  		// Calculate scaffold sequence from scaffoldPath and output it to cout
  		// NOTE: Assuming direction RIGHT!
  		// TODO: for completion include generating sequences for direction LEFT
  		SequenceStrand strand = SS_FORWARD;
  		for (auto const& edge_ptr : scaffoldPath.edges) {
  			// Determine part of the startNode that will be put into the final sequence
  			shared_ptr<Node> startNode = edge_ptr->startNode;
  			uint32_t seq_part_start, seq_part_end, seq_part_size;
  			std::string seq_part = "";
  			if (!edge_ptr->reversed) {		// Edge is not reversed
  				seq_part_start = 0;
  				seq_part_end = edge_ptr->ovl_ptr->ext_ulQBegin - edge_ptr->ovl_ptr->ext_ulTBegin;  				
  			}
  			else {							// Edge is reversed
  				seq_part_start = 0;
  				seq_part_end = edge_ptr->ovl_ptr->ext_ulTBegin - edge_ptr->ovl_ptr->ext_ulQBegin;
  			}
  			seq_part_size = seq_part_end - seq_part_start;
  			if (seq_part_size <= 0) {
  				throw std::runtime_error(std::string("SCARA BRIDGER: ERROR - invalid sequence part size: "));
  			}
  			seq_part.reserve(seq_part_size + 10);		// adding 10 just to avoid missing something by 1
  			if (strand == SS_FORWARD) {
  				// Copy relevant part of the string
  				int k=0;
  				for (; k<seq_part_size; k++) {
  					seq_part += (startNode->seq_ptr->seq_strData)[seq_part_start+k];
  				}
  				// seq_part[k] = '\0';
  			} else if (strand == SS_REVERSE) {
  				// Copy relevant part of the string
  				// If the strand is reverse, go from the end of the string and rev
  				int k=0;
  				for (; k<seq_part_size; k++) {
  					uint32_t seq_end = (startNode->seq_ptr->seq_strData).length();
  					seq_part += _bioBaseComplement((startNode->seq_ptr->seq_strData)[seq_end-k-1]);
  				}
  				// seq_part[k] = '\0';
  			} else {
  				throw std::runtime_error(std::string("SCARA BRIDGER: ERROR - invalid strand type: "));
  			}
  			// Output the sequence part to the standard output
  			cerr << "SCARA BRIDGER: Printing node " << startNode->nName << " with length " << seq_part_size << " - ";
  			cerr << seq_part.length() << "/" << startNode->seq_ptr->seq_strData.length();
  			cerr << endl;
  			cout << seq_part;

  			// Determine the strand for the next node
  			// Switch it if the overlap strand is '-' (or false)
  			if (!edge_ptr->ovl_ptr->ext_bOrientation) strand = reverseStrand(strand);
  		}

  		// At the end, add the complete last endNode
  		auto lastEndNode = scaffoldPath.edges.back()->endNode;
  		cerr << "SCARA BRIDGER: Printing node " << lastEndNode->nName << " with length ";
  		cerr << lastEndNode->seq_ptr->seq_strData.length() << "/" << lastEndNode->seq_ptr->seq_strData.length();
  		cerr << endl;
  		cout << lastEndNode->seq_ptr->seq_strData;
  		cout << endl;

  	}
  	return scaffolds.size();
  }


  // OBSOLETE
  void SBridger::Execute(void) {
  }

}
