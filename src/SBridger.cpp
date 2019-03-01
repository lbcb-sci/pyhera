#include "SBridger.h"
#include "Loader.h"
#include <vector>
// #include <string>
#include <iostream>

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
    int numPaths_MC = scara::generatePaths_MC(vPaths, mAnchorNodes, scara::MinMCPaths);
    if (scara::print_output)
  		std::cerr << "\nSCARA: Generating paths using Monte Carlo approach. Number of paths generated: " << numPaths_MC;


    return vPaths.size();
  }

  int SBridger::groupAndProcessPaths(void) {
  	int numGroups = 0;
  	
  	// Path extending to the right are reversed so that all paths extend to the left
  	for (auto const& path_ptr : vPaths) {
  		shared_ptr<Edge> firstEdge = path_ptr->edges[0];
  		Direction dir = D_LEFT;
        if (firstEdge->QES2 > firstEdge->QES1) dir = D_RIGHT;
        shared_ptr<PathInfo> pathinfo_ptr;
        if (dir == D_LEFT) pathinfo_ptr = make_shared<PathInfo>(path_ptr);
        else {
        	shared_ptr<Path> revPath = path_ptr->reversedPath();
        	pathinfo_ptr = make_shared<PathInfo>(revPath);
        }

  		vPathInfos.emplace_back(pathinfo_ptr);
  		if (scara::print_output) {
  			std::cerr << "\nPATHINFO: SNODE(" << pathinfo_ptr->startNodeName << "), ";
  			std::cerr << "ENODE(" << pathinfo_ptr->endNodeName << "), ";
  			std::cerr << "ENODE(" << pathinfo_ptr->endNodeName << "), ";
  			std::cerr << "DIRECTION(" << Direction2String(pathinfo_ptr->pathDir) << "), ";
  			std::cerr << "NODES(" << pathinfo_ptr->numNodes << "), ";
  			std::cerr << "BASES(" << pathinfo_ptr->length << "), ";
  			std::cerr << "AVG SI(" << pathinfo_ptr->avgSI << ")";
  		}
  	}


  	return numGroups;
  }

  int SBridger::generateSequences(void) {
  	// TODO: this is currently a placeholder
  	return 0;
  }

  void SBridger::Execute(void) {
  }

}
