#include <stdexcept>
#include <string>
#include <math.h>

#include "Graph.h"
#include "Overlap.h"
#include "Sequence.h"
#include "globals.h"

namespace scara {

  using namespace std;
  extern int multithreading;

  std::string NodeType2String(NodeType nType) {
  	switch (nType) {
  		case (NT_INVALID):
  			return "NT_INVALID";
  			break;
  		case (NT_ANCHOR):
  			return "NT_ANCHOR";
  			break;
  		case (NT_READ):
  			return "NT_READ";
  			break;
  		default:
  			return "UNKNOWN";
  	}
  }


  std::string Direction2String(Direction dir) {
  	switch (dir) {
  		case (D_LEFT):
  			return "LEFT";
  			break;
  		case (D_RIGHT):
  			return "RIGHT";
  			break;
  		default:
  			return "UNKNOWN";
  	}
  }

  Direction oppositeDirection(Direction dir) {
    switch (dir) {
      case (D_LEFT):
        return D_RIGHT;
        break;
      case (D_RIGHT):
        return D_LEFT;
        break;
      default:
        return D_LEFT;    // KK: Should probably be changed to raise some sort of exception
    }
  }


  Node::Node(NodeType i_nType, std::string i_nName, std::shared_ptr<Sequence> i_seq_ptr
  ) : nType(i_nType), nName(i_nName), seq_ptr(i_seq_ptr), isReverseComplement(false)
  {
  }

  Node::Node(NodeType i_nType, std::string i_nName, std::shared_ptr<Sequence> i_seq_ptr, bool i_isRC
  ) : nType(i_nType), nName(i_nName), seq_ptr(i_seq_ptr), isReverseComplement(i_isRC)
  {
  	if (isReverseComplement) {
  		nName += "_RC";
  	}
  }

  Node::Node(std::shared_ptr<Sequence> i_seq_ptr, NodeType i_nType
  ) : nType(i_nType), nName(i_seq_ptr->seq_strName), seq_ptr(i_seq_ptr), isReverseComplement(false)
  {
  }

  Node::Node(std::shared_ptr<Sequence> i_seq_ptr, NodeType i_nType, bool i_isRC
  ) : nType(i_nType), nName(i_seq_ptr->seq_strName), seq_ptr(i_seq_ptr), isReverseComplement(i_isRC)
  {
  	if (isReverseComplement) {
  		nName += "_RC";
  	}
  }

  Node::Node(NodeType i_nType, std::string i_nName, std::shared_ptr<Sequence> i_seq_ptr, std::vector<std::shared_ptr<Edge>> &i_vOutEdges
  ) : nType(i_nType), nName(i_nName), seq_ptr(i_seq_ptr), isReverseComplement(false)
  {
  	for (auto const& edge_ptr : i_vOutEdges) {
  		vOutEdges.emplace_back(std::move(edge_ptr));
  	}
  }

  Node::Node(NodeType i_nType, std::string i_nName, std::shared_ptr<Sequence> i_seq_ptr, std::vector<std::shared_ptr<Edge>> &i_vOutEdges, bool i_isRC
  ) : nType(i_nType), nName(i_nName), seq_ptr(i_seq_ptr), isReverseComplement(i_isRC)
  {
  	for (auto const& edge_ptr : i_vOutEdges) {
  		vOutEdges.emplace_back(std::move(edge_ptr));
  	}

  	if (isReverseComplement) {
  		nName += "_RC";
  	}
  }

  Edge::Edge(std::string i_startNodeName, std::shared_ptr<Node> i_startNode
  	       , std::string i_endNodeName, std::shared_ptr<Node> i_endNode, std::shared_ptr<Overlap> i_ovl_ptr
  ) : startNodeName(i_startNodeName), startNode(i_startNode), endNodeName(i_endNodeName), endNode(i_endNode), ovl_ptr(i_ovl_ptr)
  {
  	reversed = false;
  	this->calcEdgeStats();
  }

  Edge::Edge(std::shared_ptr<Overlap> i_ovl_ptr
  ) : startNodeName(i_ovl_ptr->ext_strName), endNodeName(i_ovl_ptr->ext_strTarget), ovl_ptr(i_ovl_ptr)
  {
  	reversed = false;
  	this->calcEdgeStats();
  }


  Edge::Edge(std::shared_ptr<Overlap> i_ovl_ptr, MapIdToNode &mAnchorNodes, MapIdToNode &mReadNodes
  ) : startNodeName(i_ovl_ptr->ext_strName), endNodeName(i_ovl_ptr->ext_strTarget), ovl_ptr(i_ovl_ptr)
  {
  	reversed = false;
  	this->calcEdgeStats();
  	// KK: The assumption is that all of the nodes are already loaded
  	// We are adding pointers to start and end node to each edge
  	MapIdToNode::iterator it;

  	// Look for START node in Anchor nodes
  	it = mAnchorNodes.find(startNodeName);
  	if (it != mAnchorNodes.end()) {
  		startNode = it->second;
  	}
  	else {		// Now look in Read nodes
  		it = mReadNodes.find(startNodeName);
  		if (it != mReadNodes.end()) {
  			startNode = it->second;
  		}
  		else {
  			throw std::runtime_error(std::string("Error loading graph edges. Unknown node: ") + startNodeName);
  		}
  	}

  	// Look for END node in Anchor nodes
  	it = mAnchorNodes.find(endNodeName);
  	if (it != mAnchorNodes.end()) {
  		endNode = it->second;
  	}
  	else {		// Now look in Read nodes
  		it = mReadNodes.find(endNodeName);
  		if (it != mReadNodes.end()) {
  			endNode = it->second;
  		}
  		else {
  			throw std::runtime_error(std::string("Error loading graph edges. Unknown node: ") + endNodeName);
  		}
  	}

  	// KK: the Edge will be later added to outgoing edges for both Nodes
  }

  void Edge::calcEdgeStats() {
  	if (ovl_ptr == NULL) {
  		throw std::runtime_error(std::string("Error calculating edge statistics"));
  	}

	QOH1 = ovl_ptr->ext_ulQBegin;
	QOH2 = ovl_ptr->ext_ulQLen - ovl_ptr->ext_ulQEnd;
	TOH1 = ovl_ptr->ext_ulTBegin;
	TOH2 = ovl_ptr->ext_ulTLen - ovl_ptr->ext_ulTEnd;

	QOL = ovl_ptr->ext_ulQEnd - ovl_ptr->ext_ulQBegin;
	TOL = ovl_ptr->ext_ulTEnd - ovl_ptr->ext_ulTBegin;

	SI = (float)(ovl_ptr->paf_matching_bases)/ovl_ptr->paf_overlap_length;

	float avg_ovl_len = (QOL+TOL)/2;
	OS = avg_ovl_len*SI;
	QES1 = OS + TOH1/2 - (QOH1 + TOH2)/2;
	QES2 = OS + TOH2/2 - (QOH2 + TOH1)/2;
	TES1 = OS + QOH1/2 - (QOH2 + TOH1)/2;
	TES2 = OS + QOH2/2 - (QOH1 + TOH2)/2;

	// NOTE: This seeme logical:
	// If a query extends further right or left then the target, it makes no sense to extend it in that direction
	// Therefore setting a corresponding extension score to 0
    if (QOH1 >= TOH1) {
        QES1 = 0;
    } else {
        TES1 = 0;
    }
    if (QOH2 >= TOH2) {
        QES2 = 0;
    } else {
        TES2 = 0;
    }
  }


  // Reverses the nodes in an edge
  void Edge::reverseNodes(void) {
  	std::swap(startNodeName, endNodeName);
  	std::shared_ptr<Node> tempNodePtr = startNode;
  	startNode = endNode;
  	endNode = tempNodePtr;

  	uint32_t tui32 = QOH1;
  	QOH1 = TOH1;
  	TOH1 = tui32;

  	tui32 = QOH2;
  	QOH2 = TOH2;
  	TOH2 = tui32;

  	tui32 = QOL;
  	QOL = TOL;
  	TOL = tui32;

  	// Sequence identity and overlap score remain unchanged

  	float fTemp = QES1;
  	QES1 = TES1;
  	TES1 = fTemp;

  	fTemp = QES2;
  	QES2 = TES2;
  	TES2 = fTemp;

  	reversed = !reversed;
  }



  // A function that test if an overlap (PAF line) is usable or not
  // Overall, an overlap is not usable if:
  // - one read contains the other - returns -1 
  // - length of aligned part is to short compared to overhangs - returns -2
  // - mapping quality is too low - returns -3
  // If the read is usable, the function returns 1
  int Edge::test(void) {
  	float minQOH = (QOH1 < QOH2) ? QOH1 : QOH2;          // Smaller query overhang, will be used to determine if the overlap is discarded
    float minTOH = (TOH1 < TOH2) ? TOH1 : TOH2;          // Smaller target overhang, will be used to determine if the overlap is discarded

    float minOH1 = (QOH1 < TOH1) ? QOH1 : TOH1;          // Smaller left overhang
    float minOH2 = (QOH2 < TOH2) ? QOH2 : TOH2;          // Smaller right overhang

    // Test for too short aligned length
    // In this case the overlap is discarded, but both reads are kept
    // if test_short_length:
    //     if  float(minQOH + minTOH)/avg_ovl_len > OHmax:
    //         return -2

    // New test for short overlaps (large overhangs)
    float avg_ovl_len = (QOL+TOL)/2;
    if (test_short_length) {
        if ((minOH1 + minOH2)/avg_ovl_len > scara::OHmax) {
            return -2;
        }
    }

    // Test for contained reads
    // Has to come after test for short aligned length, if the overlap is of too short a length
    // Its probably a false overlap
    if (test_contained_reads){
        if (QOH1 >= TOH1 && QOH2 >= TOH2) {
            // Target is contained within the query
            // Discarding the overlap and target read
            std::string tname = ovl_ptr->ext_strTarget;
            // reads_to_discard[tname] = 1
            return -1;
        }
        if (TOH1 >= QOH1 && TOH2 >= QOH2) {
            // Query is contained within the target
            // Discarding the overlap and query read
            std::string qname = ovl_ptr->ext_strName;
            // reads_to_discard[qname] = 1
            return -1;
        }
    }

    // Test for low quality overlap
    if (test_low_quality) {
        if (SI < scara::SImin) {
            return -3;
        }
    }

    //If there are some overlaps with zero extension score on both ends, discard those as well
    if (QES1 <= 0 && QES2 <= 0 && TES1 <= 0 && TES2 <= 0) {
        return -4;
    }

    // If the overlap is correct, return 1

  	return true;
  }


  // KK: Since the architecture of the graph has been changed, we need to create 4 edges for each overlap
  // We have 2 nodes for each sequence (FW and RC), and need to create edges in both directions but only 
  // For appropriate node combinations
  // Example:
  // 1. if we have overlap between seq1 and seq2, on the same strand (relative strand from PAF = '+')
  // Then we need to create edges for nodes Seq1 and Seq2 (both directions), and for nodes Seq1_RC and Seq2_RC (both directions)
  // 2. if we have overlap between seq1 and seq2, on different strand (relative strand from PAF = '-')
  // Then we need to create edges for nodes Seq1 and Seq2_RC (both directions), and for nodes Seq1_RC and Seq2 (both directions)
  void createEdgesFromOverlap(std::shared_ptr<Overlap> ovl_ptr, MapIdToNode& mAnchorNodes, MapIdToNode& mReadNodes, std::vector<shared_ptr<Edge>> &vEdges) {

	// startNodeName(i_ovl_ptr->ext_strName), endNodeName(i_ovl_ptr->ext_strTarget), ovl_ptr(i_ovl_ptr)

  	std::string startNodeName = ovl_ptr->ext_strName;
  	std::string startNodeName_RC = ovl_ptr->ext_strName + "_RC";
  	std::string endNodeName = ovl_ptr->ext_strTarget;
  	std::string endNodeName_RC = ovl_ptr->ext_strTarget + "_RC";
  	std::shared_ptr<Node> startNode, startNode_RC, endNode, endNode_RC;

  	// KK: The assumption is that all of the nodes are already loaded
  	// We are adding pointers to start and end node to each edge
  	MapIdToNode::iterator it;

  	// Look for START node in Anchor nodes
  	it = mAnchorNodes.find(startNodeName);
  	if (it != mAnchorNodes.end()) {
  		startNode = it->second;
  	}
  	else {		// Now look in Read nodes
  		it = mReadNodes.find(startNodeName);
  		if (it != mReadNodes.end()) {
  			startNode = it->second;
  		}
  		else {
  			throw std::runtime_error(std::string("Error loading graph edges. Unknown node: ") + startNodeName);
  		}
  	}

  	// Look for START node REVERSE COMPLEMENT in Anchor nodes
  	it = mAnchorNodes.find(startNodeName_RC);
  	if (it != mAnchorNodes.end()) {
  		startNode_RC = it->second;
  	}
  	else {		// Now look in Read nodes
  		it = mReadNodes.find(startNodeName_RC);
  		if (it != mReadNodes.end()) {
  			startNode_RC = it->second;
  		}
  		else {
  			throw std::runtime_error(std::string("Error loading graph edges. Unknown node: ") + startNodeName_RC);
  		}
  	}

  	// Look for END node in Anchor nodes
  	it = mAnchorNodes.find(endNodeName);
  	if (it != mAnchorNodes.end()) {
  		endNode = it->second;
  	}
  	else {		// Now look in Read nodes
  		it = mReadNodes.find(endNodeName);
  		if (it != mReadNodes.end()) {
  			endNode = it->second;
  		}
  		else {
  			throw std::runtime_error(std::string("Error loading graph edges. Unknown node: ") + endNodeName);
  		}
  	}

  	// Look for END node REVERSE COMPLEMENT in Anchor nodes
  	it = mAnchorNodes.find(endNodeName_RC);
  	if (it != mAnchorNodes.end()) {
  		endNode_RC = it->second;
  	}
  	else {		// Now look in Read nodes
  		it = mReadNodes.find(endNodeName_RC);
  		if (it != mReadNodes.end()) {
  			endNode_RC = it->second;
  		}
  		else {
  			throw std::runtime_error(std::string("Error loading graph edges. Unknown node: ") + endNodeName_RC);
  		}
  	}

  	// Checking relative strand
  	if (ovl_ptr->ext_bOrientation == true) {
  		// Same strand ('+')
  		// Create edge for FW strand
  		auto edge_ptr = make_shared<Edge>(startNodeName, startNode, endNodeName, endNode, ovl_ptr);
  		vEdges.emplace_back(edge_ptr);

		// Creating the second Edge in the opposite direction (still FW strand)
		auto edge_ptr2 = make_shared<Edge>(startNodeName, startNode, endNodeName, endNode, ovl_ptr);
		edge_ptr2->reverseNodes();
		vEdges.emplace_back(edge_ptr2);

		// Create edge for RC strand
		auto edge_ptr3 = make_shared<Edge>(startNodeName_RC, startNode_RC, endNodeName_RC, endNode_RC, ovl_ptr);
  		vEdges.emplace_back(edge_ptr3);

		// Creating the second Edge in the opposite direction (for RC strand)
		auto edge_ptr4 = make_shared<Edge>(startNodeName_RC, startNode_RC, endNodeName_RC, endNode_RC, ovl_ptr);
		edge_ptr4->reverseNodes();
		vEdges.emplace_back(edge_ptr4);
  	}
  	else {
  		// Different strand ('-')
  		// Create edge for FW strand
  		auto edge_ptr = make_shared<Edge>(startNodeName, startNode, endNodeName_RC, endNode_RC, ovl_ptr);
  		vEdges.emplace_back(edge_ptr);

		// Creating the second Edge in the opposite direction (still FW strand)
		auto edge_ptr2 = make_shared<Edge>(startNodeName, startNode, endNodeName_RC, endNode_RC, ovl_ptr);
		edge_ptr2->reverseNodes();
		vEdges.emplace_back(edge_ptr2);

		// Create edge for RC strand
		auto edge_ptr3 = make_shared<Edge>(startNodeName_RC, startNode_RC, endNodeName, endNode, ovl_ptr);
  		vEdges.emplace_back(edge_ptr3);

		// Creating the second Edge in the opposite direction (for RC strand)
		auto edge_ptr4 = make_shared<Edge>(startNodeName_RC, startNode_RC, endNodeName, endNode, ovl_ptr);
		edge_ptr4->reverseNodes();
		vEdges.emplace_back(edge_ptr4);
  	}

  }



  Path::Path(std::shared_ptr<Edge> edge_ptr) :  edges()
  {
    edges.emplace_back(edge_ptr);
  }

  Path::Path(void) : edges()
  {
  }

  // Appends an Edge to a Path
  // The path must be empty or the end node of the last Edge must be the same as the start node of the new Edge
  void Path::appendEdge(std::shared_ptr<Edge> edge_ptr) {
    if (edges.empty()) {
      edges.emplace_back(edge_ptr);
    }
    else {
      std::string lastEndNode = edges.back()->endNodeName;
      std::string newStartNode = edge_ptr->startNodeName;
      if (lastEndNode.compare(newStartNode) == 0) {
        edges.emplace_back(edge_ptr);
      }
      else {
        throw std::runtime_error(std::string("Unable to append an Edge to a Path. Incompatible nodes: ") + lastEndNode + " | " + newStartNode);
      }
    }
  }

  // Size of the path
  int Path::size(void) {
  	return edges.size();
  }

  std::shared_ptr<Node> Path::endNode(void) {
  	if (edges.size() == 0) {
  		return NULL;
  	}
  	else {
  		return edges.back()->endNode;
  	}
  }


  std::shared_ptr<Node> Path::startNode(void) {
  	if (edges.size() == 0) {
  		return NULL;
  	}
  	else {
  		return edges[0]->startNode;
  	}
  }

  std::shared_ptr<Edge> Path::removeLastEdge(void){
    auto last_edge = edges.back();
  	edges.pop_back();
    return last_edge;
  }

  // Return a reversed path, reverse order of edges and each edge
  shared_ptr<Path> Path::reversedPath() {
  	shared_ptr<Path> newPath = make_shared<Path>();

  	if (edges.size() > 0) {
  		for (int i = edges.size()-1; i >= 0; i--) {
  			shared_ptr<Edge> edge = edges[i];
  			shared_ptr<Edge> newEdge;
        if (edge->reversed) {
          newEdge = make_shared<Edge>(edge->endNodeName, edge->endNode                  // If the edge is already reversed, create nonreversed version
  														, edge->startNodeName, edge->startNode, edge->ovl_ptr);   // by switching start and end nodes (and do not reverse it)
        }
        else {
          newEdge = make_shared<Edge>(edge->startNodeName, edge->startNode              // if the edgee is not reversed, create a noraml version and reverse it
                              , edge->endNodeName, edge->endNode, edge->ovl_ptr);
 			    newEdge->reverseNodes();
        }
  			newPath->appendEdge(newEdge);
  		}
  	}

  	return newPath;
  }


  PathInfo::PathInfo(shared_ptr<Path> t_path_ptr) : path_ptr(t_path_ptr)
  {
  	length = avgSI = 0.0;
    length2 = 0;
    uint32_t negativeEScount = 0;
  	uint32_t SNlen, SNbegin, SNend;
  	uint32_t ENlen, ENbegin, ENend;

  	if (t_path_ptr->edges.size() > 0) {
  		numNodes = t_path_ptr->edges.size() + 1;
  		shared_ptr<Edge> firstEdge = t_path_ptr->edges.front();
  		startNodeName = firstEdge->startNodeName;
  		endNodeName = t_path_ptr->edges.back()->endNodeName;
  		pathDir = D_LEFT;
  		if (firstEdge->QES1 < firstEdge->QES2) pathDir = D_RIGHT;
  		// Add the length of the first Node (start node of the first edge)
  		if (!(firstEdge->reversed)) {
  			length += firstEdge->ovl_ptr->ext_ulQLen;
  		}
  		else {
  			length += firstEdge->ovl_ptr->ext_ulTLen;
  		}
	  	for (auto const& edge : t_path_ptr->edges) {
	  		avgSI += edge->SI;
	  		if (!(edge->reversed)) {
	  			SNlen   = edge->ovl_ptr->ext_ulQLen;
	  			SNbegin = edge->ovl_ptr->ext_ulQBegin;
	  			SNend   = edge->ovl_ptr->ext_ulQEnd;
	  			ENlen   = edge->ovl_ptr->ext_ulTLen;
	  			ENbegin = edge->ovl_ptr->ext_ulTBegin;
	  			ENend   = edge->ovl_ptr->ext_ulTEnd;
	  		}
	  		else {
	  			SNlen   = edge->ovl_ptr->ext_ulTLen;
	  			SNbegin = edge->ovl_ptr->ext_ulTBegin;
	  			SNend   = edge->ovl_ptr->ext_ulTEnd;
	  			ENlen   = edge->ovl_ptr->ext_ulQLen;
	  			ENbegin = edge->ovl_ptr->ext_ulQBegin;
	  			ENend   = edge->ovl_ptr->ext_ulQEnd;
	  		}
	  		// For each edge, add to the length of the path part of the end node that 
	  		// does not overlap with the start node (End overhang - start overhang)
	        if (pathDir == D_RIGHT) {
		  		  length += (ENlen - ENend) - (SNlen - SNend);
	          if (SNbegin <= ENbegin) {
				  negativeEScount++;
	          }
	          length2 += SNbegin - ENbegin;
	        }
	        else {
	          length += ENbegin - SNbegin;
	          if (ENbegin <= SNbegin) {
	            negativeEScount++;
	       	  }
	          length2 += (SNlen - SNend) - (ENlen - ENend);
	        }
	  	}
      length2 += ENlen;
	  	avgSI /= t_path_ptr->edges.size();
      if (negativeEScount > 0) 
        throw std::runtime_error(std::string("SCARA BRIDGER: ERROR - path with negative extensions (" + std::to_string(negativeEScount) + ") SN(" + startNodeName + ") EN(" + endNodeName + ")"));
  	}
  	else {
  		numNodes = 0;
  		startNodeName = "";
  		endNodeName = "";
  	}
  }

  PathGroup::PathGroup() : startNodeName(""), endNodeName(""), length(0.0), numPaths(0)
  {
  }

  PathGroup::PathGroup(std::string t_startNodeName, std::string t_endNodeName, double t_length
    ) : startNodeName(t_startNodeName), endNodeName(t_endNodeName), length(t_length), numPaths(0)
  {
  }
  

  PathGroup::PathGroup(shared_ptr<PathInfo> pinfo_ptr
    ) : startNodeName(pinfo_ptr->startNodeName), endNodeName(pinfo_ptr->endNodeName), length(pinfo_ptr->length), numPaths(1)
  {
    vPathInfos.emplace_back(pinfo_ptr);
  }

  bool PathGroup::addPathInfo(shared_ptr<PathInfo> pinfo_ptr) {
    if ((pinfo_ptr->startNodeName.compare(startNodeName) != 0) ||
        (pinfo_ptr->endNodeName.compare(endNodeName) != 0) ||
        (fabs(pinfo_ptr->length - length) > scara::pathGroupHalfSize)) return false;

    vPathInfos.emplace_back(pinfo_ptr);
    numPaths += 1;
    return true;

  }

}