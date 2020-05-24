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

  // Return a reverse complement node name, it the original name ends with _RC, remove _RC from the end
  // Otherwise, ad _RC at the end
  std::string getRCNodeName(const std::string nodeName) {
  	if (nodeName.rfind("_RC") != nodeName.size() - 3) return nodeName + "_RC";
  	else return nodeName.substr(0, nodeName.size()-3);
  }


  // Return an original node name, without _RC at the end
  std::string getOGNodeName(const std::string nodeName) {
  	if (nodeName.rfind("_RC") != nodeName.size() - 3) return nodeName;
  	else return nodeName.substr(0, nodeName.size()-3);
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

  Edge::Edge(std::shared_ptr<Node> i_startNode, std::shared_ptr<Node> i_endNode, std::unique_ptr<Overlap> const& i_ovl_ptr
  ) : startNode(i_startNode), endNode(i_endNode)
  {
  	this->copyDataFromOverlap(i_ovl_ptr);
  	this->calcEdgeStats();
  }

  Edge::Edge(std::unique_ptr<Overlap> const& i_ovl_ptr
  )
  {
  	this->copyDataFromOverlap(i_ovl_ptr);
  	this->calcEdgeStats();
  }


  Edge::Edge(std::unique_ptr<Overlap> const& i_ovl_ptr, MapIdToNode &mAnchorNodes, MapIdToNode &mReadNodes
  )
  {
  	this->copyDataFromOverlap(i_ovl_ptr);
  	this->calcEdgeStats();
  	// KK: The assumption is that all of the nodes are already loaded
  	// We are adding pointers to start and end node to each edge
  	MapIdToNode::iterator it;

  	// Look for START node in Anchor nodes
  	it = mAnchorNodes.find(i_ovl_ptr->ext_strName);
  	if (it != mAnchorNodes.end()) {
  		startNode = it->second;
  	}
  	else {		// Now look in Read nodes
  		it = mReadNodes.find(i_ovl_ptr->ext_strName);
  		if (it != mReadNodes.end()) {
  			startNode = it->second;
  		}
  		else {
  			throw std::runtime_error(std::string("Error loading graph edges. Unknown node: ") + getStartNodeName());
  		}
  	}

  	// Look for END node in Anchor nodes
  	it = mAnchorNodes.find(i_ovl_ptr->ext_strTarget);
  	if (it != mAnchorNodes.end()) {
  		endNode = it->second;
  	}
  	else {		// Now look in Read nodes
  		it = mReadNodes.find(i_ovl_ptr->ext_strTarget);
  		if (it != mReadNodes.end()) {
  			endNode = it->second;
  		}
  		else {
  			throw std::runtime_error(std::string("Error loading graph edges. Unknown node: ") + getEndNodeName());
  		}
  	}

  	// KK: the Edge will be later added to outgoing edges for both Nodes
  }

  // Create an empty Edge
  // KK: Do I need to initializa all the variables?
  Edge::Edge() {
  	startNode = endNode = NULL;
    
    SLen = SStart = SEnd = 0;
    ELen = EStart = EEnd = 0;
    ovl_bOrientation = true;

    paf_matching_bases = paf_overlap_length = paf_mapping_quality = 0;

    QOH1 = QOH2 = TOH1 = TOH2 = 0;
    QOL = TOL = 0;
	   SI = 0.0;

	  OS = QES1 = QES2 = TES1 = TES2 = 0.0;
  }

  std::string Edge::getStartNodeName() {
    return startNode->nName;
  }
  std::string Edge::getEndNodeName() {
    return endNode->nName;
  }

  void Edge::copyDataFromOverlap(std::unique_ptr<Overlap> const& ovl_ptr) {
  	if (ovl_ptr == NULL) {
  		throw std::runtime_error(std::string("Error calculating edge statistics"));
  	}
    // Copying data from overlap, Start node is Query, End node is Target
    SLen   = ovl_ptr->ext_ulQLen;
    SStart = ovl_ptr->ext_ulQBegin;
    SEnd   = ovl_ptr->ext_ulQEnd;
    ELen   = ovl_ptr->ext_ulTLen;
    EStart = ovl_ptr->ext_ulTBegin;
    EEnd   = ovl_ptr->ext_ulTEnd;
    ovl_bOrientation = ovl_ptr->ext_bOrientation;

    this->paf_matching_bases = ovl_ptr->paf_matching_bases;
    this->paf_overlap_length = ovl_ptr->paf_overlap_length;
    this->paf_mapping_quality = ovl_ptr->paf_mapping_quality;

  }

  void Edge::calcEdgeStats() {
  	QOH1 = SStart;
  	QOH2 = SLen - SEnd;

    if (ovl_bOrientation) {
  	  TOH1 = EStart;
  	  TOH2 = ELen - EEnd;
    } else {
      TOH2 = EStart;
      TOH1 = ELen - EEnd;
    }

  	QOL = SEnd - SStart;
  	TOL = EEnd - EStart;

  	SI = (float)(paf_matching_bases)/paf_overlap_length;

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

  std::shared_ptr<Edge> Edge::getReversedEdge() {
  	std::shared_ptr<Edge> edge_ptr2 = make_shared<Edge>();

  	edge_ptr2->startNode = this->endNode;
  	edge_ptr2->endNode = this->startNode;

  	edge_ptr2->SLen 	= this->ELen;
  	edge_ptr2->SStart 	= this->EStart;
  	edge_ptr2->SEnd 	= this->EEnd;
  	edge_ptr2->ELen 	= this->SLen;
  	edge_ptr2->EStart 	= this->SStart;
  	edge_ptr2->EEnd 	= this->SEnd;
  	edge_ptr2->ovl_bOrientation = this->ovl_bOrientation;

  	edge_ptr2->paf_matching_bases  = this->paf_matching_bases;
    edge_ptr2->paf_overlap_length  = this->paf_overlap_length;
    edge_ptr2->paf_mapping_quality = this->paf_mapping_quality;

    edge_ptr2->calcEdgeStats();

    return edge_ptr2;
  }


  // Reverses the nodes in an edge
  void Edge::reverseNodes(void) {
  	std::shared_ptr<Node> tempNodePtr = startNode;
  	startNode = endNode;
  	endNode = tempNodePtr;

  	uint32_t tui32 = SLen;
  	SLen = ELen;
  	ELen = tui32;

  	tui32 = SStart;
  	SStart = EStart;
  	EStart = tui32;

  	tui32 = SEnd;
  	SEnd = EEnd;
  	EEnd = tui32;

  	tui32 = QOH1;
  	QOH1 = TOH1;
  	TOH1 = tui32;

  	tui32 = QOH2;
  	QOH2 = TOH2;
  	TOH2 = tui32;

  	tui32 = QOL;
  	QOL = TOL;
  	TOL = tui32;

  	// Relative strand, sequence identity and overlap score remain unchanged

  	float fTemp = QES1;
  	QES1 = TES1;
  	TES1 = fTemp;

  	fTemp = QES2;
  	QES2 = TES2;
  	TES2 = fTemp;
  }

  // Change the edge data so that Query (Start node) and Target (End node) are on reverse strands
  // NOTE: relative strand remain unchanged
  void Edge::reverseStrand(void) {
  	uint32_t tui32 = this->SStart;
  	this->SStart = this->SLen - this->SEnd;
  	this->SEnd = this->SLen - tui32;

  	tui32 = this->EStart;
  	this->EStart = this->ELen - this->EEnd;
  	this->EEnd = this->ELen - tui32;

  	this->calcEdgeStats();

  }

  void Edge::reverseStrandSNode(void) {
  	this->ovl_bOrientation = !(this->ovl_bOrientation);

  	uint32_t tui32 = this->SStart;
  	this->SStart = this->SLen - this->SEnd;
  	this->SEnd = this->SLen - tui32;

  	this->calcEdgeStats();
  }


  void Edge::reverseStrandENode(void) {
  	this->ovl_bOrientation = !(this->ovl_bOrientation);

	uint32_t tui32 = this->EStart;
  	this->EStart = this->ELen - this->EEnd;
  	this->EEnd = this->ELen - tui32;

  	this->calcEdgeStats();
  }



  // A function that test if an overlap (PAF line) is usable or not
  // Overall, an overlap is not usable if:
  // - one read contains the other - returns -1 
  // - length of aligned part is to short compared to overhangs - returns -2
  // - mapping quality is too low - returns -3
  // If the read is usable, the function returns 1
  int Edge::test(void) {
  	// float minQOH = (QOH1 < QOH2) ? QOH1 : QOH2;          // Smaller query overhang, will be used to determine if the overlap is discarded
    // float minTOH = (TOH1 < TOH2) ? TOH1 : TOH2;          // Smaller target overhang, will be used to determine if the overlap is discarded

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
            // std::string tname = endNodeName;
            // reads_to_discard[tname] = 1
            return -1;
        }
        if (TOH1 >= QOH1 && TOH2 >= QOH2) {
            // Query is contained within the target
            // Discarding the overlap and query read
            // std::string qname = startNodeName;
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

  	return 1;
  }


  // KK: Since the architecture of the graph has been changed, we need to create 2 edges for each overlap
  // Since the paths are generated from all contigs (anchop nodes), they will be generated only in direction
  // RIGHT. If the nodes are connected, path should be generated between them. 
  // We have 2 nodes for each sequence (FW and RC), and need to create edges in directions RIGHT but only 
  // For appropriate node combinations
  // Example:
  // 1. if we have overlap between Seq1 and Seq2, on the same strand (relative strand from PAF = '+')
  // Then we need to create edge for nodes Seq1 and Seq2 , and for nodes Seq1_RC and Seq2_RC
  // 2. if we have overlap between Seq1 and Seq2, on different strands (relative strand from PAF = '-')
  // Then we need to create edges for nodes Seq1 and Seq2_RC, and for nodes Seq1_RC and Seq2
  void createEdgesFromOverlap(std::unique_ptr<Overlap> const& ovl_ptr, MapIdToNode& mAnchorNodes
  							, MapIdToNode& mReadNodes, std::vector<shared_ptr<Edge>> &vEdges) {

	// startNodeName(i_ovl_ptr->ext_strName), endNodeName(i_ovl_ptr->ext_strTarget), ovl_ptr(i_ovl_ptr)

  	std::string startNodeName = ovl_ptr->ext_strName;
  	std::string startNodeName_RC = getRCNodeName(ovl_ptr->ext_strName);
  	std::string endNodeName = ovl_ptr->ext_strTarget;
  	std::string endNodeName_RC = getRCNodeName(ovl_ptr->ext_strTarget);
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
  	// KK: This could be put in the previous code section, because original node and RC should be in the same node set
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

  	// If relative overlap strand is +, edges are ctreated for (SNode, Enode) and (SNodeRC, ENodeRC)
  	if (ovl_ptr->ext_bOrientation) {
		// Create edge for both nodes on original strands
		auto edge_ptr = make_shared<Edge>(startNode, endNode, ovl_ptr);
		if (edge_ptr->QOH2 > edge_ptr->TOH2) edge_ptr->reverseNodes();		// If right overhang for query is larger than for target
																			// Reverse nodes, we are extending start node only to the RIGHT
		vEdges.emplace_back(edge_ptr);

		// Creating the second Edge for both nodes on their respective RC strands
		auto edge_ptr2 = make_shared<Edge>(startNode_RC, endNode_RC, ovl_ptr);
		edge_ptr2->reverseStrand();
		if (edge_ptr2->QOH2 > edge_ptr2->TOH2) edge_ptr2->reverseNodes();		// It right overhang for query is larger than for target
																				// Reverse nodes, we are extending start node only to the RIGHT
		vEdges.emplace_back(edge_ptr2);
	}
	// If relative strand is -, edges are created for (SNode, ENodeRC) and (SNodeRC, ENode)
	else {
		// (SNode, ENodeRC)
		auto edge_ptr = make_shared<Edge>(startNode, endNode_RC, ovl_ptr);
		edge_ptr->reverseStrandENode();
		if (edge_ptr->QOH2 > edge_ptr->TOH2) edge_ptr->reverseNodes();		// If right overhang for query is larger than for target
																			// Reverse nodes, we are extending start node only to the RIGHT
		vEdges.emplace_back(edge_ptr);

		// (SNodeRC, ENode)
		auto edge_ptr2 = make_shared<Edge>(startNode_RC, endNode, ovl_ptr);
		edge_ptr2->reverseStrandSNode();
		if (edge_ptr2->QOH2 > edge_ptr2->TOH2) edge_ptr2->reverseNodes();		// It right overhang for query is larger than for target
																				// Reverse nodes, we are extending start node only to the RIGHT
		vEdges.emplace_back(edge_ptr2);
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
      std::string lastEndNode = edges.back()->getEndNodeName();
      std::string newStartNode = edge_ptr->getStartNodeName();
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
  shared_ptr<Path> Path::reversedPath(void) {
  	shared_ptr<Path> newPath = make_shared<Path>();

  	if (edges.size() > 0) {
  		for (int i = edges.size()-1; i >= 0; i--) {
  			shared_ptr<Edge> edge = edges[i];
  			shared_ptr<Edge> newEdge = edge->getReversedEdge();
  			newPath->appendEdge(newEdge);
  		}
  	}

  	return newPath;
  }

  std::string Path::toString(void){
  	size_t size = this->edges.size();
  	std::string strPath = "Nodes: " + std::to_string(size) + " | ";
  	if (size > 10) {
  		auto first = this->edges.front();
  		auto last  = this->edges.back();
  		strPath += "(" + first->getStartNodeName() + ", " + first->getEndNodeName() + ") ... ";
  		strPath += "(" + last->getStartNodeName() + ", " + last->getEndNodeName() + ")";
  	}
  	else {
  		for (auto const& t_edge_ptr : this->edges) {
  			strPath += "(" + t_edge_ptr->getStartNodeName() + ", " + t_edge_ptr->getEndNodeName() + ") ";
  		}	
  	}

	 return strPath; 
  }


  // Summary information about paths
  // Used to later group paths by start and end node
  PathInfo::PathInfo(shared_ptr<Path> t_path_ptr) : path_ptr(t_path_ptr)
  {
  	this->length = this->length2 = this->avgSI = 0.0;
  	this->path_ptr = t_path_ptr;
    uint32_t negativeEScount = 0;

  	if (t_path_ptr->edges.size() > 0) {
  		this->numNodes = t_path_ptr->edges.size() + 1;
  		shared_ptr<Edge> firstEdge = t_path_ptr->edges.front();
  		shared_ptr<Edge> lastEdge = t_path_ptr->edges.back();
  		this->startNodeName = firstEdge->getStartNodeName();
  		this->endNodeName = lastEdge->getEndNodeName();
  		this->pathDir = D_LEFT;
      	if (firstEdge->QES1 < firstEdge->QES2) this->pathDir = D_RIGHT;

  		// NOTE: All edges should be extending the Query with the Target to the right
  		// Correct values should be stored in the edge object
  		// There should be no switching strands

  		// Add the length of the first Node (start node of the first edge)
  		this->length += firstEdge->SLen;
	  	for (auto const& edge : t_path_ptr->edges) {
	  		this->avgSI += edge->SI;
	  		// For each edge, add to the length of the path part of the end node that 
	  		// does not overlap with the start node (End overhang - start overhang)
	  		this->length += (edge->ELen - edge->EEnd) - (edge->SLen - edge->SEnd);
	  		this->length2 += edge->SStart - edge->EStart;
	  		
	        if ((edge->SStart <= edge->EStart) || ((edge->ELen - edge->EEnd) <= (edge->SLen - edge->SEnd))){
				  negativeEScount++;
	        }
	    }
	    this->length2 += lastEdge->ELen;

	  	this->avgSI /= t_path_ptr->edges.size();
      	if (negativeEScount > 0) 
         	throw std::runtime_error(std::string("SCARA BRIDGER: ERROR - path with negative extensions (" 
        		+ std::to_string(negativeEScount) + ") SN(" + startNodeName + ") EN(" + endNodeName + ")"
         		+ " - " + t_path_ptr->toString()));
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

  // KK: Allow adding a path info for a reverse path!
  bool PathGroup::addPathInfo(shared_ptr<PathInfo> pinfo_ptr) {
    bool equal = false;

    // New path belong in the group in original orientation
    if ((pinfo_ptr->startNodeName.compare(this->startNodeName) == 0) &&
        (pinfo_ptr->endNodeName.compare(this->endNodeName) == 0) &&
        (fabs(pinfo_ptr->length - this->length) <= scara::pathGroupHalfSize)) equal = true;

    /*
    // New path belongs in the group in reverse orientation
    if ((pinfo_ptr->startNodeName.compare(getRCNodeName(this->endNodeName)) == 0) &&
        (pinfo_ptr->endNodeName.compare(getRCNodeName(this->startNodeName)) == 0) &&
        (fabs(pinfo_ptr->length - this->length) <= scara::pathGroupHalfSize)) equal = true;
    */

    if (equal) {
      vPathInfos.emplace_back(pinfo_ptr);
      numPaths += 1;
    }
    return equal;

  }

}