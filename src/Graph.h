#include "Types.h"
#include <string>
#include "Sequence.h"
#include "Overlap.h"

using namespace std;

namespace scara {

  enum NodeType {
    NT_INVALID = 0,
    NT_ANCHOR = 1 << 0,
    NT_READ = 1 << 1,
  };

  enum Direction {
    D_LEFT = 0,
    D_RIGHT = 1,
  };

  enum PathGenerationType {
  	PGT_INVALID = 0,
  	PGT_MAXOS = 1 << 0,
  	PGT_MAXESLEFT = 1 << 1,
  	PGT_MAXESRIGHT = 1 << 2,
  	PGT_MAXES = 1 << 3,
  	PGT_MC = 1 << 4,
  };

  std::string NodeType2String(NodeType nType);

  std::string Direction2String(Direction dir);

  Direction oppositeDirection(Direction dir);

  std::string getRCNodeName(const std::string nodeName);

  std::string getOGNodeName(const std::string nodeName);

  class Node {
  public:
    NodeType nType;
    std::string nName;
    std::shared_ptr<Sequence> seq_ptr;

    bool isReverseComplement;

    std::vector<std::shared_ptr<Edge>> vOutEdges;

    Node(NodeType nType, std::string nName, std::shared_ptr<Sequence> seq_ptr);
    Node(NodeType nType, std::string nName, std::shared_ptr<Sequence> seq_ptr, bool isRC);

    Node(std::shared_ptr<Sequence> seq, NodeType nType);
    Node(std::shared_ptr<Sequence> seq, NodeType nType, bool isRC);

    Node(NodeType nType, std::string nName, std::shared_ptr<Sequence> seq_ptr, std::vector<std::shared_ptr<Edge>> &vOutEdges);
    Node(NodeType nType, std::string nName, std::shared_ptr<Sequence> seq_ptr, std::vector<std::shared_ptr<Edge>> &vOutEdges, bool isRC);
  };

  class Edge {
  private:
  	void calcEdgeStats();

  public:
    // KK: removing names for efficiency
    // std::string startNodeName;				// Start Node is Query
    // std::string endNodeName;       // End Node is Target

    std::shared_ptr<Node> startNode;    
    std::shared_ptr<Node> endNode;

    // TODO: Refactor the code to use either start and end terminolgy or query and target

    // Basic data from PAF file, copied from Overlap
    uint32_t SLen;
    uint32_t SStart;
    uint32_t SEnd;
    uint32_t ELen;
    uint32_t EStart;
    uint32_t EEnd;
    bool ovl_bOrientation;

    uint32_t paf_matching_bases;
    uint32_t paf_overlap_length;
    uint32_t paf_mapping_quality;

    // Calculated data (Query = start node, target = end node)
    uint32_t QOH1;		// Query left overhang
    uint32_t QOH2;		// Query right overhang
    uint32_t TOH1;		// Target left overhang
    uint32_t TOH2;		// Target right overhang

    uint32_t QOL;		// Query overlap length
    uint32_t TOL;		// Target overlap length

    float SI;		// Sequence identity

    float OS;		// Overlap score
    float QES1;		// Extension score for extending Query with Target to the left
    float QES2;		// Extension score for extending Query with Target to the right
    float TES1;		// Extension score for extending Target with Query to the left
    float TES2;		// Extension score for extending Target with Query to the right


    Edge(std::shared_ptr<Node> startNode, std::shared_ptr<Node> endNode, std::unique_ptr<Overlap> const& ovl_ptr);

    Edge(std::unique_ptr<Overlap> const& ovl_ptr, MapIdToNode& mAnchorNodes, MapIdToNode& mReadNodes);

    Edge(std::unique_ptr<Overlap> const& ovl_ptr);

  	Edge();						// Create an empty Edge

    std::shared_ptr<Edge> getReversedEdge();

    std::string getStartNodeName();
    std::string getEndNodeName();

    void copyDataFromOverlap(std::unique_ptr<Overlap> const& ovl_ptr);

    void reverseNodes(void);
    void reverseStrand(void);
    void reverseStrandSNode(void);
    void reverseStrandENode(void);

    int test(void);

  };

  // KK: Since the architecture of the graph has been changed, we need to create 4 edges for each overlap
  // We have 2 nodes for each sequence (FW and RC), and need to create edges in both directions but only 
  // For appropriate node combinations
  // Example:
  // 1. if we have overlap between seq1 and seq2, on the same strand (relative strand from PAF = '+')
  // Then we need to create edges for nodes Seq1 and Seq2 (both directions), and for nodes Seq1_RC and Seq2_RC (both directions)
  // 2. if we have overlap between seq1 and seq2, on different strand (relative strand from PAF = '-')
  // Then we need to create edges for nodes Seq1 and Seq2_RC (both directions), and for nodes Seq1_RC and Seq2 (both directions)
  void createEdgesFromOverlap(std::unique_ptr<Overlap> const& ovl_ptr, MapIdToNode& mAnchorNodes
                            , MapIdToNode& mReadNodes, std::vector<shared_ptr<Edge>> &vEdges);



  class Path {
  public:
    std::vector<std::shared_ptr<Edge>> edges;

    Path(std::shared_ptr<Edge> edge_ptr);
    Path(void);

    int size(void);
    std::shared_ptr<Node> endNode(void);
    std::shared_ptr<Node> startNode(void);

    void appendEdge(std::shared_ptr<Edge> edge_ptr);
    std::shared_ptr<Edge> removeLastEdge(void);

    shared_ptr<Path> reversedPath(void);

    std::string toString(void);
  };


  /* KK:
   * Functions for algorithms on the graph
   */
  int checkPath(shared_ptr<Path> path);

  int generatePathsDeterministic(std::vector<shared_ptr<Path>> &vPaths, MapIdToNode &aNodes, PathGenerationType pgType);
  int generatePaths_MC(std::vector<shared_ptr<Path>> &vPaths, MapIdToNode &aNodes, uint32_t minNumPaths);

  int generatePathsForNode_MC(std::vector<shared_ptr<Path>> &vPaths, shared_ptr<Node> aNode, uint32_t minNumPaths, uint32_t maxNumIterations);

 
  unique_ptr<vector<shared_ptr<Edge>>> getBestNEdges(vector<shared_ptr<Edge>> &edges, uint32_t N, PathGenerationType pgType);


  class PathInfo {
  public:
  	std::string startNodeName;
  	std::string endNodeName;
  	Direction pathDir;
  	uint32_t numNodes;
  	double length;
  	uint32_t length2;
  	double avgSI;
  	shared_ptr<Path> path_ptr;

  	PathInfo(shared_ptr<Path> t_path_ptr);
  };

  class PathGroup {
  public:
  	std::string startNodeName;
  	std::string endNodeName;
  	// A representative length for the group, currently is set to the length of the first path in the group
  	// but could later be changed into something like average length of the paths in the group or some other 
  	// appropriate measure
  	double length;
  	uint32_t numPaths;
  	std::vector<shared_ptr<PathInfo>> vPathInfos;

  	PathGroup();
  	PathGroup(std::string t_startNodeName, std::string t_endNodeName, double t_length);
  	PathGroup(shared_ptr<PathInfo> pathinfo_ptr);
  	bool addPathInfo(shared_ptr<PathInfo> pinfo_ptr);
  };
}
