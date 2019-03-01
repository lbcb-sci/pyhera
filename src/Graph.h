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


  class Node {
  public:
    NodeType nType;
    std::string nName;
    std::shared_ptr<Sequence> seq_ptr;

    std::vector<std::shared_ptr<Edge>> vOutEdges;

    Node(NodeType nType, std::string nName, std::shared_ptr<Sequence> seq_ptr);

    Node(std::shared_ptr<Sequence> seq, NodeType nType);

    Node(NodeType nType, std::string nName, std::shared_ptr<Sequence> seq_ptr, std::vector<std::shared_ptr<Edge>> vOutEdges);
  };

  class Edge {
  private:
  	void calcEdgeStats();

  public:
    std::string startNodeName;
    std::shared_ptr<Node> startNode;
    std::string endNodeName;
    std::shared_ptr<Node> endNode;
    std::shared_ptr<Overlap> ovl_ptr;

    bool reversed;

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

    Edge(std::string startNodeName, std::shared_ptr<Node> startNode, std::string endNodeName, std::shared_ptr<Node> endNode, std::shared_ptr<Overlap> ovl_ptr);

    Edge(std::shared_ptr<Overlap> ovl, MapIdToNode& mAnchorNodes, MapIdToNode& mReadNodes);

    Edge(std::shared_ptr<Overlap> ovl);

    void reverseNodes(void);

    int test(void);
  };





  /*
   * KK: OBSOLETE
   *
  // KK: A class for comparing Edges according to OS score (for use with the priority queue)
  class EdgeCompareOS {
  public:
  	bool operator()(const std::shared_ptr<Edge> &E1, const shared_ptr<Edge> &E2) const
  	{
  		return (E1->OS < E2->OS);
  	};
  };

  // KK: A class for comparing Edges according to ES score to the left (for use with the priority queue)
  class EdgeCompareESright {
  public:
  	bool operator()(const std::shared_ptr<Edge> &E1, const shared_ptr<Edge> &E2) const
  	{
  		return (E1->TES1 < E2->TES1);
  	};
  };

  // KK: A class for comparing Edges according to ES score to the right (for use with the priority queue)
  class EdgeCompareESright {
  public:
  	bool operator()(const std::shared_ptr<Edge> &E1, const shared_ptr<Edge> &E2) const
  	{
  		return (E1->TES2 < E2->TES2);
  	};
  };
  */
  

  class Path {
  public:
    std::vector<std::shared_ptr<Edge>> edges;

    Path(std::shared_ptr<Edge> edge_ptr);
    Path(void);

    int size(void);
    std::shared_ptr<Node> endNode(void);
    std::shared_ptr<Node> startNode(void);

    void appendEdge(std::shared_ptr<Edge> edge_ptr);
    void removeLastEdge(void);

    shared_ptr<Path> reversedPath();
  };


  /* KK:
   * Functions for algorithms on the graph
   */
  bool checkPath(Path path);

  int generatePathsDeterministic(std::vector<shared_ptr<Path>> &vPaths, MapIdToNode &aNodes, PathGenerationType pgType);
  int generatePaths_MC(std::vector<shared_ptr<Path>> &vPaths, MapIdToNode &aNodes, int minNumPaths);
 
  shared_ptr<vector<shared_ptr<Edge>>> getBestNEdges(vector<shared_ptr<Edge>> &edges, uint32_t N, PathGenerationType pgType);


  class PathInfo {
  public:
  	std::string startNodeName;
  	std::string endNodeName;
  	Direction pathDir;
  	uint32_t numNodes;
  	double length;
  	double avgSI;
  	shared_ptr<Path> path_ptr;

  	PathInfo(shared_ptr<Path> t_path_ptr);
  };

  class PathGroup {
  public:
  	std::string startNodeName;
  	std::string endNodeName;
  	std::vector<shared_ptr<PathInfo>> vPathInfos;

  	PathGroup();
  	void AddPathInfo(shared_ptr<PathInfo> pinfo_ptr);
  };
}
