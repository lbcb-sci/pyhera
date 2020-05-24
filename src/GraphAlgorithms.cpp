#include <set>
#include <stack>
#include <queue>
#include <stdexcept>
#include <algorithm>
#include <random>
#include <iostream>

#include "Graph.h"
#include "Overlap.h"
#include "Sequence.h"
#include "globals.h"

namespace scara {

  using namespace std;
  extern int multithreading;


  // Checking if a Path is consistent
  // Checks if the path maintains direction along all edges and 
  // Checks if adjecent edges share the same node
  // returns 0 if the path is consistent
  int checkPath(shared_ptr<Path> path_ptr) {
  	if (path_ptr->edges.size() > 0) {
      shared_ptr<Edge> firstedge = path_ptr->edges.front();
      Direction dir = D_LEFT;
      if (firstedge->QES1 < firstedge->QES2) dir = D_RIGHT;
      for (uint32_t i = 1; i < path_ptr->edges.size(); i++) {
      	shared_ptr<Edge> edge = path_ptr->edges[i];
        Direction dir2 = D_LEFT;
        if (edge->QES1 < edge->QES2) dir2 = D_RIGHT;
        if (dir != dir2) return 1;
        if (firstedge->endNode != edge->startNode) return 2;
        firstedge = edge;
      }
    }

    return 0;
  }


  // Compare functions for sorting
  bool sortByOS(const shared_ptr<Edge> &lhs, const shared_ptr<Edge> &rhs) { return lhs->OS > rhs->OS; }

  bool sortByESLeft(const shared_ptr<Edge> &lhs, const shared_ptr<Edge> &rhs) { return lhs->QES1 > rhs->QES1; }

  bool sortByESRight(const shared_ptr<Edge> &lhs, const shared_ptr<Edge> &rhs) { return lhs->QES2 > rhs->QES2; }



  /*
   * A function that returns N best edges according to a specified criterion
   * Vector containing Edges is passed as reference because the sort will change it!
   */
  unique_ptr<vector<shared_ptr<Edge>>> getBestNEdges(vector<shared_ptr<Edge>> &edges, uint32_t N, PathGenerationType pgType) {
  	unique_ptr<vector<shared_ptr<Edge>>> pvEdges = make_unique<vector<shared_ptr<Edge>>>();

  	if (N > edges.size()) N = edges.size();

  	if (N > 0) {
	  	switch (pgType) {
	  		case (PGT_MAXOS):
	  			std::sort(edges.begin(), edges.end(), sortByOS);
	  			break;
	  		case (PGT_MAXESLEFT):
	  			std::sort(edges.begin(), edges.end(), sortByESLeft);
	  			break;
	  		case (PGT_MAXESRIGHT):
	  			std::sort(edges.begin(), edges.end(), sortByESRight);
	  			break;
	  		default:
	  			throw std::runtime_error(std::string("SCARA BRIDGER: ERROR - invalid path generation type!"));
	  	}
  	}
  	for (uint32_t i = 0u; i < N; i++) {
  		pvEdges->emplace_back(edges[i]);
  	}

  	return pvEdges;
  }


  /*
   * Generate paths choosing an edge with maximum overlap score or maximum extension scorein each step
   * In the first stop, for each anchor node consider all outgoiing edges
   */
  int generatePathsDeterministic(vector<shared_ptr<Path>> &vPaths, MapIdToNode &aNodes, PathGenerationType pgType){
  	int pathsGenerated = 0;

  	/* Each read can only be used once
  	 * TODO:
  	 * Currently placing read names in a set
  	 * It might be more efficient if Node pointers were used!
  	 */
  	std::set<std::string> readsUsed;
  	// uint32_t numNodes = 10;		// Number of nodes placed on the stack in each step of graph traversal
  	uint32_t numNodes = scara::NumDFSNodes;

  	if (scara::print_output)
  		std::cerr << "\nSCARA: Generating deterministic paths: ";

  	for (auto const& itANode : aNodes) {
  		if (scara::print_output)
  			std::cerr << ".";			// Printing one dot for each attempt at generating a path
  		std::string aNodeName = getOGNodeName(itANode.first);
  		std::shared_ptr<Node> aNode = itANode.second;
  		for (auto const& edge_ptr : aNode->vOutEdges) {
  			auto newPath = make_shared<Path>(edge_ptr);
  			std::stack<std::shared_ptr<Edge>> eStack;

  			// IMPORTANT: we are always extending query with the target and to the right
  			//			   using appropriate extension score

  			// KK: Control, check if extension scores are greater than 0
            if ((edge_ptr->QES2 <= 0)) continue;
            
            eStack.push(edge_ptr);
            /* For each inital node, place only its edge on the stack
             * In each step of graph traversal:
             * - Pop the last node
             * - Check if it can connect to an anchor node
             * - If it can, the path is complete
             * - If not, get a number of connected read nodes with the greatest score and place them on the stack
             * - If no reads are available, adjust the path and continue
             */ 
            while (!eStack.empty()) {
                std::shared_ptr<Edge> redge_ptr = eStack.top();                    // Pop an edge from the stack
                eStack.pop();
                if (redge_ptr == NULL) throw std::runtime_error(std::string("NULL pointer edge on the DFS stack!"));
                std::shared_ptr<Node> rnode_ptr = redge_ptr->endNode;              // And the corresponding node

                // Check if the node from the stack can continue the current path
                if ((newPath->size() > 0) && (newPath->endNode() != redge_ptr->startNode)) {
                    // If not, put the edge back on the stack
                    eStack.push(redge_ptr);
                    // And remove the last edge from the path
                    newPath->removeLastEdge();
                    // Skip to next iteration
                    continue;
                }

                // Check if the path is too long, skip this iteration and let
                // the above code eventually reduce the path
                if ((uint32_t)(newPath->size()) >= scara::HardNodeLimit) continue;

                newPath->appendEdge(redge_ptr);                           // Add edge to the path
                readsUsed.insert(rnode_ptr->nName);                       // And mark the node as traversed

                std::vector<shared_ptr<Edge>> Aedges;                                     // Edges to anchor nodes
                std::vector<shared_ptr<Edge>> Redges;                                     // Edges to read nodes

                for (auto const& edge2_ptr : rnode_ptr->vOutEdges) {
                    // KK: Control
                    if ((edge2_ptr->QES2 <= 0)) continue;

                    shared_ptr<Node> endNode = edge2_ptr->endNode;
                    if (readsUsed.find(endNode->nName) != readsUsed.end())        // Each read can only be used once
                        continue;

                    if (endNode->nType == NT_ANCHOR) {
                        if (getOGNodeName(endNode->nName).compare(aNodeName) != 0)          // We only want nodes that are different from the starting node!
                            Aedges.emplace_back(edge2_ptr);               					// NOTE: this might change, as we migh want to scaffold circulat genomes!
                    }
                    else if (endNode->nType == NT_READ)
                        Redges.emplace_back(edge2_ptr);
                    else throw std::runtime_error(std::string("SCARA BRIDGER: ERROR - invalid node type: ") + scara::NodeType2String(endNode->nType));
                }

                PathGenerationType pgType2 = PGT_INVALID;
                if (pgType == PGT_MAXOS) {
                	pgType2 = PGT_MAXOS;
                } else {
                	pgType2 = PGT_MAXESRIGHT;
                }

                unique_ptr<vector<shared_ptr<Edge>>> bestAedges = getBestNEdges(Aedges, 1, pgType2);
                unique_ptr<vector<shared_ptr<Edge>>> bestRedges = getBestNEdges(Redges, numNodes, pgType2);

                if (bestAedges->size() > 0u) {                                  // If anchor nodes have been reached find the best one
                    shared_ptr<Edge> aedge = (*bestAedges)[0];					// Since we are using a priority queue, no need to sort the elements
                    newPath->appendEdge(aedge);									// Create a path and end this instance of tree traversal
                    vPaths.emplace_back(newPath);
                    pathsGenerated++;
                    break;
                } else if (bestRedges->size() > 0u) {                             // If no anchor nodes have been found we have to continue with read nodes
                	uint32_t N = (numNodes < bestRedges->size() ? numNodes : bestRedges->size());
                    for (int i=N-1; i>=0; i--) {								// Place N best edges on the stack in reverse order, so that the best one ends on top
                    	shared_ptr<Edge> t_edge = (*bestRedges)[i];
                    	eStack.push(t_edge);
                    }                    
                } else {                                                       // Graph traversal has come to a dead end
             		// Remove the last edge from the path, and switch direction back if necessary
                	auto lastEdge_ptr = newPath->removeLastEdge();
                	std::set<std::string>::iterator it = readsUsed.find(rnode_ptr->nName);
                    readsUsed.erase(it);                         				// Remove current read node from the list of traversed ones
                }
            }
  		}
  	}

  	return pathsGenerated;
  }


  /*
   * Generate paths choosing an edge with the probability proportional to the extension score
   * Using Monte Carlo approach
   */
  int generatePaths_MC(std::vector<shared_ptr<Path>> &vPaths, MapIdToNode &aNodes, uint32_t minNumPaths) {
  	uint32_t pathsGenerated = 0;

  	/* Each read can only be used once in a path!
  	 * TODO:
  	 * Currently placing read names in a set
  	 * It might be more efficient if Node pointers were used!
  	 */
  	uint32_t maxIterations = scara::MaxMCIterations;
  	uint32_t iteration = 0;
  	uint32_t numNodes = scara::NumDFSNodes;

  	// Setting up random number generator
  	std::random_device rd;
  	std::default_random_engine generator{rd()};
  	std::uniform_int_distribution<uint32_t> dist{0, aNodes.size()-1};		// For randomly choosing an Anchor node
  	std::vector<shared_ptr<Node>> vANodes;								              // A vector for faster random Node access
  	for (auto const& itANode : aNodes) vANodes.emplace_back(itANode.second);

  	while (pathsGenerated < minNumPaths && iteration < maxIterations) {
  		iteration += 1;

  		// Randomly choose an anchor Node
  		std::shared_ptr<Node> aNode = vANodes[dist(generator)];
  		std::string aNodeName = getOGNodeName(aNode->nName);

  		if (aNode->vOutEdges.size() == 0) continue;			// Probably not necessary
  		float totalES = 0.0;
  		for (auto const& edge_ptr : aNode->vOutEdges) {
  			totalES += edge_ptr->QES2;
  		}

  		// Generating real numbers to select an edge with probability proportional to Extension score
  		std::uniform_real_distribution<float> dist2{0, totalES};
  		float rndCumES = dist2(generator);
  		float cumulativeES = 0.0;
  		std::shared_ptr<Edge> chosen_edge_ptr;
  		for (auto const& edge_ptr : aNode->vOutEdges) {
  			cumulativeES += edge_ptr->QES2;
  			if (cumulativeES > rndCumES) {
  				chosen_edge_ptr = edge_ptr;
  				break;
  			}
  		}

  		// Initialize new path and stack for graph traversal
  		auto newPath = make_shared<Path>(chosen_edge_ptr);
  		std::stack<std::shared_ptr<Edge>> eStack;

		// KK: Control, check if estension scores are greater than 0
        if ((chosen_edge_ptr->QES2 <= 0)) continue;

        eStack.push(chosen_edge_ptr);
        /* For each inital node, place only its edge on the stack
         * In each step of graph traversal:
         * - Pop the last node
         * - Check if it can connect to an anchor node
         * - If it can, the path is complete
         * - If not,randomly generate a number of connected read nodes with the probability of generation
         *   proportional to ES and place them on the stack
         * - If no reads are available, adjust the path and continue
         */
        std::set<std::string> readsUsed;		// A read cannot be used more than once within the same path
        while (!eStack.empty()) {
            std::shared_ptr<Edge> redge_ptr = eStack.top();                    // Pop an edge from the stack
            eStack.pop();
            std::shared_ptr<Node> rnode_ptr = redge_ptr->endNode;              // And the corresponding node

            // Check if the node from the stack can continue the current path
            if ((newPath->size() > 0) && (newPath->endNode() != redge_ptr->startNode)) {
                // If not, put the edge back on the stack
                eStack.push(redge_ptr);
                // And remove the last edge from the path
                newPath->removeLastEdge();
                // Skip to next step
                continue;
            }

            // Check if the path is too long skip this iteration and let
            // the above code eventually reduce the path
            if ((uint32_t)(newPath->size()) >= scara::HardNodeLimit) continue;

            newPath->appendEdge(redge_ptr);                           // Add edge to the path
            readsUsed.insert(rnode_ptr->nName);                       // And mark the node as traversed

            std::vector<shared_ptr<Edge>> Aedges;                                     // Edges to anchor nodes
            std::vector<shared_ptr<Edge>> Redges;                                     // Edges to read nodes

            for (auto const& edge2_ptr : rnode_ptr->vOutEdges) {
                // KK: Control
                if ((edge2_ptr->QES2 <= 0) && (edge2_ptr->QES1 <= 0)) continue;

                shared_ptr<Node> endNode = edge2_ptr->endNode;
                if (readsUsed.find(endNode->nName) != readsUsed.end())        // Each read can only be used once
                    continue;

                if (endNode->nType == NT_ANCHOR) {
                    if (getOGNodeName(endNode->nName).compare(aNodeName) != 0)             // We only want nodes that are different from the starting node!
                        Aedges.emplace_back(edge2_ptr);               		   // NOTE: this might change, as we migh want scaffold circulat genomes!
                }
                else if (endNode->nType == NT_READ)
                    Redges.emplace_back(edge2_ptr);
                else throw std::runtime_error(std::string("SCARA BRIDGER: ERROR - invalid node type: ") + scara::NodeType2String(endNode->nType));
            }

            PathGenerationType pgType = PGT_MAXESRIGHT;

            shared_ptr<vector<shared_ptr<Edge>>> bestAedges = getBestNEdges(Aedges, 1, pgType);
            // shared_ptr<vector<shared_ptr<Edge>>> bestRedges = getBestNEdges(Redges, numNodes, pgType2);

            if (bestAedges->size() > 0) {                                   // If anchor nodes have been reached find the best one
                shared_ptr<Edge> aedge = (*bestAedges)[0];
                newPath->appendEdge(aedge);									// Create a path and end this instance of tree traversal
                vPaths.emplace_back(newPath);
                pathsGenerated++;
                break;
            } 

            if (Redges.size() > 0) {                             			// If no anchor nodes have been found we have to continue with read nodes
            	float totalES = 0.0;
		  		for (auto const& edge_ptr : Redges) {
		  			totalES += edge_ptr->QES2;
		  		}

		  		// Generating real numbers to select an edge with probability proportional to Extension score
		  		// Do it numNodes times, it is possible to generate the same edge more than once
		  		// Place each edge on the stack
		  		std::uniform_real_distribution<float> dist2{0, totalES};
            	uint32_t N = (numNodes < Redges.size() ? numNodes : Redges.size());
                for (uint32_t i=0; i<N; i++) {
		  			float rndCumES = dist2(generator);
                	float cumulativeES = 0.0;
		  			std::shared_ptr<Edge> chosen_edge_ptr;
		  			for (auto const& edge_ptr : Redges) {
		  				cumulativeES += edge_ptr->QES2;
			  			if (cumulativeES > rndCumES) {
			  				eStack.push(edge_ptr);
			  			}
			  		}
                }                    
            } else {                                                       // Graph traversal has come to a dead end
                newPath->removeLastEdge();                                 // Remove the last edge from the path
                std::set<std::string>::iterator it = readsUsed.find(rnode_ptr->nName);
                readsUsed.erase(it);                         				// Remove current read node from the list of traversed ones
            }
        }
  	}

  	if (scara::print_output) std::cerr << "\nFinished Monte Carlo with " << iteration << " iterations!";
  	return pathsGenerated;
  }


  /*
   * Generate paths choosing an edge with the probability proportional to the extension score
   * Using Monte Carlo approach 
   * Generate paths only for a single node
   */
  int generatePathsForNode_MC(std::vector<shared_ptr<Path>> &vPaths, shared_ptr<Node> aNode, int minNumPaths, int maxNumIterations) {
  	uint32_t pathsGenerated = 0;

  	/* Each read can only be used once in a path!
  	 * TODO:
  	 * Currently placing read names in a set
  	 * It might be more efficient if Node pointers were used!
  	 */
  	uint32_t maxIterations = maxNumIterations;
  	uint32_t iteration = 0;
  	uint32_t numNodes = scara::NumDFSNodes;

  	// Setting up random number generator
  	std::random_device rd;
  	std::default_random_engine generator{rd()};

  	while (pathsGenerated < minNumPaths && iteration < maxIterations) {
  		iteration += 1;

  		// Anchor Node is set through arguments
  		std::string aNodeName = getOGNodeName(aNode->nName);

  		if (aNode->vOutEdges.size() == 0) break;			// Probably not necessary
  		float totalES = 0.0;
  		for (auto const& edge_ptr : aNode->vOutEdges) {
  			float maxES = (edge_ptr->QES1 > edge_ptr->QES2) ? edge_ptr->QES1 : edge_ptr->QES2;
  			totalES += maxES;
  		}

  		// Generating real numbers to select an initial edge with probability proportional to Extension score
  		std::uniform_real_distribution<float> dist2{0, totalES};
  		float rndCumES = dist2(generator);
  		float cumulativeES = 0.0;
  		std::shared_ptr<Edge> chosen_edge_ptr;
  		// TODO: this could be done using binary search!
  		for (auto const& edge_ptr : aNode->vOutEdges) {
  			float maxES = (edge_ptr->QES1 > edge_ptr->QES2) ? edge_ptr->QES1 : edge_ptr->QES2;
  			cumulativeES += maxES;
  			if (cumulativeES > rndCumES) {
  				chosen_edge_ptr = edge_ptr;
  				break;
  			}
  		}

  		// Initialize new path and stack for graph traversal
  		auto newPath = make_shared<Path>(chosen_edge_ptr);
  		std::stack<std::shared_ptr<Edge>> eStack;

		// KK: Control, check if estension scores are greater than 0
        if ((chosen_edge_ptr->QES2 <= 0) && (chosen_edge_ptr->QES1 <= 0)) continue;

        eStack.push(chosen_edge_ptr);
        /* For each inital node, place only its edge on the stack
         * In each step of graph traversal:
         * - Pop the last node
         * - Check if it can connect to an anchor node
         * - If it can, the path is complete
         * - If not,randomly generate a number of connected read nodes with the probability of generation
         *   proportional to ES and place them on the stack
         * - If no reads are available, adjust the path and continue
         */
        std::set<std::string> readsUsed;		// A read cannot be used more than once within the same path
        while (!eStack.empty()) {
            std::shared_ptr<Edge> redge_ptr = eStack.top();                    // Pop an edge from the stack
            eStack.pop();
            std::shared_ptr<Node> rnode_ptr = redge_ptr->endNode;              // And the corresponding node

            // Check if the node from the stack can continue the current path
            if ((newPath->size() > 0) && (newPath->endNode() != redge_ptr->startNode)) {
                // If not, put the edge back on the stack
                eStack.push(redge_ptr);
                // And remove the last edge from the path
                newPath->removeLastEdge();
                // Skip to next step
                continue;
            }

            // Check if the path is too long skip this iteration and let
            // the above code eventually reduce the path
            if ((uint32_t)(newPath->size()) >= scara::HardNodeLimit) continue;

            newPath->appendEdge(redge_ptr);                           // Add edge to the path
            readsUsed.insert(rnode_ptr->nName);                       // And mark the node as traversed

            std::vector<shared_ptr<Edge>> Aedges;                                     // Edges to anchor nodes
            std::vector<shared_ptr<Edge>> Redges;                                     // Edges to read nodes

            for (auto const& edge2_ptr : rnode_ptr->vOutEdges) {
                // KK: Control
                if ((edge2_ptr->QES2 <= 0) && (edge2_ptr->QES1 <= 0)) continue;

                shared_ptr<Node> endNode = edge2_ptr->endNode;
                if (readsUsed.find(endNode->nName) != readsUsed.end())        // Each read can only be used once
                    continue;

                if (endNode->nType == NT_ANCHOR) {
                    if (getOGNodeName(endNode->nName).compare(aNodeName) != 0)          // We only want nodes that are different from the starting node!
                        Aedges.emplace_back(edge2_ptr);               		   			// NOTE: this might change, as we migh want scaffold circulat genomes!
                }
                else if (endNode->nType == NT_READ)
                    Redges.emplace_back(edge2_ptr);
                else throw std::runtime_error(std::string("SCARA BRIDGER: ERROR - invalid node type: ") + scara::NodeType2String(endNode->nType));
            }

            PathGenerationType pgType = PGT_MAXESRIGHT;

            shared_ptr<vector<shared_ptr<Edge>>> bestAedges = getBestNEdges(Aedges, 1, pgType);
            // shared_ptr<vector<shared_ptr<Edge>>> bestRedges = getBestNEdges(Redges, numNodes, pgType2);

            if (bestAedges->size() > 0) {                                   // If anchor nodes have been reached find the best one
                shared_ptr<Edge> aedge = (*bestAedges)[0];
                newPath->appendEdge(aedge);									// Create a path and end this instance of tree traversal
                vPaths.emplace_back(newPath);
                pathsGenerated++;
                break;
            } 

            if (Redges.size() > 0) {                             			// If no anchor nodes have been found we have to continue with read nodes
            	float totalES = 0.0;
		  		for (auto const& edge_ptr : Redges) {
		  			totalES += edge_ptr->QES2;
		  		}

		  		// Generating real numbers to select an edge with probability proportional to Extension score
		  		// Do it numNodes times, it is possible to generate the same edge more than once
		  		// Place each edge on the stack
		  		std::uniform_real_distribution<float> dist2{0, totalES};
          uint32_t N = (numNodes < Redges.size() ? numNodes : Redges.size());
          for (uint32_t i=0; i<N; i++) {
		  			float rndCumES = dist2(generator);
            float cumulativeES = 0.0;
		  			std::shared_ptr<Edge> chosen_edge_ptr;
		  			for (auto const& edge_ptr : Redges) {
		  				float maxES = (edge_ptr->QES1 > edge_ptr->QES2) ? edge_ptr->QES1 : edge_ptr->QES2;
		  				cumulativeES += maxES;
			  			if (cumulativeES > rndCumES) {
			  				eStack.push(edge_ptr);
			  			}
			  		}
          }                    
        
        } else {                                                       // Graph traversal has come to a dead end
          newPath->removeLastEdge();                                 // Remove the last edge from the path
          std::set<std::string>::iterator it = readsUsed.find(rnode_ptr->nName);
          readsUsed.erase(it);                         				// Remove current read node from the list of traversed ones
        }
      }
  	}

  	if (scara::print_output) std::cerr << "\nFinished Monte Carlo with " << iteration << " iterations!";
  	return pathsGenerated;
  }
}
