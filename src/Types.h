#pragma once

#include <map>
#include <vector>
#include <memory>
#include <string>

namespace scara {
  class Sequence;
  class Overlap;
  class Node;
  class Edge;

  using MapIdToOvl = std::map<std::pair<std::string, uint32_t>, std::vector<std::shared_ptr<Overlap>>>;
  using MapIdToSeq = std::map<std::string, std::shared_ptr<Sequence>>;

  using MapIdToNode = std::map<std::string, std::shared_ptr<Node>>;

}
