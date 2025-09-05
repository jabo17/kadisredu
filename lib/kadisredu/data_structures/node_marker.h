//
// Created by jannickb on 10/19/24.
//

#pragma once

#include "kadisredu/definitions.h"
#include "kadisredu/data_structures/fast_set.h"
#include "kadisredu/data_structures/sized_vector.h"

namespace kadisredu {

struct NodeMarker {
  explicit NodeMarker(NodeID capacity)
      : current_marked_nodes(capacity), next_marked_nodes(capacity), marked_next_set(capacity) {}

  // note that we use sized_vector here as it allows us to "clear" vector in
  // constant time
  sized_vector<NodeID> current_marked_nodes;
  sized_vector<NodeID> next_marked_nodes;
  fast_set marked_next_set;

  sized_vector<NodeID>& current() { return current_marked_nodes; }

  bool mark(NodeID node) {
    if (marked_next_set.add(node)) {
      next_marked_nodes.push_back(node);
      return true;
    }
    return false;
  }

  void next() {
    using std::swap;
    swap(current_marked_nodes, next_marked_nodes);
    next_marked_nodes.clear();
    marked_next_set.clear();
  }
};

}  // namespace kadisredu