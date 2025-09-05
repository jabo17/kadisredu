//
// Created by jannickb on 11/27/24.
//

#pragma once

#include "kadisredu/data_structures/dist_dynamic_graph.h"
#include "kadisredu/definitions.h"

#include <kamping/communicator.hpp>

namespace kadisredu {

class graph_summary {
 public:
  GlobalNodeID nodes;
  NodeID max_nodes;
  NodeID min_nodes;
  GlobalEdgeID edges;
  GlobalNodeID ghosts;
  NodeID max_ghosts;
  NodeID min_ghosts;
  GlobalEdgeID cut;
  double imbalance;
  NodeID max_degree;
  NodeID min_degree;

  static constexpr NodeID MAX_DEGREE = 10;  // MAX_DEGREE+1 buckets for degree counts
  std::vector<NodeID> degree_counts;

  static std::optional<graph_summary> agg_graph_summary_on_root(const kamping::Communicator<> &comm,
                                                                dist_dynamic_graph &graph);
};
}  // namespace kadisredu
