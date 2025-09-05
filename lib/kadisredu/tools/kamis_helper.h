//
// Created by jannickb on 10/31/24.
//

#pragma once

#include "kadisredu/data_structures/dist_dynamic_graph.h"
#include "kadisredu/definitions.h"

#include <graph_access.h>
#include <mis_config.h>
#include <span>

namespace kadisredu {

class kamis_helper {
 public:
  static void build_induced_graph(const fast_set &nodes_set, std::span<const NodeID> nodes,
                                  dist_dynamic_graph &graph, sized_vector<NodeID> &map_to_sub,
                                  graph_access &neighborhood_graph);

  static MISConfig create_mis_config();
};

}  // namespace kadisredu
