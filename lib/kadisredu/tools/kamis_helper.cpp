//
// Created by jannickb on 10/31/24.
//

#include "kadisredu/tools/kamis_helper.h"

#include "kadisredu/tools/logger.h"

#include <branch_and_reduce_algorithm.h>

namespace kadisredu {

void kadisredu::kamis_helper::build_induced_graph(const fast_set& nodes_set,
                                                  std::span<const NodeID> nodes,
                                                  dist_dynamic_graph& graph,
                                                  sized_vector<NodeID>& map_to_sub,
                                                  graph_access& neighborhood_graph) {
  map_to_sub.resize(graph.num_nodes() + graph.num_ghosts());

  ::EdgeID edges_count = 0;
  ::NodeID nodes_count = 0;
  for (auto u : nodes) {
    KASSERT(nodes_set.contains(u));
    for (auto x : graph[u]) {
      if (nodes_set.contains(x)) {
        ++edges_count;
      }
    }
    map_to_sub[u] = nodes_count++;
  }

  neighborhood_graph.start_construction(nodes_count, edges_count);

  for (auto u : nodes) {
    auto new_node = neighborhood_graph.new_node();
    neighborhood_graph.setNodeWeight(new_node, graph.get_weight(u));

    for (auto x : graph[u]) {
      if (nodes_set.contains(x)) {
        auto new_edge = neighborhood_graph.new_edge(new_node, map_to_sub[x]);
        neighborhood_graph.setEdgeWeight(new_edge, 1);
      }
    }
  }

  neighborhood_graph.finish_construction();
}
MISConfig kamis_helper::create_mis_config() {
  MISConfig exact_mis_config;
  // Randomization
  exact_mis_config.seed = 0;
  // Output
  exact_mis_config.console_log = false;
  exact_mis_config.check_sorted = true;
  // ILS
  exact_mis_config.ils_iterations = 15000;
  exact_mis_config.force_cand = 4;
  exact_mis_config.sort_freenodes = true;
  // Reductions
  exact_mis_config.perform_reductions = true;
  exact_mis_config.reduction_style = MISConfig::Reduction_Style::NORMAL;

  return exact_mis_config;
}
}  // namespace kadisredu
