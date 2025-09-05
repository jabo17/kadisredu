//
// Created by jannickb on 6/18/24.
//

#include "kadisredu/algorithms/mwis/greedy.h"

namespace kadisredu::mwis {

void greedy::exclude_node(dist_dynamic_graph& graph, NodeID node) {
  // exclude
  status[node] = wis_status::EXCLUDED;
  graph.hide_node(node);
  if (graph.is_active_border_node(node)) {
    on_modified_border_node(graph, node);
  }
  // mark neighbors
  for (auto& target : graph[node]) {
    if (!graph.is_ghost(target)) {
      candidates.mark(target);
    }
  }
}
void greedy::include_node(dist_dynamic_graph& graph, NodeID node) {
  // include node
  KASSERT(status[node] == wis_status::UNSET);
  status[node] = wis_status::INCLUDED;
  graph.hide_node(node);
  if (graph.is_active_border_node(node)) {
    on_modified_border_node(graph, node);
  }
  if (!graph.is_ghost(node)) {
    wis_weight += graph.get_weight(node);
  }
  // exclude neighbors
  for (auto& target : graph[node]) {
    exclude_node(graph, target);
  }
}
bool greedy::check_for_include(const dist_dynamic_graph& graph, NodeID node) {
  KASSERT(!graph.is_ghost(node));
  if (status[node] != wis_status::UNSET) {
    return false;
  }
  auto& neigh = graph[node];
  auto rating = ratings[node];
  if (std::any_of(neigh.begin(), neigh.end(),
                  [&, global_node = node + graph.get_first_global_node()](NodeID target) {
                    KASSERT(graph.is_visible(target) && target != node);
                    return ratings[target] > rating ||
                           (ratings[target] == rating && graph.is_ghost(target) &&
                            graph.get_ghost_block(target) > comm.rank());
                  })) {
    return false;
  }
  return true;
}
}  // namespace kadisredu::mwis
