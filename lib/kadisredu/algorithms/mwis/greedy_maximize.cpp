//
// Created by jannickb on 6/18/24.
//

#include "kadisredu/algorithms/mwis/greedy_maximize.h"

namespace kadisredu::mwis {

void greedy_maximize::exclude_node(const dist_dynamic_graph& graph, NodeID node) {
  // exclude
  if (status[node] != wis_status::UNSET) {
    return;
  }
  status[node] = wis_status::EXCLUDED;
  if (graph.is_active_border_node(node)) {
    on_modified_border_node(graph, node);
  }
  // mark neighbors
  for (auto& target : graph[node]) {
    if (!graph.is_ghost(target) && status[target] == wis_status::UNSET &&
        consider_next_set.add(target)) {
      consider_next.push_back(target);
    }
  }
}
void greedy_maximize::include_node(const dist_dynamic_graph& graph, NodeID node) {
  // include node
  KASSERT(status[node] == wis_status::UNSET);
  status[node] = wis_status::INCLUDED;
  if (graph.is_active_border_node(node)) {
    on_modified_border_node(graph, node);
  }
  if (!graph.is_ghost(node)) {
    wis_weight += graph.get_weight(node);
    // KADISREDU_LOG << "included a node of weight " << graph.get_weight(node);
  }
  // exclude neighbors
  for (auto& target : graph[node]) {
    exclude_node(graph, target);
  }
}
void greedy_maximize::sync_ghost_status(const dist_dynamic_graph& graph) {
  // sync solution of ghosts
  auto ghost_status = graph.sync_ghosts_with_lazy_updates<wis_status>(
      graph.get_border(), [&](NodeID node) -> std::optional<wis_status> {
        if (status[node] == wis_status::INCLUDED) {
          return status[node];
        } else {
          return std::nullopt;
        }
      });
  for (auto& [ghost, g_status] : ghost_status) {
    KASSERT(g_status == wis_status::INCLUDED, static_cast<int>(g_status));
    status[graph.get_ghost(ghost)] = wis_status::INCLUDED;
  }
}
void greedy_maximize::init_free_nodes(const dist_dynamic_graph& graph) {
  static_assert(std::numeric_limits<GlobalNodeWeight>::is_signed);

  NodeID count_free = 0;
  graph.for_each_visible_node([&](NodeID node) {
    // free node?
    if (status[node] != wis_status::INCLUDED &&
        std::find_if(graph[node].begin(), graph[node].end(), [&](NodeID neighbor) {
          return status[neighbor] == wis_status::INCLUDED;
        }) == graph[node].end()) {
      consider_next_set.add(node);
      status[node] = wis_status::UNSET;
      ++count_free;

      // add rating (larger rating speaks in favor for include)
      // @todo do we know the correct weight of ghost (the correct rating for border vertices)?
      ratings[node] = graph.get_weight(node);  //- graph.get_weight_neigh(node);
    }
  });

  // insert candidates
  consider_next.reserve(count_free);
  consider.reserve(count_free);
  graph.for_each_visible_node([&](NodeID node) {
    if (status[node] == wis_status::UNSET) {
      consider_next.push_back(node);
    }
  });

  auto rated_ghosts = graph.sync_ghosts_with_lazy_updates<NodeWeight>(
      graph.get_border(), [&](NodeID node) -> std::optional<NodeWeight> {
        KASSERT(graph.is_original_border_node(node));
        if (status[node] == wis_status::UNSET) {
          return ratings[node];
        } else {
          return std::nullopt;
        }
      });

  for (auto [global_node, rating] : rated_ghosts) {
    auto ghost = graph.get_ghost(global_node);
    ratings[ghost] = rating;
    status[ghost] = wis_status::UNSET;
#if KASSERT_ENABLED(KASSERT_ASSERTION_LEVEL_NORMAL)
    bool free = true;
    for (auto target : graph[ghost]) {
      if (status[target] == wis_status::INCLUDED) {
        free = false;
      }
    }
    KASSERT(free, "include(s) at the border were not sent to ghost ranks");
#endif
  }
}
bool greedy_maximize::check_for_include(const dist_dynamic_graph& graph, NodeID node) {
  KASSERT(!graph.is_ghost(node));
  if (status[node] != wis_status::UNSET) {
    return false;
  }
  auto& neigh = graph[node];
  auto rating = ratings[node];
  if (std::any_of(neigh.begin(), neigh.end(),
                  [&, global_node = node + graph.get_first_global_node()](NodeID target) {
                    KASSERT(graph.is_visible(target) && target != node);
                    return status[target] == wis_status::UNSET &&
                           (ratings[target] > rating ||
                            (ratings[target] == rating && graph.is_ghost(target) &&
                             graph.get_ghost_block(target) > comm.rank()));
                  })) {
    return false;
  }
  return true;
}
}  // namespace kadisredu::mwis
