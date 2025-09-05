//
// Created by jannickb on 6/18/24.
//

#include "kadisredu/algorithms/mwis/sync_greedy_maximize.h"

#include "kadisredu/tools/logger.h"

namespace kadisredu::mwis::sync {
void sync_greedy_maximize::run(const dist_dynamic_graph& graph, const Solution& solution) {
  using namespace kamping;

  KASSERT(solution.node_status.size() == graph.num_nodes() + graph.num_ghosts());
  status = solution.node_status;
  wis_weight = solution.solution_weight;

  // sync ghost status to obtain included ghosts
  KADISREDU_RLOG << "[greedy maximize]: start";
  sync_ghost_status(graph);
  KADISREDU_RLOG << "synced ghosts";
  // sync free nodes at the border
  init_free_nodes(graph);
  KADISREDU_RLOG << "init free nodes";

  // add nodes greedily to solution
  while (true) {
    while (!consider_next.empty()) {
      std::swap(consider, consider_next);
      consider_next_set.clear();
      for (auto node : consider) {
        KASSERT(!graph.is_ghost(node));
        if (check_for_include(graph, node)) {
          // node has maximal rating in remaining neighborhood
          set_node(graph, node, wis_status::INCLUDED);
        }
      };
      consider.clear();
    }

    bool progress =
        comm.allreduce_single(send_buf((kabool)!border_updates.empty()), op(ops::logical_or<>()));
    if (!progress) {
      break;
    }

    // sync borders
    auto ghost_updates = graph.sync_ghosts(border_updates);
    ++rounds;

    border_updates.clear();
    for (auto [global_node, update] : ghost_updates) {
      // update ghost node
      auto ghost = graph.get_ghost(global_node);
      // ghost is either not set yet or was excluded due to including a
      // border node
      KASSERT(status[ghost] == wis_status::UNSET || status[ghost] == wis_status::EXCLUDED);
      // ghost is either not set or has the is also excluded at the
      // respective rank
      KASSERT(status[ghost] == wis_status::UNSET || update == status[ghost]);
      if (status[ghost] == wis_status::UNSET) {
        set_node(graph, ghost, update);
      }
    }
  }

  KASSERT(consider_next.empty(), "[greedy maximize] some non-ghosts are left in the graph.");
#if KASSERT_ENABLED(KASSERT_ASSERTION_LEVEL_NORMAL)
  graph.for_each_visible_node([&](NodeID node) { KASSERT(status[node] != wis_status::UNSET); });
#endif

  KADISREDU_RLOG << "[greedy maximize] Finished after " << rounds << " rounds (=ghosts syncs.)";
}

void sync_greedy_maximize::set_node(const dist_dynamic_graph& graph, NodeID node,
                                    wis_status new_status) {
  if (new_status == wis_status::INCLUDED) {
    include_node(graph, node);
  } else {
    KASSERT(new_status == wis_status::EXCLUDED);
    exclude_node(graph, node);
  }
}
void sync_greedy_maximize::on_modified_border_node(const dist_dynamic_graph& graph,
                                                   NodeID border_node) {
  border_updates.push_back({border_node, status[border_node]});
}
}  // namespace kadisredu::mwis::sync
