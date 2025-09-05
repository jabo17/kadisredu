//
// Created by jannickb on 6/18/24.
//

#include "kadisredu/algorithms/mwis/async_greedy_maximize.h"

#include "kadisredu/tools/logger.h"

namespace kadisredu::mwis::async {
void async_greedy_maximize::run(const dist_dynamic_graph& graph, const Solution& solution) {
  using namespace kamping;

  KASSERT(solution.node_status.size() == graph.num_nodes() + graph.num_ghosts());
  status = solution.node_status;
  wis_weight = solution.solution_weight;

  // sync ghost status to obtain included ghosts
  sync_ghost_status(graph);
  KADISREDU_RLOG << "synced ghosts";

  // sync free nodes at the border
  init_free_nodes(graph);
  KADISREDU_RLOG << "init free nodes";

  message_queue::FlushStrategy flush_strategy = message_queue::FlushStrategy::local;
  border_update_queue.local_threshold_bytes(message_queue_local_threshold_);
  border_update_queue.flush_strategy(flush_strategy);

  // add nodes greedily to solution
  bool terminated = false;
  while (!terminated) {
    while (!consider_next.empty()) {
      std::swap(consider, consider_next);
      consider_next_set.clear();
      for (NodeID node : consider) {
        KASSERT(!graph.is_ghost(node));
        if (check_for_include(graph, node)) {
          // node has maximal rating in remaining neighborhood
          set_node(graph, node, wis_status::INCLUDED);
        }
      };
      consider.clear();

      // poll incoming border updates
      poll_border_updates(graph);
    }

    // no progress; check for termination
    terminated = terminate_border_update_queue(graph);
    KASSERT(!terminated || consider_next.empty(),
            "termination must imply that there are no free nodes left");
    ++rounds;  // count tries for termination
  }

  KASSERT(consider_next.empty(), "[greedy maximize] some non-ghosts are left in the graph.");
#if KASSERT_ENABLED(KASSERT_ASSERTION_LEVEL_NORMAL)
  graph.for_each_visible_node([&](NodeID node) { KASSERT(status[node] != wis_status::UNSET); });
#endif

  KADISREDU_RLOG << "[greedy maximize] Finished after " << rounds << " rounds (=ghosts syncs.)";
}

void async_greedy_maximize::set_node(const dist_dynamic_graph& graph, NodeID node,
                                     wis_status new_status, bool poll) {
  if (new_status == wis_status::INCLUDED) {
    include_node(graph, node);
  } else {
    KASSERT(new_status == wis_status::EXCLUDED);
    exclude_node(graph, node);
  }
  if (poll) {
    poll_border_updates(graph);
  }
}
void async_greedy_maximize::on_modified_border_node(const dist_dynamic_graph& graph,
                                                    NodeID border_node) {
  auto global_node_id = graph.get_global_node_id(border_node);
  reached_ranks.clear();

  auto border_node_status = status[border_node];

  for (auto& target : graph[border_node]) {
    if (graph.is_ghost(target) && status[target] == wis_status::UNSET) {
      auto rank = static_cast<message_queue::PEID>(graph.get_ghost_block(target));
      if (reached_ranks.add(rank)) {
        border_update_queue.post_message({.node = global_node_id, .status = border_node_status},
                                         rank);
      }
    }
  }
}
bool async_greedy_maximize::terminate_border_update_queue(const dist_dynamic_graph& graph) {
  return border_update_queue.terminate(
      [&](message_queue::Envelope<BorderUpdateMessage> auto envelope) {
        for (auto [global_node, update] : envelope.message) {
          // update ghost node
          auto ghost = graph.get_ghost(global_node);
          // ghost is either not set yet or was excluded due to including a
          // border node
          KASSERT(status[ghost] == wis_status::UNSET || status[ghost] == wis_status::EXCLUDED);
          // ghost is either not set or has the is also excluded at the
          // respective rank
          KASSERT(status[ghost] == wis_status::UNSET || update == status[ghost]);
          if (status[ghost] == wis_status::UNSET) {
            set_node(graph, ghost, update, false);
          }
        }
        if (!consider_next.empty()) {
          border_update_queue.reactivate();
        }
      });
}
void async_greedy_maximize::poll_border_updates(const dist_dynamic_graph& graph) {
  border_update_queue.poll([&](message_queue::Envelope<BorderUpdateMessage> auto envelope) {
    for (auto [global_node, update] : envelope.message) {
      // update ghost node
      auto ghost = graph.get_ghost(global_node);
      // ghost is either not set yet or was excluded due to including a
      // border node
      KASSERT(status[ghost] == wis_status::UNSET || status[ghost] == wis_status::EXCLUDED);
      // ghost is either not set or has the is also excluded at the
      // respective rank
      KASSERT(status[ghost] == wis_status::UNSET || update == status[ghost]);
      if (status[ghost] == wis_status::UNSET) {
        set_node(graph, ghost, update, false);
      }
    }
  });
}
}  // namespace kadisredu::mwis::async
