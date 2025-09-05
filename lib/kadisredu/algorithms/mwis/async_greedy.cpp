//
// Created by jannickb on 6/18/24.
//

#include "kadisredu/algorithms/mwis/async_greedy.h"

#include "kadisredu/tools/logger.h"

namespace kadisredu::mwis::async {
void async_greedy::run(dist_dynamic_graph& graph) {
  using namespace kamping;

  KASSERT(status.size() >= graph.num_nodes() + graph.num_ghosts());
  KASSERT(wis_weight == 0);

  message_queue::FlushStrategy flush_strategy = message_queue::FlushStrategy::local;
  border_update_queue.local_threshold_bytes(message_queue_local_threshold_);
  border_update_queue.flush_strategy(flush_strategy);

  sized_vector<std::pair<NodeID, GlobalNodeWeight>> rated_border_nodes(graph.get_border().size());

  static_assert(std::numeric_limits<GlobalNodeWeight>::is_signed);
  graph.for_each_visible_node([&](NodeID visible_node) {
    // add rating
    ratings[visible_node] = graph.get_weight(visible_node) - graph.get_weight_neigh(visible_node);

    if (graph.is_original_border_node(visible_node)) {
      // visible node is broder node, share rating
      rated_border_nodes.push_back({visible_node, ratings[visible_node]});
    }
  });

  auto rated_ghosts = graph.sync_ghosts(rated_border_nodes);

  for (auto [global_node, rating] : rated_ghosts) {
    ratings[graph.get_ghost(global_node)] = rating;
  }

  // init round
  graph.for_each_visible_node([&](NodeID node) { candidates.mark(node); });

  // add nodes greedily to solution
  bool terminated = false;
  while (!terminated) {
    while (!candidates.next_marked_nodes.empty()) {
      candidates.next();
      for (auto node : candidates.current()) {
        KASSERT(!graph.is_ghost(node));
        if (check_for_include(graph, node)) {
          // node has maximal rating in remaining neighborhood
          set_node(graph, node, wis_status::INCLUDED);
        }
      };

      // poll incoming border updates
      poll_border_updates(graph);
    }

    // no progress; check for termination
    terminated = terminate_border_update_queue(graph);
    KASSERT(!terminated || graph.num_visible_nodes() == 0,
            "termination must imply that the local graph was fully processed by greedy algorithm");
    ++rounds;  // count tries for termination
  }

  KASSERT(graph.num_visible_nodes() == 0, "[greedy] some non-ghosts are left in the graph.");

  KADISREDU_RLOG << "[greedy] Finished after " << rounds << " rounds (=max rounds.)";
}
void async_greedy::set_node(dist_dynamic_graph& graph, NodeID node, wis_status new_status,
                            bool poll) {
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
void async_greedy::on_modified_border_node(dist_dynamic_graph& graph, NodeID border_node) {
  auto global_node_id = graph.get_global_node_id(border_node);
  reached_ranks.clear();

  auto border_node_status = status[border_node];

  for (auto& target : graph[border_node]) {
    if (graph.is_ghost(target)) {
      auto rank = static_cast<message_queue::PEID>(graph.get_ghost_block(target));
      if (reached_ranks.add(rank)) {
        border_update_queue.post_message({.node = global_node_id, .status = border_node_status},
                                         rank);
      }
    }
  }
}
bool async_greedy::terminate_border_update_queue(dist_dynamic_graph& graph) {
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
        if (!candidates.next_marked_nodes.empty()) {
          border_update_queue.reactivate();
        }
      });
}
void async_greedy::poll_border_updates(dist_dynamic_graph& graph) {
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
