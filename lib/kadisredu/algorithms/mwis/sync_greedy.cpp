//
// Created by jannickb on 6/18/24.
//

#include "kadisredu/algorithms/mwis/sync_greedy.h"

#include "kadisredu/tools/logger.h"

namespace kadisredu::mwis::sync {
void sync_greedy::run(dist_dynamic_graph& graph) {
  using namespace kamping;

  KASSERT(status.size() >= graph.num_nodes() + graph.num_ghosts());
  KASSERT(wis_weight == 0);

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
  auto n = comm.allreduce_single(send_buf(graph.num_visible_nodes()), op(ops::plus<>()));
  while (n > 0) {
    while (!candidates.next_marked_nodes.empty()) {
      candidates.next();
      for (auto node : candidates.current()) {
        KASSERT(!graph.is_ghost(node));
        if (check_for_include(graph, node)) {
          // node has maximal rating in remaining neighborhood
          set_node(graph, node, wis_status::INCLUDED);
        }
      };
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
    // is still work to do?
    GlobalNodeID n_new =
        comm.allreduce_single(send_buf(graph.num_visible_nodes()), op(ops::plus<>()));

    // KASSERT(!comm.is_root() || n_new < n, "[greedy]: no global progress was made");
    n = n_new;
  }

  KASSERT(graph.num_visible_nodes() == 0, "[greedy] some non-ghosts are left in the graph.");

  KADISREDU_RLOG << "[greedy] Finished after " << rounds << " rounds (=ghosts syncs.)";
}
void sync_greedy::set_node(dist_dynamic_graph& graph, NodeID node, wis_status new_status) {
  if (new_status == wis_status::INCLUDED) {
    include_node(graph, node);
  } else {
    KASSERT(new_status == wis_status::EXCLUDED);
    exclude_node(graph, node);
  }
}
void sync_greedy::on_modified_border_node(dist_dynamic_graph& graph, NodeID border_node) {
  border_updates.push_back({border_node, status[border_node]});
}
}  // namespace kadisredu::mwis::sync
