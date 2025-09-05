#include "kadisredu/algorithms/mwis/peeler.h"

#include "branch_and_reduce_algorithm.h"
#include "kadisredu/algorithms/mwis/mwis_def.h"
#include "kadisredu/algorithms/mwis/reducer.h"
#include "kadisredu/data_structures/maxNodeHeap.h"
#include "kadisredu/tools/kamis_helper.h"
#include "kadisredu/tools/logger.h"
#include "kassert/kassert.hpp"
#include "message-queue/buffered_queue.hpp"

#include <kamping/measurements/measurement_aggregation_definitions.hpp>

namespace kadisredu::mwis {
bool peeler::reduce(weighted_reduce_ptr algo) {
  auto &status = algo->status_;
  auto &graph = status.graph;
  auto &solution = status.solution;
  auto &node_status = solution.node_status;

  auto get_rating = [&graph](NodeID node) -> NodeWeight {
    // KASSERT(graph.get_weight_neigh(node) - graph.get_weight(node) > 0, "missed reduction");
    /*KASSERT((graph.get_weight_neigh(node) - graph.get_weight(node)) <=
                (graph.get_weight_neigh(node) - graph.get_weight(node)) *
                    static_cast<NodeWeight>(graph[node].ghost_deg() + 1),
            "overflow");*/
    return (graph.get_weight_neigh(node) - graph.get_weight(node)); /* *
           static_cast<NodeWeight>(graph[node].ghost_deg() + 1);*/
  };

  /*auto get_rating = [&graph](NodeID node) -> NodeWeight {
    KASSERT(graph.get_weight_neigh(node) - graph.get_weight(node) > 0, "missed reduction");
    NodeWeight neigh = 0;
    for (auto target : graph[node]) {
      if (graph.is_ghost(target)) {
        neigh += graph.get_weight(target) / 2;
      } else {
        neigh += graph.get_weight(target);
      }
    }

    return (neigh - graph.get_weight(node)) * static_cast<NodeWeight>(graph[node].ghost_deg() + 1);
  };*/

  /*auto get_rating = [&graph](NodeID node) -> NodeWeight { return graph[node].ghost_deg(); };*/

  // do not poll before initialization of pq_; otherwise we loose the marked nodes if the poll is
  // successful
  if (graph.num_visible_nodes() > 0) {
    if (pq_.empty()) {
      pq_.reserve(graph.num_visible_nodes());
      for_each_unset_marked_node(status.solution, marker,
                                 [&](GlobalNodeID node) { pq_.insert(node, get_rating(node)); });
      KASSERT(pq_.size() == graph.num_visible_nodes());
    } else {
      for_each_unset_marked_node(status.solution, marker, [&](GlobalNodeID node) {
        KASSERT(pq_.contains(node));
        pq_.changeKey(node, get_rating(node));
      });
    }
  }

  const bool sync_mode =
      graph.num_visible_nodes() == 0 || !init || peeled == static_cast<NodeID>(sync_threshold);
  // const bool sync_mode = graph.num_visible_nodes() == 0 || !init;
  if (!init) {
    init = true;
  }
  if (!sync_mode) {
    // try to circumvent a peeling step
    if (algo->poll_updates_during_reduce_phase()) {
      return true;
    }
  } else {
    //  Peel only if global progress can be made
    if (!algo->exchange_updates()) {
      // Check for local reductions applications before peeling
      return true;
    } else {
      // No global progress was made
      using namespace kamping;
      GlobalNodeID global_n =
          algo->comm_.allreduce_single(send_buf(graph.num_visible_nodes()), op(ops::plus<>()));

      if (global_n > 0) {
        sync_threshold = 1.1 * sync_threshold;
        peeled = 0;

        algo->reactivate_reduction_phase();
      } else {
        return false;
      }
      if (graph.num_visible_nodes() == 0) {
        // local rank finished, but it needs to wait for the others;
        return true;
      }
    }
  }

  constexpr NodeID solve_subproblem_threshold = 0;
  // status.reduce_timer->start("subproblem");
  if (graph.num_visible_nodes() < solve_subproblem_threshold) {
    bool fail = false;

    if (graph.num_ghosts() != 0) {
      for (NodeID node = 0; node < graph.num_nodes(); ++node) {
        if (graph.is_visible(node) && graph[node].ghost_deg() > 0) {
          fail = true;
          break;
        }
      }
    }
    if (!fail) {
      KADISREDU_LOG << "solving subproblem of size " << graph.num_visible_nodes();
      auto &nodes_set = algo->set_1;
      auto &nodes = algo->buffers[0];
      auto &map_to_local = algo->buffers[1];
      graph_access local_graph;

      nodes.clear();
      nodes_set.clear();
      graph.for_each_visible_node([&](NodeID node) {
        nodes.push_back(node);
        nodes_set.add(node);
      });

      kamis_helper::build_induced_graph(nodes_set, nodes, graph, map_to_local, local_graph);
      MISConfig cfg = kamis_helper::create_mis_config();
      cfg.time_limit = static_cast<double>(solve_subproblem_threshold) / 10.0;

      branch_and_reduce_algorithm solver(local_graph, cfg, true);
      solver.run_branch_reduce();

      solver.apply_branch_reduce_solution(local_graph);

      for (NodeID node = 0; node < local_graph.number_of_nodes(); ++node) {
        if (local_graph.getPartitionIndex(node) == 1) {
          algo->set(nodes[node], wis_status::INCLUDED);
        }
      }
      KASSERT(graph.num_visible_nodes() == 0, "solution by KaMIS was not maximal");

      return true;
    }
  }
  /*status.reduce_timer->stop_and_add({kamping::measurements::GlobalAggregationMode::min,
                                     kamping::measurements::GlobalAggregationMode::max});*/

  KASSERT(!pq_.empty());
  while (!graph.is_visible(pq_.maxElement())) {
    pq_.deleteMax();
  }

  if (algo->poll_updates_during_reduce_phase()) {
    return true;
  }

  KASSERT(!pq_.empty());
  algo->set(pq_.maxElement(), wis_status::EXCLUDED);
  ++peeled;

  pq_.deleteMax();

  return true;
}
}  // namespace kadisredu::mwis
