//
// Created by jannickb on 6/25/24.
//

#include "kadisredu/algorithms/mwis/sync_reducer.h"

#include "kadisredu/algorithms/mwis/debug.h"
#include "kadisredu/tools/logger.h"

#include <kamping/checking_casts.hpp>
#include <kamping/communicator.hpp>

namespace kadisredu::mwis {

void SyncReducer::sync_border_status_updates() {
  auto &detached_interface_vertices = set_1;
  detached_interface_vertices.clear();
  auto &graph = status_.graph;
  auto &solution = status_.solution;
  auto &node_status = solution.node_status;

  auto recv_border_updates = graph.sync_ghosts(border_updates, [&](NodeID node) {
    if (node_status[node] == wis_status::DETACHED) {
      if (std::all_of(graph[node].begin(), graph[node].end(),
                      [&](NodeID target) { return node_status[target] == wis_status::EXCLUDED; })) {
        node_status[node] = wis_status::INCLUDED;
        solution.solution_weight += graph.get_weight(node);
        return false;
      } else {
        detached_interface_vertices.add(node);
      }
    }
    return true;
  });
  ++status_.alltoallv;
  detached_since_last_sync = false;
  border_updates.clear();

  status_.reduce_timer->start("border updates");
  auto old_n = graph.num_visible_nodes();
  status_.recv_border_updates += recv_border_updates.size();
  for (auto [global_node, update] : recv_border_updates) {
    bool redundant = true;
    if (global_node >= graph.get_first_global_node() &&
        global_node < graph.get_first_global_node() + graph.num_nodes()) {
      // received update for detached node
      auto local_node = global_node - graph.get_first_global_node();
      KASSERT(node_status[local_node] == wis_status::DETACHED ||
                  node_status[local_node] == update.new_status,
              static_cast<int>(node_status[local_node])
                  << " " << static_cast<int>(update.new_status));
      if (node_status[local_node] == wis_status::DETACHED) {
        node_status[local_node] = update.new_status;
        if (update.new_status == wis_status::INCLUDED) {
          solution.solution_weight += graph.get_weight(local_node);
        };
        redundant = false;
      }
    } else {
      // update ghost node
      auto ghost = graph.get_ghost(global_node);
      // ghost is either not set yet or was excluded due to including a
      // border node
      KASSERT(
          node_status[ghost] == wis_status::UNSET || node_status[ghost] == wis_status::EXCLUDED,
          "The vertex received already (falsely) in the last reduction round a status: included or "
          "folded or detached.");

      if (node_status[ghost] == wis_status::UNSET) {
        // ghost is not reduced yet
        if (update.new_status == wis_status::DETACHED) {
          auto &detached_neigh = buffers[BUFFERS - 1];
          detached_neigh.clear();
          for (auto target_it = graph[ghost].end();
               graph[ghost].deg() + detached_neigh.size() < MAX_DETACH_DEGREE &&
               target_it != graph[ghost].hidden_end();
               target_it++) {
            if (detached_interface_vertices.contains(*target_it)) {
              detached_neigh.push_back(*target_it);
            }
          };
          // ghost was detached to this rank
          // reduce ghost
          for (auto &reduction : status_.reductions) {
            bool reduced = reduction->reduce_detached_ghost(this, ghost, detached_neigh);

            if (reduced) {
              break;  // jump back to first reduction
            }
          }
          KASSERT(node_status[ghost] != wis_status::UNSET || graph[ghost].deg() == 0,
                  "failed to reduce detached interface vertex\n"
                      << debug::get_neighborhood_string(solution, graph, ghost));
          redundant = false;
        } else if (update.new_status == wis_status::EXCLUDED ||
                   update.new_status == wis_status::INCLUDED) {
          set(ghost, update.new_status);
          redundant = false;
        }
      } else if (node_status[ghost] == wis_status::EXCLUDED &&
                 update.new_status == wis_status::INCLUDED) {
        // ghost was reduced and recv a conflicting status -> tie break solution
        KASSERT(std::find_if(graph[ghost].end(), graph[ghost].hidden_end(), [&](NodeID target) {
                  return node_status[target] == wis_status::INCLUDED;
                }) != graph[ghost].hidden_end());
        auto included_node = *std::find_if(
            graph[ghost].end(), graph[ghost].hidden_end(),
            [&](NodeID target) { return node_status[target] == wis_status::INCLUDED; });
        KASSERT(graph.get_weight(ghost) == graph.get_weight(included_node));
        if (!local_vertex_wins_tie(included_node, ghost)) {
          // KADISREDU_LOG << "flipped decision";
          solution.solution_weight -= graph.get_weight(included_node);
          node_status[included_node] = wis_status::EXCLUDED;
          node_status[ghost] = wis_status::INCLUDED;
        }
        redundant = false;
      } else if (update.new_status == wis_status::DETACHED) {
        KASSERT(node_status[ghost] == wis_status::EXCLUDED);
      }
    }
    if (redundant) {
      ++status_.recv_redundant_border_updates;
    }
  }
  status_.reduce_timer->stop_and_add({kamping::measurements::GlobalAggregationMode::min,
                                      kamping::measurements::GlobalAggregationMode::max});

  progress = progress || old_n - graph.num_visible_nodes() > 0;
}

std::pair<bool, bool> SyncReducer::reduce(double time_limit) {
  using namespace kamping;
  status_.reduce_timer->start("reduce-algorithm");

  KADISREDU_RLOG << "[reducer] " << "timelimit [s]: " << time_limit;
  status_.time_limit = time_limit;
  status_.t.restart();

  //GlobalNodeID global_visible_nodes;

  // mark all nodes
  status_.graph.for_each_visible_node([&](NodeID node) { mark_node(node); });

  // has something changed so that another reduction might be applicable?
  progress = false;
  detached_since_last_sync = false;
  timeout = false;

  auto all_finished = [this]() {
    // all are finished if no local reducer made progress
    ++status_.allreduce;
    return comm_.allreduce_single(
        send_buf((kabool)(!progress && border_updates.empty() && weight_modified_nodes.empty()) ||
                 time_limit_reached()),
        op(ops::logical_and<>()));
  };

  do {
    progress = false;

    status_.reduce_timer->start("reductions");

    do {
      for (std::size_t i = 0; i < status_.reductions.size(); ++i) {
        active_reduction_index = i;
        auto &reduction = *status_.reductions[i];
        // consider next marked nodes buffer
        reduction.marker.next();
        // run reductions
        {
          status_.reduce_timer->start(to_string(reduction.get_reduction_type()));
          progress = reduction.reduce(this);
          status_.reduce_timer->stop_and_add({kamping::measurements::GlobalAggregationMode::min,
                                              kamping::measurements::GlobalAggregationMode::max});
        }

        if (progress) {
          break;  // jump back to first reduction
        }
      }
      // LOCALLY EXHAUSTIVE -> reduce until no progress can be made locally anymore
    } while (progress && LOCALLY_EXHAUSTIVE);

    status_.reduce_timer->stop_and_add({kamping::measurements::GlobalAggregationMode::min,
                                        kamping::measurements::GlobalAggregationMode::max});

    check_and_sync_borders();

    if (timeout) {
      KADISREDU_RLOG << "timelimit reached";
      break;
    }

    // as long as at least one process made progress and time is left
  } while (!all_finished());

  KASSERT(border_updates.empty(), "border updates ignored");

  // clean up isolated visible ghosts
  if (status_.graph.num_visible_ghosts() > 0) {
    status_.graph.for_each_visible_ghost([&](auto visible_ghost) {
      if (status_.graph[visible_ghost].deg() == 0) {
        status_.graph.hide_node(visible_ghost);
        status_.modified_queue.push_back(visible_ghost);
      }
    });
  }
#ifndef NDEBUG
  // check ghost interface
  // all visible border nodes are sent to their adjacent visible ghost partitions
  // we expect that all received ghosts are visible border nodes at each partition
  auto &graph = status_.graph;
  sized_vector<std::pair<NodeID, GlobalNodeWeight>> updates(
      std::min(graph.get_border().size(), static_cast<std::size_t>(graph.num_visible_nodes())));
  std::for_each(
      graph.get_border().begin(), graph.get_border().end(),
      [&graph, &node_status = status_.solution.node_status, &updates](NodeID border_node) {
        if (graph.is_visible(border_node)) {
          KASSERT(node_status[border_node] == wis_status::UNSET);
          updates.push_back({border_node, graph.get_weight(border_node)});
        }
      });
  auto recv_ghosts = graph.sync_ghosts(updates);
  for (auto recv_ghost : recv_ghosts) {
    auto local_ghost = graph.get_ghost(recv_ghost.global_id);
    if (status_.solution.node_status[local_ghost] != wis_status::UNSET) {
      std::cout << static_cast<std::size_t>(status_.solution.node_status[local_ghost]) << std::endl;
    }
    KASSERT(graph.is_visible(local_ghost), "ghost is hidden");
    KASSERT(graph[local_ghost].deg() > 0, "ghost should be connected to local graph");
    KASSERT(graph.get_weight(local_ghost) >= recv_ghost.update, "ghost weight was underestimated");
  }
  KASSERT(recv_ghosts.size() == graph.num_visible_ghosts(),
          "amount of received visible ghosts does not match local visible ghosts.");
#endif
  KASSERT(border_updates.empty());

  status_.reduce_timer->stop_and_add({kamping::measurements::GlobalAggregationMode::min,
                                      kamping::measurements::GlobalAggregationMode::max});

  KADISREDU_RLOG << "[reducer] " << "Finished with " << status_.alltoallv << " ghost syncs, "
                 << status_.allreduce << " blocking update checks!";

  return std::make_pair(
      comm_.allreduce_single(send_buf(status_.graph.num_visible_nodes()), op(ops::plus<>())) == 0,
      timeout);
}

void SyncReducer::apply_and_restore() {
  auto &graph = status_.graph;
  auto &solution = status_.solution;
  auto &node_status = solution.node_status;

  while (!status_.modified_queue.empty()) {
    auto modified_node = status_.modified_queue.back();
    status_.modified_queue.pop_back();

    if (modified_node == MODIFIED_EDGE) {
      auto &modifying_reduction =
          *status_.reductions[status_.reduction_map[static_cast<std::size_t>(
              status_.folded_queue.back())]];
      modifying_reduction.apply(this, modified_node);

      status_.folded_queue.pop_back();

    } else if (modified_node == BARRIER) {
      auto updates = graph.sync_ghosts(border_updates);
      border_updates.clear();

      for (auto [global_node, update] : updates) {
        if (global_node >= graph.get_first_global_node() &&
            global_node < graph.get_first_global_node() + graph.num_nodes()) {
          // received update for detached node
          auto local_node = global_node - graph.get_first_global_node();
          KASSERT(node_status[local_node] == wis_status::DETACHED ||
                  node_status[local_node] == update.new_status);
          if (node_status[local_node] == wis_status::DETACHED) {
            node_status[local_node] = update.new_status;
            if (update.new_status == wis_status::INCLUDED) {
              solution.solution_weight += graph.get_weight(local_node);
            };
          }
        }
      }
    } else {
      // If modified_node is neither in nor out of the solution,
      //   then it is at this state an isolated ghost
      // The actual status is stored at the respective block.
      KASSERT(node_status[modified_node] != wis_status::UNSET ||
              graph[modified_node].deg() == 0 && graph.is_ghost(modified_node));

      if (node_status[modified_node] == wis_status::FOLDED) {
        // reductions made modifications that are beyond excluding and including
        // nodes
        //  -> it has its own special handling
        KASSERT(!status_.folded_queue.empty());
        KASSERT(static_cast<std::size_t>(status_.folded_queue.back()) <
                status_.reduction_map.size());

        auto &modifying_reduction =
            *status_.reductions[status_.reduction_map[static_cast<std::size_t>(
                status_.folded_queue.back())]];
        modifying_reduction.apply(this, modified_node);

        status_.folded_queue.pop_back();
      } else {
        auto restored = graph.restore_last();
        KASSERT(restored == modified_node, "restored the wrong modified node");
        if (node_status[modified_node] == wis_status::DETACHED) {
          KASSERT(std::any_of(
              graph[modified_node].begin(), graph[modified_node].end(),
              [&](NodeID target) { return node_status[target] == wis_status::INCLUDED; }));
          node_status[modified_node] = wis_status::EXCLUDED;
        }
      }
    }
  }
  KASSERT(status_.sol_offset == 0);
  KASSERT(border_updates.empty());
}

void SyncReducer::sync_border_weight_updates() {
  status_.reduce_timer->start("ghost-weight-sync");

  KASSERT(border_weight_updates.empty());

  for (auto border_node : weight_modified_nodes) {
    if (status_.solution.node_status[border_node] == wis_status::UNSET ||
        status_.solution.node_status[border_node] == wis_status::DETACHED) {
      KASSERT(status_.graph.get_weight(border_node) > 0, "The border node should be excluded.");
      border_weight_updates.push_back({border_node, status_.graph.get_weight(border_node)});
    }
  }
  weight_modified.clear();
  weight_modified_nodes.clear();

  auto ghost_weight_updates = status_.graph.sync_ghosts(border_weight_updates);
  ++status_.alltoallv;
  border_weight_updates.clear();

  auto old_n = status_.graph.num_visible_nodes();
  status_.recv_border_updates += ghost_weight_updates.size();
  for (auto [global_node, update] : ghost_weight_updates) {
    // update ghost node
    auto ghost = status_.graph.get_ghost(global_node);
    if (status_.graph.is_visible(ghost)) {
      KASSERT(update > 0);
      KASSERT(status_.graph.get_weight(ghost) > update);
      // update weight
      status_.graph.set_weight(ghost, update);
      // mark neighbors (local border vertices)
      mark_neigh(ghost);
    } else {
      ++status_.recv_redundant_border_updates;
    }
  }

  status_.reduce_timer->stop_and_add({kamping::measurements::GlobalAggregationMode::min,
                                      kamping::measurements::GlobalAggregationMode::max});

  progress = progress || old_n - status_.graph.num_visible_nodes() > 0;
}
bool SyncReducer::check_and_sync_borders() {
  using namespace kamping;
  // does any process wants to synchronize borders
  status_.reduce_timer->start("ghost-sync-phase");
  bool sync_border_now = false;
  bool barrier = false;
  std::tie(sync_border_now, barrier, timeout) = check_sync_borders();
  if (sync_border_now) {
    int old_active_reduction_index = active_reduction_index;
    active_reduction_index = status_.reductions.size();
    if (barrier) {
      status_.modified_queue.push_back(BARRIER);
    }
    sync_border_weight_updates();
    sync_border_status_updates();
    active_reduction_index = old_active_reduction_index;
    //auto global_visible_nodes =
    //    comm_.reduce_single(send_buf(status_.graph.num_visible_nodes()), op(ops::plus<>()));
    //	KADISREDU_RLOG << "[reducer] "
    //               << "synced borders; visible nodes: " << global_visible_nodes.value();
  }
  status_.reduce_timer->stop_and_add({kamping::measurements::GlobalAggregationMode::min,
                                      kamping::measurements::GlobalAggregationMode::max});

  return sync_border_now;
}

void SyncReducer::on_border_status_change(NodeID border_node, wis_status new_status) {
  detached_since_last_sync = detached_since_last_sync || new_status == wis_status::DETACHED;
  border_updates.push_back({border_node, {new_status}});
};
void SyncReducer::on_border_weight_shift(NodeID border_node, NodeWeight new_weight) {
  if (weight_modified.add(border_node)) {
    weight_modified_nodes.push_back(border_node);
  }
};

std::tuple<bool, bool, bool> SyncReducer::check_sync_borders() {
  using namespace kamping;
  struct reducer_state_message {
    int no_progress_pes_counter = 0;
    bool barrier = false;
    bool timeout = false;
  };
  auto agg_state = comm_.allreduce_single(
      // vote for synchronisation if local reducer is out of work
      send_buf(
          reducer_state_message{(int)(!progress), detached_since_last_sync, time_limit_reached()}),
      op(
          [&](auto &p1, auto &p2) {
            return reducer_state_message{p1.no_progress_pes_counter + p2.no_progress_pes_counter,
                                         p1.barrier || p2.barrier, p1.timeout || p2.timeout};
          },
          ops::commutative));
  ++status_.allreduce;
  return std::make_tuple(agg_state.timeout || LOCALLY_EXHAUSTIVE ||
                             agg_state.no_progress_pes_counter >= comm_.size() / 2,
                         agg_state.barrier, agg_state.timeout);
};

}  // namespace kadisredu::mwis
