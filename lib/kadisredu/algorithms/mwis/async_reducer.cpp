//
// Created by jannickb on 6/25/24.
//

#include "kadisredu/algorithms/mwis/async_reducer.h"

#include "kadisredu/algorithms/mwis/debug.h"
#include "kadisredu/algorithms/mwis/mwis_def.h"
#include "kadisredu/tools/logger.h"
#include "kassert/kassert.hpp"

#include <algorithm>
#include <kamping/checking_casts.hpp>
#include <kamping/communicator.hpp>
#include <message-queue/buffered_queue.hpp>

namespace kadisredu::mwis {

std::pair<bool, bool> AsyncReducer::reduce(double time_limit) {
  using namespace kamping;
  status_.reduce_timer->start("reduce-algorithm");
  KADISREDU_RLOG << "[reducer] " << "timelimit [s]: " << time_limit;
  status_.time_limit = time_limit;
  status_.t.restart();

  // has something changed so that another reduction might be applicable?
  bool progress = false;

  // configure message-queue
  message_queue::FlushStrategy flush_strategy = message_queue::FlushStrategy::local;
  const size_t local_threshold = status_.cfg.message_queue_local_threshold;
  border_update_queue.local_threshold_bytes(local_threshold);
  border_update_queue.flush_strategy(flush_strategy);

  // mark all nodes
  status_.graph.for_each_visible_node([&](NodeID node) { mark_node(node); });

  do {
    progress = false;

    for (std::size_t i = 0; i < status_.reductions.size(); ++i) {
      active_reduction_index = i;
      if (time_limit_reached()) {
        KADISREDU_RLOG << "[reducer]" << "timeout" << time_limit;
        // set progress to true, so that queue will terminate
        progress = false;
        break;
      }

      auto &reduction = *status_.reductions[i];
      // consider next marked nodes buffer
      reduction.marker.next();
      // run reductions
      {
        status_.reduce_timer->start(kadisredu::mwis::to_string(reduction.get_reduction_type()));
        progress = reduction.reduce(this);
        status_.reduce_timer->stop_and_add({kamping::measurements::GlobalAggregationMode::max});
      }
      progress = poll_updates_during_reduce_phase() || progress;  // @todo capture progress
      if (progress) break;                                        // jump back to first reduction
    }

    if (!progress) {
      auto finished_globally = exchange_updates();
      if (finished_globally) {
        break;  // time limit reached or no global progress
      }
    }
  } while (true);

  KADISREDU_RLOG << "finished with reductions";

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

  status_.reduce_timer->stop_and_add({kamping::measurements::GlobalAggregationMode::min,
                                      kamping::measurements::GlobalAggregationMode::max});

  KADISREDU_RLOG << "[reducer] " << "Finished!";

  return std::make_pair(
      comm_.allreduce_single(send_buf(status_.graph.num_visible_nodes() == 0),
                             op(ops::logical_and<>())),
      comm_.allreduce_single(send_buf(time_limit_reached()), op(ops::logical_or<>())));
}

bool AsyncReducer::exchange_updates() {
  bool progress = false;
  return border_update_queue.terminate(
      [&](message_queue::Envelope<AsyncReducer::BorderUpdateMessage> auto envelope) {
        progress = on_recv_border_update(std::move(envelope));
      },
      [&]() {
        if (progress && !time_limit_reached()) {
          border_update_queue.reactivate();
        }
      });
}

void AsyncReducer::apply_and_restore() {
  auto &graph = status_.graph;
  auto &solution = status_.solution;
  auto &node_status = solution.node_status;
#ifndef NDEBUG
  graph.for_each_visible_node_inc_ghosts([&](GlobalNodeID node) {
    auto s = node_status[node];
    KASSERT(graph.is_ghost(node) || s == wis_status::INCLUDED || s == wis_status::EXCLUDED,
            static_cast<int>(s));
    KASSERT(!graph.is_ghost(node) || s == wis_status::INCLUDED || s == wis_status::EXCLUDED ||
                s == wis_status::UNSET,
            static_cast<int>(s));
  });
#endif

  this->border_update_queue.reactivate();

  while (!status_.modified_queue.empty()) {
    auto modified_node = status_.modified_queue.back();
    status_.modified_queue.pop_back();

    if (modified_node == MODIFIED_EDGE) {
      auto &modifying_reduction =
          *status_.reductions[status_.reduction_map[static_cast<std::size_t>(
              status_.folded_queue.back())]];
      modifying_reduction.apply(this, modified_node);

      status_.folded_queue.pop_back();

    } else if (node_status[modified_node] == wis_status::FOLDED) {
      // reductions made modifications that are beyond excluding and including
      // nodes
      //  -> it has its own special handling
      KASSERT(!status_.folded_queue.empty());
      KASSERT(static_cast<std::size_t>(status_.folded_queue.back()) < status_.reduction_map.size());

      // ensure that neighbors of modified_node get attached before unfolding
      for (std::size_t i = 0; i < graph[modified_node].deg(); i++) {
        auto target = graph[modified_node][i];
        if (node_status[target] == wis_status::DETACHED) {
          auto on_message =
              [&](message_queue::Envelope<AsyncReducer::BorderUpdateMessage> auto &&envelope) {
                for (auto [global_node, new_status, new_weight
#ifndef NDEBUG
                           ,
                           debug
#endif
                ] : envelope.message) {
                  if (global_node >= graph.get_first_global_node() &&
                      global_node < graph.get_first_global_node() + graph.num_nodes()) {
                    // received update for detached node
                    auto local_node = global_node - graph.get_first_global_node();
                    KASSERT((node_status[local_node] == wis_status::DETACHED &&
                             (new_status == wis_status::INCLUDED ||
                              new_status == wis_status::EXCLUDED)) ||
                                node_status[local_node] == new_status,
                            std::to_string(static_cast<int>(node_status[local_node])) + " " +
                                std::to_string(static_cast<int>(new_status)));
                    if (node_status[local_node] == wis_status::DETACHED) {
                      node_status[local_node] = new_status;
                      if (new_status == wis_status::INCLUDED) {
                        solution.solution_weight += graph.get_weight(local_node);
                      };
                    }
                    if (local_node == target) {
                      border_update_queue.reactivate();
                    }
                  } else {
                    message_queue::atomic_debug(fmt::format("Interesting {} {} {}", global_node,
                                                            static_cast<int>(new_status),
                                                            new_weight));
                  }
                }
              };
          border_update_queue.terminate(on_message);
          border_update_queue.reactivate();
          KASSERT(node_status[target] == wis_status::EXCLUDED ||
                  node_status[target] == wis_status::INCLUDED);
        }
      }

      auto &modifying_reduction =
          *status_.reductions[status_.reduction_map[static_cast<std::size_t>(
              status_.folded_queue.back())]];
      modifying_reduction.apply(this, modified_node);

      status_.folded_queue.pop_back();
    } else {
      auto restored = graph.restore_last();
      KASSERT(restored == modified_node, "restored the wrong modified node");
      if (MAX_DETACH_DEGREE == 1 && node_status[modified_node] == wis_status::DETACHED) {
        auto target = graph[modified_node][0];
        if (node_status[target] == wis_status::INCLUDED) {
          node_status[modified_node] = wis_status::EXCLUDED;
        }
      }
    }
  }

  this->border_update_queue.terminate(
      [&](message_queue::Envelope<AsyncReducer::BorderUpdateMessage> auto &&envelope) {
        for (auto [global_node, new_status, new_weight
#ifndef NDEBUG
                   ,
                   debug
#endif
        ] : envelope.message) {
          if (global_node >= graph.get_first_global_node() &&
              global_node < graph.get_first_global_node() + graph.num_nodes()) {
            // received update for detached node
            auto local_node = global_node - graph.get_first_global_node();
            KASSERT(node_status[local_node] == wis_status::DETACHED ||
                    node_status[local_node] == new_status);
            if (node_status[local_node] == wis_status::DETACHED) {
              node_status[local_node] = new_status;
              if (new_status == wis_status::INCLUDED) {
                solution.solution_weight += graph.get_weight(local_node);
              };
            }
          }
        }
      });
}

void AsyncReducer::on_border_status_change(NodeID border_node, wis_status new_status) {
  auto &graph = status_.graph;
  auto global_node_id = status_.graph.get_global_node_id(border_node);
  reached_ranks.clear();

  KASSERT(new_status == status_.solution.node_status[border_node]);
  auto border_node_weight = graph.get_weight(border_node);
  if (new_status == wis_status::DETACHED) {
    border_node_weight = static_cast<NodeWeight>(graph.get_global_node_id(graph[border_node][0]));
  }

  if (graph.is_ghost(border_node)) {
    border_update_queue.post_message(
        {.node = global_node_id, .status = new_status, .weight = border_node_weight},
        graph.get_ghost_block(border_node));
  } else {
    for (auto &target : graph[border_node]) {
      if (graph.is_ghost(target)) {
        auto rank = static_cast<message_queue::PEID>(graph.get_ghost_block(target));
        if (reached_ranks.add(rank)) {
          border_update_queue.post_message(
              {.node = global_node_id, .status = new_status, .weight = border_node_weight}, rank);
        }
      }
    }
  }
}

void AsyncReducer::on_border_weight_shift(NodeID border_node, NodeWeight new_weight) {
  auto &graph = status_.graph;
  auto global_node_id = status_.graph.get_global_node_id(border_node);
  reached_ranks.clear();

  KASSERT(wis_status::UNSET == status_.solution.node_status[border_node]);
  KASSERT(new_weight == graph.get_weight(border_node));

  for (auto &target : graph[border_node]) {
    if (graph.is_ghost(target)) {
      auto rank = static_cast<message_queue::PEID>(graph.get_ghost_block(target));
      if (reached_ranks.add(rank)) {
        border_update_queue.post_message(
            {.node = global_node_id, .status = wis_status::UNSET, .weight = new_weight}, rank);
      }
    }
  }
}

bool AsyncReducer::on_recv_border_update(
    message_queue::Envelope<AsyncReducer::BorderUpdateMessage> auto envelope) {
  auto &graph = status_.graph;
  auto &solution = status_.solution;
  auto &node_status = solution.node_status;
#if KASSERT_ENABLED(KASSERT_ASSERTION_LEVEL_NORMAL)
  KASSERT(!is_polling_border_updates);
  is_polling_border_updates = true;
#endif
  int old_active_reduction_index = active_reduction_index;
  active_reduction_index = status_.reductions.size();
  bool progress = false;
  // message_queue::atomic_debug(fmt::format("Received message from R{}", envelope.sender));
  for (auto [global_node, new_status, new_weight
#ifndef NDEBUG
             ,
             debug
#endif
  ] : envelope.message) {
    ++status_.recv_border_updates;
    bool redundant = true;
    if (global_node >= graph.get_first_global_node() &&
        global_node < graph.get_first_global_node() + graph.num_nodes()) {
      // received update for detached node
      auto local_node = global_node - graph.get_first_global_node();
      KASSERT(node_status[local_node] == wis_status::DETACHED ||
              node_status[local_node] == new_status);
      if (node_status[local_node] == wis_status::DETACHED) {
        node_status[local_node] = new_status;
        if (new_status == wis_status::INCLUDED) {
          solution.solution_weight += graph.get_weight(local_node);
        };
        redundant = false;
      }
    } else {
      // update ghost node
      auto ghost = graph.get_ghost(global_node);
      // ghost is either not set yet or was excluded due to including a border node
      KASSERT(
          node_status[ghost] == wis_status::UNSET || node_status[ghost] == wis_status::EXCLUDED,
          "The vertex received already (falsely) in the last reduction round a status: included or "
          "folded or detached."
              << global_node << " " << static_cast<int>(node_status[ghost]) << " "
              << static_cast<int>(new_status));

      if (node_status[ghost] == wis_status::UNSET) {
        if (new_status == wis_status::DETACHED) {
          progress = progress || graph[ghost].deg() > 0;
          KASSERT(MAX_DETACH_DEGREE == 1, "higher detach degree not supported yet");
          auto detached_neighbor =
              static_cast<NodeID>(new_weight) - status_.graph.get_first_global_node();
          auto &detached_neigh = buffers[BUFFERS - 1];
          detached_neigh.clear();
          if (node_status[detached_neighbor] == wis_status::DETACHED) {
            KASSERT(status_.graph[detached_neighbor][0] == ghost);
            detached_neigh.push_back(detached_neighbor);
          }
          // reduce ghost
          for (auto &reduction : status_.reductions) {
            bool reduced = reduction->reduce_detached_ghost(this, ghost, detached_neigh);

            if (reduced) {
              break;  // jump back to first reduction
            }
          }
          KASSERT(status_.solution.node_status[ghost] != wis_status::UNSET ||
                      status_.graph[ghost].deg() == 0,
                  "failed to reduce detached interface vertex\n"
                      << debug::get_neighborhood_string(solution, graph, ghost));
          redundant = false;
        } else if (new_status == wis_status::EXCLUDED || new_status == wis_status::INCLUDED) {
          progress = true;
          graph.set_weight(ghost, new_weight);
          set(ghost, new_status, false);
          redundant = false;
        } else {
          progress = true;
          set_new_weight(ghost, new_weight);
          redundant = false;
        }
      } else if (node_status[ghost] == wis_status::EXCLUDED && new_status == wis_status::INCLUDED) {
        auto included_node_it = std::find_if(
            graph[ghost].end(), graph[ghost].hidden_end(),
            [&](NodeID target) { return node_status[target] == wis_status::INCLUDED; });
        KASSERT(included_node_it != graph[ghost].hidden_end());
        KASSERT(graph.get_weight(ghost) == graph.get_weight(*included_node_it));
        if (!local_vertex_wins_tie(*included_node_it, ghost)) {
          solution.solution_weight -= graph.get_weight(*included_node_it);
          node_status[*included_node_it] = wis_status::EXCLUDED;
          node_status[ghost] = wis_status::INCLUDED;
        }
        redundant = false;
      } else if (MAX_DETACH_DEGREE > 1 && node_status[ghost] == wis_status::EXCLUDED &&
                 new_status == wis_status::DETACHED) {
        // Ghost was detached to this block before it had a chance to recv the implicit exclude
        // In the special case were the max. detach degree = 1, we do not need to send now an extra
        // exclude message. This is because, the other rank will receive the include message of
        // ghost's neighbor and can reconstruct the solution with the included neighbor.
        on_border_status_change(ghost, wis_status::EXCLUDED);
      }
    }
    if (redundant) {
      ++status_.recv_redundant_border_updates;
    }
  }
  active_reduction_index = old_active_reduction_index;
#if KASSERT_ENABLED(KASSERT_ASSERTION_LEVEL_NORMAL)
  is_polling_border_updates = false;
#endif
  return progress;
};

bool AsyncReducer::poll_updates_during_reduce_phase() {
  bool progress = false;
  border_update_queue.poll(
      [&](message_queue::Envelope<AsyncReducer::BorderUpdateMessage> auto envelope) {
        progress |= on_recv_border_update(std::move(envelope));
      });
  return progress;
}

void AsyncReducer::reactivate_reduction_phase() { border_update_queue.reactivate(); };

void AsyncReducer::printing_cleaner(BorderUpdateBuffer &buf, message_queue::PEID receiver) {
  auto &graph = status_.graph;
  auto &solution = status_.solution;
  auto &node_status = solution.node_status;
  auto &revoked = set_2;  // ensure not mixing up with reductions and there usage of set_2

  revoked.clear();
  std::for_each(buf.begin(), buf.end(), [&](auto &msg) {
    // filter out detached vertices that can now be included
    if (msg.node >= graph.get_first_global_node() &&
        msg.node < graph.get_first_global_node() + graph.num_nodes()) {
      auto local_node = msg.node - graph.get_first_global_node();
      if (msg.status == wis_status::DETACHED && node_status[local_node] == wis_status::DETACHED) {
        // weight message or detach message
        // try to revoke detach status
        if (std::all_of(graph[local_node].begin(), graph[local_node].end(), [&](NodeID target) {
              return node_status[target] == wis_status::EXCLUDED;
            })) {
          node_status[local_node] = wis_status::INCLUDED;
          solution.solution_weight += graph.get_weight(local_node);
          revoked.add(local_node);
        }
      }
    }
  });
  auto new_end = std::remove_if(buf.begin(), buf.end(), [&](auto &msg) {
    if (msg.node >= graph.get_first_global_node() &&
        msg.node < graph.get_first_global_node() + graph.num_nodes()) {
      auto local_node = msg.node - graph.get_first_global_node();

      // filter out revoked detach messages and previous weight updates
      if ((msg.status == wis_status::UNSET || msg.status == wis_status::DETACHED) &&
          revoked.contains(local_node)) {
        return true;
      }

      // filter out outdated weight message
      if (msg.status == wis_status::UNSET && msg.weight > graph.get_weight(local_node)) {
        // weight update message
        // but node is going to receive exclude (with weight), include (with weight), or smaller
        // weight
        return true;
      }
    }
    return false;
  });
  buf.erase(new_end, buf.end());
}
}  // namespace kadisredu::mwis
