//
// Created by jannickb on 6/18/24.
//

#pragma once

#include "kadisredu/algorithms/mwis/mwis_def.h"
#include "kadisredu/data_structures/dist_dynamic_graph.h"
#include "kadisredu/data_structures/fast_set.h"
#include "kadisredu/definitions.h"

#include <algorithm>

namespace kadisredu::mwis {

class greedy_maximize {
 public:
  greedy_maximize(MPI_Comm comm, NodeID num_nodes, NodeID num_ghosts)
      : comm(comm),
        consider_next_set(num_nodes),
        status(num_nodes + num_ghosts, wis_status::UNSET),
        ratings(num_nodes + num_ghosts, std::numeric_limits<NodeWeight>::min()),
        wis_weight(0) {}

  virtual ~greedy_maximize() = default;

  virtual void run(const dist_dynamic_graph& graph, const Solution& solution) = 0;

  Solution take_solution() { return {wis_weight, std::move(status)}; }

  [[nodiscard]] std::size_t get_rounds() const { return rounds; };

  void check_local_maximal_independent_set(const dist_dynamic_graph& graph) {
    KASSERT(status.size() == graph.num_nodes() + graph.num_ghosts());
    for (NodeID node = 0; node < graph.num_nodes(); ++node) {
      if (status[node] == wis_status::INCLUDED) {
#if KASSERT_ENABLED(KASSERT_ASSERTION_LEVEL_NORMAL)
        for (auto target = graph[node].begin(); target != graph[node].hidden_end(); ++target) {
          KASSERT(status[*target] != wis_status::INCLUDED, "greedy solution is not independent");
        }
#endif
      } else {
        bool free = true;
        for (auto target = graph[node].begin(); target != graph[node].hidden_end(); ++target) {
          if (status[*target] == wis_status::INCLUDED) {
            free = false;
            break;
          }
        }
        KASSERT(!free, "greedy solution locally not maximal");
      }
    }
    for (NodeID node = 0; node < graph.num_nodes(); ++node) {
      KASSERT(status[node] == wis_status::INCLUDED || status[node] == wis_status::EXCLUDED);
    }
    for (NodeID node = graph.num_nodes(); node < graph.num_nodes() + graph.num_ghosts(); ++node) {
      KASSERT(status[node] == wis_status::INCLUDED || status[node] == wis_status::EXCLUDED ||
              status[node] == wis_status::UNSET);
    }
  }

 protected:
  virtual void sync_ghost_status(const dist_dynamic_graph& graph);
  virtual void init_free_nodes(const dist_dynamic_graph& graph);

  void exclude_node(const dist_dynamic_graph& graph, NodeID node);

  void include_node(const dist_dynamic_graph& graph, NodeID node);

  virtual void on_modified_border_node(const dist_dynamic_graph& graph, NodeID border_node) = 0;

  bool check_for_include(const dist_dynamic_graph& graph, NodeID node);

  kamping::Communicator<> comm;
  std::vector<NodeID> consider;
  std::vector<NodeID> consider_next;
  fast_set consider_next_set;
  std::vector<wis_status> status;
  std::vector<NodeWeight> ratings;
  GlobalNodeWeight wis_weight;

  // stats
  std::size_t rounds = 0;
};

}  // namespace kadisredu::mwis
