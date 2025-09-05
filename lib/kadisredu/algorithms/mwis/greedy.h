//
// Created by jannickb on 6/18/24.
//

#pragma once

#include "kadisredu/algorithms/mwis/mwis_def.h"
#include "kadisredu/data_structures/dist_dynamic_graph.h"
#include "kadisredu/data_structures/node_marker.h"
#include "kadisredu/definitions.h"

#include <kamping/communicator.hpp>
#include <vector>

namespace kadisredu::mwis {

class greedy {
 public:
  greedy(MPI_Comm comm, NodeID num_nodes, NodeID num_ghosts)
      : comm(comm),
        candidates(num_nodes),
        status(num_nodes + num_ghosts, wis_status::UNSET),
        wis_weight(0),
        ratings(num_nodes + num_ghosts) {}

  virtual ~greedy() = default;

  virtual void run(dist_dynamic_graph& graph) = 0;

  Solution take_solution() { return {wis_weight, std::move(status)}; }

  [[nodiscard]] std::size_t get_rounds() const { return rounds; };

  void check_local_maximal_independent_set(dist_dynamic_graph& graph) {
    KASSERT(status.size() == graph.num_nodes() + graph.num_ghosts());
    for (NodeID node = 0; node < graph.num_nodes(); ++node) {
      if (status[node] == wis_status::INCLUDED) {
        for (auto target = graph[node].begin(); target != graph[node].hidden_end(); ++target) {
          KASSERT(status[*target] != wis_status::INCLUDED, "greedy solution is not independent");
        }
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
  void exclude_node(dist_dynamic_graph& graph, NodeID node);

  void include_node(dist_dynamic_graph& graph, NodeID node);

  virtual void on_modified_border_node(dist_dynamic_graph& graph, NodeID border_node) = 0;

  bool check_for_include(const dist_dynamic_graph& graph, NodeID node);

  kamping::Communicator<> comm;
  NodeMarker candidates;
  std::vector<wis_status> status;
  GlobalNodeWeight wis_weight;
  std::vector<GlobalNodeWeight> ratings;

  // stats
  std::size_t rounds = 0;
};

}  // namespace kadisredu::mwis
