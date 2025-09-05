//
// Created by jannickb on 6/18/24.
//

#pragma once

#include "kadisredu/algorithms/mwis/greedy_maximize.h"
#include "kadisredu/data_structures/dist_dynamic_graph.h"

#include <message-queue/buffered_queue.hpp>

namespace kadisredu::mwis::sync {

class sync_greedy_maximize : public greedy_maximize {
 public:
  sync_greedy_maximize(MPI_Comm comm, NodeID num_nodes, NodeID num_ghosts)
      : greedy_maximize(comm, num_nodes, num_ghosts), border_updates(num_nodes) {}

  ~sync_greedy_maximize() override = default;

  void run(const dist_dynamic_graph& graph, const Solution& solution) override;

 private:
  void set_node(const dist_dynamic_graph& graph, NodeID node, wis_status new_status);

  void on_modified_border_node(const dist_dynamic_graph& graph, NodeID border_node) override;

  // manage border updates
  sized_vector<std::pair<NodeID, wis_status>> border_updates;
};

}  // namespace kadisredu::mwis::sync
