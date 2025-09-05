//
// Created by jannickb on 6/18/24.
//

#pragma once

#include "kadisredu/algorithms/mwis/greedy.h"
#include "kadisredu/algorithms/mwis/mwis_def.h"
#include "kadisredu/data_structures/dist_dynamic_graph.h"
#include "kadisredu/data_structures/fast_set.h"
#include "kadisredu/definitions.h"

namespace kadisredu::mwis::sync {

class sync_greedy : public greedy {
 public:
  sync_greedy(MPI_Comm comm, NodeID num_nodes, NodeID num_ghosts)
      : greedy(comm, num_nodes, num_ghosts), border_updates(num_nodes) {}

  ~sync_greedy() override = default;

  void run(dist_dynamic_graph& graph) override;

 private:
  void set_node(dist_dynamic_graph& graph, NodeID node, wis_status new_status);

  void on_modified_border_node(dist_dynamic_graph& graph, NodeID border_node) override;

  sized_vector<std::pair<NodeID, wis_status>> border_updates;

};

}  // namespace kadisredu::mwis::sync
