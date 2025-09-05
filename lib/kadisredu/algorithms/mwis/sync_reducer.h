//
// Created by jannickb on 6/24/24.
//

#pragma once

#include "kadisredu/algorithms/mwis/reducer.h"

namespace kadisredu::mwis {

class SyncReducer : public Reducer {
 public:
  SyncReducer(MPI_Comm comm, NodeID capacity)
      : Reducer(comm, capacity),
        border_updates(capacity),
        border_weight_updates(capacity),
        weight_modified(capacity),
        weight_modified_nodes(capacity),
        detached_since_last_sync(false),
        progress(false),
        timeout(false) {}

  [[nodiscard]] ReducerType get_type() const override { return ReducerType::SYNC; }

  std::pair<bool, bool> reduce(double time_limit) override;
  void apply_and_restore() override;

 protected:
  struct BorderUpdate {
    wis_status new_status;
  };
  sized_vector<std::pair<NodeID, BorderUpdate>> border_updates;

  fast_set weight_modified;
  sized_vector<NodeID> weight_modified_nodes;
  sized_vector<std::pair<NodeID, NodeWeight>> border_weight_updates;

  const bool LOCALLY_EXHAUSTIVE = true;
  bool detached_since_last_sync;
  bool progress;
  bool timeout;

  void sync_border_status_updates();
  void sync_border_weight_updates();

  std::tuple<bool, bool, bool> check_sync_borders();

  /*
   * @brief return true if borders were updated
   */
  bool check_and_sync_borders();

  bool exchange_updates() override { 
	using namespace kamping;
	check_and_sync_borders();
	// return false if progress was made
	return comm_.allreduce_single(send_buf((kabool)!progress),
                                  op(ops::logical_and<>()));
  };

  void on_border_status_change(NodeID border_node, wis_status new_status) override;
  void on_border_weight_shift(NodeID border_node, NodeWeight new_weight) override;
};

}  // namespace kadisredu::mwis
