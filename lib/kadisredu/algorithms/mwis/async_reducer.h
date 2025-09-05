//
// Created by jannickb on 6/24/24.
//

#pragma once

#include "kadisredu/algorithms/mwis/reducer.h"

#include <fmt/format.h>
#include <message-queue/buffered_queue.hpp>

namespace kadisredu::mwis {

class AsyncReducer : public Reducer {
 public:
  AsyncReducer(MPI_Comm comm, NodeID capacity)
      : Reducer(comm, capacity),
        reached_ranks(comm_.size()),
        border_update_queue(comm, 8, message_queue::ReceiveMode::poll,
                            message_queue::aggregation::AppendMerger{},
                            message_queue::aggregation::NoSplitter{},
                            [this](BorderUpdateBuffer& buf, message_queue::PEID receiver) {
                              this->printing_cleaner(buf, receiver);
                            }) {}

  [[nodiscard]] ReducerType get_type() const override { return ReducerType::ASYNC; }

  std::pair<bool, bool> reduce(double time_limit) override;
  void apply_and_restore() override;

 protected:
#if KASSERT_ENABLED(KASSERT_ASSERTION_LEVEL_NORMAL)
  bool is_polling_border_updates = false;
#endif
  struct BorderUpdateMessage {
    GlobalNodeID node;
    wis_status status;
    NodeWeight weight;
#ifndef NDEBUG
    NodeWeight debug = 0;
#endif
  };
  using BorderUpdateBuffer = std::vector<BorderUpdateMessage>;

  void printing_cleaner(BorderUpdateBuffer& buf, message_queue::PEID receiver);

  using BorderUpdateQueue = message_queue::BufferedMessageQueue<
      BorderUpdateMessage, BorderUpdateMessage, BorderUpdateBuffer, BorderUpdateBuffer,
      message_queue::aggregation::AppendMerger, message_queue::aggregation::NoSplitter,
      std::function<void(BorderUpdateBuffer&, message_queue::PEID receiver)>>;

  bool poll_updates_during_reduce_phase() override;
  bool exchange_updates() override;
  void reactivate_reduction_phase() override;

  bool on_recv_border_update(
      message_queue::Envelope<AsyncReducer::BorderUpdateMessage> auto envelope);
  void on_border_status_change(NodeID border_node, wis_status new_status) override;
  void on_border_weight_shift(NodeID border_node, NodeWeight new_weight) override;

  // manage border updates
  BorderUpdateQueue border_update_queue;
  fast_set reached_ranks;
};

}  // namespace kadisredu::mwis
