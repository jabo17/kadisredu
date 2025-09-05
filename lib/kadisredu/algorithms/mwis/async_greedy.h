//
// Created by jannickb on 6/18/24.
//

#pragma once

#include "kadisredu/algorithms/mwis/greedy.h"
#include "kadisredu/data_structures/dist_dynamic_graph.h"

#include <message-queue/buffered_queue.hpp>

namespace kadisredu::mwis::async {

class async_greedy : public greedy {
 public:
  async_greedy(MPI_Comm comm, NodeID num_nodes, NodeID num_ghosts,
               unsigned message_queue_local_threshold)
      : greedy(comm, num_nodes, num_ghosts),
        border_update_queue(comm, 8, message_queue::ReceiveMode::poll,
                            message_queue::aggregation::AppendMerger{},
                            message_queue::aggregation::NoSplitter{}, printing_cleaner),
        reached_ranks(this->comm.size()),
        message_queue_local_threshold_(message_queue_local_threshold) {}

  ~async_greedy() override = default;

  void run(dist_dynamic_graph& graph) override;

 private:
  struct BorderUpdateMessage {
    GlobalNodeID node;
    wis_status status;
  };
  using BorderUpdateBuffer = std::vector<BorderUpdateMessage>;

  static void printing_cleaner(BorderUpdateBuffer& buf, message_queue::PEID receiver) {
    // message_queue::atomic_debug(fmt::format("Preparing buffer {} to {}.", buf.size(), receiver));
  }

  using BorderUpdateQueue = message_queue::BufferedMessageQueue<
      BorderUpdateMessage, BorderUpdateMessage, BorderUpdateBuffer, BorderUpdateBuffer,
      message_queue::aggregation::AppendMerger, message_queue::aggregation::NoSplitter,
      std::function<void(BorderUpdateBuffer&, message_queue::PEID receiver)>>;

  void set_node(dist_dynamic_graph& graph, NodeID node, wis_status new_status, bool poll = true);

  void on_modified_border_node(dist_dynamic_graph& graph, NodeID border_node) override;

  void poll_border_updates(dist_dynamic_graph& graph);
  bool terminate_border_update_queue(dist_dynamic_graph& graph);

  // manage border updates
  BorderUpdateQueue border_update_queue;
  fast_set reached_ranks;
  unsigned message_queue_local_threshold_;
};

}  // namespace kadisredu::mwis::async
