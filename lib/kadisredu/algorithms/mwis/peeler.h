#pragma once

#include "kadisredu/algorithms/mwis/weighted_reductions.h"
#include "kadisredu/data_structures/maxNodeHeap.h"

namespace kadisredu::mwis {

class peeler : public general_weighted_reduction {
 public:
  peeler(NodeID num_nodes, NodeID num_ghosts) : general_weighted_reduction(num_nodes, num_ghosts) {}
  ~peeler() override = default;

  ReductionType get_reduction_type() final { return ReductionType::peeler; }

  bool reduce(weighted_reduce_ptr algo) override;
  reduction_ptr clone() override { return std::make_unique<peeler>(*this); }

 protected:
  maxNodeHeap pq_;
  bool init = false;
  NodeID peeled = 0;
  double sync_threshold = 1.0;
};

}  // namespace kadisredu::mwis
