//
// Created by jannickb on 12/31/24.
//

#pragma once

#include "kadisredu/algorithms/mwis/mwis_def.h"
#include "kadisredu/algorithms/mwis/peeler.h"
#include "kadisredu/algorithms/mwis/reduce_measurement.h"
#include "kadisredu/algorithms/mwis/weighted_reductions.h"
#include "kadisredu/data_structures/dist_dynamic_graph.h"
#include "kadisredu/data_structures/sized_vector.h"
#include "kadisredu/graphutils/graph_summary.h"
#include "kadisredu/tools/logger.h"
#include "kadisredu/tools/random_functions.h"
#include "kadisredu/tools/timer.h"

#include <algorithm>
#include <kamping/collectives/reduce.hpp>
#include <kamping/communicator.hpp>
#include <kamping/measurements/timer.hpp>
#include <limits>
#include <type_traits>

namespace kadisredu::mwis {

class Reducer {
 public:
  friend neighborhood_removal;
  friend domination;
  friend clique;
  friend simplicial_weight_transfer;
  friend degree_one;
  friend v_shape;
  friend basic_single_edge;
  friend extended_single_edge;
  friend generalized_neighborhood_removal;
  friend generalized_neighborhood_folding;
  friend neighborhood_folding;
  friend peeler;
  friend extended_domination;

  static constexpr NodeID BARRIER = std::numeric_limits<NodeID>::max();
  static constexpr NodeID MODIFIED_EDGE = std::numeric_limits<NodeID>::max() - 1;

  struct KernelStatistic {
    graph_summary kernel_graph;

    GlobalNodeWeight offset;

    std::size_t allreduce;
    std::size_t alltoallv;
    std::size_t max_recv_border_updates;
    std::size_t sum_recv_border_updates;
    std::size_t sum_recv_redundant_border_updates;
    double max_recv_redundant_border_updates_rel;

    struct ReductionStatistic {
      ReductionType type;
      NodeID min_reduced_local_nodes;
      NodeID max_reduced_local_nodes;
      NodeID min_reduced_ghosts;
      NodeID max_reduced_ghosts;
    };

    std::vector<ReductionStatistic> reduction_stats;
  };

  /**
    config with hyperparameters for the reduce algorithm
  */
  struct Config {
    ReductionStyle style;
    bool optimize_candidates_order;
    std::size_t max_allowed_subproblem_size;

    // async only
    unsigned message_queue_local_threshold;

    Config() = default;
    explicit Config(const reduction_context& ctx)
        : style(ctx.reduction_style),
          optimize_candidates_order(ctx.optimize_candidates_order),
          max_allowed_subproblem_size(ctx.max_allowed_subproblem_size),
          message_queue_local_threshold(ctx.message_queue_local_threshold) {}
  };

  class ReduceStatus {
   public:
    // friend Reducer;
    friend neighborhood_removal;
    friend domination;
    friend clique;
    friend simplicial_weight_transfer;
    friend degree_one;
    friend v_shape;
    friend basic_single_edge;
    friend extended_single_edge;
    friend generalized_neighborhood_removal;
    friend generalized_neighborhood_folding;
    friend neighborhood_folding;

    ReduceStatus() = default;

    ReduceStatus(Reducer::Config cfg, dist_dynamic_graph g, kamping::measurements::Timer<>* t);

    Solution solution;
    // sol_offset is used for weight shifts; if vertices are included or excluded, the weight is set
    // in solution.solution_weight
    GlobalNodeWeight sol_offset = 0;

    [[nodiscard]] GlobalNodeWeight get_solution_weight() const {
      return solution.solution_weight + sol_offset;
    }

    std::vector<BlockID> block_rank;

    dist_dynamic_graph graph;
    Reducer::Config cfg{};
    kamping::measurements::Timer<>* reduce_timer{nullptr};

    double time_limit = 0;  ////< current reduce time limit
    timer t;                ////< elapsed time in current reduce call

    // protected:
    std::vector<std::size_t> reduction_map;
    std::vector<ReductionPtr<general_weighted_reduction>> reductions;

    std::vector<NodeID> modified_queue;
    sized_vector<ReductionType> folded_queue;

    std::vector<ReduceMeasurement> reduction_measurements;

    // general communication stats
    std::size_t allreduce = 0;
    std::size_t alltoallv = 0;
    std::size_t recv_border_updates = 0;
    std::size_t recv_redundant_border_updates = 0;
  };

  static_assert(std::is_move_constructible<ReduceStatus>::value,
                "ReduceStatus must have a move constructor!");

  using weighted_reduce_ptr = Reducer*;

  const unsigned BUFFERS = 3;

  Reducer(MPI_Comm comm, NodeID capacity)
      : comm_(comm),
        set_1(capacity),
        set_2(capacity),
        buffers(BUFFERS, sized_vector<NodeID>(capacity)) {}

  [[nodiscard]] virtual ReducerType get_type() const = 0;

  static Reducer::ReduceStatus make_status(Reducer::Config cfg, dist_dynamic_graph graph,
                                           kamping::measurements::Timer<>* t) {
    ReduceStatus status{cfg, std::move(graph), t};
    init_reduction_style(status);
    return status;
  }

  void set_status(Reducer::ReduceStatus status) { status_ = std::move(status); }

  [[nodiscard]] const Reducer::ReduceStatus& get_status() const { return status_; }

  Reducer::ReduceStatus& get_status() { return status_; }

  void swap_status(Reducer::ReduceStatus& status) {
    using std::swap;
    swap(status, status_);
  }

  virtual std::pair<bool, bool> reduce(double time_limit) = 0;

  virtual void apply_and_restore() = 0;

  std::optional<Reducer::KernelStatistic> aggregate_reduce_statistic();

 protected:
  static void init_reduction_style(ReduceStatus& status);

  // Reduce Operations
  void set(NodeID node, wis_status new_node_status, bool poll = true);
  virtual bool poll_updates_during_reduce_phase() { return false; };
  /*
   * Exchange border updates
   * Return true if all border updates were exchanged and no global progress was made.
   */
  virtual bool exchange_updates() { return false; };

  /*
   * Reactivates the reduction phase if the reduction phase terminated globally.
   */
  virtual void reactivate_reduction_phase() {};

  void exclude_node(NodeID node);
  void fold_node(NodeID node);
  void include_node(NodeID node);
  void detach_node(NodeID node);
  void set_new_weight(NodeID node, NodeWeight new_weight);
  void include_node_with_bulk_hide(NodeID node, fast_set& neighbors, bool poll = true);
  void bulk_exclude(std::span<const NodeID> exclude, fast_set& exclude_set, bool poll = true);

  // helpers
  bool time_limit_reached() { return status_.t.elapsed_seconds() >= status_.time_limit; }
  void mark_node(NodeID node);
  void mark_neigh(NodeID node);
  virtual void on_border_status_change(NodeID border_node, wis_status new_status) = 0;
  virtual void on_border_weight_shift(NodeID border_node, NodeWeight new_weight) = 0;
  bool local_vertex_wins_tie(NodeID local, NodeID ghost) {
    auto& graph = status_.graph;
    KASSERT(graph.get_weight(local) == graph.get_weight(ghost));
    return status_.block_rank[comm_.rank()] < status_.block_rank[graph.get_ghost_block(ghost)];
  }

  const NodeID MAX_DETACH_DEGREE = 1;

  int active_reduction_index = 0;
  kamping::Communicator<> comm_;
  // current reduce status
  ReduceStatus status_;
  // buffers that can be used across different reduce status instances
  fast_set set_1;
  fast_set set_2;
  std::vector<sized_vector<NodeID>> buffers;
};

}  // namespace kadisredu::mwis
