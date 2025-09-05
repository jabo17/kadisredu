//
// Created by borowitz on 30.05.2 .
//

#pragma once

#include "kadisredu/algorithms/mwis/async_reducer.h"
#include "kadisredu/algorithms/mwis/mwis_def.h"
#include "kadisredu/algorithms/mwis/reducer.h"
#include "kadisredu/algorithms/mwis/sync_reducer.h"
#include "kadisredu/data_structures/dist_dynamic_graph.h"
#include "kadisredu/definitions.h"
#include "kadisredu/graphutils/redistribution.h"
#include "kadisredu/tools/json_printer.h"
#include "kadisredu/tools/timer.h"

#include <kagen/kagen.h>
#include <nlohmann/json.hpp>

namespace kadisredu::mwis {

template <typename... Args>
std::unique_ptr<Reducer> create_reducer(ReducerType choice, Args &&...args) {
  switch (choice) {
    case ReducerType::ASYNC:
      return std::make_unique<AsyncReducer>(std::forward<Args>(args)...);
    case ReducerType::SYNC:
      return std::make_unique<SyncReducer>(std::forward<Args>(args)...);
  }
  __builtin_unreachable();
}

class kadisredu_mwis {
 public:
  kadisredu_mwis(MPI_Comm comm, application_context _ctx)
      : comm_(comm),
        app_ctx_(std::move(_ctx)),
        t_(kamping::measurements::timer()),
        current_global_offset_on_root_(0),
        restore_data_(),
        fully_reduced_(false),
        timeout_(false) {}

  /// @brief Reduces the input graph @p graph and then applies a solver @p solver to the
  /// reduced graph @p solver.
  /// @param input graph
  /// @param solver as a binary function that takes the global reduction offset and reduced
  /// kagen graph as input. @p solver is assumed to return a Solution for the reduced graph.
  /// @return the solution for graph
  template <typename Solver>
  Solution reduce_and_solve_kagen_graph(kagen::Graph graph, Solver &&solver) {
    return impl_reduce_and_solve(std::move(graph), [&](GlobalNodeWeight global_offset,
                                                       const dist_dynamic_graph &reduced_graph) {
      auto [reduced_kagen_graph, _] = reduced_graph.build_visible_graph();
      return solver(global_offset, std::move(reduced_kagen_graph));
    });
  };

  /// @brief Reduces the input graph @p graph and then applies a greedy algorithm to the reduced
  /// graph. If no reduction phase is configured, the greedy solver is called
  /// @param input graph
  /// @return the solution for graph
  Solution reduce_and_greedy(kagen::Graph graph);

  nlohmann::json take_json_dump() {
    json_dump_["context"] = JsonPrinter::to_json(app_ctx_);
    json_dump_["context"]["np"] = comm_.size();
    return std::move(json_dump_);
  }

 private:
  kamping::Communicator<> comm_;
  const application_context app_ctx_;
  timer local_t_;
  kamping::measurements::Timer<> &t_;

  std::vector<std::unique_ptr<Reducer>> reduce_steps_;
  dist_dynamic_graph current_graph_;
  GlobalNodeWeight current_global_offset_on_root_;
  std::vector<GlobalNodeWeight> global_offsets_on_root_;

  std::vector<redistribution::RedistributionRestoreData> restore_data_;

  bool fully_reduced_;
  bool timeout_;

  nlohmann::json json_dump_;

  template <typename Solver>
  Solution impl_reduce_and_solve(kagen::Graph &&graph, Solver &&solver) {
    KASSERT(app_ctx_.red_ctx.size() > 0, "Expected at least one reduction context.");
    local_t_.restart();
    reduce(std::move(graph));
    KASSERT(!reduce_steps_.empty());

    Solution kernel_solution;
    // TODO: how to handle timeout_==true case
    if (!fully_reduced_) {
      kernel_solution = solver(current_global_offset_on_root_, get_reduced_graph_from_status());
    }

    auto solution = reconstruct(kernel_solution);
    set_best_solution_weight_stat_on_root(solution);
    return solution;
  };

  Solution greedy(GlobalNodeWeight global_offset_on_root, dist_dynamic_graph graph);

  void perform_reduce_step(const reduction_context &red_ctx, const std::string &red_phase_ident);

  dist_dynamic_graph &get_reduced_graph_from_status();

  std::pair<kagen::Graph, std::vector<NodeID>> build_reduced_kagen_graph();

  void reduce(kagen::Graph &&graph);

  static std::string reduce_phase_ident(std::size_t step) { return "R" + std::to_string(step); }

  void build_initial_distributed_dynamic_graph(kagen::Graph &&graph,
                                               const partitioner_context &part_ctx);

  void build_next_distributed_dynamic_graph(const partitioner_context &part_ctx);

  /// @brief Partitions @p graph with dKaMinPar and then build a dynamic graph that is distributed
  /// according to the partitioning.
  /// @returns a pair holding the redistribution data to undo the redistribution and the
  /// redistributed distributed dynamic graph
  [[nodiscard]] std::pair<redistribution::RedistributionRestoreData, dist_dynamic_graph>
  partition_and_rearrange(kagen::Graph &&graph, const partitioner_context &part_ctx) const;

  // @brief Reconstruct a solution for the reduced graph.
  // @param kernel_solution provide a solution for the reduced graph
  // @return solution for the reduced graph. The solution is restricted to the number of local
  // vertices.
  Solution reconstruct(const Solution &kernel_solution);

  /// @brief Reconstruct an MWIS for the graph that was reduced given a solution for the reduced
  /// graph. The independent set is at least maximal and optimal if @p kernel_solution is optimal.
  /// @return the reconstructed solution where the status of each local vertex is set together with
  /// the local solution weight.
  Solution reconstruct_last_reduce_step(const Solution &kernel_solution,
                                        const reduction_context &red_ctx,
                                        const std::string &red_phase_ident);

  /// todo: this should be part of the reducer
  /// @brief Sets solution for local visible vertices given the kernel solution and computes the
  /// local solution weight. Expects for the \c i-th visible vertex that its solution status is the
  /// @p kernel_solution.node_status[i].
  static void set_local_solution(const Solution &kernel_solution, Reducer::ReduceStatus &status);

  /// @brief Sets the solution status for ghosts.
  static void sync_ghost_solution(Reducer::ReduceStatus &status);

  void rearrange_solution(Solution &solution) const;

  void agg_reduce_stats_on_root(Reducer &reducer, const std::string &red_phase_ident);

  void agg_graph_summary_on_root(dist_dynamic_graph &graph, const std::string &red_phase_ident);

  std::optional<GlobalNodeWeight> get_global_solution_weight_on_root(const Solution &solution);

  void set_best_solution_weight_stat_on_root(const Solution &solution);

  void set_greedy_sol_stat_on_root(const Solution &solution);

  void print_current_best_solution_on_root(GlobalNodeWeight global_offset_on_root,
                                           const Solution &local_solution);
};

}  // namespace kadisredu::mwis
