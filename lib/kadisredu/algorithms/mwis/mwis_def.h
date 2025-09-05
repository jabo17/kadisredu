//
// Created by jannickb on 6/18/24.
//

#pragma once

#include "kadisredu/definitions.h"

namespace kadisredu::mwis {

enum class wis_status : int { INCLUDED = 0, EXCLUDED = 1, UNSET = 2, FOLDED = 3, DETACHED = 4 };

struct Solution {
  Solution() = default;
  Solution(NodeID nodes, GlobalNodeID ghosts)
      : node_status(ghosts + nodes, wis_status::UNSET), solution_weight(0) {}

  Solution(GlobalNodeWeight sol_weight, std::vector<wis_status> status)
      : node_status(std::move(status)), solution_weight(sol_weight) {}

  std::vector<wis_status> node_status;
  GlobalNodeWeight solution_weight = 0;
  bool optimal = false;
};

enum class ReductionStyle { PRELIMINARY, NORMAL, FULL, REDUCE_AND_PEEL };
static std::string to_string(const ReductionStyle& style) {
  switch (style) {
    case ReductionStyle::PRELIMINARY:
      return "preliminary";
    case ReductionStyle::NORMAL:
      return "normal";
    case ReductionStyle::FULL:
      return "full";
    case ReductionStyle::REDUCE_AND_PEEL:
      return "reduce-and-peel";
    default:
      return "unknown";
  }
}

enum class ReducerType { SYNC, ASYNC };
static std::string to_string(const ReducerType& comm_type) {
  switch (comm_type) {
    case ReducerType::SYNC:
      return "sync";
    case ReducerType::ASYNC:
      return "async";
    default:
      return "unknown";
  }
}

enum class GreedyAlgo { SYNC, ASYNC };
static std::string to_string(const GreedyAlgo& comm_type) {
  switch (comm_type) {
    case GreedyAlgo::SYNC:
      return "sync";
    case GreedyAlgo::ASYNC:
      return "async";
    default:
      return "unknown";
  }
}

struct partitioner_context {
  bool disable_partitioning = false;

  // If the global number of vertices is smaller than minimum_vertices_per_rank * comm.size(), do
  // not partition the graph and proceed on a single process.
  NodeID minimum_vertices_per_rank = 10;
  double epsilon = 0.03;

  bool operator==(const partitioner_context&) const = default;
};

struct reduction_context {
  ReducerType reducer_type = ReducerType::SYNC;
  ReductionStyle reduction_style = ReductionStyle::FULL;

  partitioner_context initial_partitioning;
  bool maximize_solution_greedily = true;

  unsigned message_queue_local_threshold = 16'000;

  std::size_t distributed_degree_one_threshold = 0;  //(not used by async approach)
  bool optimize_candidates_order = false;            // true;
  std::size_t max_allowed_subproblem_size = 100;     // std::numeric_limits<std::size_t>::max();
  bool operator==(const reduction_context&) const = default;
};

struct greedy_solver_context {
  partitioner_context initial_partitioning;
  GreedyAlgo greedy_algo = GreedyAlgo::SYNC;
  unsigned message_queue_local_threshold = 16'000;
};

/**
 * @struct application_context
 */
struct application_context {
  enum class preset { NONE, G, RG, PRG, RPRG, aG, aPRG, aRG, aRPRG };

  // general options
  int seed = 0;                  // random seed used by the algorithm
  std::string json_output_path;  // json output file with statistics and context
  unsigned time_limit = 3600;    // seconds
  bool quiet = false;            // suppress all standard output
  bool experiment = true;        // enables output of extra statistics
  bool print_kernel = false;
  std::string kernel_filename;
  bool print_independent_set = false;
  std::string independent_set_filename;
  bool warmup_mpi = false;

  // first reduction phase
  std::vector<reduction_context> red_ctx;

  greedy_solver_context greedy_ctx;

  // KaGen options
  bool check_input_graph = false;
  std::string kagen_option_string;

  void load_preset_context(preset p, BlockID comm_size) {
    // assumes default values are set
    if (p == preset::G || p == preset::RG || p == preset::RPRG || p == preset::PRG) {
      // sync
      greedy_ctx.initial_partitioning.disable_partitioning = true;
    } else {
      // async
      greedy_ctx = {.initial_partitioning = {.disable_partitioning = true},
                    .greedy_algo = GreedyAlgo::ASYNC};
    }

    if (p == preset::RG) {
      red_ctx.push_back({.initial_partitioning = {.disable_partitioning = true}});
    } else if (p == preset::PRG) {
      red_ctx.push_back({.initial_partitioning = {.disable_partitioning = false}});
    } else if (p == preset::RPRG) {
      red_ctx.push_back({.initial_partitioning = {.disable_partitioning = true}});
      red_ctx.push_back({.initial_partitioning = {.disable_partitioning = false}});
    } else if (p == preset::aRG) {
      red_ctx.push_back({.reducer_type = ReducerType::ASYNC,
                         .initial_partitioning = {.disable_partitioning = true}});
    } else if (p == preset::aPRG) {
      red_ctx.push_back({.reducer_type = ReducerType::ASYNC,
                         .initial_partitioning = {.disable_partitioning = false}});
    } else if (p == preset::aRPRG) {
      red_ctx.push_back({.reducer_type = ReducerType::ASYNC,
                         .initial_partitioning = {.disable_partitioning = true}});
      red_ctx.push_back({.reducer_type = ReducerType::ASYNC,
                         .initial_partitioning = {.disable_partitioning = false}});
    }
  };
};

}  // namespace kadisredu::mwis
