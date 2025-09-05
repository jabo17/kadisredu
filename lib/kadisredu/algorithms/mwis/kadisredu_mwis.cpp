//
// Created by borowitz on 30.05.25.
//

#include "kadisredu_mwis.h"

#include "kadisredu/algorithms/mwis/debug.h"
#include "kadisredu/algorithms/mwis/mwis_def.h"
#include "kadisredu/algorithms/mwis/reducer.h"
#include "kadisredu/algorithms/mwis/solver_factories.h"
#include "kadisredu/data_structures/dist_dynamic_graph.h"
#include "kadisredu/definitions.h"
#include "kadisredu/graphutils/dist_dynamic_graph_factories.h"
#include "kadisredu/graphutils/graph_summary.h"
#include "kadisredu/graphutils/redistribution.h"
#include "kadisredu/tools/graph_io.h"
#include "kadisredu/tools/json_printer.h"
#include "kadisredu/tools/kaminpar_helper.h"
#include "kadisredu/tools/logger.h"
#include "kadisredu/tools/timer.h"
#include "nlohmann/json_fwd.hpp"

namespace kadisredu::mwis {

Solution kadisredu_mwis::reduce_and_greedy(kagen::Graph graph) {
  auto &greedy_ctx = app_ctx_.greedy_ctx;
  if (app_ctx_.red_ctx.size() == 0) {
    // Build first distributed dynamic graph from a kagen graph.
    // This is the first distributed dynamic graph because no reduction were applied before.
    comm_.barrier();
    t_.start("greedy");
    build_initial_distributed_dynamic_graph(std::move(graph), greedy_ctx.initial_partitioning);
    auto result = greedy(0, std::move(current_graph_));
    set_greedy_sol_stat_on_root(result);
    set_best_solution_weight_stat_on_root(result);
    t_.stop();  // stop greedy
    return result;
  }
  return impl_reduce_and_solve(
      std::move(graph),
      [&](GlobalNodeWeight global_offset,
          const dist_dynamic_graph &current_graph) -> decltype(auto) {
        // Build next distributed dynamic graph
        t_.start("greedy");
        KASSERT(&current_graph == &get_reduced_graph_from_status());
        build_next_distributed_dynamic_graph(greedy_ctx.initial_partitioning);
        auto result = greedy(global_offset, std::move(current_graph_));
        set_greedy_sol_stat_on_root(result);
        t_.stop();
        return result;
      });
}
Solution kadisredu_mwis::greedy(GlobalNodeWeight global_offset_on_root, dist_dynamic_graph graph) {
  agg_graph_summary_on_root(graph, "greedy");
  auto &greedy_ctx = app_ctx_.greedy_ctx;
  // build and apply greedy solver
  auto greedy_solver = create_greedy_algo(greedy_ctx, comm_.mpi_communicator(), graph.num_nodes(),
                                          graph.num_ghosts());

  KADISREDU_RLOG << "[greedy]: Running " << to_string(greedy_ctx.greedy_algo)
                 << " greedy algorithm";

  greedy_solver->run(graph);
  auto greedy_solution = std::move(greedy_solver->take_solution());
  // restrict greedy solution to local solution
  greedy_solution.node_status.resize(graph.num_nodes());

  print_current_best_solution_on_root(global_offset_on_root, greedy_solution);

  KASSERT(greedy_ctx.initial_partitioning.disable_partitioning || !restore_data_.empty(),
          "Expected RedistributionRestoreData when graph was redistributed before greedily "
          "solving the graph.");
  if (!greedy_ctx.initial_partitioning.disable_partitioning) {
    rearrange_solution(greedy_solution);
  }
  return greedy_solution;
}
void kadisredu_mwis::perform_reduce_step(const reduction_context &red_ctx,
                                         const std::string &red_phase_ident) {
  using namespace kamping;
  agg_graph_summary_on_root(current_graph_, red_phase_ident);
  KADISREDU_RLOG << "[" << red_phase_ident << "]: " << "Running " << to_string(red_ctx.reducer_type)
                 << " reduce algorithm with style " << to_string(red_ctx.reduction_style);
  t_.start("reduce-phase");
  // Setup reducer
  Reducer::Config config(red_ctx);
  auto reducer = create_reducer(red_ctx.reducer_type, comm_.mpi_communicator(),
                                current_graph_.num_nodes() + current_graph_.num_ghosts());
  reducer->set_status(Reducer::make_status(config, std::move(current_graph_), &t_));

  // Run reductions
  std::tie(fully_reduced_, timeout_) =
      reducer->reduce(app_ctx_.time_limit - local_t_.elapsed_seconds());

  auto global_offset =
      comm_.reduce_single(send_buf(reducer->get_status().get_solution_weight()), op(ops::plus<>()));
  if (comm_.is_root()) {
    KASSERT(global_offset.has_value());
    global_offsets_on_root_.push_back(global_offset.value());
    current_global_offset_on_root_ += global_offset.value();
  }

  t_.stop();  // stop reduce-phase
  agg_reduce_stats_on_root(*reducer, red_phase_ident);
  reduce_steps_.push_back(std::move(reducer));
  KADISREDU_RLOG << "Reducer finished!";
}
dist_dynamic_graph &kadisredu_mwis::get_reduced_graph_from_status() {
  KASSERT(!reduce_steps_.empty());
  dist_dynamic_graph &graph = reduce_steps_.back()->get_status().graph;
  NodeID n = graph.num_visible_nodes();
  using namespace kamping;
  graph.comm.allreduce_inplace(send_recv_buf(n), op(ops::plus<>()));
  return graph;
}
std::pair<kagen::Graph, std::vector<NodeID>> kadisredu_mwis::build_reduced_kagen_graph() {
  t_.start("build-reduced-kagen-graph");
  auto result = get_reduced_graph_from_status().build_visible_graph();
  t_.stop();
  return result;
}
void kadisredu_mwis::reduce(kagen::Graph &&graph) {
  t_.start("reduce");
  const auto &red_ctx = app_ctx_.red_ctx;
  KASSERT(red_ctx.size() > 0, "Expected at least one reduction phase.");

  // TODO: preprocessing
  // idea: test and apply some simple reductions directly on the kagen graph
  // -> this may reduce the size of the dynamic graph and further communication of vtx. weights
  // (when the dyn graph is build)

  // rearrange or build directly
  t_.start(reduce_phase_ident(0));
  build_initial_distributed_dynamic_graph(std::move(graph), red_ctx[0].initial_partitioning);
  perform_reduce_step(red_ctx[0], reduce_phase_ident(0));
  t_.stop();

  // further reduction steps
  for (std::size_t i = 1; i < red_ctx.size() && !fully_reduced_ && !timeout_; ++i) {
    t_.start(reduce_phase_ident(i));
    build_next_distributed_dynamic_graph(red_ctx[i].initial_partitioning);
    // Apply next reduction context in the next phase
    perform_reduce_step(red_ctx[i], reduce_phase_ident(i));
    t_.stop();
  }

  t_.stop();  // stop reduce
  if (app_ctx_.print_kernel) {
    t_.start("print-kernel");
    auto [kernel, _] = get_reduced_graph_from_status().build_visible_graph();
    write_metis_graph(comm_.mpi_communicator(), kernel, app_ctx_.kernel_filename);
    t_.stop();
  }
}
void kadisredu_mwis::build_initial_distributed_dynamic_graph(kagen::Graph &&graph,
                                                             const partitioner_context &part_ctx) {
  t_.start("build-distributed-dynamic-graph");
  if (part_ctx.disable_partitioning) {
    current_graph_ = build_dist_dynamic_graph(comm_, graph, t_);
  } else {
    auto [rest_data, dyn_graph] = partition_and_rearrange(std::move(graph), part_ctx);
    current_graph_ = std::move(dyn_graph);
    restore_data_.push_back(std::move(rest_data));
  }
  t_.stop();
}
void kadisredu_mwis::build_next_distributed_dynamic_graph(const partitioner_context &part_ctx) {
  t_.start("build-distributed-dynamic-graph");
  if (part_ctx.disable_partitioning) {
    // Build reduced graph by copying visible part of dynamic graph
    auto &graph = get_reduced_graph_from_status();
    std::vector<NodeID> new_to_old(graph.num_visible_nodes() + graph.num_visible_ghosts());
    current_graph_ = dist_dynamic_graph(graph, new_to_old);
  } else {
    // Build kagen graph of current reduced graph for the partitioner.
    auto [g, old_to_new] = build_reduced_kagen_graph();
    // Obtain restore data to undo the redistribution for the reconstructed solution
    // and set the redistributed dynamic graph as the current graph
    restore_data_.emplace_back();
    std::tie(restore_data_.back(), current_graph_) =
        partition_and_rearrange(std::move(g), part_ctx);
  }
  t_.stop_and_append();
}
std::pair<redistribution::RedistributionRestoreData, dist_dynamic_graph>
kadisredu_mwis::partition_and_rearrange(kagen::Graph &&graph,
                                        const partitioner_context &part_ctx) const {
  // compute partition into comm.size() blocks with dKaMinPar
  comm_.barrier();
  KADISREDU_RLOG << "Partition graph";
  t_.start("partition");
  std::vector<BlockID> partition(graph.NumberOfLocalVertices(), comm_.rank());
  kaminpar_helper::partition_graph(graph, comm_.mpi_communicator(), app_ctx_.seed, part_ctx,
                                   partition);
  t_.stop();
  comm_.barrier();
  // rearrange graph according to partition
  KADISREDU_RLOG << "Redistributing graph according to partitioning";
  t_.start("redistribute-and-build-dynamic-graph");
  auto result = redistribution::redistribute_and_build_dist_dynamic_graph(
      comm_.mpi_communicator(), std::move(graph), partition);
  t_.stop();

  return result;
}
Solution kadisredu_mwis::reconstruct(const Solution &kernel_solution) {
  t_.start("reconstruct");
  const int reduce_steps = static_cast<int>(reduce_steps_.size());
  KASSERT(reduce_steps > 0);

  Solution solution = reconstruct_last_reduce_step(
      kernel_solution, app_ctx_.red_ctx[reduce_steps - 1], reduce_phase_ident(reduce_steps - 1));
  for (int i = reduce_steps - 2; i >= 0; i--) {
    solution = reconstruct_last_reduce_step(solution, app_ctx_.red_ctx[i], reduce_phase_ident(i));
  }
  t_.stop();  // stop reconstruct
  return std::move(solution);
}
Solution kadisredu_mwis::reconstruct_last_reduce_step(const Solution &kernel_solution,
                                                      const reduction_context &red_ctx,
                                                      const std::string &red_phase_ident) {
  t_.start(red_phase_ident);
  auto reducer = std::move(reduce_steps_.back());
  reduce_steps_.pop_back();
  auto &status = reducer->get_status();

  t_.start("reducer");
  set_local_solution(kernel_solution, status);
  sync_ghost_solution(status);
  if (comm_.is_root()) {
    KASSERT(current_global_offset_on_root_ == std::reduce(global_offsets_on_root_.begin(),
                                                          global_offsets_on_root_.end(),
                                                          GlobalNodeWeight{0}));
    current_global_offset_on_root_ -= global_offsets_on_root_.back();
    global_offsets_on_root_.pop_back();
  }

  reducer->apply_and_restore();
  t_.stop();
  debug::check_solution_weight(status.solution, status.graph);

  // take solution, graph, and cfg from status before freeing reducer
  Solution solution = std::move(status.solution);
  auto graph = std::move(status.graph);
  auto cfg = std::move(status.cfg);
  auto comm_type = reducer->get_type();
  reducer.reset();

  // maximize solution
  if (red_ctx.maximize_solution_greedily) {
    t_.start("maximize");
    auto maximizer = create_greedy_maximizer(comm_type, cfg, comm_.mpi_communicator(),
                                             graph.num_nodes(), graph.num_ghosts());
    maximizer->run(graph, solution);
    solution = std::move(maximizer->take_solution());
    print_current_best_solution_on_root(current_global_offset_on_root_, solution);
    t_.stop();  // stop maximize
  }
  solution.node_status.resize(graph.num_nodes());

  // rearrange solution according to prior vertex distribution
  if (!red_ctx.initial_partitioning.disable_partitioning) {
    KASSERT(!restore_data_.empty());
    t_.start("reconstruct");
    rearrange_solution(solution);
    t_.stop();  // stop reconstruct
  }
  t_.stop();  // stop reduce_phase_ident

  return solution;
}
void kadisredu_mwis::set_local_solution(const Solution &kernel_solution,
                                        Reducer::ReduceStatus &status) {
  auto &graph = status.graph;
  auto &sol = status.solution;
  auto &node_status = sol.node_status;
  NodeID v = 0;
  KASSERT(graph.num_visible_nodes() == kernel_solution.node_status.size());
  graph.for_each_visible_node(
      [&](NodeID former_v) { node_status[former_v] = kernel_solution.node_status[v++]; });
  sol.solution_weight += kernel_solution.solution_weight;
}
void kadisredu_mwis::sync_ghost_solution(Reducer::ReduceStatus &status) {
  auto &graph = status.graph;
  auto &sol = status.solution;

  sol.node_status.resize(graph.num_nodes() + graph.num_ghosts());

  auto ghost_status = graph.sync_ghosts_with_lazy_updates<wis_status>(
      std::span(graph.get_border()), [&](NodeID v) -> std::optional<wis_status> {
        if (graph.is_visible(v)) {
          KASSERT(debug::is_excluded_or_included(sol.node_status[v]),
                  "Interface vertex is neither included nor excluded.");
          return sol.node_status[v];
        }
        return std::nullopt;
      });
  for (auto &[global_id, g_status] : ghost_status) {
    KASSERT(status.graph.is_visible(graph.get_ghost(global_id)),
            "Received ghost status for a ghost that is hidden while the corresponding interface "
            "vertex at the adjacent PE is visible.");
    sol.node_status[graph.get_ghost(global_id)] = g_status;
  }
}
void kadisredu_mwis::rearrange_solution(Solution &solution) const {
  KASSERT(!restore_data_.empty());

  using namespace kamping;

  auto &rest_data = restore_data_.back();
  auto &new_to_old = rest_data.new_to_old;
  auto &old_weights = rest_data.old_weights;

  // Obtain solution status for vertices before redistribution.
  auto node_status = comm_.alltoallv(send_buf(solution.node_status),
                                     send_counts(rest_data.reverse_partition_counts),
                                     recv_counts(rest_data.partition_counts));
  KASSERT(node_status.size() == new_to_old.size());

  solution.solution_weight = 0;
  solution.node_status.resize(node_status.size());

  // reverse local permutation and then set node status
  for (NodeID v = 0; v < new_to_old.size(); ++v) {
    NodeID old_v = new_to_old[v];
    solution.node_status[old_v] = node_status[v];
    if (node_status[v] == wis_status::INCLUDED) {
      solution.solution_weight += old_weights[old_v];
    }
  }
}
void kadisredu_mwis::agg_reduce_stats_on_root(Reducer &reducer,
                                              const std::string &red_phase_ident) {
  auto stats = reducer.aggregate_reduce_statistic();
  if (comm_.is_root()) {
    KASSERT(stats.has_value());
    auto &s = stats.value();

    KASSERT(stats.has_value());

    KADISREDU_RLOG << std::setw(15) << "new kernel" << std::setw(15) << "overall" << std::setw(15)
                   << "min" << std::setw(15) << "max";
    KADISREDU_RLOG << std::setw(15) << "nodes" << std::setw(15) << s.kernel_graph.nodes
                   << std::setw(15) << s.kernel_graph.min_nodes << std::setw(15)
                   << s.kernel_graph.max_nodes;
    KADISREDU_RLOG << std::setw(15) << "ghosts" << std::setw(15) << s.kernel_graph.ghosts
                   << std::setw(15) << s.kernel_graph.min_ghosts << std::setw(15)
                   << s.kernel_graph.max_ghosts;
    KADISREDU_RLOG << std::setw(15) << "undirected edges" << std::setw(15) << s.kernel_graph.edges;
    KADISREDU_RLOG << std::setw(15) << "undirected cut-edges" << std::setw(15)
                   << s.kernel_graph.cut;
    KADISREDU_RLOG << std::setw(15) << "red. offset" << std::setw(15) << s.offset;

    json_dump_[red_phase_ident]["reducer"] = JsonPrinter::get_json_reduce_statistic(stats.value());
  }
}
void kadisredu_mwis::agg_graph_summary_on_root(dist_dynamic_graph &graph,
                                               const std::string &red_phase_ident) {
  if (app_ctx_.experiment) {
    auto summary = graph_summary::agg_graph_summary_on_root(comm_, graph);
    if (comm_.is_root()) {
      KASSERT(summary.has_value());
      json_dump_[red_phase_ident]["graph"] = JsonPrinter::to_json(summary.value());
    }
  }
}
std::optional<GlobalNodeWeight> kadisredu_mwis::get_global_solution_weight_on_root(
    const Solution &solution) {
  using namespace kamping;

  auto sol_weight = comm_.reduce_single(send_buf(solution.solution_weight), op(ops::plus<>()));
  KASSERT(!comm_.is_root() || sol_weight.has_value());
  return sol_weight;
}
void kadisredu_mwis::set_best_solution_weight_stat_on_root(const Solution &solution) {
  using namespace kamping;

  auto sol_weight = get_global_solution_weight_on_root(solution);
  if (comm_.is_root()) {
    KADISREDU_RLOG << "Best WIS weight: " << std::to_string(sol_weight.value());
    json_dump_["solution_weight"] = sol_weight.value();
  }
}
void kadisredu_mwis::set_greedy_sol_stat_on_root(const Solution &solution) {
  using namespace kamping;

  auto sol_weight = get_global_solution_weight_on_root(solution);
  if (comm_.is_root()) {
    json_dump_["greedy"]["solution_weight"] = sol_weight.value();
  }
}
void kadisredu_mwis::print_current_best_solution_on_root(GlobalNodeWeight global_offset_on_root,
                                                         const Solution &local_solution) {
  auto sol_weight = get_global_solution_weight_on_root(local_solution);
  if (comm_.is_root()) {
    KASSERT(sol_weight.has_value());
    KADISREDU_RLOG << "Current best solution weight "
                   << (sol_weight.value() + global_offset_on_root);
  }
}
}  // namespace kadisredu::mwis
