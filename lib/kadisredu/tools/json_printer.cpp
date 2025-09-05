//
// Created by jannickb on 11/27/24.
//

#include "kadisredu/tools/json_printer.h"

#include <utils/output.hpp>

namespace kadisredu {
nlohmann::json JsonPrinter::to_json(const mwis::partitioner_context &ctx) {
  using namespace nlohmann;
  return {
      {"disable_partitioning", ctx.disable_partitioning},
      {"epsilon", ctx.epsilon},
      {"minimum_vertices_per_rank", ctx.minimum_vertices_per_rank},
  };
}

nlohmann::json JsonPrinter::to_json(const mwis::greedy_solver_context &ctx) {
  using namespace nlohmann;
  return {
      {"initial_partitioning", to_json(ctx.initial_partitioning)},
      {"greedy_algo", to_string(ctx.greedy_algo)},
      {"message_queue_local_threshold", ctx.message_queue_local_threshold},
  };
}

nlohmann::json JsonPrinter::to_json(const mwis::reduction_context &ctx) {
  using namespace nlohmann;
  return {
      {"initial_partitioning", to_json(ctx.initial_partitioning)},
      {"reduction_style", ctx.reduction_style},
      {"reducer_type", ctx.reducer_type},
      {"max_allowed_subproblem_size", ctx.max_allowed_subproblem_size},
      {"optimize_candidates_order", ctx.optimize_candidates_order},
      {"maximize_solution_greedily", ctx.maximize_solution_greedily},
      {"message_queue_local_threshold", ctx.message_queue_local_threshold},
  };
}
nlohmann::json JsonPrinter::to_json(const mwis::application_context &ctx) {
  using namespace nlohmann;
  return {{"seed", ctx.seed},
          {"time_limit", ctx.time_limit},
          {"red_ctx",
           [&]() {
             std::vector<nlohmann::json> vec(ctx.red_ctx.size());
             std::ranges::transform(ctx.red_ctx, vec.begin(),
                                    [&](auto &red_ctx) { return to_json(red_ctx); });
             return vec;
           }()},
          {"greedy_solver_context", to_json(ctx.greedy_ctx)},
          {"check_input_graph", ctx.check_input_graph},
          {"kagen_option_string", ctx.kagen_option_string},
          {"experiment", ctx.experiment}};
}
nlohmann::json JsonPrinter::get_json_reduce_statistic(mwis::Reducer::KernelStatistic &stats) {
  using namespace nlohmann;
  return {{"kernel_graph", to_json(stats.kernel_graph)},
          {"offset", stats.offset},
          {"alltoallv", stats.alltoallv},
          {"allreduce", stats.allreduce},
          {"sum_recv_border_updates", stats.sum_recv_border_updates},
          {"max_recv_border_updates", stats.max_recv_border_updates},
          {"sum_recv_redundant_border_updates", stats.sum_recv_redundant_border_updates},
          {"max_recv_redundant_border_updates_rel", stats.max_recv_redundant_border_updates_rel},
          {"reduction_map",
           [&stats]() {
             std::vector<std::string> vec(stats.reduction_stats.size());
             std::ranges::transform(stats.reduction_stats, vec.begin(),
                                    [&](auto &stat) { return mwis::to_string(stat.type); });
             return vec;
           }()},
          {"reduction_statistics", [&stats]() {
             json red_stats;
             for (auto &stat : stats.reduction_stats) {
               red_stats[mwis::to_string(stat.type)] = {
                   {"min_reduced_local_nodes", stat.min_reduced_local_nodes},
                   {"max_reduced_local_nodes", stat.max_reduced_local_nodes},
                   {"min_reduced_ghosts", stat.min_reduced_ghosts},
                   {"max_reduced_ghosts", stat.max_reduced_ghosts}};
             }
             return red_stats;
           }()}};
}
nlohmann::json JsonPrinter::get_json_from_stream(std::stringstream &stream) {
  using namespace nlohmann;
  json j;
  stream >> j;
  return j;
}
void JsonPrinter::pretty_print_json(nlohmann::json &dump, std::ostream &out) {
  out << std::setw(4) << dump << std::endl;
}
nlohmann::json JsonPrinter::to_json(const graph_summary &summary) {
  return {{"nodes", summary.nodes},
          {"max_nodes", summary.max_nodes},
          {"min_nodes", summary.min_nodes},
          {"imbalance", summary.imbalance},
          {"edges", summary.edges},
          {"ghosts", summary.ghosts},
          {"max_ghosts", summary.max_ghosts},
          {"min_ghosts", summary.min_ghosts},
          {"cut", summary.cut},
          {"max_degree", summary.max_degree},
          {"min_degree", summary.min_degree},
          {"degree_dist", [&summary]() {
             nlohmann::json j{};
             for (NodeID degree = 0; degree < summary.degree_counts.size(); degree++) {
               j[degree] = static_cast<double>(summary.degree_counts[degree]) /
                           static_cast<double>(summary.nodes);
             }
             return j;
           }()}};
}
}  // namespace kadisredu
