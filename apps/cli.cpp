//
// Created by jannickb on 11/27/24.
//

#include "cli.h"

#include "kadisredu/tools/logger.h"

#include <CLI/CLI.hpp>
#include <CLI/Option.hpp>
#include <CLI/Validators.hpp>
#include <kadisredu/algorithms/mwis/mwis_def.h>
#include <toml++/toml.hpp>

namespace kadisredu::cli {

int parse_config(const std::string &filename, mwis::application_context &ctx) {
  using namespace kadisredu;
  using namespace kadisredu::mwis;
  using ReductionStyle = mwis::ReductionStyle;
  using ReducerType = mwis::ReducerType;
  using GreedyAlgo = mwis::GreedyAlgo;

  toml::table tbl;
  try {
    tbl = toml::parse_file(filename);
  } catch (const toml::parse_error &err) {
    std::cerr << "Parsing toml file" << filename << " failed: " << std::endl << err << std::endl;
    return 1;
  }

  std::map<std::string, ReductionStyle> red_style_map{
      {"normal", ReductionStyle::NORMAL},
      {"preliminary", ReductionStyle::PRELIMINARY},
      {"full", ReductionStyle::FULL},
      {"reduce-and-peel", ReductionStyle::REDUCE_AND_PEEL}};
  std::map<std::string, ReducerType> reducer_type_map{
      {"sync", ReducerType::SYNC},
      {"async", ReducerType::ASYNC},
  };
  std::map<std::string, GreedyAlgo> greedy_algo_map{
      {"sync", GreedyAlgo::SYNC},
      {"async", GreedyAlgo::ASYNC},
  };

  auto parse_partitioning_table = [](const toml::table &partitioning_table,
                                     partitioner_context &part_ctx) {
    if (auto entry = partitioning_table["disable_partitioning"].template value<bool>()) {
      part_ctx.disable_partitioning = *entry;
    }
    if (auto entry = partitioning_table["epsilon"].template value<double>()) {
      part_ctx.epsilon = *entry;
    }
    if (auto entry = partitioning_table["minimum_vertices_per_rank"].template value<int64_t>()) {
      part_ctx.minimum_vertices_per_rank = *entry;
    }
  };

  if (toml::array *reduce_phases = tbl["reduce_phase"].as_array()) {
    auto &red_ctx = ctx.red_ctx;
    red_ctx.resize(reduce_phases->size());
    auto current_red_ctx = red_ctx.begin();

    reduce_phases->for_each([&](auto &&reduce_phase) {
      if constexpr (toml::is_table<decltype(reduce_phase)>) {
        if (auto style = reduce_phase["style"].template value<std::string>()) {
          current_red_ctx->reduction_style = red_style_map[*style];
        }
        if (auto type = reduce_phase["type"].template value<std::string>()) {
          current_red_ctx->reducer_type = reducer_type_map[*type];
        }
        if (auto entry =
                reduce_phase["max_allowed_subproblem_size"].template value<std::size_t>()) {
          current_red_ctx->max_allowed_subproblem_size = *entry;
        }
        if (auto entry = reduce_phase["optimize_candidates_order"].template value<bool>()) {
          current_red_ctx->optimize_candidates_order = *entry;
        }
        if (auto entry = reduce_phase["message_queue_local_threshold"].template value<int64_t>()) {
          current_red_ctx->message_queue_local_threshold = *entry;
        }
        if (toml::table *part_table = reduce_phase["initial_partitioning"].as_table()) {
          auto &part_ctx = current_red_ctx->initial_partitioning;
          parse_partitioning_table(*part_table, part_ctx);
        }
      }
      ++current_red_ctx;
    });
  }

  if (auto maybe_greedy_table = tbl["greedy"].as_table()) {
    auto &greedy_ctx = ctx.greedy_ctx;
    auto &greedy_table = *maybe_greedy_table;
    if (auto entry = greedy_table["algo"].value<std::string>()) {
      greedy_ctx.greedy_algo = greedy_algo_map[*entry];
    }
    if (auto entry = greedy_table["message_queue_local_threshold"].value<int64_t>()) {
      greedy_ctx.message_queue_local_threshold = *entry;
    }
    if (toml::table *part_table = greedy_table["initial_partitioning"].as_table()) {
      auto &part_ctx = greedy_ctx.initial_partitioning;
      parse_partitioning_table(*part_table, part_ctx);
    }
  }

  return -1;
}

int parse_parameters(int argc, char *argv[], mwis::application_context &ctx, BlockID comm_size) {
  CLI::App app{"Distributed MWIS solver using data reductions", "KaDisRedu-MWIS"};
  app.add_option("--kagen_option_string", ctx.kagen_option_string,
                 "graph file or kagen option string")
      ->required();

  std::string config;
  app.add_option("--config", config, "load a toml config");

  auto print_kernel_option =
      app.add_option("--kernel", ctx.kernel_filename, "Writes kernel to this file (METIS)");

  auto print_independent_set_option = app.add_option(
      "--independent_set", ctx.independent_set_filename,
      "Writes independent set into this file. The i-th vertex received 0 (excluded) or "
      "1 (included) in i-th line.");

  app.add_option("--time_limit", ctx.time_limit, "choose an overall time limit (seconds)");
  app.add_option("-j,--json_output_path", ctx.json_output_path,
                 "json output file (statistics, time measurements)")
      ->required();
  app.add_option("--seed", ctx.seed, "random seed");
  app.add_flag("--warmup_mpi", ctx.warmup_mpi, "warm-up MPI");

  unsigned message_queue_local_threshold;
  auto opt_message_queue_local_threshold =
      app.add_option("--message_queue_local_threshold", message_queue_local_threshold,
                     "Set message_queue_local_threshold (bytes) for async.");

  app.parse_complete_callback([&]() {
    if (!config.empty()) {
      parse_config(config, ctx);
    }

    if (print_kernel_option->count() > 0) {
      ctx.print_kernel = true;
    }

    if (print_independent_set_option->count() > 0) {
      ctx.print_independent_set = true;
    }

    if (opt_message_queue_local_threshold->count() > 0) {
      ctx.greedy_ctx.message_queue_local_threshold = message_queue_local_threshold;
      for (auto &red_ctx : ctx.red_ctx) {
        red_ctx.message_queue_local_threshold = message_queue_local_threshold;
      }
    }
  });

  CLI11_PARSE(app, argc, argv);

  return -1;
}

}  // namespace kadisredu::cli
