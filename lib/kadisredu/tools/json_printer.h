//
// Created by jannickb on 11/27/24.
//

#pragma once

#include "kadisredu/algorithms/mwis/mwis_def.h"
#include "kadisredu/algorithms/mwis/reducer.h"
#include "kadisredu/graphutils/graph_summary.h"

#include <nlohmann/json.hpp>
#include <ostream>
#include <sstream>

namespace kadisredu {

class JsonPrinter {
 public:
  static nlohmann::json to_json(const mwis::partitioner_context &ctx);

  static nlohmann::json to_json(const mwis::greedy_solver_context &ctx);

  static nlohmann::json to_json(const mwis::reduction_context &ctx);

  static nlohmann::json to_json(const mwis::application_context &ctx);

  static nlohmann::json get_json_reduce_statistic(mwis::Reducer::KernelStatistic &stats);

  static nlohmann::json get_json_from_stream(std::stringstream &stream);

  static nlohmann::json to_json(const graph_summary &summary);

  static void pretty_print_json(nlohmann::json &dump, std::ostream &out);
};

}  // namespace kadisredu
