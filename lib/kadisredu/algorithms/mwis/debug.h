//
// Created by borowitz on 30.05.25.
//

#pragma once
#include "kadisredu/algorithms/mwis/mwis_def.h"
#include "kadisredu/data_structures/dist_dynamic_graph.h"

#include <kassert/kassert.hpp>

namespace kadisredu::mwis::debug {

inline void check_solution_weight(const Solution &solution, const dist_dynamic_graph &graph) {
#if KASSERT_ENABLED(KASSERT_ASSERTION_LEVEL_NORMAL)
  GlobalNodeWeight check = 0;
  graph.for_each_visible_node([&](NodeID node) {
    if (solution.node_status[node] == wis_status::INCLUDED) {
      check += graph.get_weight(node);
    }
  });
  KASSERT(check == solution.solution_weight);
#endif
}

inline bool is_excluded_or_included(const wis_status &status) {
  return status == wis_status::INCLUDED || status == wis_status::EXCLUDED;
}

inline std::string get_neighborhood_string(Solution &sol, dist_dynamic_graph &graph, NodeID node) {
  std::stringstream s;
  s << "neighborhood of " << node << ":" << std::endl;
  s << node << " " << graph.get_weight(node) << " " << graph.is_visible(node) << " "
    << static_cast<int>(sol.node_status[node]) << std::endl;
  for (auto target : graph[node]) {
    s << target << " " << graph.get_weight(target) << " " << graph.is_visible(target) << " "
      << static_cast<int>(sol.node_status[target]) << std::endl;
  }
  return s.str();
};

}  // namespace kadisredu::mwis::debug
