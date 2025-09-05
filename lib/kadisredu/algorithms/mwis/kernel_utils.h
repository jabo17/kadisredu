//
// Created by jannickb on 11/21/24.
//

#pragma once

#include "mwis_def.h"
#include "kadisredu/data_structures/sized_vector.h"
#include "kadisredu/data_structures/dist_dynamic_graph.h"

#include <algorithm>

namespace kadisredu::mwis {

void inline sort_visible_candidates_by_degree(sized_vector<NodeID>& candidates,
                                              dist_dynamic_graph& graph) {
  // filter out non-visible candidates
  candidates.pop_back(candidates.end() -
                      std::remove_if(candidates.begin(), candidates.end(),
                                     [&](auto& v) { return !graph.is_visible(v); }));
  // sort remaining candidates ascending by their degree
  std::sort(candidates.begin(), candidates.end(),
            [&](auto& first, auto& second) { return graph[first].deg() < graph[second].deg(); });
}

}