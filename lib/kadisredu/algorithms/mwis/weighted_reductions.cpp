//
// Created by jannickb on 6/24/24.
//

#include "kadisredu/algorithms/mwis/weighted_reductions.h"

#include "kadisredu/algorithms/mwis/debug.h"
#include "kadisredu/algorithms/mwis/kernel_utils.h"
#include "kadisredu/algorithms/mwis/reducer.h"
#include "kadisredu/tools/kamis_helper.h"
#include "kadisredu/tools/logger.h"
#include "kadisredu/tools/utils.h"

#include <branch_and_reduce_algorithm.h>
#include <graph_access.h>

namespace kadisredu::mwis {

bool neighborhood_folding::reduce(weighted_reduce_ptr algo) {
  KADISREDU_MEASURE_TIME_METHOD(reduce_timer);

  auto& status = algo->status_;
  auto& graph = status.graph;
  auto old_visible_nodes = graph.num_visible_nodes();

  if (status.cfg.optimize_candidates_order) {
    sort_visible_candidates_by_degree(marker.current(), graph);
  };

  for_each_unset_marked_node(status.solution, marker, [&](GlobalNodeID node) {
    KASSERT(!graph.is_ghost(node));

    auto weight_neigh = graph.get_weight_neigh(node);
    auto weight_v = graph.get_weight(node);

    if (weight_v >= weight_neigh) {
      algo->set(node, wis_status::INCLUDED);
      return;
    }

    if (weight_v >= weight_neigh) {
      return;
    }

    GlobalNodeWeight min_weight = graph.get_weight(
        *std::min_element(graph[node].begin(), graph[node].end(), [&](auto first, auto second) {
          return graph.get_weight(first) < graph.get_weight(second);
        }));

    if (weight_v <= weight_neigh - min_weight) {
      return;
    }

    auto is_border_node = [&](GlobalNodeID v) { return graph.is_active_border_node(v); };

    if (is_border_node(node)) {
      return;
    }

    if (std::any_of(graph[node].begin(), graph[node].end(), is_border_node)) {
      return;
    }

    auto& neighbors = algo->set_1;
    neighbors.clear();
    utils::add_to_set(neighbors, graph[node].begin(), graph[node].end());

    // check if neighbors form an independent set
    if (std::any_of(graph[node].begin(), graph[node].end(), [&](auto v) {
          return !utils::disjoint(neighbors, graph[v].begin(), graph[v].end());
        })) {
      return;
    }

    fold(algo, node, neighbors, weight_neigh);
  });

  return old_visible_nodes - graph.num_visible_nodes() > 0;
}

bool generalized_neighborhood_folding::reduce(weighted_reduce_ptr algo) {
  KADISREDU_MEASURE_TIME_METHOD(reduce_timer);

  auto& status = algo->status_;
  auto& graph = status.graph;
  auto old_visible_nodes = graph.num_visible_nodes();

  if (status.cfg.optimize_candidates_order) {
    sort_visible_candidates_by_degree(marker.current(), graph);
  }

  for_each_unset_marked_node(status.solution, marker, [&](GlobalNodeID node) {
    KASSERT(!graph.is_ghost(node));

    if (graph[node].deg() > status.cfg.max_allowed_subproblem_size) {
      return;
    }

    auto weight_neigh = graph.get_weight_neigh(node);
    auto weight_v = graph.get_weight(node);

    if (weight_v >= weight_neigh) {
      algo->set(node, wis_status::INCLUDED);
      return;
    }

    auto is_border_node = [&](NodeID v) { return graph.is_active_border_node(v); };

    if (is_border_node(node)) {
      return;
    }

    if (std::any_of(graph[node].begin(), graph[node].end(), is_border_node)) {
      return;
    }

    // find MWIS of G[N(v)]
    graph_access induced_subgraph;
    auto& mapping = algo->buffers[1];
    auto& neighbors = algo->buffers[2];
    auto& neighborhood_v = algo->set_1;
    auto& MWIS_set = algo->set_2;

    neighborhood_v.clear();
    utils::add_to_set(neighborhood_v, graph[node].begin(), graph[node].end());
    kamis_helper::build_induced_graph(neighborhood_v, std::span{graph[node]}, graph, mapping,
                                      induced_subgraph);

    MISConfig exact_mis_config = kamis_helper::create_mis_config();
    exact_mis_config.time_limit = graph[node].deg() / 10;
    branch_and_reduce_algorithm exact_solver(induced_subgraph, exact_mis_config, true);

    cout_handler::disable_cout();
    auto solved = exact_solver.run_branch_reduce();
    cout_handler::enable_cout();

    if (!solved) {
      return;
    }

    exact_mis_config.time_limit = std::max(2.0, graph[node].deg() / 10.0);

    GlobalNodeWeight MWIS_neigh = exact_solver.get_current_is_weight();
    if (weight_v >= MWIS_neigh) {
      algo->set(node, wis_status::INCLUDED);
      return;
    }

    exact_solver.apply_branch_reduce_solution(induced_subgraph);
    MWIS_set.clear();

    GlobalNodeWeight min_MWIS_neighbor_weight = MWIS_neigh;

    forall_nodes(induced_subgraph, v) {
      if (induced_subgraph.getPartitionIndex(v) == 1) {
        auto neighbor = graph[node][v];
        MWIS_set.add(neighbor);

        if (graph.get_weight(neighbor) < min_MWIS_neighbor_weight) {
          min_MWIS_neighbor_weight = status.graph.get_weight(neighbor);
        }
      }
    }
    endfor

        if (graph.get_weight(node) < MWIS_neigh - min_MWIS_neighbor_weight) {
      return;
    }

    bool check_failed = false;
    neighbors.clear();
    for (auto neighbor : graph[node]) {
      neighbors.push_back(neighbor);
    }
    for (auto neighbor : graph[node]) {
      if (!MWIS_set.contains(neighbor)) {
        continue;
      }

      using std::swap;
      swap(*std::find(neighbors.begin(), neighbors.end(), neighbor), neighbors.back());
      KASSERT(neighbors.back() == neighbor);
      neighbors.pop_back();
      KASSERT(neighbors.empty() || neighbors.back() != neighbor);
      neighborhood_v.remove(neighbor);

      kamis_helper::build_induced_graph(neighborhood_v, std::span{neighbors}, graph, mapping,
                                        induced_subgraph);

      branch_and_reduce_algorithm exact_solver2(induced_subgraph, exact_mis_config, true);

      cout_handler::disable_cout();
      solved = exact_solver2.run_branch_reduce();
      cout_handler::enable_cout();

      if (!solved) {
        KADISREDU_LOG << "exact_solver timeout for subgraph";
      }

      check_failed = !solved || exact_solver2.get_current_is_weight() >= graph.get_weight(node);

      neighbors.push_back(neighbor);
      neighborhood_v.add(neighbor);

      if (check_failed) {
        break;
      }
    }

    if (!check_failed) {
      fold(algo, node, MWIS_set, MWIS_neigh);
    }

    // @todo try to exclude some neighbors
  });

  return old_visible_nodes - graph.num_visible_nodes() > 0;
}
void generalized_neighborhood_folding::fold(weighted_reduce_ptr algo, NodeID main, fast_set& mwis,
                                            GlobalNodeWeight mwis_weight) {
  auto& status = algo->status_;
  auto& graph = status.graph;

  restore_data data;
  data.main_weight = graph.get_weight(main);
  data.MWIS_weight = mwis_weight;

  auto& nodes = data.nodes;
  nodes.main = main;

  data.main_neighbor_list = status.graph[main];

  for (auto neighbor : data.main_neighbor_list) {
    if (mwis.contains(neighbor))
      nodes.MWIS.push_back(neighbor);
    else
      algo->set(neighbor, wis_status::EXCLUDED, false);
  }

  for (unsigned i = nodes.MWIS.size() - 1; i >= 1; i--) {
    graph.hide_node(nodes.MWIS[i]);
    status.reduction_measurements[algo->active_reduction_index].reduced_local_nodes += 1;
    status.solution.node_status[nodes.MWIS[i]] = wis_status::FOLDED;
  }

  algo->set(nodes.MWIS[0], wis_status::FOLDED, false);

  data.main_neighbor_list = status.graph[main];

  status.sol_offset += data.main_weight;

  graph.set_weight(nodes.main, mwis_weight - data.main_weight);

  std::vector<NodeID> new_neighbors;
  auto& neighbors = mwis;
  neighbors.clear();
  neighbors.add(main);

  for (NodeID MWIS_node : nodes.MWIS) {
    std::vector<NodeID> node_vec;

    for (auto neighbor : status.graph[MWIS_node]) {
      if (neighbors.add(neighbor)) {
        new_neighbors.push_back(neighbor);

        status.graph.restore_and_replace_edge(neighbor, nodes.main);

        node_vec.push_back(neighbor);
      }
    }
    data.MWIS_node_vecs.push_back(std::move(node_vec));
  }

  status.graph.build_neighborhood(main, std::move(new_neighbors));
  applications.push_back(std::move(data));
  status.folded_queue.push_back(get_reduction_type());

  algo->mark_node(nodes.main);
  algo->mark_neigh(nodes.main);
}
void generalized_neighborhood_folding::apply(weighted_reduce_ptr algo, NodeID modified_node) {
  auto& status = algo->status_;
  auto nodes = applications.back().nodes;
  auto MWIS_weight = applications.back().MWIS_weight;
  auto main_status = status.solution.node_status[nodes.main];
  restore(algo, modified_node);

  if (main_status == wis_status::INCLUDED) {
    status.solution.node_status[nodes.main] = wis_status::EXCLUDED;

    for (auto node : nodes.MWIS) {
      status.solution.node_status[node] = wis_status::INCLUDED;
    }
  } else {
    status.solution.node_status[nodes.main] = wis_status::INCLUDED;

    for (auto node : nodes.MWIS) {
      status.solution.node_status[node] = wis_status::EXCLUDED;
    }
  }
  status.solution.solution_weight += status.graph.get_weight(nodes.main);
}

void generalized_neighborhood_folding::restore(weighted_reduce_ptr algo, NodeID modified_node) {
  auto& status = algo->status_;
  auto& graph = status.graph;
  auto& data = applications.back();

  graph.hide_node(data.nodes.main);
  graph.set_neighborhood(data.nodes.main, std::move(data.main_neighbor_list));
  KASSERT(graph[data.nodes.main].deg() == 0);
  KASSERT(!graph.is_visible(data.nodes.main));
  auto last_restored = status.graph.restore_last();
  KASSERT(last_restored == data.nodes.main);
  KASSERT(graph.is_visible(data.nodes.main));

  for (size_t i = 0; i < data.nodes.MWIS.size(); i++) {
    for (auto neighbor : data.MWIS_node_vecs[i]) {
      graph.replace_last_hidden_edge(neighbor, data.nodes.MWIS[i]);
    }

    status.solution.node_status[data.nodes.MWIS[i]] = wis_status::UNSET;
    last_restored = graph.restore_last();
    KASSERT(last_restored == data.nodes.MWIS[i]);
  }

  graph.set_weight(data.nodes.main, data.main_weight);
  status.sol_offset -= data.main_weight;

  applications.pop_back();
}

bool neighborhood_removal::reduce(weighted_reduce_ptr algo) {
  KADISREDU_MEASURE_TIME_METHOD(reduce_timer);
  auto& status = algo->status_;
  auto& graph = status.graph;
  auto old_visible_nodes = graph.num_visible_nodes();

  for_each_unset_marked_node(status.solution, marker, [&](NodeID node) {
    KASSERT(!graph.is_ghost(node));

    auto weight_neigh = graph.get_weight_neigh(node);
    auto weight_v = graph.get_weight(node);

    if (weight_v >= weight_neigh) {
      algo->set(node, wis_status::INCLUDED);
    }
  });

  return old_visible_nodes - graph.num_visible_nodes() > 0;
}

bool generalized_neighborhood_removal::reduce(weighted_reduce_ptr algo) {
  auto& status = algo->status_;
  auto& graph = status.graph;
  auto old_visible_nodes = graph.num_visible_nodes();

  if (status.cfg.optimize_candidates_order) {
    sort_visible_candidates_by_degree(marker.current(), graph);
  }

  for_each_unset_marked_node(status.solution, marker, [&](NodeID v) {
    KASSERT(!graph.is_ghost(v));

    if (graph[v].deg() > status.cfg.max_allowed_subproblem_size) {
      return;
    }

    auto weight_neigh = graph.get_weight_neigh(v);
    auto weight_v = graph.get_weight(v);

    if (weight_v >= weight_neigh) {
      algo->set(v, wis_status::INCLUDED);
      return;
    }

    // ensure that v has at least maximum weight
    for (auto u : graph[v]) {
      if (weight_v < graph.get_weight(u)) {
        return;
      }
    }

    return;

    // find MWIS of G[N(v)]
    graph_access induced_subgraph;
    auto& reverse_mapping = algo->buffers[0];
    auto& mapping = algo->buffers[1];
    auto& neighborhood_v = algo->set_1;
    neighborhood_v.clear();

    utils::add_to_set(neighborhood_v, graph[v].begin(), graph[v].end());
    kamis_helper::build_induced_graph(neighborhood_v, std::span{graph[v]}, graph, mapping,
                                      induced_subgraph);

    MISConfig exact_mis_config = kamis_helper::create_mis_config();
    exact_mis_config.time_limit = graph[v].deg() / 10;
    branch_and_reduce_algorithm exact_solver(induced_subgraph, exact_mis_config, true);

    cout_handler::disable_cout();
    auto solved = exact_solver.run_branch_reduce();
    cout_handler::enable_cout();

    if (!solved) {
      return;
    }

    if (weight_v >= exact_solver.get_current_is_weight()) {
      algo->set(v, wis_status::INCLUDED);
      return;
    }
  });

  return old_visible_nodes - graph.num_visible_nodes() > 0;
}

bool degree_one::reduce(weighted_reduce_ptr algo) {
  KADISREDU_MEASURE_TIME_METHOD(reduce_timer);
  using namespace kamping;

  auto& status = algo->status_;
  auto& graph = status.graph;
  auto old_visible_nodes = graph.num_visible_nodes();
  NodeID detached = 0;

  // exhaustively reduce inner degree-one vertices and border vertices that can be included
  while (!marker.current_marked_nodes.empty()) {
    KASSERT(marker.next_marked_nodes.empty());
    for_each_unset_marked_node(status.solution, marker, [&](NodeID node) {
      KASSERT(!graph.is_ghost(node));

      if (graph[node].deg() == 0) {
        algo->set(node, wis_status::INCLUDED);
      } else if (graph[node].deg() == 1) {
        auto neighbor = graph[node][0];
        auto weight_neigh = graph.get_weight(neighbor);
        auto weight_v = graph.get_weight(node);

        // degree one cases:
        if (weight_v >= weight_neigh) {
          algo->set(node, wis_status::INCLUDED);
        } else if (weight_v < weight_neigh && !graph.is_ghost(neighbor)) {
          this->fold(algo, node);
        } else {
          KASSERT(graph.is_ghost(neighbor));
          algo->set(node, wis_status::DETACHED);
          ++detached;
        }
      }
    });
    marker.next();
  }

  return old_visible_nodes - graph.num_visible_nodes() + detached > 0;
}
void degree_one::fold(weighted_reduce_ptr algo, NodeID degree_one_node) {
  auto& status = algo->status_;
  auto& graph = status.graph;
  KASSERT(graph[degree_one_node].deg() == 1);
  auto neighbor = graph[degree_one_node][0];
  auto weight_neigh = graph.get_weight(neighbor);
  auto weight_v = graph.get_weight(degree_one_node);

  KASSERT(weight_neigh > weight_v);
  algo->set_new_weight(neighbor, weight_neigh - weight_v);
  status.sol_offset += weight_v;
  applications.push_back({.offset = weight_v, .node = degree_one_node});
  algo->set(degree_one_node, wis_status::FOLDED, false);
  status.folded_queue.push_back(get_reduction_type());
}

void degree_one::restore(weighted_reduce_ptr algo, NodeID modified_node) {
  auto& status = algo->status_;
  auto& graph = status.graph;
  auto& rest_data = applications.back();

  auto restored_node = graph.restore_last();
  KASSERT(restored_node == modified_node);
  auto neighbor = status.graph[modified_node][0];
  graph.set_weight(neighbor, graph.get_weight(neighbor) + rest_data.offset);

  status.sol_offset -= rest_data.offset;
  applications.pop_back();
}
void degree_one::apply(weighted_reduce_ptr algo, NodeID modified_node) {
  auto& status = algo->status_;
  KASSERT(modified_node < status.graph.num_nodes() + status.graph.num_ghosts());
  KASSERT(status.solution.node_status[modified_node] == wis_status::FOLDED);
  KASSERT(status.graph[modified_node].deg() == 1, "not a degree-one fold");
  auto& rest_data = applications.back();
  KASSERT(rest_data.node == modified_node);
  auto& graph = status.graph;

  // some inner folded node
  auto neighbor = graph[modified_node][0];
  KASSERT(!status.graph.is_ghost(neighbor));
  KASSERT(status.solution.node_status[neighbor] == wis_status::INCLUDED ||
              status.solution.node_status[neighbor] == wis_status::EXCLUDED,
          debug::get_neighborhood_string(status.solution, graph, neighbor));
  bool include_modified_node = status.solution.node_status[neighbor] != wis_status::INCLUDED;

  if (!include_modified_node || !graph.is_ghost(modified_node)) {
    // only add offset if the included vertex is a local node
    status.solution.solution_weight += rest_data.offset;
  }

  restore(algo, modified_node);

  if (include_modified_node) {
    status.solution.node_status[modified_node] = wis_status::INCLUDED;
  } else {
    status.solution.node_status[modified_node] = wis_status::EXCLUDED;
  }

  if (graph.is_ghost(modified_node)) {
    algo->on_border_status_change(modified_node, status.solution.node_status[modified_node]);
  }
}
bool degree_one::reduce_detached_ghost(weighted_reduce_ptr algo, NodeID detached_ghost,
                                       sized_vector<NodeID>& detached_neigh) {
  auto& status = algo->status_;
  auto& graph = status.graph;
  auto& solution = status.solution;
  auto& node_status = solution.node_status;

  KASSERT(status.solution.node_status[detached_ghost] == wis_status::UNSET);

  KASSERT(algo->MAX_DETACH_DEGREE == 1, "higher degree are not yet supported");
  KASSERT(detached_neigh.size() <= algo->MAX_DETACH_DEGREE);

  if (graph[detached_ghost].deg() == 0) {
    if (!detached_neigh.empty()) {
      auto detached_neighbor = detached_neigh[0];
      KASSERT(graph[detached_neighbor][0] == detached_ghost);
      // tie breaking
      if (graph.get_weight(detached_neighbor) > graph.get_weight(detached_ghost) ||
          (graph.get_weight(detached_neighbor) == graph.get_weight(detached_ghost) &&
           algo->local_vertex_wins_tie(detached_neighbor, detached_ghost))) {
        algo->set(detached_ghost, wis_status::EXCLUDED, false);
        node_status[detached_neighbor] = wis_status::INCLUDED;
        solution.solution_weight += graph.get_weight(detached_neighbor);
      } else {
        algo->set(detached_ghost, wis_status::INCLUDED, false);
        node_status[detached_neighbor] = wis_status::EXCLUDED;
      }
    } else {
      algo->set(detached_ghost, wis_status::INCLUDED, false);
      algo->on_border_status_change(detached_ghost, wis_status::INCLUDED);
    }
    return true;
  } else if (graph[detached_ghost].deg() == 1) {
    KASSERT(detached_neigh.empty());
    auto neighbor = graph[detached_ghost][0];
    KASSERT(status.solution.node_status[neighbor] != wis_status::DETACHED);
    if (graph.get_weight(detached_ghost) >= graph.get_weight(neighbor)) {
      // exclude neighbor first, so that detached ghost receives update
      algo->set(detached_ghost, wis_status::INCLUDED, false);
      // the interface counterpart of the detached_ghost is hidden and requires some when an update
      algo->on_border_status_change(detached_ghost, wis_status::INCLUDED);
    } else {
      // fold
      fold(algo, detached_ghost);
    }
    return true;
  }
  return false;
}

bool v_shape::reduce(weighted_reduce_ptr algo) {
  KADISREDU_MEASURE_TIME_METHOD(reduce_timer);
  using std::swap;

  auto& status = algo->status_;
  auto& graph = status.graph;
  auto old_visible_nodes = graph.num_visible_nodes();

  for_each_unset_marked_node(status.solution, marker, [&](NodeID node) {
    KASSERT(!graph.is_ghost(node));

    if (graph[node].deg() == 2) {
      auto x = graph[node][0];
      auto y = graph[node][1];

      if (graph.get_weight(node) >= graph.get_weight(x) + graph.get_weight(y)) {
        algo->set(node, wis_status::INCLUDED);
        return;
      }

      if (graph[x].deg() > graph[y].deg()) {
        swap(x, y);
      }

      if (graph[x].deg() <= 1) {
        return;  // reduced by degree one reduction
      }

      if (graph.is_active_border_node(node) || graph.is_active_border_node(x) ||
          graph.is_active_border_node(y)) {
        return;  // reduction must act entirely local
      }
      if (std::find(graph[x].begin(), graph[x].end(), y) != graph[x].end()) {
        // triangle
        if (graph.get_weight(node) >= std::max(graph.get_weight(x), graph.get_weight(y))) {
          algo->set(node, wis_status::INCLUDED);
        }
        return;
      }

      if (graph.get_weight(x) > graph.get_weight(y)) {
        swap(x, y);
      }
      auto weight_node = graph.get_weight(node);
      // weight_x <= weight_y
      auto weight_x = graph.get_weight(x);
      auto weight_y = graph.get_weight(y);

      // v-shape (2 cases)
      if (weight_node >= weight_x) {
        if (weight_node < weight_y) {
          // fold mid
          fold_mid(algo, {.main = node, .rest = {x, y}});
        } else {
          KASSERT(weight_node < weight_x + weight_y);
          // fold max
          fold_max(algo, {.main = node, .rest = {x, y}});
        }
      }
    }
  });
  return old_visible_nodes - graph.num_visible_nodes() > 0;
}

void v_shape::restore(weighted_reduce_ptr algo, NodeID modified_node) {
  auto& status = algo->status_;
  auto& graph = status.graph;
  auto& rest = applications.back();

  if (rest.fold_case == 1) {  // mid
    KASSERT(!graph.is_visible(rest.nodes.main));
    KASSERT(graph.is_visible(rest.nodes.rest[0]));
    KASSERT(graph.is_visible(rest.nodes.rest[1]));
    graph.hide_node(rest.nodes.rest[0]);
    graph.set_neighborhood(rest.nodes.rest[0], std::move(rest.main_neighbor_list));
    KASSERT(!graph.is_visible(rest.nodes.rest[0]));
    auto last_restored = status.graph.restore_last();
    KASSERT(last_restored == rest.nodes.rest[0]);
    KASSERT(graph.is_visible(rest.nodes.rest[0]));

    for (auto neighbor : rest.node_vecs[1]) {
      status.graph.remove_hidden_edge(neighbor, rest.nodes.rest[0]);
    }

    status.solution.node_status[rest.nodes.rest[1]] = wis_status::UNSET;

    last_restored = graph.restore_last();
    KASSERT(last_restored == rest.nodes.main);

    graph.set_weight(rest.nodes.rest[1], graph.get_weight(rest.nodes.rest[1]) + rest.main_weight);
    status.sol_offset -= rest.main_weight;
  } else if (rest.fold_case == 2) {  // max
    KASSERT(graph.is_visible(rest.nodes.main));
    KASSERT(!graph.is_visible(rest.nodes.rest[0]));
    KASSERT(!graph.is_visible(rest.nodes.rest[1]));
    graph.hide_node(rest.nodes.main);
    graph.set_neighborhood(rest.nodes.main, std::move(rest.main_neighbor_list));
    KASSERT(graph[rest.nodes.main].deg() == 0);
    KASSERT(!graph.is_visible(rest.nodes.main));
    auto last_restored = status.graph.restore_last();
    KASSERT(last_restored == rest.nodes.main);
    KASSERT(graph.is_visible(rest.nodes.main));

    for (size_t i = 0; i < 2; i++) {
      for (auto neighbor : rest.node_vecs[i]) {
        graph.replace_last_hidden_edge(neighbor, rest.nodes.rest[i]);
      }

      status.solution.node_status[rest.nodes.rest[i]] = wis_status::UNSET;
      last_restored = graph.restore_last();
      KASSERT(last_restored == rest.nodes.rest[i]);
    }

    graph.set_weight(rest.nodes.main, rest.main_weight);
    status.sol_offset -= rest.main_weight;
  }

  applications.pop_back();
}
void v_shape::apply(weighted_reduce_ptr algo, NodeID modified_node) {
  auto& status = algo->status_;
  auto& graph = status.graph;
  auto& rest = applications.back();

  auto nodes = rest.nodes;
  auto main_status = status.solution.node_status[nodes.main];
  auto x_status = status.solution.node_status[nodes.rest[0]];
  auto y_status = status.solution.node_status[nodes.rest[1]];

  int fold_case = rest.fold_case;
  restore(algo, modified_node);

  if (fold_case == 1) {  // mid
    if (x_status == wis_status::INCLUDED) {
      // KASSERT(y_status == wis_status::INCLUDED, static_cast<int>(y_status));
      status.solution.node_status[nodes.main] = wis_status::EXCLUDED;
      status.solution.node_status[nodes.rest[0]] = wis_status::INCLUDED;
      status.solution.node_status[nodes.rest[1]] = wis_status::INCLUDED;
      if (y_status == wis_status::EXCLUDED) {
        // possible if current solution not maximal (eg. because kernel solution not optimal)
        status.solution.solution_weight += graph.get_weight(nodes.rest[1]);
      } else {
        status.solution.solution_weight += graph.get_weight(nodes.main);
      }
    } else if (y_status == wis_status::INCLUDED) {
      status.solution.node_status[nodes.main] = wis_status::EXCLUDED;
      status.solution.node_status[nodes.rest[0]] = wis_status::EXCLUDED;
      status.solution.node_status[nodes.rest[1]] = wis_status::INCLUDED;
      status.solution.solution_weight += graph.get_weight(nodes.main);
    } else {
      KASSERT(x_status == wis_status::EXCLUDED);
      KASSERT(y_status == wis_status::EXCLUDED);
      status.solution.node_status[nodes.main] = wis_status::INCLUDED;
      status.solution.node_status[nodes.rest[0]] = wis_status::EXCLUDED;
      status.solution.node_status[nodes.rest[1]] = wis_status::EXCLUDED;
      status.solution.solution_weight += graph.get_weight(nodes.main);
    }

  } else if (fold_case == 2) {  // max
    if (main_status == wis_status::INCLUDED) {
      status.solution.node_status[nodes.main] = wis_status::EXCLUDED;
      status.solution.node_status[nodes.rest[0]] = wis_status::INCLUDED;
      status.solution.node_status[nodes.rest[1]] = wis_status::INCLUDED;
    } else {
      KASSERT(main_status == wis_status::EXCLUDED);
      status.solution.node_status[nodes.main] = wis_status::INCLUDED;
      status.solution.node_status[nodes.rest[0]] = wis_status::EXCLUDED;
      status.solution.node_status[nodes.rest[1]] = wis_status::EXCLUDED;
    }
    status.solution.solution_weight += graph.get_weight(nodes.main);
  }
}
void v_shape::fold_max(weighted_reduce_ptr algo, const fold_nodes& nodes) {
  auto& status = algo->status_;
  auto& graph = status.graph;

  graph.hide_node(nodes.rest[1]);
  status.reduction_measurements[algo->active_reduction_index].reduced_local_nodes += 1;
  status.solution.node_status[nodes.rest[1]] = wis_status::FOLDED;
  algo->set(nodes.rest[0], wis_status::FOLDED, false);

  applications.push_back({
      .nodes = nodes,
      .main_weight = graph.get_weight(nodes.main),
      .fold_case = 2,
      .main_neighbor_list = std::move(status.graph[nodes.main])  // neighbors of main are hidden
  });
  auto& data = applications.back();

  status.sol_offset += data.main_weight;
  graph.set_weight(nodes.main, graph.get_weight(nodes.rest[1]) + graph.get_weight(nodes.rest[0]) -
                                   data.main_weight);

  std::vector<NodeID> new_neighbors;
  auto& neighbors = algo->set_1;
  neighbors.clear();
  neighbors.add(nodes.main);

  for (std::size_t i = 0; i < 2; i++) {
    for (auto neighbor : status.graph[nodes.rest[i]]) {
      if (neighbors.add(neighbor)) {
        new_neighbors.push_back(neighbor);

        status.graph.restore_and_replace_edge(neighbor, nodes.main);

        data.node_vecs[i].push_back(neighbor);
      }
    }
  }

  status.graph.build_neighborhood(nodes.main, std::move(new_neighbors));
  status.folded_queue.push_back(get_reduction_type());

  algo->mark_node(nodes.main);
  algo->mark_neigh(nodes.main);
}
void v_shape::fold_mid(weighted_reduce_ptr algo, const fold_nodes& nodes) {
  auto& status = algo->status_;
  auto& graph = status.graph;

  algo->set(nodes.main, wis_status::FOLDED, false);

  applications.push_back({.nodes = nodes,
                          .main_weight = graph.get_weight(nodes.main),
                          .fold_case = 1,
                          .main_neighbor_list = std::move(status.graph[nodes.rest[0]])});
  auto& data = applications.back();

  status.sol_offset += data.main_weight;
  graph.set_weight(nodes.rest[1], graph.get_weight(nodes.rest[1]) - data.main_weight);

  std::vector<NodeID> new_neighbors;
  auto& neighbors = algo->set_1;
  neighbors.clear();
  for (auto neighbor : data.main_neighbor_list) {
    neighbors.add(neighbor);
    new_neighbors.push_back(neighbor);
  }
  for (auto neighbor : status.graph[nodes.rest[1]]) {
    if (neighbors.add(neighbor)) {
      new_neighbors.push_back(neighbor);
      status.graph.add_edge(neighbor, nodes.rest[0]);
      data.node_vecs[1].push_back(neighbor);
    }
  }

  status.graph.build_neighborhood(nodes.rest[0], std::move(new_neighbors));
  status.folded_queue.push_back(get_reduction_type());

  algo->mark_node(nodes.rest[0]);
  algo->mark_node(nodes.rest[1]);
  algo->mark_neigh(nodes.rest[0]);
}

bool extended_domination::reduce(weighted_reduce_ptr algo) {
  KADISREDU_MEASURE_TIME_METHOD(reduce_timer);

  auto& status = algo->status_;
  auto& graph = status.graph;
  auto& neighbors = algo->set_1;
  auto old_visible_nodes = graph.num_visible_nodes();

  for_each_unset_marked_node(status.solution, marker, [&](NodeID v) {
    KASSERT(!graph.is_ghost(v));

    auto weight_neigh_v = graph.get_weight_neigh(v);
    auto weight_v = graph.get_weight(v);
    auto deg_v = graph[v].deg();

    // check for neighborhood removal
    if (weight_v >= weight_neigh_v) {
      algo->set(v, wis_status::INCLUDED);
      return;
    }

    neighbors.clear();
    for (auto target : graph[v]) {
      neighbors.add(target);
    }
    neighbors.add(v);

    // search for dominated
    for (auto dominated : graph[v]) {
      // dominated cannot be a ghost
      // dominated has no smaller weight
      // dominated has no larger degree (subset condition)
      if (graph.is_ghost(dominated) || graph[dominated].deg() > deg_v) {
        continue;
      }

      // deg(v) >= deg(dominated)
      // w(v) <= w(dominated)

      bool subset = true;
      for (auto target : graph[dominated]) {
        if (!neighbors.contains(target)) {
          subset = false;
          break;
        }
      }

      if (subset) {
        if (graph.get_weight(dominated) >= weight_v) {
          KASSERT(graph.get_weight(dominated) >= weight_v);
          algo->set(v, wis_status::EXCLUDED);
          break;
        } else {
          KASSERT(graph.get_weight(dominated) < weight_v);
          fold(algo, v, dominated);
          break;
        }
      }
    }
  });

  return old_visible_nodes - graph.num_visible_nodes() > 0;
}

void extended_domination::fold(weighted_reduce_ptr algo, NodeID v, NodeID u) {
  auto& status = algo->status_;
  auto& graph = status.graph;

  graph.hide_edge(v, u);
  graph.hide_edge(u, v);

  algo->set_new_weight(v, graph.get_weight(v) - graph.get_weight(u));

  algo->mark_node(u);
  algo->mark_neigh(u);

  applications.push_back(v);
  status.modified_queue.push_back(algo->MODIFIED_EDGE);
  status.folded_queue.push_back(get_reduction_type());
}

void extended_domination::restore(weighted_reduce_ptr algo, NodeID modified_node) {
  KASSERT(modified_node == algo->MODIFIED_EDGE);
  auto& status = algo->status_;
  NodeID v = applications.back();
  applications.pop_back();
  auto& graph = status.graph;

  graph.restore_last_hidden_edge(v);
  NodeID u = graph[v][graph[v].deg() - 1];
  graph.restore_last_hidden_edge(u);
  graph.set_weight(v, graph.get_weight(v) + graph.get_weight(u));
}

void extended_domination::apply(weighted_reduce_ptr algo, NodeID modified_node) {
  KASSERT(modified_node == algo->MODIFIED_EDGE);
  auto& status = algo->status_;
  auto& node_status = status.solution.node_status;
  NodeID v = applications.back();
  auto& graph = status.graph;

  restore(algo, modified_node);

  NodeID u = graph[v][graph[v].deg() - 1];
  if (node_status[v] == wis_status::INCLUDED) {
    node_status[u] = wis_status::EXCLUDED;
  }
}

bool domination::reduce(weighted_reduce_ptr algo) {
  KADISREDU_MEASURE_TIME_METHOD(reduce_timer);

  auto& status = algo->status_;
  auto& graph = status.graph;
  auto& neighbors = algo->set_1;
  auto old_visible_nodes = graph.num_visible_nodes();

  for_each_unset_marked_node(status.solution, marker, [&](NodeID v) {
    KASSERT(!graph.is_ghost(v));

    auto weight_neigh_v = graph.get_weight_neigh(v);
    auto weight_v = graph.get_weight(v);
    auto deg_v = graph[v].deg();

    // check for neighborhood removal
    if (weight_v >= weight_neigh_v) {
      algo->set(v, wis_status::INCLUDED);
      return;
    }

    neighbors.clear();
    for (auto target : graph[v]) {
      neighbors.add(target);
    }
    neighbors.add(v);

    // search for dominated
    for (auto dominated : graph[v]) {
      // dominated cannot be a ghost
      // dominated has no smaller weight
      // dominated has no larger degree (subset condition)
      if (graph.is_ghost(dominated) || graph[dominated].deg() > deg_v ||
          graph.get_weight(dominated) < weight_v) {
        continue;
      }

      // deg(v) >= deg(dominated)
      // w(v) <= w(dominated)

      bool subset = true;
      for (auto target : graph[dominated]) {
        if (!neighbors.contains(target)) {
          subset = false;
          break;
        }
      }

      if (subset) {
        KASSERT(graph.get_weight(dominated) >= weight_v);
        algo->set(v, wis_status::EXCLUDED);
        break;
      }
    }
  });

  return old_visible_nodes - graph.num_visible_nodes() > 0;
}

bool basic_single_edge::reduce(weighted_reduce_ptr algo) {
  KADISREDU_MEASURE_TIME_METHOD(reduce_timer);

  auto& status = algo->status_;
  auto& graph = status.graph;
  auto& neighbors = algo->set_1;
  auto old_visible_nodes = graph.num_visible_nodes();

  for_each_unset_marked_node(status.solution, marker, [&](NodeID v) {
    KASSERT(!graph.is_ghost(v));

    auto weight_neigh_v = graph.get_weight_neigh(v);
    auto weight_v = graph.get_weight(v);
    auto deg_v = graph[v].deg();

    // check for neighborhood removal
    if (weight_v >= weight_neigh_v) {
      algo->set(v, wis_status::INCLUDED);
      return;
    }

    neighbors.clear();
    utils::add_to_set(neighbors, graph[v].begin(), graph[v].end());

    // search for dominated
    for (auto u : graph[v]) {
      // dominated cannot be a ghost
      // dominated has no smaller weight
      // dominated has no larger degree (subset condition)
      GlobalNodeWeight weight_u = graph.get_weight(u);
      if (graph.is_ghost(u) || weight_u < weight_v) {
        continue;
      }

      GlobalNodeWeight upper = weight_v;  ////< upper bound for MWIS of G[N(u)\N(v)]
      for (auto target : graph[u]) {
        if (!neighbors.contains(target)) {
          upper += graph.get_weight(target);
          if (upper > weight_u) {
            break;  ////< reduction not applicable
          }
        }
      }

      if (upper <= weight_u) {
        algo->set(v, wis_status::EXCLUDED);
        break;
      }
    }
  });

  return old_visible_nodes - graph.num_visible_nodes() > 0;
}
bool extended_single_edge::reduce(weighted_reduce_ptr algo) {
  KADISREDU_MEASURE_TIME_METHOD(reduce_timer);

  auto& status = algo->status_;
  auto& graph = status.graph;
  auto& neighbors = algo->set_1;
  auto& intersection = algo->buffers[0];
  auto& neighbor_buf_v = algo->buffers[1];
  auto old_visible_nodes = graph.num_visible_nodes();

  for_each_unset_marked_node(status.solution, marker, [&](NodeID v) {
    KASSERT(!graph.is_ghost(v));

    auto weight_neigh_v = graph.get_weight_neigh(v);
    auto weight_v = graph.get_weight(v);
    auto deg_v = graph[v].deg();

    // check for neighborhood removal
    if (weight_v >= weight_neigh_v) {
      algo->set(v, wis_status::INCLUDED);
      return;
    }

    neighbors.clear();
    neighbor_buf_v.clear();
    for (auto target : graph[v]) {
      neighbors.add(target);
      neighbor_buf_v.push_back(target);
    }

    // search for dominated
    for (auto u : neighbor_buf_v) {
      if (!graph.is_visible(u)) {
        continue;  ////< we already excluded u in some intersection
      }
      // dominated cannot be a ghost
      // dominated has no smaller weight
      // dominated has no larger degree (subset condition)
      GlobalNodeWeight weight_u = graph.get_weight(u);
      if (graph.is_ghost(u) || weight_neigh_v - weight_u > weight_v) {
        continue;
      }

      intersection.clear();
      for (auto target : graph[u]) {
        if (neighbors.contains(target) && !graph.is_ghost(target)) {
          intersection.push_back(target);
          weight_neigh_v -= graph.get_weight(target);
        }
      }

      for (auto target : intersection) {
        algo->set(target, wis_status::EXCLUDED, false);
      }

      if (weight_v >= weight_neigh_v) {
        algo->set(v, wis_status::INCLUDED, true);
        return;
      }
    }
  });

  return old_visible_nodes - graph.num_visible_nodes() > 0;
}

bool clique::reduce(weighted_reduce_ptr algo) {
  KADISREDU_MEASURE_TIME_METHOD(reduce_timer);

  auto& status = algo->status_;
  auto& graph = status.graph;
  auto& neighbors = algo->set_1;
  auto old_visible_nodes = graph.num_visible_nodes();

  if (status.cfg.optimize_candidates_order) {
    sort_visible_candidates_by_degree(marker.current(), graph);
  }

  for_each_unset_marked_node(status.solution, marker, [&](NodeID v) {
    KASSERT(!graph.is_ghost(v));
    KASSERT(graph.is_visible(v));

    auto weight_neigh_v = graph.get_weight_neigh(v);
    auto weight_v = graph.get_weight(v);

    // check for neighborhood removal
    if (weight_v >= weight_neigh_v) {
      algo->set(v, wis_status::INCLUDED);
      return;
    }

    auto deg_v = graph[v].deg();
    auto max_weight_v = graph.get_max_weight_neigh(v);
    auto max_weight = std::max(weight_v, max_weight_v);

    if (max_weight_v > weight_v) {
      return;
    }

    // check degree and ghosts
    unsigned ghost_counter = 0;
    for (auto target : graph[v]) {
      if (graph[target].deg() < deg_v) {
        return;
      }
      if (graph.is_ghost(target)) {
        if (++ghost_counter > 1) {
          return;
        }
      }
    }

    neighbors.clear();
    utils::add_to_set(neighbors, graph[v].begin(), graph[v].end());
    neighbors.add(v);

    // check intersection of neighborhoods
    for (auto target : graph[v]) {
      auto intersection_counter =
          utils::intersection(neighbors, graph[target].begin(), graph[target].end());
      if (intersection_counter < deg_v) {
        return;  // no clique
      }
    }

    algo->include_node_with_bulk_hide(v, neighbors);
    return;
  });

  return old_visible_nodes - graph.num_visible_nodes() > 0;
}

bool simplicial_weight_transfer::reduce(weighted_reduce_ptr algo) {
  KADISREDU_MEASURE_TIME_METHOD(reduce_timer);
  auto& status = algo->status_;
  auto& graph = status.graph;
  auto& neighbors = algo->set_1;
  auto old_visible_nodes = graph.num_visible_nodes();

  if (status.cfg.optimize_candidates_order) {
    sort_visible_candidates_by_degree(marker.current(), graph);
  }

  for_each_unset_marked_node(status.solution, marker, [&](NodeID v) {
    KASSERT(!graph.is_ghost(v));

    auto weight_neigh_v = graph.get_weight_neigh(v);
    auto weight_v = graph.get_weight(v);

    // check for neighborhood removal
    if (weight_v >= weight_neigh_v) {
      algo->set(v, wis_status::INCLUDED);
      return;
    }

    auto deg_v = graph[v].deg();
    auto max_weight_v = graph.get_max_weight_neigh(v);
    auto max_weight = std::max(weight_v, max_weight_v);

    // check degree and ghosts
    unsigned ghost_counter = 0;
    NodeWeight ghost_weight = 0;
    NodeID ghost;
    for (auto target : graph[v]) {
      if (graph[target].deg() < deg_v) {
        return;
      }
      if (graph.is_ghost(target)) {
        if (ghost_counter > 0) {  // TODO one ghost
          return;
        }
        ghost_weight = graph.get_weight(target);
        ghost = target;
        ++ghost_counter;
      }
    }

    neighbors.clear();
    utils::add_to_set(neighbors, graph[v].begin(), graph[v].end());
    neighbors.add(v);

    NodeID max_simplicial_node = v;
    GlobalNodeWeight max_simplicial_weight = weight_v;

    // check intersection of neighborhoods
    for (auto target : graph[v]) {
      auto intersection_counter =
          utils::intersection(neighbors, graph[target].begin(), graph[target].end());
      if (intersection_counter < deg_v) {
        return;  // no clique
      } else if (!graph.is_ghost(target) && graph[target].deg() == deg_v) {
        // another simplicial non-ghost vertex
        KASSERT(intersection_counter == deg_v);
        if (graph.get_weight(target) > max_simplicial_weight) {
          max_simplicial_node = target;
          max_simplicial_weight = graph.get_weight(target);
        }
      }
    }

    if (max_simplicial_weight == max_weight) {
      algo->include_node_with_bulk_hide(max_simplicial_node, neighbors);
      return;
    }

    if (ghost_counter > 0) {
      return;
    }

    this->fold(algo, max_simplicial_node);
  });

  return old_visible_nodes - graph.num_visible_nodes() > 0;
}

void simplicial_weight_transfer::fold(weighted_reduce_ptr algo, NodeID max_simplicial_node) {
  auto& status = algo->status_;
  auto& graph = status.graph;

  auto& exclude = algo->buffers[0];
  exclude.clear();
  auto& exclude_set = algo->set_1;
  exclude_set.clear();

  auto max_simplicial_weight = status.graph.get_weight(max_simplicial_node);

  std::for_each(graph[max_simplicial_node].begin(), graph[max_simplicial_node].end(),
                [&](auto node) {
                  if (graph.get_weight(node) <= max_simplicial_weight) {
                    exclude.push_back(node);
                    exclude_set.add(node);
                  }
                });

  algo->bulk_exclude(exclude, exclude_set, false);
  algo->set(max_simplicial_node, wis_status::FOLDED, false);

  status.sol_offset += max_simplicial_weight;
  for (auto target : status.graph[max_simplicial_node]) {
    KASSERT(!graph.is_ghost(target));
    KASSERT(graph.get_weight(target) > max_simplicial_weight);
    algo->set_new_weight(target, graph.get_weight(target) - max_simplicial_weight);
  }
  applications.push_back(
      {.offset = max_simplicial_weight /*, .non_simplicials = std::move(non_simplicials)*/});
  status.folded_queue.push_back(get_reduction_type());
}

void simplicial_weight_transfer::restore(weighted_reduce_ptr algo, NodeID modified_node) {
  auto& status = algo->status_;
  auto& graph = status.graph;
  auto& rest_data = applications.back();

  for (auto target : graph[modified_node]) {
    // if (!graph.is_ghost(target)) {
    graph.set_weight(target, graph.get_weight(target) + rest_data.offset);
    //}
  }
  auto restored_node = graph.restore_last();
  KASSERT(restored_node == modified_node,
          "last restored node is not the current folded node that should be restored");
  status.sol_offset -= rest_data.offset;
  applications.pop_back();
}
void simplicial_weight_transfer::apply(weighted_reduce_ptr algo, NodeID modified_node) {
  auto& status = algo->status_;
  KASSERT(status.solution.node_status[modified_node] == wis_status::FOLDED);
  auto& rest_data = applications.back();
  auto& graph = status.graph;

  bool set_simplicial = true;
  for (auto target : graph[modified_node]) {
    KASSERT(status.solution.node_status[target] == wis_status::INCLUDED ||
            status.solution.node_status[target] == wis_status::EXCLUDED);
    if (status.solution.node_status[target] == wis_status::INCLUDED) {
      set_simplicial = false;
      break;
    }
  }
  status.solution.solution_weight += rest_data.offset;

  restore(algo, modified_node);

  if (set_simplicial) {
    status.solution.node_status[modified_node] = wis_status::INCLUDED;
  } else {
    status.solution.node_status[modified_node] = wis_status::EXCLUDED;
  }
}

}  // namespace kadisredu::mwis
