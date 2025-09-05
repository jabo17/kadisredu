//
// Created by jannickb on 11/27/24.
//

#include "kadisredu/graphutils/graph_summary.h"

#include <kamping/collectives/reduce.hpp>
#include <kamping/communicator.hpp>

namespace kadisredu {

std::optional<graph_summary> graph_summary::agg_graph_summary_on_root(
    const kamping::Communicator<> &comm, dist_dynamic_graph &graph) {
  using namespace kamping;

  struct summary_message {
    GlobalNodeID nodes;
    NodeID max_nodes;
    NodeID min_nodes;
    GlobalEdgeID edges;
    GlobalNodeID ghosts;
    NodeID max_ghosts;
    NodeID min_ghosts;
    GlobalEdgeID cut;
    NodeID max_degree;
    NodeID min_degree;
  };

  summary_message local = {
      .nodes = graph.num_visible_nodes(),
      .max_nodes = graph.num_visible_nodes(),
      .min_nodes = graph.num_visible_nodes(),
      .edges = 0,
      .ghosts = graph.num_visible_ghosts(),
      .max_ghosts = graph.num_visible_ghosts(),
      .min_ghosts = graph.num_visible_ghosts(),
      .cut = 0,
      .max_degree = 0,
      .min_degree = std::numeric_limits<NodeID>::max(),
  };

  graph.for_each_visible_node([&g = graph, &local](NodeID node) {
    local.edges += g[node].deg();
    local.max_degree = std::max(local.max_degree, g[node].deg());
    local.min_degree = std::min(local.min_degree, g[node].deg());
  });
  graph.for_each_visible_ghost([&g = graph, &local](NodeID node) { local.cut += g[node].deg(); });

  auto res = comm.reduce_single(send_buf(local),
                                op(
                                    [](auto &d1, auto &d2) {
                                      return summary_message{
                                          .nodes = d1.nodes + d2.nodes,
                                          .max_nodes = std::max(d1.max_nodes, d2.max_nodes),
                                          .min_nodes = std::min(d1.min_nodes, d2.min_nodes),
                                          .edges = d1.edges + d2.edges,
                                          .ghosts = d1.ghosts + d2.ghosts,
                                          .max_ghosts = std::max(d1.max_ghosts, d2.max_ghosts),
                                          .min_ghosts = std::min(d1.min_ghosts, d2.min_ghosts),
                                          .cut = d1.cut + d2.cut,
                                          .max_degree = std::max(d1.max_degree, d2.max_degree),
                                          .min_degree = std::min(d1.min_degree, d2.min_degree)};
                                    },
                                    ops::commutative));

  std::vector<NodeID> degree_counts(MAX_DEGREE + 1, 0);
  graph.for_each_visible_node([&g = graph, &degree_counts](NodeID node) {
    if (g[node].deg() <= MAX_DEGREE) {
      degree_counts[g[node].deg()]++;
    }
  });
  auto global_degree_counts = comm.reduce(send_buf(degree_counts), op(ops::plus<>()));

  if (comm.is_root()) {
    KASSERT(res.has_value());
    auto &data = res.value();

    return graph_summary{
        .nodes = data.nodes,
        .max_nodes = data.max_nodes,
        .min_nodes = data.min_nodes,
        .edges = data.edges / 2,
        .ghosts = data.ghosts,
        .max_ghosts = data.max_ghosts,
        .min_ghosts = data.min_ghosts,
        .cut = data.cut / 2,
        .imbalance = static_cast<double>(data.max_nodes) /
                     (static_cast<double>(data.nodes) / static_cast<double>(graph.blocks())),
        .max_degree = data.max_degree,
        .min_degree = data.min_degree,
        .degree_counts = std::move(global_degree_counts)};
  }
  return std::nullopt;
}

}  // namespace kadisredu
